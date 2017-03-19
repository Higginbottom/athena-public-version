//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//#define DEBUG_SOLUTION // prints steps in root finder
#define STIFF_TRICK    // guard against H/C crashes
//#define XSTAR_HC       // use XSTAR tables instead of analytic Blondin

#include <math.h>
#define MAXIT 100
#define TOL 1e-10  			//tolerance for root finding

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

// User code
#ifdef XSTAR_HC
  #include "../user/xstar_table.h"  // XSTAR TABLES
#endif

// HP15 UNITS
/***********************************/
#define T0 114087.95981974274
#define T_C 14000000.0
#define xi_disk 125.89254117941674
#define L_x 3.3e37
#define tdyn_over_tth 1.0
/***********************************/

#define k_B 1.38e-16;
#define m_p 1.67e-24;

#define Gc_fac 8.890467707253008e-36 //used in H/C models below
//#define T_x 1.16e8  // 10 keV 
#define T_x 5.6e7  // 4.826 keV used in HP15 

double blondin_h(double T,double xi) 
{
	/* Cooling rate only */
	double Gc,Gx;
	/* Compton */
	Gc = Gc_fac * xi * (T_x - 4.*T);
	/* photo heating */
	Gx = 1.5e-21 * pow(xi,0.25) / sqrt(T) * (1.- T/T_x);

	return (Gc + Gx);
}

double blondin_c(double T,double xi) 
{
	/* Cooling rate only */
	double Ll,Lb;
	double T_L = 1.3e5;
	// brems
	Lb = 3.3e-27 * sqrt(T);
	/* line cooling */
	Ll = 1.7e-18 * exp(-T_L/T) / xi / sqrt(T) + 1e-24;
	return (Lb + Ll);
}

double blondin_hc(double T,double xi) 
{
	/* Net cooling rate */
	double Gc,Gx,Ll,Lb;
	double T_L = 1.3e5;
	/* Compton */
	Gc  = Gc_fac * xi * (T_x - 4.*T);
	/* photo heating */
	Gx = 1.5e-21 * pow(xi,0.25) / sqrt(T) * (1.- T/T_x);
	/* bremss, this takes into account the prefactors for case C */
	Lb = 3.3e-27 * sqrt(T);
	/* line cooling */
	Ll = 1.7e-18 * exp(-T_L/T) / xi / sqrt(T) + 1e-24;
	return (Lb + Ll - Gc - Gx);
}


using namespace std;

static Real Gamma, gm1, GM, rho_0,R_IC;


void HeatingAndCooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

// this function is also called for restarts so put 
// any pertinent restart functionality here
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  //Enroll H/C source term
    Gamma = pin->GetReal("hydro","gamma");
    gm1 = (pin->GetReal("hydro","gamma") - 1.0); 
    GM = pin->GetOrAddReal("problem","GM",1.); // activates gravity
    rho_0 = pin->GetOrAddReal("problem","rho_0",1.);
    R_IC = pin->GetOrAddReal("problem","R_IC",1.);
	

	
  EnrollUserExplicitSourceFunction(HeatingAndCooling);
  
  return;
}

Real solve4pressure(const Real dt, const Real d, const Real p, const Real r, const int i, const Real theta, const int j);

int bracket_root(const Real dt, const Real d, const Real p_old, const Real r, Real* pL, Real* pR, int);

Real polish_root(const Real dt, const Real d, const Real p_old, const Real r, const Real pL, const Real pR, const int i);

Real zbrent(const Real dt, const Real d, const Real p_old, const Real r, const Real pL, const Real pR, const int i);

Real root_bisection(const Real d, const Real r);

// Global Variables to store BC


static Real hc_norm;

#ifdef XSTAR_HC  // use tables generated from XSTAR
  XSTAR_LUT hc_xstar("/Users/atwaters/Documents/Athena++/nov-athena-master/src/user/Blondin_125by200.tab","bicubic");
  XSTAR_LUT h_xstar("/Users/atwaters/Documents/Athena++/nov-athena-master/src/user/Blondin_125by200_h.tab","bicubic");
  
  // c7 SED
  //XSTAR_LUT hc_xstar("/Users/atwaters/Documents/Athena++/nov-athena-master/src/user/hc_datafiles/c7/tim/sed_c7_net.dat","bicubic");
  //XSTAR_LUT h_xstar("/Users/atwaters/Documents/Athena++/nov-athena-master/src/user/hc_datafiles/c7/tim/sed_c7_heating.dat","bicubic");
#endif 

Real tabularHC(const Real d, const Real p, const Real r)
{ 
//  Real T = T0*(Gamma*p/d);
	
	Real n_h=d/1.67e-24/1.41;
	Real n_e=n_h*1.21;
  
  Real T = 0.6*1.67e-24*p/1.38e-16/d;
  Real xi = L_x/r/r/n_h;
  Real hc_rate;
  
//    if (d>1e-23)
//	  {
//		    printf ("p=%e d=%e xi=%e T=%e\n",p,d,xi,T);
//		}
  
 // cout << "r = " << r << "T = " << T << "xi = " << xi << endl;
#ifdef XSTAR_HC
  hc_rate = tdyn_over_tth*d*d*hc_xstar.get_rate(T,xi)*hc_norm;  // tabular 
#else
  hc_rate = n_e*n_h*blondin_hc(T,xi); // analytic Blondin
//   printf ("T=%e T0=%e Gamma=%e gm1=%e p=%e d=%e xi=%e blondin_hc=%e\n",T,T0,Gamma,gm1,p,d,xi,blondin_hc(T,xi));
#endif
  
//      if (d>8e-10)
//	  	  {
//			  	    printf ("d=%e xi=%e T=%e hc=%e\n",d,xi,T,hc_rate/n_e/n_h);
//							}

  
  return hc_rate; 
}


void MeshBlock::ProblemGenerator(ParameterInput *pin)
{  
// Read problem parameters

  Real r,theta,x,z;
  Real d_atmsph,p_atmsph,p_disk;
  Real deg = 180./PI;
  Real e_disk,e_atmsph;
  Real c_s;
  


  
#ifdef XSTAR_HC  
  hc_xstar.initialize_table();
  h_xstar.initialize_table();
  hc_norm = 1./h_xstar.get_rate(T0,xi_disk);
#else
  hc_norm = 1./blondin_h(T0,xi_disk);
#endif 
  e_disk=3./2.*rho_0*1.38e-16*T0/0.6/1.67e-24;
  
    c_s=sqrt(Gamma*gm1*e_disk/rho_0);
  
printf ("e_disk=%e speed of sound at disk=%e time scale=%e\n",e_disk,c_s,R_IC/c_s/c_s);
  
  // Midplane values
  for (int i=is; i<=ie+2; i++) 
  {
    r = pcoord->x1v(i);
    
    // density profile
    phydro->u(IDN,ks,je,i) = rho_0*pow(r/R_IC,-2.);
  
    // set keplerian rotation
    phydro->u(IM3,ks,je,i) = phydro->u(IDN,ks,je,i)*pow(GM/r,0.5);
    
    // pressure corresponds to making temp T0 in code units
//    p_disk = phydro->u(IDN,ks,je,i)/Gamma; // T0 = Gamma*p_disk/d 
    //p_disk = root_bisection(phydro->u(IDN,ks,je,i),r); // equivalent: solve for T0=T_eq given xi_disk
    e_disk=3./2.*phydro->u(IDN,ks,je,i)*1.38e-16*T0/0.6/1.67e-24;
	
	
	
    // total energy
    phydro->u(IEN,ks,je,i) = e_disk + 0.5*SQR(phydro->u(IM3,ks,je,i))/phydro->u(IDN,ks,je,i);
    
  }
  
  // Initialize all other active zones 
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=(je-1); j++) {
  for (int i=is; i<=ie+2; i++) {
    d_atmsph = 2.67475745e-18; // 5e-7 still no floors hit, 1e-7 hits floors
//    p_atmsph = (T_C/T0)*d_atmsph/Gamma; 
	
	e_atmsph=3./2.*d_atmsph*1.38e-16*T_C/0.6/1.67e-24;
    
    phydro->u(IDN,k,j,i) = d_atmsph;
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;   
    phydro->u(IEN,k,j,i) = e_atmsph + 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
  }}}
  
  
  // Print out ICs
  Real d,vt,vp,p,xi,T;
  for (int k=ks; k<=ke; ++k) {
  for (int j=je; j<=je; ++j) {  
    for (int i=is; i<=ie; ++i) {
      r = pcoord->x1v(i);
    
      d = phydro->u(IDN,k,j,i);
      vt = phydro->u(IM2,k,j,i)/phydro->u(IDN,k,j,i);
      vp = phydro->u(IM3,k,j,i)/phydro->u(IDN,k,j,i);
      p = gm1*(phydro->u(IEN,k,j,i) - 0.5*phydro->u(IDN,k,j,i)*(SQR(vt) + SQR(vp)));
	  T = 0.6*1.67e-24*p/1.38e-16/d;
      xi = L_x*1.67e-24*1.41/r/r/d;
	        printf("ALONG DISK: i = %2d r(R_IC) = %.4f  (d,M,E) = (%.4e,%.4e,%.4e)  (xi,T) = (%.4e,%.4e) \n",
			      i,r,phydro->u(IDN,k,j,i),phydro->u(IM3,k,j,i),phydro->u(IEN,k,j,i),xi,T);
      }
      
  }}
  
  for (int k=ks; k<=ke; ++k) 
  for (int j=js; j<=js; ++j)  
  for (int i=is; i<=ie; ++i)
  {
      r = pcoord->x1v(i);
      theta = pcoord->x2v(j);
      d = phydro->u(IDN,k,j,i);
	  
      x = r*sin(theta);
      z = r*cos(theta);
      p = gm1*(phydro->u(IEN,k,j,i) - 0.5*phydro->u(IDN,k,j,i)*(SQR(vt) + SQR(vp)));
	  
	  T = 0.6*1.67e-24*p/1.38e-16/d;
      xi = L_x*1.67e-24*1.41/r/r/d;
      printf("ALONG POLE: i = %2d x = %.4f z = %.4f (d,M,E) = (%.4e,%.4e,%.4e) (xi,T) = (%.4e,%.4e) \n",
      i,x,z,phydro->u(IDN,k,j,i),phydro->u(IM3,k,j,i),phydro->u(IEN,k,j,i),xi,T);     
  }
  
  for (int k=ks; k<=ke; ++k) 
  for (int j=js; j<=je; ++j)  
  for (int i=is; i<=is; ++i)
  {
      r = pcoord->x1v(i);
      theta = pcoord->x2v(j);
      d = phydro->u(IDN,k,j,i);
	  
      x = r*sin(theta);
      z = r*cos(theta);
      p = gm1*(phydro->u(IEN,k,j,i) - 0.5*phydro->u(IDN,k,j,i)*(SQR(vt) + SQR(vp)));
	  
	  T = 0.6*1.67e-24*p/1.38e-16/d;
      xi = L_x*1.67e-24*1.41/r/r/d;
	  printf("ALONG R_IN: j = %2d r = %.4f theta = %.4f (d,M,E) = (%.4e,%.4e,%.4e) (xi,T) = (%.4e,%.4e) \n",
	  j,r,theta*deg,phydro->u(IDN,k,j,i),phydro->u(IM3,k,j,i),phydro->u(IEN,k,j,i),xi,T);     
  }
  
  for (int k=ks; k<=ke; ++k) 
  for (int j=js; j<=je; ++j)  
  for (int i=ie; i<=ie; ++i)
  {
      r = pcoord->x1v(i);
      theta = pcoord->x2v(j);
      d = phydro->u(IDN,k,j,i);
	  
      x = r*sin(theta);
      z = r*cos(theta);
      p = gm1*(phydro->u(IEN,k,j,i) - 0.5*phydro->u(IDN,k,j,i)*(SQR(vt) + SQR(vp)));
	  
	  T = 0.6*1.67e-24*p/1.38e-16/d;
      xi = L_x*1.67e-24*1.41/r/r/d;
	  printf("ALONG R_OUT: j = %2d r = %.4f theta = %.4f (d,M,E) = (%.4e,%.4e,%.4e) (xi,T) = (%.4e,%.4e) \n",
	  j,r,theta*deg,phydro->u(IDN,k,j,i),phydro->u(IM3,k,j,i),phydro->u(IEN,k,j,i),xi,T);     
  }

   cout << "[ProblemGenerator]: hc_norm = " << hc_norm << endl;
   
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
void MeshBlock::UserWorkInLoop(void)
{
  Real kinetic_energy,p_disk,e_disk;
  
  for (int k = ks; k <=ke; ++k)
    for (int i = is; i <=ie; ++i)
      {
		  
// 		 if (phydro->u(IDN,k,je,i)>8e-10)
 //			 printf("old density=%e %e new density=%e\n",phydro->u(IDN,k,je,i),phydro->u(IDN,k,je-1,i),rho_0*pow(pcoord->x1v(i)/R_IC,-2.));
		  
        phydro->u(IDN,k,je,i) = rho_0*pow(pcoord->x1v(i)/R_IC,-2.);

		phydro->u(IM1,k,je,i) = 0.;
		
		
		
        phydro->u(IM3,k,je,i) = phydro->u(IDN,k,je,i)*pow(GM/pcoord->x1v(i),0.5);
		
		
		
		        kinetic_energy = 0.5*(SQR(phydro->u(IM2,k,je,i)) + SQR(phydro->u(IM3,k,je,i)))/phydro->u(IDN,k,je,i);
		
	    p_disk = phydro->u(IDN,k,je,i)/Gamma; // T0 = Gamma*p_disk/d 
		
			    e_disk=3./2.*phydro->u(IDN,k,je,i)*1.38e-16*T0/0.6/1.67e-24;
		
		
		      phydro->u(IEN,k,je,i) = e_disk + kinetic_energy; 
      }
   
  /*  Verify that disk stays at 10^4 K 
  Real Ts,Te;
  Ts = phydro->u(IEN,ks,je,is) - 0.5*( SQR(phydro->u(IM1,ks,je,is)) + SQR(phydro->u(IM2,ks,je,is)) + SQR(phydro->u(IM3,ks,je,is)) )/phydro->u(IDN,ks,je,is);
  Te = phydro->u(IEN,ks,je,ie) - 0.5*( SQR(phydro->u(IM1,ks,je,ie)) + SQR(phydro->u(IM2,ks,je,ie)) + SQR(phydro->u(IM3,ks,je,ie)) )/phydro->u(IDN,ks,je,ie);
  Ts = T0*Gamma*(Ts*gm1)/phydro->u(IDN,ks,je,is); //
  Te = T0*Gamma*(Te*gm1)/phydro->u(IDN,ks,je,ie);
  printf("is: den = %E M2 = %E M3 = %E E = %E e/d = %E T = %E\n",phydro->u(IDN,ks,je,is),phydro->u(IM2,ks,je,is),phydro->u(IM3,ks,je,is),phydro->u(IEN,ks,je,is),e_bc[is]/phydro->u(IDN,ks,je,is),Ts);
  printf("ie: den = %E M2 = %E M3 = %E E = %E e/d = %E T = %E\n",phydro->u(IDN,ks,je,ie),phydro->u(IM2,ks,je,ie),phydro->u(IM3,ks,je,ie),phydro->u(IEN,ks,je,ie),e_bc[ie]/phydro->u(IDN,ks,je,ie),Te);
 */

  return;
}


void HeatingAndCooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real dens,pres,radius,theta,RHS;
  Coordinates *pco = pmb->pcoord;
  Real xi,T,HC,dtrhoL;
  Real kinetic_energy;
  Real internal_energy;
  Real sign,pres_floor;
  Real RHSmax;
  Real new_temp;
  Real v,p_root;
  
  
     if (pmb->pmy_mesh->ncycle % 100 == 0) 
		 	    cout << "cycle = " << pmb->pmy_mesh->ncycle << " dt = " << dt << endl;
  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
    //Real phi=pmb->pcoord->x2v(j) ;
    for (int i=pmb->is; i<=pmb->ie; ++i) {
      // update energy equation 
      dens = prim(IDN,k,j,i);
      pres = prim(IEN,k,j,i);
      radius = pco->x1v(i);
      theta = pco->x2v(j);

/* Implicit method */
  p_root = solve4pressure(dt,dens,pres,radius,i,theta,j);
#ifdef DEBUG_SOLUTION
  cout << i << ": p_old = " << pres << "; p_new = " << p_root << endl;
#endif
  RHS = dt*tabularHC(dens,p_root,radius);
  dtrhoL = RHS;

  cons(IEN,k,j,i) -= RHS;
  
  
  kinetic_energy = 0.5*(SQR(prim(IM2,k,j,i)) + SQR(prim(IM3,k,j,i)))*dens;
  
  internal_energy=cons(IEN,k,j,i)-kinetic_energy;
	  
	  new_temp=internal_energy*(2./3.)*0.6*m_p;
	  new_temp=new_temp/dens/k_B;
  
//  if (j==pmb->je-10 )
//	  printf ("r=%e dens=%e old pressure=%e root=%e old temp=%e new temp=%e\n",radius,dens,pres,p_root,0.6*1.67e-24*pres/1.38e-16/dens,new_temp);
  
  
#ifdef DEBUG_SOLUTION
  cout << "Etot = " << cons(IEN,k,j,i) << "RHS = " << RHS << endl;
#endif     
     //if (pmb->pmy_mesh->ncycle % 10000 == 0) //10000 == 0)
     if (j==pmb->js || j==pmb->je)
      if (i==pmb->is || i==pmb->ie)
      {
        xi = xi_disk/(cons(IDN,k,j,i)*SQR(radius));
        T = (Gamma*pres/dens);
        HC = tabularHC(dens,pres,radius)/tdyn_over_tth/SQR(dens);
#ifdef VERBOSE
        (j==pmb->js) ? printf("POLE:") : printf("DISK:");
        printf("i = %2d r = %E den = %E vel = %E pres = %E dtrhoL = %E RHS = %E xi = %E T = %E HC = %E\n",i,radius,prim(IDN,k,j,i),prim(IM1,k,j,i),prim(IEN,k,j,i),dtrhoL,RHS,xi,T,HC);
#endif
	    
      }
    } // i
  }} // k j
  
}

Real solve4pressure(const Real dt, const Real d, const Real p, const Real r, const int i, const Real theta, const int j)
{
   //first bracket the root
   Real nroots,pL,pR,p_root;
   // just 2 brackets suffices most of the time
   
//    printf ("Got to solve4pressure, d=%e, p=%e, r=%e, dt=%e tdyn_over_tth=%e\n",d,p,r,dt,tdyn_over_tth);
   
   
   nroots = bracket_root(dt,d,p,r,&pL,&pR,2);
   if (nroots == 0)  // try 100 brackets
   {
     //cout << "i = " << i << " nroots = " << nroots << " with 2 brackets." << endl;
     //cout << "Attempting 100 brackets..." << endl;
     nroots = bracket_root(dt,d,p,r,&pL,&pR,100);
    }
   if (nroots == 0)  // try 1000 brackets
   {
     //cout << "i = " << i << " nroots = " << nroots << " with 100 brackets." << endl;
     //cout << "Attempting 1,000 brackets..." << endl;
     nroots = bracket_root(dt,d,p,r,&pL,&pR,1000);
     //cout << "Result: nroots = " << nroots << endl;
    }
   if (nroots == 0)
   {
     cout << "i = " << i << " nroots = " << nroots << " with 1,000 brackets." << endl;
     cout << "Attempting 10,000 brackets..." << endl;
     nroots = bracket_root(dt,d,p,r,&pL,&pR,10000);
     cout << "Result: nroots = " << nroots << endl;
    }
   if (nroots == 0)
   {
     cout << "i = " << i << " nroots = " << nroots << " with 10,000 brackets." << endl;
     cout << "Attempting 100,000 brackets..." << endl;
     nroots = bracket_root(dt,d,p,r,&pL,&pR,100000);
     cout << "Result: nroots = " << nroots << endl;
     stringstream msg;
     msg << "Implicit solve for pressure failed!" << endl;
     msg << "Setting p_new = p_old at r[" << i << "] = " << r << "theta[" << j << "] = " << theta << endl;
//     if (nroots == 0)
       //p_root = p; // equation is very stiff
//       throw runtime_error(msg.str().c_str());
    }
  
#ifdef STIFF_TRICK  
   if (nroots == -1)
     p_root = p; // equation is very stiff
   else  //polish the root
#endif
     p_root = zbrent(dt,d,p,r,pL,pR,i);		// Brent's method
     //p_root = polish_root(dt,d,p,r,pL,pR,i); // bisection method
   
   return p_root;
}

int bracket_root(const Real dt, const Real d, const Real p_old, const Real r, Real* pL, Real* pR, int Nbrak)
{
  // define bracket range
  Real p_min,p_max;
  
#ifdef XSTAR_HC 
  p_min = (hc_xstar.T_min/T0)*d/Gamma + TOL;
  p_max = (hc_xstar.T_max/T0)*d/Gamma - TOL;
  

  if (Nbrak > 100)
  {
    p_min = 1e-1*p_min;
    p_max = 1e1*p_max;
  }
  
  if (Nbrak > 1000)
  {
    p_min = 1e-2*p_min;
    p_max = p_max;
  }
#else
  p_min = 1e-1*p_old;
  p_max = 1e1*p_old;
//  if (d>8e-10)
//	  printf ("p_min=%e p_max=%e NBrak=%i\n",p_min,p_max,Nbrak);
  
  

  if (Nbrak > 100)
  {
    p_min = 5e-2*p_old;
    p_max = 2e1*p_old;
  }
  if (Nbrak > 1000)
  { 
    p_min = 1e-2*p_min;
    p_max = p_max;
  }
  if (Nbrak > 100000)
  { 
    p_min = 1e-2*p_old;
    p_max = p_old;
  }
#endif  
  
  //cout << "[bracket_root]: p_min = " << p_min << " p_max = " << p_max << endl;
  
  // define values used by algorithm 
  Real pl,pr,dp,fl,fr;
  int n_ctr = 0; // current root tally
  int Nmax = 2;  // maximum number of roots allowed
  dp = (p_max-p_min)/Nbrak;
  Real ps_L[2], ps_R[2]; 
  Real p_correction;
  
  // algorithm: search all brackets for roots from left to right
  pl = p_min; pr = p_min;
  p_correction = d*dt*gm1*tabularHC(d,pl,r);
//  if (d>8e-10)
//	  printf ("p_correction=%e\n",p_correction);
  
  fl = pl - p_old + p_correction;
#ifdef DEBUG_SOLUTION
  printf("r = %E den = %E pl = %E H/C = %E\n",r,d,pl,tabularHC(d,pl,r));
  cout << "[bracket root]: p_old = " << p_old << " correction = " << d*dt*gm1*tabularHC(d,pl,r) << endl;
#endif

#ifdef STIFF_TRICK
  if (fabs(p_correction) > 0.2*p_old) n_ctr = -1;
  else
#endif
  for (int j=1; j<Nbrak; j++)
  {
      pr += dp;
      fr = pr - p_old + d*dt*gm1*tabularHC(d,pr,r);
      if (fl*fr <= 0.0)
      {
		ps_L[++n_ctr-1] = pr-dp;
		ps_R[n_ctr - 1] = pr;
		if (n_ctr == Nmax) break; 
  	  }
      fl = fr;
  }
  
  // record result
  *pL = ps_L[0];
  *pR = ps_R[0];

  #ifdef DEBUG_SOLUTION
  cout << "nroots found: " << n_ctr << endl;
  int n;
  for (n = 0; n < n_ctr; n++) 
    cout << "bracket " << n << ": [pL,pR] = [" << ps_L[n] << ", " << ps_R[n] << "]" << endl;
  #endif
  
  return n_ctr;

}  
  
Real polish_root(const Real dt, const Real d, const Real p_old, const Real r, const Real pL, const Real pR, const int i)
{
  //bisection method
  double p_root,pr,pl,dp,fl,fr;

  dp = pR - pL;
  pl = pL;
  pr = pR;
  while (fabs(dp) > TOL) 
  {
	 p_root = 0.5*(pl+pr);
	 fl = pl - p_old + d*dt*gm1*tabularHC(d,pl,r);
	 fr = p_root - p_old + d*dt*gm1*tabularHC(d,p_root,r);
	 if (fl*fr < 0.)
	   pr  = p_root;
	 else
	   pl  = p_root;
     dp = pr-pl;
#ifdef DEBUG_SOLUTION
     if (i == 2) cout << "[inside bisection]  [" << pl << ", " << pr << "]" << endl;
#endif
  }
  return p_root;
}


#include <math.h>      
#define ITMAX 100     
   
//float zbrent(Real (*func)(Real), Real x1, Real x2, Real tol)  
Real zbrent(const Real dt, const Real dens, const Real p_old, const Real rad, const Real pL, const Real pR, const int i)
{   
    int iter;   
    Real a=pL,b=pR,c=pR,d,e,min1,min2;   
    //Real fa=(*func)(a)
    Real fa = a - p_old + dt*gm1*tabularHC(dens,a,rad);
    Real fb = b - p_old + dt*gm1*tabularHC(dens,b,rad);
    Real fc,p,q,r,s,tol1,xm;   
   
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))   
        cout << "Root must be bracketed in zbrent" << endl;   
    fc=fb;   
    for (iter=1;iter<=ITMAX;iter++) {   
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {   
            c=a;   
            fc=fa;   
            e=d=b-a;   
        }   
        if (fabs(fc) < fabs(fb)) {   
            a=b;   
            b=c;   
            c=a;   
            fa=fb;   
            fb=fc;   
            fc=fa;   
        }   
        tol1=2.0*TOL*fabs(b)+0.5*TOL;   
        xm=0.5*(c-b);   
        if (fabs(xm) <= tol1 || fb == 0.0) return b;   
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {   
            s=fb/fa;   
            if (a == c) {   
                p=2.0*xm*s;   
                q=1.0-s;   
            } else {   
                q=fa/fc;   
                r=fb/fc;   
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));   
                q=(q-1.0)*(r-1.0)*(s-1.0);   
            }   
            if (p > 0.0) q = -q;   
            p=fabs(p);   
            min1=3.0*xm*q-fabs(tol1*q);   
            min2=fabs(e*q);   
            if (2.0*p < (min1 < min2 ? min1 : min2)) {   
                e=d;   
                d=p/q;   
            } else {   
                d=xm;   
                e=d;   
            }   
        } else {   
            d=xm;   
            e=d;   
        }   
        a=b;   
        fa=fb;   
        if (fabs(d) > tol1)   
            b += d;   
        else   
            b += SIGN(xm)*fabs(tol1);   
        //fb=(*func)(b);   
        fb = b - p_old + dt*gm1*tabularHC(dens,b,rad);
#ifdef DEBUG_SOLUTION
     if (i == 2) cout << "[inside zbrent]  [" << a << ", " << b << "]" << endl;
#endif
    }   
    cout << "Maximum number of iterations exceeded in zbrent" << endl;   
    return 0.0;   
}   
#undef ITMAX        

Real root_bisection(const Real d, const Real r)
{
  // we seek the temperature, so we solve for pressure,
  // i.e. x is pressure
  double x1,x2,dx,root,product;
  int iter = 0;
  root = 0.;
   //loop through all brackets and find the roots
    x1 = 1e-8;
    x2 = 1e8;
    dx = x2-x1;
    while (fabs(dx) > TOL) 
    {
	  root = 0.5*(x1+x2);
	  product = tabularHC(d,x1,r)*tabularHC(d,root,r); 
	  if (product < 0.)
	  {
	    x2  = root;
	    dx = x2-x1;
	  }
	  else
	  {
	    x1  = root;
	    dx = x2-x1;
	  }
    }
    //if (root>TINY_NUMBER) 
    return root;
  
}
