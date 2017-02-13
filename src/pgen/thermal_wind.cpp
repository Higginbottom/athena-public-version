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

using namespace std;

// Global arrays to store midplane BC
static Real* d_bc;
static Real* M3_bc;
static Real* e_bc;

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{  
// Read problem parameters
  Real Gamma = peos->GetGamma();
  Real gm1 = (Gamma - 1.0); 
  Real GM = pin->GetOrAddReal("problem","GM",1.); // activates gravity
  Real c_s = pin->GetOrAddReal("problem","c_s",1.); 
  Real n_g = pin->GetOrAddReal("problem","n_g",1.); 
  Real font_case = pin->GetOrAddInteger("problem","font_case",0); 

  // Initial conditions for pressure 
  Real p0 = SQR(c_s)*n_g; // reference pressure
  Real e0 = p0/gm1; // reference energy
 
  // Allocate memory for disk BC
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  d_bc  = new Real[nx1];
  M3_bc = new Real[nx1];
  e_bc  = new Real[nx1];


  Real r;
  
  //Font cases: density prescriptions along the disk
  for (int i=is; i<=ie; i++) 
  {
    r = pcoord->x1v(i);
    
    if (font_case == 0)
      phydro->u(IDN,ks,je,i) = n_g*pow(r,-2.); // eq (9), alpha = 2
    else if (font_case == 1)
      phydro->u(IDN,ks,je,i) = n_g*pow(r,-1.5); // eq (9), alpha = 3/2
    else if (font_case == 2)
      phydro->u(IDN,ks,je,i) = n_g*pow(r,-2.5); // eq (9), alpha = 5/2
	else if (font_case == 3)
      phydro->u(IDN,ks,je,i) = n_g*pow(2.0/(pow(r,7.5)+pow(r,12.5)),0.2); // "Photo evaporative disc model eq (10)"
  
    // set keplerian rotation
    phydro->u(IM3,ks,je,i) = phydro->u(IDN,ks,je,i)*pow(GM/r,0.5);
    
    // set total energy
    phydro->u(IEN,ks,je,i) = p0*pow(phydro->u(IDN,ks,je,i),Gamma)/gm1 + 0.5*SQR(phydro->u(IM3,ks,je,i))/phydro->u(IDN,ks,je,i); 
    
    // Save the boundary state
    d_bc[i] = phydro->u(IDN,ks,je,i);
    M3_bc[i] = phydro->u(IM3,ks,je,i);
    e_bc[i] = p0*pow(phydro->u(IDN,ks,je,i),Gamma)/gm1;
  }
  
  // Initialize all other active zones 
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=(je-1); j++) {
  for (int i=is; i<=ie; i++) {
    r = pcoord->x1v(i);
    phydro->u(IDN,k,j,i) = 1e-4*n_g;
    phydro->u(IEN,k,j,i) = p0*pow(phydro->u(IDN,k,j,i),Gamma)/gm1;
    phydro->u(IM1,k,j,i) = 0.0;
    phydro->u(IM2,k,j,i) = 0.0;
    phydro->u(IM3,k,j,i) = 0.0;   
  }}}


  // Print out ICs along the midplane
  Real d,v,p,xi,T;
  for (int k=ks; k<=ke; ++k) {
  for (int j=je; j<=je; ++j) {  
    for (int i=is; i<=ie; ++i) 
    {
      r = pcoord->x1v(i);

      d = phydro->u(IDN,k,j,i);
      v = phydro->u(IM3,k,j,i)/phydro->u(IDN,k,j,i);
      p = gm1*(phydro->u(IEN,k,j,i) - 0.5*phydro->u(IDN,k,j,i)*SQR(v));

      printf("ALONG DISK: i = %2d Radius = %.4f  (d,v,p) = (%.4e,%.4e,%.4e) \n",i,r,d,v,p);
    }     
  }}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
void MeshBlock::UserWorkInLoop(void)
{
  Real kinetic_energy;
  
  // hold initial density/temperature along midplane fixed and enforce Keplerian rotation
  for (int k = ks; k <=ke; ++k)
    for (int i = is; i <=ie; ++i)
      {
        phydro->u(IDN,k,je,i) = d_bc[i];
        phydro->u(IM1,k,je,i) = 0.;
        phydro->u(IM3,k,je,i) = M3_bc[i];
        kinetic_energy = 0.5*(SQR(phydro->u(IM2,k,je,i)) + SQR(phydro->u(IM3,k,je,i)))/d_bc[i];
        phydro->u(IEN,k,je,i) = e_bc[i] + kinetic_energy; 
      }

  return;
}



