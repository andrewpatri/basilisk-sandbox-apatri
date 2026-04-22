/* This file is set to rappresent the experimental setup from https://www.sciencedirect.com/science/article/pii/S2352710224001268?via%3Dihub*/

// define section
#define NO_ADVECTION_DIV 1 // no divergene term in advection
#define SOLVE_TEMPERATURE 1 // solve also temperature equation ( entalpic balance)
#define RADIATION_INTERFACE 0.9 // valore dell'emissività
#define MOLAR_DIFFUSION 1 // use molar diffusion instead of mass
#define FICK_CORRECTED 1 // enable fick correction for multicomponent
#define MASS_DIFFUSION_ENTHALPY 1 // enable enthalpic contribution to mass  flux 

// ifndef section, here i would implement the possibility to change the q source, but need to be refined since it is a power 
// at this stage T_ENV so i can test the geometry 
#ifndef T_ENV
# define T_ENV 773
#endif
 /* ifndef a_q
 # define a_q 20
#endif */
/* ifndef b_q
 # define b_q 20
#endif */

// include section

// in 2d no axi.h since no viable axial simmetry is found and also then considering the gravity in 3d ok like in the photos of wood unburnt burnt etc
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h" // properties
#include "two-phase.h"// for VOF and defines fraction and colors
#include "shrinking.h" // shriniking models ATTENTION choose the good one remember to define
#include "multicomponent-varprop.h" // for all the define for the properties
#include "darcy.h" // to take into account flow resistance due to porosity
#include "view.h"

// dati base
double Uin = 0.; //no velocity in exp
double tend = 600.; //simulation time.  If need to compute temperature for its test condition before  the insertion to measure the heat flux

// Boundary condition
u.n[left]    = dirichlet (Uin);
u.t[left]    = dirichlet (0.);
p[left]      = neumann (0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

u.n[top]      = dirichlet (0.);
u.t[top]      = dirichlet (0.);
p[top]        = neumann (0.);
pf[top]       = neumann (0.);
psi[top]      = dirichlet (0.);
  		


u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[bottom]      = neumann (0.);
u.t[bottom]      = neumann (0.);
p[bottom]        = neumann (0.);
pf[bottom]       = neumann (0.);
psi[bottom]      = neumann (0.);

int maxlevel = 7; int minlevel = 2; // risoluzione minima e massima 128 o 4 celle epr lato
double H0 = 2e-2; // initially
double solid_mass0 = 0., moisture0 = 0.; // massa della fase solida iniziale, contenuto di umidità iniziale

