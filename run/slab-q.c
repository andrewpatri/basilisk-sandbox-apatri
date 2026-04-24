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
#ifndef a_q
# define a_q 0.28
#endif
#ifndef b_q
# define b_q 1.79
#endif 

// include section

#include "axi.h"
#include "navier-stokes/centered-phasechange.h"
#include "opensmoke-properties.h" // properties
#include "two-phase.h"// for VOF and defines fraction and colors
#include "shrinking.h" // shriniking models ATTENTION choose the good one remember to define
#include "multicomponent-varprop.h" // for all the define for the properties MAYBE CHANGE TO USE MODIFIED DATA
#include "darcy.h" // to take into account flow resistance due to porosity
#include "view.h"
#include "superquadric.h"
// dati base
double Uin = 0.; //no velocity in exp
double tend = 600.; //simulation time.  If need to compute temperature for its test condition before  the insertion to measure the heat flux

// Boundary condition
u.n[right]    = neumann (0.);
u.t[right]    = neumann (0.);
p[right]      = dirichlet (0.);
pf[right]     = dirichlet (0.);
psi[right]    = dirichlet (0.);

u.n[top]      = neumann (0.);
u.t[top]      = neumann (0.);
p[top]        = dirichlet (0.);
pf[top]       = dirichlet (0.);
psi[top]      = dirichlet (0.);

u.n[left]    = dirichlet (0.);
u.t[left]    = dirichlet (0.);
p[left]      = neumann(0.);
pf[left]     = neumann (0.);
psi[left]    = dirichlet (0.);

/*u.n[bottom]      = dirichlet (0.);
u.t[bottom]      = dirichlet (0.);
p[bottom]        = neumann (0.);
pf[bottom]       = neumann(0.);
psi[bottom]      = neumann (0.);
*/

#define q_time(a,b,t)(a*pow(t,b))
// q sorg 
double q_sorg; 
//# define t_time(T,dt)(T*dT) 
int maxlevel = 7; int minlevel = 2; // risoluzione minima e massima 128 o 4 celle epr lato
double H0 = 2e-2; // initially
double solid_mass0 = 0., moisture0 = 0.; // massa della fase solida iniziale, contenuto di umidità iniziale
//double Temperatura_daupdate = T_ENV;

int main() {
  
  lambdaS = 0.1987; // on the pubblication another value used but it's coming from optimization of their parameters
  lambdaSmodel = L_HUANG;
  TS0 = 300.; TG0 = T_ENV; // the change this to 300 K
  rhoS = 720;  // kg/m3
  eps0 = 0.39;

  //dummy properties
  rho1 = 1., rho2 = 1.;
  mu1 = 1., mu2 = 1.;

  zeta_policy = ZETA_CONST;


L0 = H0*4; // first try just the block of wood
  

origin (0, 0);

  DT = 1e-1;

  shift_prod = true;
  kinfolder = "biomass/dummy-solid";
  init_grid(1 << maxlevel);
  run();
}

#define rectangle(x,y,H0)()
//#define circle(x,y,R)(sq(R) - sq(x) - sq(y))
/*//#define unione(a,b) ((a) > (b) ? (a) : (b))

// Quadrato 1: centro (0., 0.), lato H0
#define phi1(x,y,H0)(min(H0/2-fabs(x),-fabs(y)+H0/2))

// Quadrato 2: centro (0, H0), lato 0.2
#define phi2(x,y,H0)(max( fabs(x), fabs(y - H0)-H0))*/
event init(i=0) {

  // fraction(f, circle(x,y,H0));
  fraction (f, superquadric(x, y, 20, H0, 2*H0));
  // fraction (f,phi1(x,y,H0));
  gas_start[OpenSMOKE_IndexOfSpecies ("N2")] = 1.;
  gas_start[OpenSMOKE_IndexOfSpecies ("TAR")] = 0.;
  gas_start[OpenSMOKE_IndexOfSpecies ("H2O")] = 0.;
  
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("BIOMASS")]  = 1;
  sol_start[OpenSMOKE_IndexOfSolidSpecies ("CHAR")]  = 0;

  foreach()
    porosity[] = eps0*f[];
	
  solid_mass0 = 0.;
  foreach (reduction(+:solid_mass0))
    solid_mass0 += (f[]-porosity[])*rhoS*dv(); //Note: (1-e) = (1-ef)!= (1-e)f
	
  fprintf(stderr, "DEBUG = %g\n", solid_mass0);
  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    if (jj == OpenSMOKE_IndexOfSpecies ("N2")) { // change when adding also 02
      YG[right] = dirichlet (1.);
    } else {
      YG[right] = dirichlet (0.);
    }
  }

  foreach()
    u.x[] = f[] > F_ERR ? 0. : Uin;
// Temperature

    TG[right] = dirichlet(TG0);// to use when then setting T0 in gas phase
    TG[left] = neumann (0.);
    TG[top] =  neumann(0.);
  //  TG[bottom] = neumann(0.);
 
    //q_sorg = q_time(a_q,b_q,t);
   //fprintf(stderr, "DEBUG q_sorginizio =%g\n", q_sorg);
}
/*event T_update(t += 0.5){
	Temperatura_daupdate = Temperatura_daupdate*1.01;
	TG[right] = Temperatura_daupdate;
}*/

/*event movie(t += 1){
clear();
cells();
view();
squares("T", spread=-1,linear = true, min =TS0, max = TG0*1.2);
draw_vof("f");

save("T.mp4");
}*/
event output (t += 1) {
  //fprintf (stderr, "%g\n", t);

  char name[80];
  sprintf(name, "OutputData_T1-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  //log mass profile
  double solid_mass = 0.;
  foreach (reduction(+:solid_mass))
    solid_mass += (f[]-porosity[])*rhoS*dv();
 //fprintf (stderr, "DEBUG solid_mass = %g\n", solid_mass);
 
 
/ T on the surface
double T_surf = 0.;
T_surf = TInt[x,0];
fprintf (stderr, "DEBUG Tsurfi= %g\n", T_surf);
//average temperature of the surface
  double Tsurf_avg = 0.; 
  int count = 0;
  foreach(reduction(+:Tsurf_avg)reduction(+:count)) {
  	 if (f[] > F_ERR && f[] < 1.-F_ERR) {
   	Tsurf_avg += TInt[];
    	  count++;
   	}
     }
    Tsurf_avg /= count;
fprintf (stderr, "DEBUG Tsurf= %g\n", Tsurf_avg);


  double Tcore  = interpolate (T, 0., 0.);
fprintf (stderr, "DEBUG Tcore= %g\n", Tcore);

  double Th2    = interpolate (T, H0/2, 0);
fprintf (stderr, "DEBUG Th2= %g\n", Th2);
 
 q_sorg = q_time(a_q,b_q,t);
 fprintf(stderr, "DEBUG q_sorginizio =%g\n", q_sorg);

  fprintf (fp, "%g %g\n", 
            t, solid_mass/solid_mass0);//,// Tcore, Th2, Tsurf_avg, q_sorg); 
            // radius/(D0/2.)  r/r0);

  fflush(fp);
}

#if TREE
event adapt (i++) {
  scalar inert = YGList_G[OpenSMOKE_IndexOfSpecies ("N2")];
  adapt_wavelet_leave_interface ({T, u.x, u.y, inert}, {f},
    (double[]){1.e0, 1.e-1, 1.e-1, 1e-1}, maxlevel, minlevel, 2);
}
#endif

event stop (t = tend);
