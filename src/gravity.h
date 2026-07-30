coord G = {0.,0.,0.}, Z = {0.,0.,0.};

#include "curvature.h"

event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

event acceleration (i++)
{
  coord G1;
  foreach_dimension()
    G1.x = G.x;

  scalar phig[];
  position (f, phig, G1, Z, add = false);

#if TREE
  for (scalar f in {interfaces})
    set_prolongation (f, p.prolongation);
#endif

  scalar rhovar[];
  for (scalar f in {interfaces})
    foreach()
      rhovar[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);

#if TREE
  set_prolongation (rhovar, p.prolongation);
#endif

  face vector av = a, sth[];
  foreach_face() {
    sth.x[] = 1.;
    if (f[] != f[-1] && fm.x[] > 0.) {
      double phif =
        (phig[] < nodata && phig[-1] < nodata) ?
        (phig[] + phig[-1])/2. :
        phig[] < nodata ? phig[] :
        phig[-1] < nodata ? phig[-1] :
        0.;
      sth.x[] = (phif == 0.) ? 1. : 0.;

      av.x[] -= alpha.x[]/(fm.x[] + SEPS)*phif*(rhovar[] - rhovar[-1])/Delta;
    }
  }

  /**
  Far from interfacial faces, we use the coordinates of the centroids
  of the cells intead of the interface centroids. */

  foreach_face() {
    coord o = {x,y,z};
    double phiof = 0.;
    foreach_dimension()
      phiof += (o.x - Z.x)*G1.x;
    phiof *= sth.x[];

    av.x[] -= alpha.x[]/(fm.x[] + SEPS)*phiof*(rhovar[] - rhovar[-1])/Delta;
  }

#if TREE
  for (scalar f in {interfaces})
    set_prolongation (f, fraction_refine);
#endif
}