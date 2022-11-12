#include "ib.h"


#if NDIMS == 2

int ib_reset_particle_increments(ib_t *ib){
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    ib->particles[n]->dux = 0.;
    ib->particles[n]->duy = 0.;
    ib->particles[n]->dvz = 0.;
    ib->particles[n]->dx  = 0.;
    ib->particles[n]->dy  = 0.;
  }
  return 0;
}

#else // NDIMS == 3

int ib_reset_particle_increments(ib_t *ib){
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    ib->particles[n]->dux = 0.;
    ib->particles[n]->duy = 0.;
    ib->particles[n]->duz = 0.;
    ib->particles[n]->dvx = 0.;
    ib->particles[n]->dvy = 0.;
    ib->particles[n]->dvz = 0.;
    ib->particles[n]->dx  = 0.;
    ib->particles[n]->dy  = 0.;
    ib->particles[n]->dz  = 0.;
  }
  return 0;
}

#endif // NDIMS
