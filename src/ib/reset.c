#include "ib.h"



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

