#include "common.h"
#include "ib.h"



int ib_finalise(ib_t *ib){
  // particles
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    common_free(ib->particles[n]);
  }
  common_free(ib->particles);
  // Euler variables
  common_free(ib->dux);
  common_free(ib->duy);
  common_free(ib->duz);
  // main structure
  common_free(ib);
  return 0;
}

