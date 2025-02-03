#include <string.h>
#include "ib.h"
#include "ib_solver.h"

static int reset_particles(
    ib_t * ib
){
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    p->dx = 0.;
    p->dy = 0.;
    p->dux = 0.;
    p->duy = 0.;
    p->dvz = 0.;
  }
  return 0;
}

static int reset_arrays(
    ib_t * ib
){
  array_t ibfx = ib->ibfx;
  array_t ibfy = ib->ibfy;
  memset(ibfx.data, 0, ibfx.datasize);
  memset(ibfy.data, 0, ibfy.datasize);
  return 0;
}

int ib_reset_variables(
    ib_t * ib
){
  if(0 != reset_particles(ib)){
    return 1;
  }
  if(0 != reset_arrays(ib)){
    return 1;
  }
  return 0;
}

