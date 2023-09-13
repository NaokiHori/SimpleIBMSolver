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
#if NDIMS == 3
    p->dz = 0.;
#endif
    p->dux = 0.;
    p->duy = 0.;
#if NDIMS == 3
    p->duz = 0.;
#endif
#if NDIMS == 3
    p->dvx = 0.;
    p->dvy = 0.;
#endif
    p->dvz = 0.;
  }
  return 0;
}

static int reset_arrays(
    ib_t * ib
){
  array_t ibfx = ib->ibfx;
  array_t ibfy = ib->ibfy;
#if NDIMS == 3
  array_t ibfz = ib->ibfz;
#endif
  memset(ibfx.data, 0, ibfx.datasize);
  memset(ibfy.data, 0, ibfy.datasize);
#if NDIMS == 3
  memset(ibfz.data, 0, ibfz.datasize);
#endif
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

