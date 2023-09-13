#include "domain.h"
#include "ib.h"
#include "ib_solver.h"

static int update_particle_velocities(
    ib_t * ib
){
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    p->ux += p->dux;
    p->uy += p->duy;
#if NDIMS == 3
    p->uz += p->duz;
#endif
#if NDIMS == 3
    p->vx += p->dvx;
    p->vy += p->dvy;
#endif
    p->vz += p->dvz;
  }
  return 0;
}

static int update_particle_positions(
    const domain_t * domain,
    ib_t * ib
){
  const double ly = domain->lengths[1];
#if NDIMS == 3
  const double lz = domain->lengths[2];
#endif
  const size_t nitems = ib->nitems;
  particle_t * particles = ib->particles;
  for(size_t n = 0; n < nitems; n++){
    particle_t * p = particles + n;
    p->x += p->dx;
    p->y += p->dy;
#if NDIMS == 3
    p->z += p->dz;
#endif
    // correct periodicity
    if(p->y < 0.) p->y += ly;
    if(p->y > ly) p->y -= ly;
#if NDIMS == 3
    if(p->z < 0.) p->z += lz;
    if(p->z > lz) p->z -= lz;
#endif
  }
  return 0;
}

int ib_update_particles(
    const domain_t * domain,
    ib_t * ib
){
  update_particle_velocities(ib);
  update_particle_positions(domain, ib);
  return 0;
}

