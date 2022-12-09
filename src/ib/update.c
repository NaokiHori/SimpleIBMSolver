#include "common.h" // M_PI
#include "fluid.h"
#include "ib.h"
#include "arrays/fluid.h"
#include "arrays/ib.h"



static int update_particle_velocities(ib_t *ib){
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = ib->particles[n];
    p->ux += p->dux;
    p->uy += p->duy;
    p->uz += p->duz;
    p->vx += p->dvx;
    p->vy += p->dvy;
    p->vz += p->dvz;
  }
  return 0;
}

static int update_particle_positions(const domain_t *domain, ib_t *ib){
  const double ly = domain->lengths[1];
  const double lz = domain->lengths[2];
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = ib->particles[n];
    p->x += p->dx;
    p->y += p->dy;
    p->z += p->dz;
    // correct periodicity
    if(p->y < 0.) p->y += ly;
    if(p->y > ly) p->y -= ly;
    if(p->z < 0.) p->z += lz;
    if(p->z > lz) p->z -= lz;
  }
  return 0;
}

int ib_update_momentum_fleid(const domain_t *domain, fluid_t *fluid, const ib_t *ib){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double *ux = fluid->ux;
  double *uy = fluid->uy;
  double *uz = fluid->uz;
  const double *dux = ib->dux;
  const double *duy = ib->duy;
  const double *duz = ib->duz;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        UX(i, j, k) += 0.5*(DUX(i-1, j  , k  )+DUX(i  , j  , k  ));
      }
    }
  }
  fluid_update_boundaries_ux(domain, ux);
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UY(i, j, k) += 0.5*(DUY(i  , j-1, k  )+DUY(i  , j  , k  ));
      }
    }
  }
  fluid_update_boundaries_uy(domain, uy);
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        UZ(i, j, k) += 0.5*(DUZ(i  , j  , k-1)+DUZ(i  , j  , k  ));
      }
    }
  }
  fluid_update_boundaries_uz(domain, uz);
  return 0;
}


int ib_update_particles(const domain_t *domain, ib_t *ib){
  update_particle_velocities(ib);
  update_particle_positions(domain, ib);
  return 0;
}

