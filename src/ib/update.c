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
    p->vz += p->dvz;
  }
  return 0;
}

static int update_particle_positions(const domain_t *domain, ib_t *ib){
  const double ly = domain->lengths[1];
  const int n_particles = ib->n_particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = ib->particles[n];
    p->x += p->dx;
    p->y += p->dy;
    // correct periodicity
    if(p->y < 0.) p->y += ly;
    if(p->y > ly) p->y -= ly;
  }
  return 0;
}

int ib_update_momentum_fleid(const domain_t *domain, fluid_t *fluid, const ib_t *ib){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  double *ux = fluid->ux;
  double *uy = fluid->uy;
  const double *dux = ib->dux;
  const double *duy = ib->duy;
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      UX(i, j) += 0.5*(DUX(i-1, j  )+DUX(i  , j  ));
    }
  }
  fluid_update_boundaries_ux(domain, ux);
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      UY(i, j) += 0.5*(DUY(i  , j-1)+DUY(i  , j  ));
    }
  }
  fluid_update_boundaries_uy(domain, uy);
  return 0;
}


int ib_update_particles(const domain_t *domain, ib_t *ib){
  update_particle_velocities(ib);
  update_particle_positions(domain, ib);
  return 0;
}

