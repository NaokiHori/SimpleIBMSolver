#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "ib.h"
#include "arrays/fluid.h"



static int compute_inertia(const domain_t *domain, const int cnstep, const fluid_t *fluid, ib_t *ib){
  // \int ux dV
  // \int uy dV
  // \int ( - ry ux + rx uy ) dV
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ioffs = domain->offsets[0];
  const int joffs = domain->offsets[1];
  const double ly = domain->lengths[1];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double xoffs = dx * ioffs;
  const double yoffs = dy * joffs;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const int n_particles = ib->n_particles;
  particle_t **particles = ib->particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    // constant parameters
    const double pr  = p->r;
    const double px  = p->x + p->dx;
    const double py  = p->y + p->dy;
    const double pm  = ib_compute_mass(1., pr);
    const double pim = ib_compute_moment_of_inertia(1., pr);
    // buffers (for simplicity)
    double iux = 0.;
    double iuy = 0.;
    double ivz = 0.;
    for(int periodic = -1; periodic <= 1; periodic++){
      double py_ = py + ly * periodic;
      int imin, imax, jmin, jmax;
      ib_decide_loop_size(1, isize, dx, pr, px -xoffs, &imin, &imax);
      ib_decide_loop_size(1, jsize, dy, pr, py_-yoffs, &jmin, &jmax);
      for(int j = jmin; j <= jmax; j++){
        double y = yoffs + 0.5 * (2 * j - 1) * dy;
        for(int i = imin; i <= imax; i++){
          double x = xoffs + 0.5 * (2 * i - 1) * dx;
          double w = ib_v_weight(fmin(dx, dy), pr, px, py_, x, y);
          double valx = w * 0.5*(UX(i  , j  )+UX(i+1, j  ))*(dx*dy);
          double valy = w * 0.5*(UY(i  , j  )+UY(i  , j+1))*(dx*dy);
          iux += valx/pm;
          ivz += -(y-py_)*valx/pim;
          iuy += valy/pm;
          ivz += +(x-px )*valy/pim;
        }
      }
    }
    // assign results
    p->iux[cnstep] = iux;
    p->iuy[cnstep] = iuy;
    p->ivz[cnstep] = ivz;
  }
  return 0;
}

static int synchronise_information(const int cnstep, ib_t *ib){
  const int n_particles = ib->n_particles;
  particle_t **particles = ib->particles;
  // message buffer
  double *buf = common_calloc(3*n_particles, sizeof(double));
  // pack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    buf[3*n+0] = p->iux[cnstep];
    buf[3*n+1] = p->iuy[cnstep];
    buf[3*n+2] = p->ivz[cnstep];
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, 3*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // unpack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->iux[cnstep] = buf[3*n+0];
    p->iuy[cnstep] = buf[3*n+1];
    p->ivz[cnstep] = buf[3*n+2];
  }
  common_free(buf);
  return 0;
}


int ib_compute_inertia(const domain_t *domain, const int cnstep, const fluid_t *fluid, ib_t *ib){
  // update for each particle LOCALLY
  compute_inertia(domain, cnstep, fluid, ib);
  // communicate updated information
  synchronise_information(cnstep, ib);
  return 0;
}

