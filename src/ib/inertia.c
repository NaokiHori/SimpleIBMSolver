#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "ib.h"
#include "arrays/fluid.h"



static int compute_inertia(const domain_t *domain, const int cnstep, const fluid_t *fluid, ib_t *ib){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const int ioffs = domain->offsets[0];
  const int joffs = domain->offsets[1];
  const int koffs = domain->offsets[2];
  const double ly = domain->lengths[1];
  const double lz = domain->lengths[2];
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double xoffs = dx * ioffs;
  const double yoffs = dy * joffs;
  const double zoffs = dz * koffs;
  const double *ux = fluid->ux;
  const double *uy = fluid->uy;
  const double *uz = fluid->uz;
  const int n_particles = ib->n_particles;
  particle_t **particles = ib->particles;
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    // constant parameters
    const double pr  = p->r;
    const double px  = p->x + p->dx;
    const double py  = p->y + p->dy;
    const double pz  = p->z + p->dz;
    const double pm  = ib_compute_mass(1., pr);
    const double pim = ib_compute_moment_of_inertia(1., pr);
    // buffers (for simplicity)
    double iux = 0.;
    double iuy = 0.;
    double iuz = 0.;
    double ivx = 0.;
    double ivy = 0.;
    double ivz = 0.;
    for(int zper = -1; zper <= 1; zper++){
      double pz_ = pz + lz * zper;
      for(int yper = -1; yper <= 1; yper++){
        double py_ = py + ly * yper;
        int imin, imax, jmin, jmax, kmin, kmax;
        ib_decide_loop_size(1, isize, dx, pr, px -xoffs, &imin, &imax);
        ib_decide_loop_size(1, jsize, dy, pr, py_-yoffs, &jmin, &jmax);
        ib_decide_loop_size(1, ksize, dz, pr, pz_-zoffs, &kmin, &kmax);
        for(int k = kmin; k <= kmax; k++){
          double z = zoffs + 0.5 * (2 * k - 1) * dz;
          for(int j = jmin; j <= jmax; j++){
            double y = yoffs + 0.5 * (2 * j - 1) * dy;
            for(int i = imin; i <= imax; i++){
              double x = xoffs + 0.5 * (2 * i - 1) * dx;
              double w = ib_v_weight(fmin(dx, fmin(dy, dz)), pr, px, py_, pz_, x, y, z);
              double valx = w * 0.5*(UX(i  , j  , k  )+UX(i+1, j  , k  ))*(dx*dy*dz);
              double valy = w * 0.5*(UY(i  , j  , k  )+UY(i  , j+1, k  ))*(dx*dy*dz);
              double valz = w * 0.5*(UZ(i  , j  , k  )+UZ(i  , j  , k+1))*(dx*dy*dz);
              iux += valx/pm;
              iuy += valy/pm;
              iuz += valz/pm;
              ivx += -(z-pz_)*valy/pim;
              ivx += +(y-py )*valz/pim;
              ivy += -(x-px )*valz/pim;
              ivy += +(z-pz )*valx/pim;
              ivz += -(y-py_)*valx/pim;
              ivz += +(x-px )*valy/pim;
            }
          }
        }
      }
    }
    // assign results
    p->iux[cnstep] = iux;
    p->iuy[cnstep] = iuy;
    p->iuz[cnstep] = iuz;
    p->ivx[cnstep] = ivx;
    p->ivy[cnstep] = ivy;
    p->ivz[cnstep] = ivz;
  }
  return 0;
}

static int synchronise_information(const int cnstep, ib_t *ib){
  const int n_particles = ib->n_particles;
  particle_t **particles = ib->particles;
  // message buffer
  double *buf = common_calloc(6*n_particles, sizeof(double));
  // pack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    buf[6*n+0] = p->iux[cnstep];
    buf[6*n+1] = p->iuy[cnstep];
    buf[6*n+2] = p->iuz[cnstep];
    buf[6*n+3] = p->ivx[cnstep];
    buf[6*n+4] = p->ivy[cnstep];
    buf[6*n+5] = p->ivz[cnstep];
  }
  // sum up all
  MPI_Allreduce(MPI_IN_PLACE, buf, 6*n_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // unpack
  for(int n = 0; n < n_particles; n++){
    particle_t *p = particles[n];
    p->iux[cnstep] = buf[6*n+0];
    p->iuy[cnstep] = buf[6*n+1];
    p->iuz[cnstep] = buf[6*n+2];
    p->ivx[cnstep] = buf[6*n+3];
    p->ivy[cnstep] = buf[6*n+4];
    p->ivz[cnstep] = buf[6*n+5];
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

