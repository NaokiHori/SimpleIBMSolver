#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "common.h"
#include "fileio.h"


#if NDIMS == 2

typedef struct circle_t_ {
  double x, y;
  double r;
} circle_t;

#else // NDIMS == 3

typedef struct circle_t_ {
  double x, y, z;
  double r;
} circle_t;

#endif // NDIMS

static double gen_random(const double min, const double max){
  return (max-min)*rand()/RAND_MAX+min;
}

#if NDIMS == 2

static int check_stats(double *vfrac, const double lx, const double ly, const int n_particles, const double *rs){
  *vfrac = 0.;
  for(int n = 0; n < n_particles; n++){
    *vfrac += M_PI*pow(rs[n], 2.);
  }
  *vfrac /= lx*ly;
  return 0;
}

#else // NDIMS == 3

static int check_stats(double *vfrac, const double lx, const double ly, const double lz, const int n_particles, const double *rs){
  *vfrac = 0.;
  for(int n = 0; n < n_particles; n++){
    *vfrac += 4./3.*M_PI*pow(rs[n], 3.);
  }
  *vfrac /= lx*ly*lz;
  return 0;
}

#endif // NDIMS

#if NDIMS == 2

static int output(const int n_particles, const double *rs, const double *xs, const double *ys, const double *uxs, const double *uys, const double *vzs){
  fileio_w_0d_serial("..", "n_particles",  NPYIO_INT,    sizeof(int),    &n_particles);
  fileio_w_1d_serial("..", "p_rs",  NPYIO_DOUBLE, sizeof(double), n_particles,  rs);
  fileio_w_1d_serial("..", "p_xs",  NPYIO_DOUBLE, sizeof(double), n_particles,  xs);
  fileio_w_1d_serial("..", "p_ys",  NPYIO_DOUBLE, sizeof(double), n_particles,  ys);
  fileio_w_1d_serial("..", "p_uxs", NPYIO_DOUBLE, sizeof(double), n_particles, uxs);
  fileio_w_1d_serial("..", "p_uys", NPYIO_DOUBLE, sizeof(double), n_particles, uys);
  fileio_w_1d_serial("..", "p_vzs", NPYIO_DOUBLE, sizeof(double), n_particles, vzs);
  return 0;
}

#else // NDIMS == 3

static int output(const int n_particles, const double *rs, const double *xs, const double *ys, const double *zs, const double *uxs, const double *uys, const double *uzs, const double *vxs, const double *vys, const double *vzs){
  fileio_w_0d_serial("..", "n_particles",  NPYIO_INT,    sizeof(int),    &n_particles);
  fileio_w_1d_serial("..", "p_rs",  NPYIO_DOUBLE, sizeof(double), n_particles,  rs);
  fileio_w_1d_serial("..", "p_xs",  NPYIO_DOUBLE, sizeof(double), n_particles,  xs);
  fileio_w_1d_serial("..", "p_ys",  NPYIO_DOUBLE, sizeof(double), n_particles,  ys);
  fileio_w_1d_serial("..", "p_zs",  NPYIO_DOUBLE, sizeof(double), n_particles,  zs);
  fileio_w_1d_serial("..", "p_uxs", NPYIO_DOUBLE, sizeof(double), n_particles, uxs);
  fileio_w_1d_serial("..", "p_uys", NPYIO_DOUBLE, sizeof(double), n_particles, uys);
  fileio_w_1d_serial("..", "p_uzs", NPYIO_DOUBLE, sizeof(double), n_particles, uzs);
  fileio_w_1d_serial("..", "p_vxs", NPYIO_DOUBLE, sizeof(double), n_particles, vxs);
  fileio_w_1d_serial("..", "p_vys", NPYIO_DOUBLE, sizeof(double), n_particles, vys);
  fileio_w_1d_serial("..", "p_vzs", NPYIO_DOUBLE, sizeof(double), n_particles, vzs);
  return 0;
}

#endif // NDIMS

#if NDIMS == 2

int main(void){
  srand(time(NULL));
  const double lx = 1.;
  const double ly = 2.;
  const int n_particles = 32;
  double *rs   = common_calloc(n_particles, sizeof(double));
  double *xs   = common_calloc(n_particles, sizeof(double));
  double *ys   = common_calloc(n_particles, sizeof(double));
  double *uxs  = common_calloc(n_particles, sizeof(double));
  double *uys  = common_calloc(n_particles, sizeof(double));
  double *vzs  = common_calloc(n_particles, sizeof(double));
  // radius
  for(int n = 0; n < n_particles; n++){
    rs[n] = gen_random(0.10, 0.10);
  }
  // position, no overlaps
  for(int n0 = 0; n0 < n_particles; n0++){
regen:
    {
      double x0 = gen_random(rs[n0], lx-rs[n0]);
      double y0 = gen_random(    0.,        ly);
      for(int n1 = 0; n1 < n0; n1++){
        // correct periodicity
        double yoffset = 0.;
        {
          double minval = DBL_MAX;
          for(int periodic = -1; periodic <= 1; periodic++){
            double val = fabs((ys[n1]+ly*periodic)-y0);
            if(val < minval){
              minval = val;
              yoffset = ly*periodic;
            }
          }
        }
        // check collision
        {
          double r0 = rs[n0];
          double x1 = xs[n1];
          double y1 = ys[n1]+yoffset;
          double r1 = rs[n1];
          double nx = x1-x0;
          double ny = y1-y0;
          double norm = fmax(hypot(nx, ny), DBL_EPSILON);
          nx /= norm;
          ny /= norm;
          double dist = norm-r0-r1;
          if(dist < 0.){
            goto regen;
          }
        }
      }
      xs[n0] = x0;
      ys[n0] = y0;
      printf("particle %*d @ % .3f % .3f\n", 5, n0, x0, y0);
    }
  }
  {
    double vfrac;
    check_stats(&vfrac, lx, ly, n_particles, rs);
    printf("volume fraction: % .1e\n", vfrac);
  }
  // velocities, still
  memset(uxs, 0, sizeof(double)*n_particles);
  memset(uys, 0, sizeof(double)*n_particles);
  memset(vzs, 0, sizeof(double)*n_particles);
  // output
  output(n_particles, rs, xs, ys, uxs, uys, vzs);
  // clean-up buffers
  common_free( rs);
  common_free( xs);
  common_free( ys);
  common_free(uxs);
  common_free(uys);
  common_free(vzs);
  return 0;
}

#else // NDIMS == 3

int main(void){
  srand(time(NULL));
  const double lx = 1.;
  const double ly = 1.;
  const double lz = 1.;
  const int n_particles = 64;
  double *rs   = common_calloc(n_particles, sizeof(double));
  double *xs   = common_calloc(n_particles, sizeof(double));
  double *ys   = common_calloc(n_particles, sizeof(double));
  double *zs   = common_calloc(n_particles, sizeof(double));
  double *uxs  = common_calloc(n_particles, sizeof(double));
  double *uys  = common_calloc(n_particles, sizeof(double));
  double *uzs  = common_calloc(n_particles, sizeof(double));
  double *vxs  = common_calloc(n_particles, sizeof(double));
  double *vys  = common_calloc(n_particles, sizeof(double));
  double *vzs  = common_calloc(n_particles, sizeof(double));
  // radius
  for(int n = 0; n < n_particles; n++){
    rs[n] = gen_random(0.10, 0.10);
  }
  // position, no overlaps
  for(int n0 = 0; n0 < n_particles; n0++){
regen:
    {
      double x0 = gen_random(rs[n0], lx-rs[n0]);
      double y0 = gen_random(    0.,        ly);
      double z0 = gen_random(    0.,        lz);
      for(int n1 = 0; n1 < n0; n1++){
        // correct periodicity
        double yoffset = 0.;
        {
          double minval = DBL_MAX;
          for(int periodic = -1; periodic <= 1; periodic++){
            double val = fabs((ys[n1]+ly*periodic)-y0);
            if(val < minval){
              minval = val;
              yoffset = ly*periodic;
            }
          }
        }
        // correct periodicity
        double zoffset = 0.;
        {
          double minval = DBL_MAX;
          for(int periodic = -1; periodic <= 1; periodic++){
            double val = fabs((zs[n1]+lz*periodic)-z0);
            if(val < minval){
              minval = val;
              zoffset = lz*periodic;
            }
          }
        }
        // check collision
        {
          double r0 = rs[n0];
          double x1 = xs[n1];
          double y1 = ys[n1]+yoffset;
          double z1 = zs[n1]+zoffset;
          double r1 = rs[n1];
          double nx = x1-x0;
          double ny = y1-y0;
          double nz = z1-z0;
          double norm = sqrt(
              + pow(nx, 2.)
              + pow(ny, 2.)
              + pow(nz, 2.)
          );
          norm = fmax(norm, DBL_EPSILON);
          nx /= norm;
          ny /= norm;
          nz /= norm;
          double dist = norm-r0-r1;
          if(dist < 0.){
            goto regen;
          }
        }
      }
      xs[n0] = x0;
      ys[n0] = y0;
      zs[n0] = z0;
      printf("particle %*d @ % .3f % .3f % .3f\n", 5, n0, x0, y0, z0);
    }
  }
  {
    double vfrac;
    check_stats(&vfrac, lx, ly, lz, n_particles, rs);
    printf("volume fraction: % .1e\n", vfrac);
  }
  // velocities, still
  memset(uxs, 0, sizeof(double)*n_particles);
  memset(uys, 0, sizeof(double)*n_particles);
  memset(uzs, 0, sizeof(double)*n_particles);
  memset(vxs, 0, sizeof(double)*n_particles);
  memset(vys, 0, sizeof(double)*n_particles);
  memset(vzs, 0, sizeof(double)*n_particles);
  // output
  output(n_particles, rs, xs, ys, zs, uxs, uys, uzs, vxs, vys, vzs);
  // clean-up buffers
  common_free( rs);
  common_free( xs);
  common_free( ys);
  common_free( zs);
  common_free(uxs);
  common_free(uys);
  common_free(uzs);
  common_free(vxs);
  common_free(vys);
  common_free(vzs);
  return 0;
}

#endif
