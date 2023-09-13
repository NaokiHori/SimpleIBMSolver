#include <stdio.h>
#include <math.h>
#include "memory.h"
#include "fileio.h"
#include "ib.h"
#include "ib_solver.h"
#include "array_macros/ib/ibfx.h"
#include "array_macros/ib/ibfy.h"
#if NDIMS == 3
#include "array_macros/ib/ibfz.h"
#endif

static double compute_m(
    const double d,
    const double r
){
  const double pi = 3.14159265358979324;
#if NDIMS == 2
  const double vol = pi * pow(r, 2.);
#else
  const double vol = 4. / 3. * pi * pow(r, 3.);
#endif
  return d * vol;
}

static double compute_mi(
    const double d,
    const double r
){
  const double m = compute_m(d, r);
#if NDIMS == 2
  return 0.5 * m * pow(r, 2.);
#else
  return 0.4 * m * pow(r, 2.);
#endif
}

static void report(
    const sdecomp_info_t * info,
    const ib_t * ib
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(info, &myrank);
  if(root == myrank){
    printf("IB\n");
    printf("\tnitems: %zu\n", ib->nitems);
    fflush(stdout);
  }
}

int ib_init(
    const char dirname_ic[],
    const domain_t * domain,
    ib_t * ib
){
  // load number of particles
  if(0 != fileio.r_serial(dirname_ic, "p_nitems", 0, NULL, fileio.npy_size_t, sizeof(size_t), &ib->nitems)){
    return 1;
  }
  // allocate buffers
  const size_t nitems = ib->nitems;
  ib->particles = memory_calloc(nitems, sizeof(particle_t));
  if(0 != array.prepare(domain, IBFX_NADDS, sizeof(double), &ib->ibfx)) return 1;
  if(0 != array.prepare(domain, IBFY_NADDS, sizeof(double), &ib->ibfy)) return 1;
#if NDIMS == 3
  if(0 != array.prepare(domain, IBFZ_NADDS, sizeof(double), &ib->ibfz)) return 1;
#endif
  // convert SoA to AoS
  double * rs = memory_calloc(nitems, sizeof(double));
  double * ds = memory_calloc(nitems, sizeof(double));
  double * xs = memory_calloc(nitems, sizeof(double));
  double * ys = memory_calloc(nitems, sizeof(double));
#if NDIMS == 3
  double * zs = memory_calloc(nitems, sizeof(double));
#endif
  double * uxs = memory_calloc(nitems, sizeof(double));
  double * uys = memory_calloc(nitems, sizeof(double));
#if NDIMS == 3
  double * uzs = memory_calloc(nitems, sizeof(double));
#endif
#if NDIMS == 3
  double * vxs = memory_calloc(nitems, sizeof(double));
  double * vys = memory_calloc(nitems, sizeof(double));
#endif
  double * vzs = memory_calloc(nitems, sizeof(double));
  if(0 != fileio.r_serial(dirname_ic, "p_rs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), rs)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname_ic, "p_ds", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), ds)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname_ic, "p_xs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), xs)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname_ic, "p_ys", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), ys)){
    return 1;
  }
#if NDIMS == 3
  if(0 != fileio.r_serial(dirname_ic, "p_zs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), zs)){
    return 1;
  }
#endif
  if(0 != fileio.r_serial(dirname_ic, "p_uxs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), uxs)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname_ic, "p_uys", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), uys)){
    return 1;
  }
#if NDIMS == 3
  if(0 != fileio.r_serial(dirname_ic, "p_uzs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), uzs)){
    return 1;
  }
#endif
#if NDIMS == 3
  if(0 != fileio.r_serial(dirname_ic, "p_vxs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), vxs)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname_ic, "p_vys", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), vys)){
    return 1;
  }
#endif
  if(0 != fileio.r_serial(dirname_ic, "p_vzs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), vzs)){
    return 1;
  }
  particle_t * ps = ib->particles;
  for(size_t n = 0; n < nitems; n++){
    ps[n].r  = rs[n];
    ps[n].d  = ds[n];
    ps[n].m  = compute_m (ps[n].d, ps[n].r);
    ps[n].mi = compute_mi(ps[n].d, ps[n].r);
    ps[n].x  = xs[n];
    ps[n].y  = ys[n];
#if NDIMS == 3
    ps[n].z  = zs[n];
#endif
    ps[n].ux = uxs[n];
    ps[n].uy = uys[n];
#if NDIMS == 3
    ps[n].uz = uzs[n];
#endif
#if NDIMS == 3
    ps[n].vx = vxs[n];
    ps[n].vy = vys[n];
#endif
    ps[n].vz = vzs[n];
  }
  memory_free(rs);
  memory_free(xs);
  memory_free(ys);
#if NDIMS == 3
  memory_free(zs);
#endif
  memory_free(uxs);
  memory_free(uys);
#if NDIMS == 3
  memory_free(uzs);
#endif
#if NDIMS == 3
  memory_free(vxs);
  memory_free(vys);
#endif
  memory_free(vzs);
  report(domain->info, ib);
  return 0;
}
