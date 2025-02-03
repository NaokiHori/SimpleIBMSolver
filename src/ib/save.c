#include <stdbool.h>
#include "memory.h"
#include "array.h"
#include "fileio.h"
#include "domain.h"
#include "ib.h"
#include "ib_solver.h"

static const bool dump_arrays = false;

int ib_save(
    const char dirname[],
    const domain_t * domain,
    const ib_t * ib
){
  if(dump_arrays){
    array.dump(domain, dirname, "ibfx", fileio.npy_double, &ib->ibfx);
    array.dump(domain, dirname, "ibfy", fileio.npy_double, &ib->ibfy);
    array.dump(domain, dirname, "ibfz", fileio.npy_double, &ib->ibfz);
  }
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root != myrank){
    return 0;
  }
  const size_t nitems = ib->nitems;
  // convert AoS to SoA
  double * rs = memory_calloc(nitems, sizeof(double));
  double * ds = memory_calloc(nitems, sizeof(double));
  double * xs = memory_calloc(nitems, sizeof(double));
  double * ys = memory_calloc(nitems, sizeof(double));
  double * zs = memory_calloc(nitems, sizeof(double));
  double * uxs = memory_calloc(nitems, sizeof(double));
  double * uys = memory_calloc(nitems, sizeof(double));
  double * uzs = memory_calloc(nitems, sizeof(double));
  double * vxs = memory_calloc(nitems, sizeof(double));
  double * vys = memory_calloc(nitems, sizeof(double));
  double * vzs = memory_calloc(nitems, sizeof(double));
  particle_t * ps = ib->particles;
  for(size_t n = 0; n < nitems; n++){
    rs[n] = ps[n].r;
    ds[n] = ps[n].d;
    xs[n] = ps[n].x;
    ys[n] = ps[n].y;
    zs[n] = ps[n].z;
    uxs[n] = ps[n].ux;
    uys[n] = ps[n].uy;
    uzs[n] = ps[n].uz;
    vxs[n] = ps[n].vx;
    vys[n] = ps[n].vy;
    vzs[n] = ps[n].vz;
  }
  if(0 != fileio.w_serial(dirname, "p_nitems", 0, NULL, fileio.npy_size_t, sizeof(size_t), &nitems)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_rs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), rs)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_ds", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), ds)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_xs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), xs)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_ys", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), ys)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_zs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), zs)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_uxs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), uxs)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_uys", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), uys)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_uzs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), uzs)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_vxs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), vxs)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_vys", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), vys)){
    return 1;
  }
  if(0 != fileio.w_serial(dirname, "p_vzs", 1, (size_t [1]){nitems}, fileio.npy_double, sizeof(double), vzs)){
    return 1;
  }
  memory_free(rs);
  memory_free(ds);
  memory_free(xs);
  memory_free(ys);
  memory_free(zs);
  memory_free(uxs);
  memory_free(uys);
  memory_free(uzs);
  memory_free(vxs);
  memory_free(vys);
  memory_free(vzs);
  return 0;
}

