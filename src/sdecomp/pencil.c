#include <assert.h>
#include "sdecomp.h"


/*** process distributions ***/

typedef enum {
  type_nprocs,
  type_myrank
} type_t;

static int kernel_get_process_config_2d(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const type_t type){
#define SDECOMP_NDIMS 2
  // compute
  //   type == type_nprocs: number of processors
  // or
  //   type == type_myrank: my location
  // of the given decomposition "sdecomp"
  //   in the given dimension "dim"
  const int ndims = sdecomp->ndims;
  const MPI_Comm comm_cart = sdecomp->comm_cart;
  // get configuration for x1 pencil
  int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
  int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
  int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
  MPI_Cart_get(comm_cart, ndims, dims, periods, coords);
  // x1 pencil
  int x1pencil_nprocs[SDECOMP_NDIMS];
  int x1pencil_myrank[SDECOMP_NDIMS];
  x1pencil_nprocs[0] = dims[0];
  x1pencil_nprocs[1] = dims[1];
  x1pencil_myrank[0] = coords[0];
  x1pencil_myrank[1] = coords[1];
  sdecomp_free(dims);
  sdecomp_free(periods);
  sdecomp_free(coords);
  if(pencil == SDECOMP_X1PENCIL){
    return type == type_nprocs ? x1pencil_nprocs[dir] : x1pencil_myrank[dir];
  }
  // y1 pencil, determined by x1 pencil
  int y1pencil_nprocs[SDECOMP_NDIMS];
  int y1pencil_myrank[SDECOMP_NDIMS];
  y1pencil_nprocs[0] = x1pencil_nprocs[1];
  y1pencil_nprocs[1] = x1pencil_nprocs[0];
  y1pencil_myrank[0] = x1pencil_myrank[1];
  y1pencil_myrank[1] = x1pencil_myrank[0];
  if(pencil == SDECOMP_Y1PENCIL){
    return type == type_nprocs ? y1pencil_nprocs[dir] : y1pencil_myrank[dir];
  }
  return 0;
#undef SDECOMP_NDIMS
}

static int kernel_get_process_config_3d(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const type_t type){
#define SDECOMP_NDIMS 3
  // compute
  //   type == type_nprocs: number of processors
  // or
  //   type == type_myrank: my location
  // of the given decomposition "sdecomp"
  //   in the given dimension "dim"
  const int ndims = sdecomp->ndims;
  const MPI_Comm comm_cart = sdecomp->comm_cart;
  // get configuration for x1 pencil
  int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
  int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
  int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
  MPI_Cart_get(comm_cart, ndims, dims, periods, coords);
  // x1 pencil
  int x1pencil_nprocs[SDECOMP_NDIMS];
  int x1pencil_myrank[SDECOMP_NDIMS];
  x1pencil_nprocs[0] = dims[0];
  x1pencil_nprocs[1] = dims[1];
  x1pencil_nprocs[2] = dims[2];
  x1pencil_myrank[0] = coords[0];
  x1pencil_myrank[1] = coords[1];
  x1pencil_myrank[2] = coords[2];
  sdecomp_free(dims);
  sdecomp_free(periods);
  sdecomp_free(coords);
  if(pencil == SDECOMP_X1PENCIL){
    return type == type_nprocs ? x1pencil_nprocs[dir] : x1pencil_myrank[dir];
  }
  // y1 pencil, determined by x1 pencil
  int y1pencil_nprocs[SDECOMP_NDIMS];
  int y1pencil_myrank[SDECOMP_NDIMS];
  y1pencil_nprocs[0] = x1pencil_nprocs[1];
  y1pencil_nprocs[1] = x1pencil_nprocs[0];
  y1pencil_nprocs[2] = x1pencil_nprocs[2];
  y1pencil_myrank[0] = x1pencil_myrank[1];
  y1pencil_myrank[1] = x1pencil_myrank[0];
  y1pencil_myrank[2] = x1pencil_myrank[2];
  if(pencil == SDECOMP_Y1PENCIL){
    return type == type_nprocs ? y1pencil_nprocs[dir] : y1pencil_myrank[dir];
  }
  // z1 pencil, determined by y1 pencil
  int z1pencil_nprocs[SDECOMP_NDIMS];
  int z1pencil_myrank[SDECOMP_NDIMS];
  z1pencil_nprocs[0] = y1pencil_nprocs[0];
  z1pencil_nprocs[1] = y1pencil_nprocs[2];
  z1pencil_nprocs[2] = y1pencil_nprocs[1];
  z1pencil_myrank[0] = y1pencil_myrank[0];
  z1pencil_myrank[1] = y1pencil_myrank[2];
  z1pencil_myrank[2] = y1pencil_myrank[1];
  if(pencil == SDECOMP_Z1PENCIL){
    return type == type_nprocs ? z1pencil_nprocs[dir] : z1pencil_myrank[dir];
  }
  // x2 pencil, determined by z1 pencil
  int x2pencil_nprocs[SDECOMP_NDIMS];
  int x2pencil_myrank[SDECOMP_NDIMS];
  x2pencil_nprocs[0] = z1pencil_nprocs[2];
  x2pencil_nprocs[1] = z1pencil_nprocs[1];
  x2pencil_nprocs[2] = z1pencil_nprocs[0];
  x2pencil_myrank[0] = z1pencil_myrank[2];
  x2pencil_myrank[1] = z1pencil_myrank[1];
  x2pencil_myrank[2] = z1pencil_myrank[0];
  if(pencil == SDECOMP_X2PENCIL){
    return type == type_nprocs ? x2pencil_nprocs[dir] : x2pencil_myrank[dir];
  }
  // y2 pencil, determined by x2 pencil
  int y2pencil_nprocs[SDECOMP_NDIMS];
  int y2pencil_myrank[SDECOMP_NDIMS];
  y2pencil_nprocs[0] = x2pencil_nprocs[1];
  y2pencil_nprocs[1] = x2pencil_nprocs[0];
  y2pencil_nprocs[2] = x2pencil_nprocs[2];
  y2pencil_myrank[0] = x2pencil_myrank[1];
  y2pencil_myrank[1] = x2pencil_myrank[0];
  y2pencil_myrank[2] = x2pencil_myrank[2];
  if(pencil == SDECOMP_Y2PENCIL){
    return type == type_nprocs ? y2pencil_nprocs[dir] : y2pencil_myrank[dir];
  }
  // z2 pencil, determined by y2 pencil
  int z2pencil_nprocs[SDECOMP_NDIMS];
  int z2pencil_myrank[SDECOMP_NDIMS];
  z2pencil_nprocs[0] = y2pencil_nprocs[0];
  z2pencil_nprocs[1] = y2pencil_nprocs[2];
  z2pencil_nprocs[2] = y2pencil_nprocs[1];
  z2pencil_myrank[0] = y2pencil_myrank[0];
  z2pencil_myrank[1] = y2pencil_myrank[2];
  z2pencil_myrank[2] = y2pencil_myrank[1];
  if(pencil == SDECOMP_Z2PENCIL){
    return type == type_nprocs ? z2pencil_nprocs[dir] : z2pencil_myrank[dir];
  }
  return 0;
#undef SDECOMP_NDIMS
}

/**
 * @brief get number of processes in the given dimension
 * @param[in] sdecomp : struct containing information of process distribution
 * @param[in] pencil  : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in] dir     : direction which I am interested in
 * @return            : number of processes in the given dimension
 */
int sdecomp_get_nprocs(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir){
  const int ndims = sdecomp->ndims;
  assert(ndims == 2 || ndims == 3);
  assert(0 <= (int)dir && (int)dir <= ndims-1);
  if(ndims == 2){
    return kernel_get_process_config_2d(sdecomp, pencil, dir, type_nprocs);
  }else{
    return kernel_get_process_config_3d(sdecomp, pencil, dir, type_nprocs);
  }
}

/**
 * @brief get my process position in the whole domain
 * @param[in] sdecomp : struct containing information of process distribution
 * @param[in] pencil  : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in] dir     : direction which I am interested in
 * @return            : my position in the given dimension
 */
int sdecomp_get_myrank(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir){
  const int ndims = sdecomp->ndims;
  assert(ndims == 2 || ndims == 3);
  assert(0 <= (int)dir && (int)dir <= ndims-1);
  if(ndims == 2){
    return kernel_get_process_config_2d(sdecomp, pencil, dir, type_myrank);
  }else{
    return kernel_get_process_config_3d(sdecomp, pencil, dir, type_myrank);
  }
}

/*** pencil definitions ***/

int sdecomp_kernel_get_mysize(const int num_total, const int nprocs, const int myrank){
  /* number of grid points of the process */
  // example: num_total: 11, nprocs: 3 (3 processes in total)
  // myrank = 0 -> num_local = 3
  // myrank = 1 -> num_local = 4
  // myrank = 2 -> num_local = 4
  // -> sum of "num_local"s is "num_total"
  int num_local = (num_total+myrank)/nprocs;
  return num_local;
}

int sdecomp_kernel_get_offset(const int num_total, const int nprocs, const int myrank){
  /* sum up the number of grid points to the process */
  // NOTE: equivalent to MPI_Exscan, but this is intuitive and shorter
  int offset = 0;
  for(int i = 0; i < myrank; i++){
    offset += sdecomp_kernel_get_mysize(num_total, nprocs, i);
  }
  return offset;
}

/**
 * @brief wrapper function to get pencil size or offset
 * @param[in] sdecomp : struct containing information of process distribution
 * @param[in] pencil  : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in] dir     : direction which I am interested in
 * @param[in] num     : number of total grid points in the dimension
 * @param[in] func    : function pointer to kernel functions
 * @return            : number of local grid points in the dimension
 */
static int get_pencil_config(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const int num, int (*func)(const int, const int, const int)){
  const int ndims = sdecomp->ndims;
  assert(ndims == 2 || ndims == 3);
  assert(0 <= (int)dir && (int)dir <= ndims-1);
  if(ndims == 2){
    assert(pencil == SDECOMP_X1PENCIL || pencil == SDECOMP_Y1PENCIL);
  }else{
    assert(
        pencil == SDECOMP_X1PENCIL || pencil == SDECOMP_Y1PENCIL || pencil == SDECOMP_Z1PENCIL
     || pencil == SDECOMP_X2PENCIL || pencil == SDECOMP_Y2PENCIL || pencil == SDECOMP_Z2PENCIL
    );
  }
  const int nprocs = sdecomp_get_nprocs(sdecomp, pencil, (sdecomp_dir_t)dir);
  const int myrank = sdecomp_get_myrank(sdecomp, pencil, (sdecomp_dir_t)dir);
  return func(num, nprocs, myrank);
}

/**
 * @brief get number of grid points of the given pencil in the given dimension
 * @param[in] sdecomp : struct containing information of process distribution
 * @param[in] pencil  : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in] dir     : direction which I am interested in
 * @param[in] num     : number of total grid points in the dimension
 * @return            : number of local grid points in the dimension
 */
int sdecomp_get_pencil_mysize(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const int num){
  return get_pencil_config(sdecomp, pencil, dir, num, sdecomp_kernel_get_mysize);
}

/**
 * @brief get offset of grid points of the given pencil in the given dimension
 * @param[in] sdecomp : struct containing information of process distribution
 * @param[in] pencil  : type of pencil (e.g., SDECOMP_X1PENCIL)
 * @param[in] dir     : direction which I am interested in
 * @param[in] num     : number of total grid points in the dimension
 * @return            : offset in the give dimension
 */
int sdecomp_get_pencil_offset(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const int num){
  return get_pencil_config(sdecomp, pencil, dir, num, sdecomp_kernel_get_offset);
}

