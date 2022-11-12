#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include "sdecomp.h"


/**
 * @brief construct a structure sdecomp_t
 * @param[in] comm_default : MPI communicator which contains all processes
 *                             participating in the decomposition
 *                             (normally MPI_COMM_WORLD)
 * @param[in] ndims        : number of dimensions of the target domain
 * @param[in] dims         : number of processes in each dimension
 * @param[in] periods      : periodicities in each dimension
 * @return                 : a pointer to sdecomp_t, which is allocated and initialised
 */
sdecomp_t *sdecomp_init(const MPI_Comm comm_default, const int ndims, const int *dims, const int *periods){
  /* sanitise argument "ndims" */
  if(ndims != 2 && ndims != 3){
    printf("ndims should be one of:\n");
    printf("  2 (two-dimensional domain)\n");
    printf("  3 (three-dimensional domain)\n");
    printf("%d is given\n", ndims);
    MPI_Abort(comm_default, 0);
  }
  /* sanitise argument "dims" */
  {
    // check decompose automatically or not
    bool auto_decomp = true;
    for(int dim = 0; dim < ndims; dim++){
      // if at least one of dims[dim] is non-zero,
      // user specifies how the domain should be decomposed
      if(dims[dim] != 0){
        auto_decomp = false;
        break;
      }
    }
    if(!auto_decomp){
      if(dims[0] != 1){
        printf("Process distribution is specified by the user:\n");
        if(ndims == 2){
          printf("  dims[0]: %d, dims[1]: %d\n", dims[0], dims[1]);
        }else{
          printf("  dims[0]: %d, dims[1]: %d, dims[2]: %d\n", dims[0], dims[1], dims[2]);
        }
        printf("When specified, dims[0] should be 1 since we do not decompose in x\n");
        printf("Please zero-initialise \"dims\" if you want to distribute process automatically\n");
        MPI_Abort(comm_default, 0);
      }
      // compare number of processes
      int nprocs;
      MPI_Comm_size(comm_default, &nprocs);
      int nprocs_user = 1;
      for(int dim = 0; dim < ndims; dim++){
        nprocs_user *= dims[dim];
      }
      if(nprocs != nprocs_user){
        printf("Number of processes in comm_default: %d\n", nprocs);
        if(ndims == 2){
          printf("dims[0] x dims[1]: %d\n", nprocs_user);
        }else{
          printf("dims[0] x dims[1] x dims[2]: %d\n", nprocs_user);
        }
        printf("They should be identical\n");
        MPI_Abort(comm_default, 0);
      }
    }
  }
  /*
   * create a new structure storing
   *   1. number of dmensions
   *   2. Cartesian info (x1 pencil)
   */
  sdecomp_t *sdecomp = sdecomp_calloc(1, sizeof(sdecomp_t));
  sdecomp->ndims = ndims;
  {
    // number of total processes participating in this decomposition
    int nprocs;
    MPI_Comm_size(comm_default, &nprocs);
    // number of processes in each dimension
    int *dims_ = sdecomp_calloc((size_t)ndims, sizeof(int));
    // check decompose automatically or not
    bool auto_decomp = true;
    for(int dim = 0; dim < ndims; dim++){
      // if at least one of dims[dim] is non-zero,
      // user specifies how the domain should be decomposed
      if(dims[dim] != 0){
        auto_decomp = false;
        break;
      }
    }
    if(auto_decomp){
      for(int dim = 0; dim < ndims; dim++){
        // force 1st dimension NOT decomposed (assign 1)
        dims_[dim] = dim == 0 ? 1 : 0;
      }
      MPI_Dims_create(nprocs, ndims, dims_);
    }else{
      memcpy(dims_, dims, (size_t)ndims * sizeof(int));
    }
    // MPI rank in comm_default is not important for me
    const int reorder = 1;
    // create communicator
    MPI_Cart_create(comm_default, ndims, dims_, periods, reorder, &(sdecomp->comm_cart));
    // clean-up heaps
    sdecomp_free(dims_);
  }
  // for debug use
#if defined(SDECOMP_DEBUG)
  {
    // check process distribution
    {
      int *dims    = sdecomp_calloc(ndims, sizeof(int));
      int *periods = sdecomp_calloc(ndims, sizeof(int));
      int *coords  = sdecomp_calloc(ndims, sizeof(int));
      MPI_Cart_get(sdecomp->comm_cart, ndims, dims, periods, coords);
      int myrank;
      MPI_Comm_rank(sdecomp->comm_cart, &myrank);
      if(myrank == 0){
        printf("Number of processes: ");
        for(int dim = 0; dim < ndims; dim++){
          printf("%d%c", dims[dim], dim == ndims-1 ? '\n' : ' ');
        }
        printf("Periodicity: ");
        for(int dim = 0; dim < ndims; dim++){
          printf("%s%c", periods[dim] ? "true" : "false", dim == ndims-1 ? '\n' : ' ');
        }
      }
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    // check position of each process in Cartesian communicator
    {
      int nprocs, myrank;
      MPI_Comm_size(sdecomp->comm_cart, &nprocs);
      MPI_Comm_rank(sdecomp->comm_cart, &myrank);
      int *coords = sdecomp_calloc(ndims, sizeof(int));
      for(int rank = 0; rank < nprocs; rank++){
        MPI_Cart_coords(sdecomp->comm_cart, rank, ndims, coords);
        if(myrank == 0){
          printf("Rank %d is at (", rank);
          for(int dim = 0; dim < ndims; dim++){
            printf("%d%s", coords[dim], dim == ndims-1 ? ")\n" : ", ");
          }
        }
      }
      sdecomp_free(coords);
    }
  }
#endif // SDECOMP_DEBUG
  return sdecomp;
}

