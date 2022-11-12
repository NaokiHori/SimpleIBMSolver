#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "sdecomp.h"


struct sdecomp_transpose_t_ {
  // essential variables used by MPI_Alltoallw
  MPI_Comm comm_2d;
  int *sendcounts, *recvcounts;
  int *sdispls, *rdispls;
  MPI_Datatype *sendtypes, *recvtypes;
  // auxiliary variables used by tests
  sdecomp_pencil_t pencil_bef, pencil_aft;
  int *gsizes;
  size_t size_elem;
};

/**
 * @brief initialise transpose plan for 2d
 * @param[in] sdecomp      : struct contains information of process distribution
 * @param[in] pencil       : type of pencil (SDECOMP_X1PENCIL or SDECOMP_Y1PENCIL)
 * @param[in] gsizes       : global array size in each dimension
 * @param[in] size_elem    : size of each element, e.g., sizeof(double)
 * @param[in] mpi_datatype : corresponding MPI_Datatype of each element, e.g., MPI_DOUBLE
 * @return                 : a pointer to the created plan (struct)
 */
static sdecomp_transpose_t *sdecomp_transpose_init_2d(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const int *gsizes, const size_t size_elem, const MPI_Datatype mpi_datatype){
#define SDECOMP_SDECOMP_NDIMS 2
  assert(
      pencil == SDECOMP_X1PENCIL ||
      pencil == SDECOMP_Y1PENCIL
  );
  const MPI_Comm comm_cart = sdecomp->comm_cart;
  // argument "gsizes" does NOT consider memory contiguous directions
  // convert them from inner memory to outer memory
  // from x1 (0) pencil : (gsizes[0], gsizes[1])
  // from y1 (1) pencil : (gsizes[1], gsizes[0])
  int sizes[SDECOMP_SDECOMP_NDIMS] = {0};
  switch(pencil){
    case SDECOMP_X1PENCIL:
      {
        sizes[0] = gsizes[0];
        sizes[1] = gsizes[1];
        break;
      }
    case SDECOMP_Y1PENCIL:
      {
        sizes[0] = gsizes[1];
        sizes[1] = gsizes[0];
        break;
      }
    case SDECOMP_Z1PENCIL:
    case SDECOMP_X2PENCIL:
    case SDECOMP_Y2PENCIL:
    case SDECOMP_Z2PENCIL:
      {
        assert(
            0 == 1
            /* SDECOMP_Z1PENCIL is not applicable in 2D */
            /* SDECOMP_X2PENCIL is not applicable in 2D */
            /* SDECOMP_Y2PENCIL is not applicable in 2D */
            /* SDECOMP_Z2PENCIL is not applicable in 2D */
        );
        break;
      }
  }
  // NOTE: from now on, since the memory contiguity is taken into account,
  //   the memory is contiguous in x and sparse in y,
  //   which are irrespective to the physical domain
  // two-dimensional communicator, identical for two-dimensional domain
  MPI_Comm comm_2d;
  MPI_Comm_dup(comm_cart, &comm_2d);
  // create plan
  sdecomp_transpose_t *plan = NULL;
  {
    // number of total process and my position ("rank0"-th pencil)
    int nprocs, rank0;
    {
      const int ndims = 2; // 2d communicator
      int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
      MPI_Cart_get(comm_2d, ndims, dims, periods, coords);
      nprocs = dims[1];   // interested in y direction
      rank0  = coords[1]; // interested in y direction
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    int *sendcounts         = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *recvcounts         = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *sdispls            = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *rdispls            = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    MPI_Datatype *sendtypes = sdecomp_calloc((size_t)nprocs, sizeof(MPI_Datatype));
    MPI_Datatype *recvtypes = sdecomp_calloc((size_t)nprocs, sizeof(MPI_Datatype));
    // consider communication between "rank0"-th and "rank1"-th pencils
    for(int rank1 = 0; rank1 < nprocs; rank1++){
      MPI_Datatype temptype;
      // send
      {
        int chunk_isize = sdecomp_kernel_get_mysize(sizes[0], nprocs, rank1);
        int chunk_jsize = sdecomp_kernel_get_mysize(sizes[1], nprocs, rank0);
        /* datatype to be sent: contiguous in y dimension */
        MPI_Type_create_hvector(
            /* int count             */ chunk_jsize,
            /* int blocklength       */ 1,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem * (MPI_Aint)sizes[0],
            /* MPI_Datatype oldtype  */ mpi_datatype,
            /* MPI_Datatype *newtype */ &temptype
        );
        /* datatype to be sent: datatype previously defined is repeated in x direction */
        MPI_Type_create_hvector(
            /* int count             */ chunk_isize,
            /* int blocklength       */ 1,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem,
            /* MPI_Datatype oldtype  */ temptype,
            /* MPI_Datatype *newtype */ &(sendtypes[rank1])
        );
        /* commit send datatype */
        MPI_Type_commit(&(sendtypes[rank1]));
        // check size of the committed datatype
        {
          int size;
          MPI_Type_size(sendtypes[rank1], &size);
          assert(size == (int)size_elem * chunk_isize * chunk_jsize);
        }
        /* number of elements is 1, since we send single "sendtypes[rank1]" */
        sendcounts[rank1] = 1;
        /* offset to the pointer to be sent */
        sdispls[rank1] = (int)size_elem * sdecomp_kernel_get_offset(sizes[0], nprocs, rank1);
      }
      // recv
      {
        int chunk_isize = sdecomp_kernel_get_mysize(sizes[0], nprocs, rank0);
        int chunk_jsize = sdecomp_kernel_get_mysize(sizes[1], nprocs, rank1);
        /* datatype to be received: simple since data being sent is already aligned in y direction */
        MPI_Type_create_hvector(
            /* int count             */ chunk_isize,
            /* int blocklength       */ chunk_jsize,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem * (MPI_Aint)sizes[1],
            /* MPI_Datatype oldtype  */ mpi_datatype,
            /* MPI_Datatype *newtype */ &(recvtypes[rank1])
        );
        /* commit recv datatype */
        MPI_Type_commit(&(recvtypes[rank1]));
        // check size of the committed datatype
        {
          int size;
          MPI_Type_size(recvtypes[rank1], &size);
          assert(size == (int)size_elem * chunk_jsize * chunk_isize);
        }
        /* number of elements is 1, since we receive single "sendtypes[rank1]" */
        recvcounts[rank1] = 1;
        /* offset to the pointer to be received */
        rdispls[rank1] = (int)size_elem * sdecomp_kernel_get_offset(sizes[1], nprocs, rank1);
      }
    }
    // assign all to the struct
    plan = sdecomp_calloc(1, sizeof(sdecomp_transpose_t));
    plan->comm_2d    = comm_2d;
    plan->sendcounts = sendcounts;
    plan->recvcounts = recvcounts;
    plan->sdispls    = sdispls;
    plan->rdispls    = rdispls;
    plan->sendtypes  = sendtypes;
    plan->recvtypes  = recvtypes;
    plan->pencil_bef = pencil;
    plan->pencil_aft = (sdecomp_pencil_t)( ((int)pencil+1) % 2);
    plan->gsizes     = sdecomp_calloc(SDECOMP_SDECOMP_NDIMS, sizeof(int));
    for(int dim = 0; dim < SDECOMP_SDECOMP_NDIMS; dim++){
      plan->gsizes[dim] = gsizes[dim];
    }
    plan->size_elem = size_elem;
  }
  return plan;
#undef SDECOMP_SDECOMP_NDIMS
}

/**
 * @brief initialise transpose plan for 3d, forward transposes
 * @param[in] sdecomp      : struct contains information of process distribution
 * @param[in] pencil       : type of pencil (one of six types)
 * @param[in] gsizes       : global array size in each dimension
 * @param[in] size_elem    : size of each element, e.g., sizeof(double)
 * @param[in] mpi_datatype : corresponding MPI_Datatype of each element, e.g., MPI_DOUBLE
 * @return                 : a pointer to the created plan (struct)
 */
static sdecomp_transpose_t *sdecomp_transpose_fwrd_init_3d(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const int *gsizes, const size_t size_elem, const MPI_Datatype mpi_datatype){
#define SDECOMP_SDECOMP_NDIMS 3
  assert(
      pencil == SDECOMP_X1PENCIL ||
      pencil == SDECOMP_Y1PENCIL ||
      pencil == SDECOMP_Z1PENCIL ||
      pencil == SDECOMP_Z2PENCIL ||
      pencil == SDECOMP_Y2PENCIL ||
      pencil == SDECOMP_X2PENCIL
  );
  const MPI_Comm comm_cart = sdecomp->comm_cart;
  // argument "gsizes" does NOT consider memory contiguous directions
  // convert them from inner memory to outer memory
  // from x1 (0) pencil or x2 (3) pencil : (gsizes[0], gsizes[1], gsizes[2])
  // from y1 (1) pencil or y2 (4) pencil : (gsizes[1], gsizes[2], gsizes[0])
  // from z1 (2) pencil or z2 (5) pencil : (gsizes[2], gsizes[0], gsizes[1])
  // -> condition based on mod(pencil, 3)
  int sizes[SDECOMP_SDECOMP_NDIMS] = {0};
  switch((int)pencil % 3){
    case 0:
      {
        sizes[0] = gsizes[0];
        sizes[1] = gsizes[1];
        sizes[2] = gsizes[2];
        break;
      }
    case 1:
      {
        sizes[0] = gsizes[1];
        sizes[1] = gsizes[2];
        sizes[2] = gsizes[0];
        break;
      }
    case 2:
      {
        sizes[0] = gsizes[2];
        sizes[1] = gsizes[0];
        sizes[2] = gsizes[1];
        break;
      }
    default:
      assert(0 == 1 /* should not be here */);
      break;
  }
  // NOTE: from now on, since the memory contiguity is taken into account,
  //   the memory is contiguous in x, middle in y, and the most sparse in z
  //   which are irrespective to the physical domain
  // create 2d communicator for each z (in terms of memory contiguity) chunk,
  //   and compute domain size in the unchanged dimension
  MPI_Comm comm_2d;
  {
    // we split the default communicator comm_cart and create a new one comm_2d,
    //   in which only processes participating in the all-to-all communication belong
    // comm_cart is x1 pencil ->
    //   we need to consider the unchanged dimsnion on x1 pencil
    // from x1 (0) pencil : z is unchanged -> for x1pencil, z
    // from y1 (1) pencil : x is unchanged -> for x1pencil, y
    // from z1 (2) pencil : y is unchanged -> for x1pencil, z
    // from x2 (3) pencil : z is unchanged -> for x1pencil, y
    // from y2 (4) pencil : x is unchanged -> for x1pencil, z
    // from z2 (5) pencil : y is unchanged -> for x1pencil, y
    // -> unchanged dimension in x1pencil differs for odd and even pencils
    const int unchanged_dim = (int)pencil % 2 == 0 ? 2 : 1;
    int remain_dims[SDECOMP_SDECOMP_NDIMS] = {1, 1, 1};
    remain_dims[unchanged_dim] = 0;
    MPI_Cart_sub(comm_cart, remain_dims, &comm_2d);
    // to compute mysize (number of grid points) in the unchanged dimension,
    //   we need to know the number of processes and my position in the dimension
    int nprocs, myrank;
    {
      const int ndims = 3;
      int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
      MPI_Cart_get(comm_cart, ndims, dims, periods, coords);
      nprocs = dims[unchanged_dim];
      myrank = coords[unchanged_dim];
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    // correct size in the unchanged dimension,
    //   since all-to-all does not care the global size
    sizes[2] = sdecomp_kernel_get_mysize(sizes[2], nprocs, myrank);
  }
  // create plan
  sdecomp_transpose_t *plan = NULL;
  {
    // number of processes which participate in the given 2d communicator "comm_2d"
    int nprocs;
    {
      const int ndims = 2; // 2d communicator
      int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
      MPI_Cart_get(comm_2d, ndims, dims, periods, coords);
      nprocs = dims[1]; // interested in y direction
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    // my position, "rank0"-th pencil
    int rank0;
    {
      const int ndims = 2; // 2d communicator
      int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
      MPI_Cart_get(comm_2d, ndims, dims, periods, coords);
      rank0 = coords[1]; // interested in y direction
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    int *sendcounts         = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *recvcounts         = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *sdispls            = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *rdispls            = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    MPI_Datatype *sendtypes = sdecomp_calloc((size_t)nprocs, sizeof(MPI_Datatype));
    MPI_Datatype *recvtypes = sdecomp_calloc((size_t)nprocs, sizeof(MPI_Datatype));
    // consider communication between "rank0"-th and "rank1"-th pencils
    for(int rank1 = 0; rank1 < nprocs; rank1++){
      MPI_Datatype temptype;
      // send
      {
        int chunk_isize = sdecomp_kernel_get_mysize(sizes[0], nprocs, rank1);
        int chunk_jsize = sdecomp_kernel_get_mysize(sizes[1], nprocs, rank0);
        /* datatype to be sent: contiguous in y and z dimensions */
        MPI_Type_create_hvector(
            /* int count             */ chunk_jsize * sizes[2],
            /* int blocklength       */ 1,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem * (MPI_Aint)sizes[0],
            /* MPI_Datatype oldtype  */ mpi_datatype,
            /* MPI_Datatype *newtype */ &temptype
        );
        /* datatype to be sent: datatype previously defined is repeated in x direction */
        MPI_Type_create_hvector(
            /* int count             */ chunk_isize,
            /* int blocklength       */ 1,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem,
            /* MPI_Datatype oldtype  */ temptype,
            /* MPI_Datatype *newtype */ &(sendtypes[rank1])
        );
        /* commit send datatype */
        MPI_Type_commit(&(sendtypes[rank1]));
        // check size of the committed datatype
        {
          int size;
          MPI_Type_size(sendtypes[rank1], &size);
          assert(size == (int)size_elem * chunk_isize * chunk_jsize * sizes[2]);
        }
        /* number of elements is 1, since we send single "sendtypes[rank1]" */
        sendcounts[rank1] = 1;
        /* offset to the pointer to be sent */
        sdispls[rank1] = (int)size_elem * sdecomp_kernel_get_offset(sizes[0], nprocs, rank1);
      }
      // recv
      {
        int chunk_isize = sdecomp_kernel_get_mysize(sizes[0], nprocs, rank0);
        int chunk_jsize = sdecomp_kernel_get_mysize(sizes[1], nprocs, rank1);
        /* datatype to be received: simple since data being sent is already aligned in y direction */
        MPI_Type_create_hvector(
            /* int count             */ sizes[2] * chunk_isize,
            /* int blocklength       */ chunk_jsize,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem * (MPI_Aint)sizes[1],
            /* MPI_Datatype oldtype  */ mpi_datatype,
            /* MPI_Datatype *newtype */ &(recvtypes[rank1])
        );
        /* commit recv datatype */
        MPI_Type_commit(&(recvtypes[rank1]));
        // check size of the committed datatype
        {
          int size;
          MPI_Type_size(recvtypes[rank1], &size);
          assert(size == (int)size_elem * chunk_jsize * sizes[2] * chunk_isize);
        }
        /* number of elements is 1, since we receive single "sendtypes[rank1]" */
        recvcounts[rank1] = 1;
        /* offset to the pointer to be received */
        rdispls[rank1] = (int)size_elem * sdecomp_kernel_get_offset(sizes[1], nprocs, rank1);
      }
    }
    // assign all to the struct
    plan = sdecomp_calloc(1, sizeof(sdecomp_transpose_t));
    plan->comm_2d    = comm_2d;
    plan->sendcounts = sendcounts;
    plan->recvcounts = recvcounts;
    plan->sdispls    = sdispls;
    plan->rdispls    = rdispls;
    plan->sendtypes  = sendtypes;
    plan->recvtypes  = recvtypes;
    plan->pencil_bef = pencil;
    plan->pencil_aft = (sdecomp_pencil_t)(((int)pencil+1) % 6);
    plan->gsizes     = sdecomp_calloc(SDECOMP_SDECOMP_NDIMS, sizeof(int));
    for(int dim = 0; dim < SDECOMP_SDECOMP_NDIMS; dim++){
      plan->gsizes[dim] = gsizes[dim];
    }
    plan->size_elem = size_elem;
  }
  return plan;
#undef SDECOMP_SDECOMP_NDIMS
}

/**
 * @brief initialise transpose plan for 3d, backward transposes
 * @param[in] sdecomp      : struct contains information of process distribution
 * @param[in] pencil       : type of pencil (one of six types)
 * @param[in] gsizes       : global array size in each dimension
 * @param[in] size_elem    : size of each element, e.g., sizeof(double)
 * @param[in] mpi_datatype : corresponding MPI_Datatype of each element, e.g., MPI_DOUBLE
 * @return                 : a pointer to the created plan (struct)
 */
static sdecomp_transpose_t *sdecomp_transpose_bwrd_init_3d(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const int *gsizes, const size_t size_elem, const MPI_Datatype mpi_datatype){
#define SDECOMP_SDECOMP_NDIMS 3
  assert(
      pencil == SDECOMP_X1PENCIL ||
      pencil == SDECOMP_Y1PENCIL ||
      pencil == SDECOMP_Z1PENCIL ||
      pencil == SDECOMP_Z2PENCIL ||
      pencil == SDECOMP_Y2PENCIL ||
      pencil == SDECOMP_X2PENCIL
  );
  const MPI_Comm comm_cart = sdecomp->comm_cart;
  // argument "gsizes" does NOT consider memory contiguous directions
  // convert them from inner memory to outer memory
  // from x1 (0) pencil or x2 (3) pencil : (gsizes[0], gsizes[1], gsizes[2])
  // from y1 (1) pencil or y2 (4) pencil : (gsizes[1], gsizes[2], gsizes[0])
  // from z1 (2) pencil or z2 (5) pencil : (gsizes[2], gsizes[0], gsizes[1])
  // -> condition based on mod(pencil, 3)
  int sizes[SDECOMP_SDECOMP_NDIMS] = {0};
  switch((int)pencil % 3){
    case 0:
      {
        sizes[0] = gsizes[0];
        sizes[1] = gsizes[1];
        sizes[2] = gsizes[2];
        break;
      }
    case 1:
      {
        sizes[0] = gsizes[1];
        sizes[1] = gsizes[2];
        sizes[2] = gsizes[0];
        break;
      }
    case 2:
      {
        sizes[0] = gsizes[2];
        sizes[1] = gsizes[0];
        sizes[2] = gsizes[1];
        break;
      }
    default:
      assert(0 == 1 /* should not be here */);
      break;
  }
  // NOTE: from now on, since the memory contiguity is taken into account,
  //   the memory is contiguous in x, middle in y, and the most sparse in z
  //   which are irrespective to the physical domain
  // create 2d communicator for each z (in terms of memory contiguity) chunk,
  //   and compute domain size in the unchanged dimension
  MPI_Comm comm_2d;
  {
    // we split the default communicator comm_cart and create a new one comm_2d,
    //   in which only processes participating in the all-to-all communication belong
    // comm_cart is x1 pencil ->
    //   we need to consider the unchanged dimsnion on x1 pencil
    // from x1 (0) pencil : y is unchanged -> for x1pencil, y
    // from y1 (1) pencil : x is unchanged -> for x1pencil, z
    // from z1 (2) pencil : z is unchanged -> for x1pencil, y
    // from x2 (3) pencil : y is unchanged -> for x1pencil, z
    // from y2 (4) pencil : x is unchanged -> for x1pencil, y
    // from z2 (5) pencil : z is unchanged -> for x1pencil, z
    // -> unchanged dimension in x1pencil differs for odd and even pencils
    const int unchanged_dim = (int)pencil % 2 == 0 ? 1 : 2;
    int remain_dims[SDECOMP_SDECOMP_NDIMS] = {1, 1, 1};
    remain_dims[unchanged_dim] = 0;
    MPI_Cart_sub(comm_cart, remain_dims, &comm_2d);
    // to compute mysize (number of grid points) in the unchanged dimension,
    //   we need to know the number of processes and my position in the dimension
    int nprocs, myrank;
    {
      const int ndims = 3;
      int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
      MPI_Cart_get(comm_cart, ndims, dims, periods, coords);
      nprocs = dims[unchanged_dim];
      myrank = coords[unchanged_dim];
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    // correct size in the unchanged dimension,
    //   since all-to-all does not care the global size
    sizes[1] = sdecomp_kernel_get_mysize(sizes[1], nprocs, myrank);
  }
  // create plan
  sdecomp_transpose_t *plan = NULL;
  {
    // number of processes which participate in the given 2d communicator "comm_2d"
    int nprocs;
    {
      const int ndims = 2; // 2d communicator
      int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
      MPI_Cart_get(comm_2d, ndims, dims, periods, coords);
      nprocs = dims[1]; // interested in y direction
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    // my position, "rank0"-th pencil
    int rank0;
    {
      const int ndims = 2; // 2d communicator
      int *dims    = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *periods = sdecomp_calloc((size_t)ndims, sizeof(int));
      int *coords  = sdecomp_calloc((size_t)ndims, sizeof(int));
      MPI_Cart_get(comm_2d, ndims, dims, periods, coords);
      rank0 = coords[1]; // interested in y direction
      sdecomp_free(dims);
      sdecomp_free(periods);
      sdecomp_free(coords);
    }
    int *sendcounts         = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *recvcounts         = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *sdispls            = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    int *rdispls            = sdecomp_calloc((size_t)nprocs, sizeof(         int));
    MPI_Datatype *sendtypes = sdecomp_calloc((size_t)nprocs, sizeof(MPI_Datatype));
    MPI_Datatype *recvtypes = sdecomp_calloc((size_t)nprocs, sizeof(MPI_Datatype));
    // consider communication between "rank0"-th and "rank1"-th pencils
    for(int rank1 = 0; rank1 < nprocs; rank1++){
      MPI_Datatype temptype;
      // send
      {
        int chunk_isize = sdecomp_kernel_get_mysize(sizes[0], nprocs, rank1);
        int chunk_ksize = sdecomp_kernel_get_mysize(sizes[2], nprocs, rank0);
        /* datatype to be sent: contiguous in z dimension */
        MPI_Type_create_hvector(
            /* int count             */ chunk_ksize,
            /* int blocklength       */ 1,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem * (MPI_Aint)sizes[0] * (MPI_Aint)sizes[1],
            /* MPI_Datatype oldtype  */ mpi_datatype,
            /* MPI_Datatype *newtype */ &temptype
        );
        /* datatype to be sent: previously defined type is repeated in x dimension */
        MPI_Type_create_hvector(
            /* int count             */ chunk_isize,
            /* int blocklength       */ 1,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem,
            /* MPI_Datatype oldtype  */ temptype,
            /* MPI_Datatype *newtype */ &(temptype)
        );
        /* datatype to be sent: previously defined type is repeated in y dimension */
        MPI_Type_create_hvector(
            /* int count             */ sizes[1],
            /* int blocklength       */ 1,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem * (MPI_Aint)sizes[0],
            /* MPI_Datatype oldtype  */ temptype,
            /* MPI_Datatype *newtype */ &(sendtypes[rank1])
        );
        /* commit send datatype */
        MPI_Type_commit(&(sendtypes[rank1]));
        // check size of the committed datatype
        {
          int size;
          MPI_Type_size(sendtypes[rank1], &size);
          assert(size == (int)size_elem * chunk_isize * sizes[1] * chunk_ksize);
        }
        /* number of elements is 1, since we send single "sendtypes[rank1]" */
        sendcounts[rank1] = 1;
        /* offset to the pointer to be sent */
        sdispls[rank1] = (int)size_elem * sdecomp_kernel_get_offset(sizes[0], nprocs, rank1);
      }
      // recv
      {
        int chunk_isize = sdecomp_kernel_get_mysize(sizes[0], nprocs, rank0);
        int chunk_ksize = sdecomp_kernel_get_mysize(sizes[2], nprocs, rank1);
        /* datatype to be received: simple since data being sent is already aligned in z direction */
        MPI_Type_create_hvector(
            /* int count             */ chunk_isize * sizes[1],
            /* int blocklength       */ chunk_ksize,
            /* MPI_Aint stride       */ (MPI_Aint)size_elem * (MPI_Aint)sizes[2],
            /* MPI_Datatype oldtype  */ mpi_datatype,
            /* MPI_Datatype *newtype */ &(recvtypes[rank1])
        );
        /* commit recv datatype */
        MPI_Type_commit(&(recvtypes[rank1]));
        // check size of the committed datatype
        {
          int size;
          MPI_Type_size(recvtypes[rank1], &size);
          assert(size == (int)size_elem * chunk_ksize * chunk_isize * sizes[1]);
        }
        /* number of elements is 1, since we receive single "sendtypes[rank1]" */
        recvcounts[rank1] = 1;
        /* offset to the pointer to be received */
        rdispls[rank1] = (int)size_elem * sdecomp_kernel_get_offset(sizes[2], nprocs, rank1);
      }
    }
    // assign all to the struct
    plan = sdecomp_calloc(1, sizeof(sdecomp_transpose_t));
    plan->comm_2d    = comm_2d;
    plan->sendcounts = sendcounts;
    plan->recvcounts = recvcounts;
    plan->sdispls    = sdispls;
    plan->rdispls    = rdispls;
    plan->sendtypes  = sendtypes;
    plan->recvtypes  = recvtypes;
    plan->pencil_bef = pencil;
    plan->pencil_aft = (sdecomp_pencil_t)(((int)pencil+5) % 6);
    plan->gsizes     = sdecomp_calloc(SDECOMP_SDECOMP_NDIMS, sizeof(int));
    for(int dim = 0; dim < SDECOMP_SDECOMP_NDIMS; dim++){
      plan->gsizes[dim] = gsizes[dim];
    }
    plan->size_elem = size_elem;
  }
  return plan;
#undef SDECOMP_SDECOMP_NDIMS
}

/**
 * @brief initialise forward transpose plan
 * @param[in] sdecomp      : struct contains information of process distribution
 * @param[in] pencil       : type of pencil (SDECOMP_X1PENCIL or SDECOMP_Y1PENCIL)
 * @param[in] gsizes       : global array size in each dimension
 * @param[in] size_elem    : size of each element, e.g., sizeof(double)
 * @param[in] mpi_datatype : corresponding MPI_Datatype of each element, e.g., MPI_DOUBLE
 * @return                 : a pointer to the created plan (struct)
 */
sdecomp_transpose_t *sdecomp_transpose_fwrd_init(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const int *gsizes, const size_t size_elem, const MPI_Datatype mpi_datatype){
  const int ndims = sdecomp->ndims;
  if(ndims == 2){
    return sdecomp_transpose_init_2d(sdecomp, pencil, gsizes, size_elem, mpi_datatype);
  }else{
    return sdecomp_transpose_fwrd_init_3d(sdecomp, pencil, gsizes, size_elem, mpi_datatype);
  }
}

/**
 * @brief initialise backward transpose plan
 * @param[in] sdecomp      : struct contains information of process distribution
 * @param[in] pencil       : type of pencil (SDECOMP_X1PENCIL or SDECOMP_Y1PENCIL)
 * @param[in] gsizes       : global array size in each dimension
 * @param[in] size_elem    : size of each element, e.g., sizeof(double)
 * @param[in] mpi_datatype : corresponding MPI_Datatype of each element, e.g., MPI_DOUBLE
 * @return                 : a pointer to the created plan (struct)
 */
sdecomp_transpose_t *sdecomp_transpose_bwrd_init(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const int *gsizes, const size_t size_elem, const MPI_Datatype mpi_datatype){
  const int ndims = sdecomp->ndims;
  if(ndims == 2){
    return sdecomp_transpose_init_2d(sdecomp, pencil, gsizes, size_elem, mpi_datatype);
  }else{
    return sdecomp_transpose_bwrd_init_3d(sdecomp, pencil, gsizes, size_elem, mpi_datatype);
  }
}

/**
 * @brief execute transpose
 * @param[in]  plan    : transpose plan initialised by constructor
 *                         (sdecomp_transpose_fwrd_init or sdecomp_transpose_bwrd_init)
 * @param[in]  sendbuf : pointer to the input  buffer
 * @param[out] recvbuf : pointer to the output buffer
 * @return             : error code
 */
int sdecomp_transpose_execute(sdecomp_transpose_t *plan, const void * restrict sendbuf, void * restrict recvbuf){
  MPI_Alltoallw(
      sendbuf, plan->sendcounts, plan->sdispls, plan->sendtypes,
      recvbuf, plan->recvcounts, plan->rdispls, plan->recvtypes,
      plan->comm_2d
  );
  return 0;
}

/**
 * @brief finalise transpose plan
 * @param[in,out] plan : transpose plan to be cleaned-up
 * @return             : error code
 */
int sdecomp_transpose_finalise(sdecomp_transpose_t *plan){
  int nprocs;
  MPI_Comm_size(plan->comm_2d, &nprocs);
  MPI_Comm_free(&(plan->comm_2d));
  sdecomp_free(plan->sendcounts);
  sdecomp_free(plan->recvcounts);
  sdecomp_free(plan->sdispls);
  sdecomp_free(plan->rdispls);
  for(int n = 0; n < nprocs; n++){
    MPI_Type_free(&(plan->sendtypes[n]));
    MPI_Type_free(&(plan->recvtypes[n]));
  }
  sdecomp_free(plan->sendtypes);
  sdecomp_free(plan->recvtypes);
  sdecomp_free(plan->gsizes);
  sdecomp_free(plan);
  return 0;
}

/***** TEST *****/

/*
 * Extremely dirty test implementations
 *
 * Concept:
 *
 *     Say this is the x pencil
 *
 *     +---------------+
 *     | 12 13   14 15 |
 *     | 08 09   10 11 |
 *     +---------------+
 *     | 04 05   06 07 |
 *     | 00 01   02 03 |
 *     +---------------+
 *
 *     which are in memory of rank 0: 00 01 02 03 04 05 06 07
 *
 *     y1 pencil after being transposed should be
 *
 *     +-------+-------+
 *     | 12 13 | 14 15 |
 *     | 08 09 | 10 11 |
 *     |       |       |
 *     | 04 05 | 06 07 |
 *     | 00 01 | 02 03 |
 *     +-------+-------+
 *
 *     which are in memory of rank 0: 00 04 08 12 01 05 09 13
 *
 *     The contents are know a priori using sizes and offsets.
 *
 * These things are tested in this source.
 */

#define MSG

#define TAIL \
  sdecomp_free(buf_bef); \
  sdecomp_free(buf_aft);

#define HEAD \
  const int mysizes_bef[] = { \
    sdecomp_get_pencil_mysize(sdecomp, pencil_bef, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_mysize(sdecomp, pencil_bef, SDECOMP_YDIR, glsizes[1]) \
  }; \
  const int offsets_bef[] = { \
    sdecomp_get_pencil_offset(sdecomp, pencil_bef, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_offset(sdecomp, pencil_bef, SDECOMP_YDIR, glsizes[1]) \
  }; \
  const int mysizes_aft[] = { \
    sdecomp_get_pencil_mysize(sdecomp, pencil_aft, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_mysize(sdecomp, pencil_aft, SDECOMP_YDIR, glsizes[1]) \
  }; \
  const int offsets_aft[] = { \
    sdecomp_get_pencil_offset(sdecomp, pencil_aft, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_offset(sdecomp, pencil_aft, SDECOMP_YDIR, glsizes[1]) \
  }; \
  unsigned char *buf_bef = sdecomp_calloc((size_t)mysizes_bef[0] * (size_t)mysizes_bef[1], size_elem/sizeof(unsigned char)); \
  unsigned char *buf_aft = sdecomp_calloc((size_t)mysizes_aft[0] * (size_t)mysizes_aft[1], size_elem/sizeof(unsigned char));

// local values are independent of the pencil directions
#define VALUE(befaft) \
  unsigned char value = (unsigned char) ( \
    + (j + offsets_##befaft[1]) * glsizes[0] \
    + (i + offsets_##befaft[0]) \
  );

int sdecomp_test_transpose_2d_x1_to_y1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_X1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Y1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  HEAD
  // assign values to bef
  {
    for(int j = 0; j < mysizes_bef[1]; j++){
      for(int i = 0; i < mysizes_bef[0]; i++){
        int index = j * mysizes_bef[0] + i;
        VALUE(bef)
        for(size_t n = 0; n < size_elem; n++){
          memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int i = 0; i < mysizes_aft[0]; i++){
      for(int j = 0; j < mysizes_aft[1]; j++){
        int index = i * mysizes_aft[1] + j;
        VALUE(aft)
        for(size_t n = 0; n < size_elem; n++){
          assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
        }
      }
    }
  }
  TAIL
  MSG
  return 0;
}

int sdecomp_test_transpose_2d_y1_to_x1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Y1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_X1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  HEAD
  // assign values to bef
  {
    for(int i = 0; i < mysizes_bef[0]; i++){
      for(int j = 0; j < mysizes_bef[1]; j++){
        int index = i * mysizes_bef[1] + j;
        VALUE(bef)
        for(size_t n = 0; n < size_elem; n++){
          memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int j = 0; j < mysizes_aft[1]; j++){
      for(int i = 0; i < mysizes_aft[0]; i++){
        int index = j * mysizes_aft[0] + i;
        VALUE(aft)
        for(size_t n = 0; n < size_elem; n++){
          assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
        }
      }
    }
  }
  TAIL
  MSG
  return 0;
}

#undef HEAD
#undef VALUE

#define HEAD \
  const int mysizes_bef[] = { \
    sdecomp_get_pencil_mysize(sdecomp, pencil_bef, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_mysize(sdecomp, pencil_bef, SDECOMP_YDIR, glsizes[1]), \
    sdecomp_get_pencil_mysize(sdecomp, pencil_bef, SDECOMP_ZDIR, glsizes[2]) \
  }; \
  const int offsets_bef[] = { \
    sdecomp_get_pencil_offset(sdecomp, pencil_bef, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_offset(sdecomp, pencil_bef, SDECOMP_YDIR, glsizes[1]), \
    sdecomp_get_pencil_offset(sdecomp, pencil_bef, SDECOMP_ZDIR, glsizes[2]) \
  }; \
  const int mysizes_aft[] = { \
    sdecomp_get_pencil_mysize(sdecomp, pencil_aft, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_mysize(sdecomp, pencil_aft, SDECOMP_YDIR, glsizes[1]), \
    sdecomp_get_pencil_mysize(sdecomp, pencil_aft, SDECOMP_ZDIR, glsizes[2]) \
  }; \
  const int offsets_aft[] = { \
    sdecomp_get_pencil_offset(sdecomp, pencil_aft, SDECOMP_XDIR, glsizes[0]), \
    sdecomp_get_pencil_offset(sdecomp, pencil_aft, SDECOMP_YDIR, glsizes[1]), \
    sdecomp_get_pencil_offset(sdecomp, pencil_aft, SDECOMP_ZDIR, glsizes[2]) \
  }; \
  unsigned char *buf_bef = sdecomp_calloc((size_t)mysizes_bef[0] * (size_t)mysizes_bef[1] * (size_t)mysizes_bef[2], size_elem/sizeof(unsigned char)); \
  unsigned char *buf_aft = sdecomp_calloc((size_t)mysizes_aft[0] * (size_t)mysizes_aft[1] * (size_t)mysizes_aft[2], size_elem/sizeof(unsigned char));

#define VALUE(befaft) \
  unsigned char value = (unsigned char) ( \
    + (k + offsets_##befaft[2]) * glsizes[1] * glsizes[0] \
    + (j + offsets_##befaft[1]) * glsizes[0] \
    + (i + offsets_##befaft[0]) \
  );

static int sdecomp_test_transpose_3d_x_to_y(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan, const int *glsizes, const size_t size_elem, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft){
  HEAD
  // assign values to bef
  {
    for(int k = 0; k < mysizes_bef[2]; k++){
      for(int j = 0; j < mysizes_bef[1]; j++){
        for(int i = 0; i < mysizes_bef[0]; i++){
          int index = k * mysizes_bef[1] * mysizes_bef[0] + j * mysizes_bef[0] + i;
          VALUE(bef)
          for(size_t n = 0; n < size_elem; n++){
            memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
          }
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int i = 0; i < mysizes_aft[0]; i++){
      for(int k = 0; k < mysizes_aft[2]; k++){
        for(int j = 0; j < mysizes_aft[1]; j++){
          int index = i * mysizes_aft[2] * mysizes_aft[1] + k * mysizes_aft[1] + j;
          VALUE(aft)
          for(size_t n = 0; n < size_elem; n++){
            assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
          }
        }
      }
    }
  }
  TAIL
  return 0;
}

static int sdecomp_test_transpose_3d_y_to_x(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan, const int *glsizes, const size_t size_elem, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft){
  HEAD
  // assign values to bef
  {
    for(int i = 0; i < mysizes_bef[0]; i++){
      for(int k = 0; k < mysizes_bef[2]; k++){
        for(int j = 0; j < mysizes_bef[1]; j++){
          int index = i * mysizes_bef[2] * mysizes_bef[1] + k * mysizes_bef[1] + j;
          VALUE(bef)
          for(size_t n = 0; n < size_elem; n++){
            memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
          }
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int k = 0; k < mysizes_aft[2]; k++){
      for(int j = 0; j < mysizes_aft[1]; j++){
        for(int i = 0; i < mysizes_aft[0]; i++){
          int index = k * mysizes_aft[1] * mysizes_aft[0] + j * mysizes_aft[0] + i;
          VALUE(aft)
          for(size_t n = 0; n < size_elem; n++){
            assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
          }
        }
      }
    }
  }
  TAIL
  return 0;
}

static int sdecomp_test_transpose_3d_y_to_z(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan, const int *glsizes, const size_t size_elem, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft){
  HEAD
  // assign values to bef
  {
    for(int i = 0; i < mysizes_bef[0]; i++){
      for(int k = 0; k < mysizes_bef[2]; k++){
        for(int j = 0; j < mysizes_bef[1]; j++){
          int index = i * mysizes_bef[2] * mysizes_bef[1] + k * mysizes_bef[1] + j;
          VALUE(bef)
          for(size_t n = 0; n < size_elem; n++){
            memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
          }
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int j = 0; j < mysizes_aft[1]; j++){
      for(int i = 0; i < mysizes_aft[0]; i++){
        for(int k = 0; k < mysizes_aft[2]; k++){
          int index = j * mysizes_aft[0] * mysizes_aft[2] + i * mysizes_aft[2] + k;
          VALUE(aft)
          for(size_t n = 0; n < size_elem; n++){
            assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
          }
        }
      }
    }
  }
  TAIL
  return 0;
}

static int sdecomp_test_transpose_3d_z_to_y(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan, const int *glsizes, const size_t size_elem, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft){
  HEAD
  // assign values to bef
  {
    for(int j = 0; j < mysizes_bef[1]; j++){
      for(int i = 0; i < mysizes_bef[0]; i++){
        for(int k = 0; k < mysizes_bef[2]; k++){
          int index = j * mysizes_bef[0] * mysizes_bef[2] + i * mysizes_bef[2] + k;
          VALUE(bef)
          for(size_t n = 0; n < size_elem; n++){
            memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
          }
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int i = 0; i < mysizes_aft[0]; i++){
      for(int k = 0; k < mysizes_aft[2]; k++){
        for(int j = 0; j < mysizes_aft[1]; j++){
          int index = i * mysizes_aft[2] * mysizes_aft[1] + k * mysizes_aft[1] + j;
          VALUE(aft)
          for(size_t n = 0; n < size_elem; n++){
            assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
          }
        }
      }
    }
  }
  TAIL
  return 0;
}

static int sdecomp_test_transpose_3d_z_to_x(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan, const int *glsizes, const size_t size_elem, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft){
  HEAD
  // assign values to bef
  {
    for(int j = 0; j < mysizes_bef[1]; j++){
      for(int i = 0; i < mysizes_bef[0]; i++){
        for(int k = 0; k < mysizes_bef[2]; k++){
          int index = j * mysizes_bef[0] * mysizes_bef[2] + i * mysizes_bef[2] + k;
          VALUE(bef)
          for(size_t n = 0; n < size_elem; n++){
            memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
          }
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int k = 0; k < mysizes_aft[2]; k++){
      for(int j = 0; j < mysizes_aft[1]; j++){
        for(int i = 0; i < mysizes_aft[0]; i++){
          int index = k * mysizes_aft[1] * mysizes_aft[0] + j * mysizes_aft[0] + i;
          VALUE(aft)
          for(size_t n = 0; n < size_elem; n++){
            assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
          }
        }
      }
    }
  }
  TAIL
  return 0;
}

static int sdecomp_test_transpose_3d_x_to_z(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan, const int *glsizes, const size_t size_elem, const sdecomp_pencil_t pencil_bef, const sdecomp_pencil_t pencil_aft){
  HEAD
  // assign values to bef
  {
    for(int k = 0; k < mysizes_bef[2]; k++){
      for(int j = 0; j < mysizes_bef[1]; j++){
        for(int i = 0; i < mysizes_bef[0]; i++){
          int index = k * mysizes_bef[1] * mysizes_bef[0] + j * mysizes_bef[0] + i;
          VALUE(bef)
          for(size_t n = 0; n < size_elem; n++){
            memcpy(&buf_bef[size_elem * (size_t)index + n], &value, sizeof(unsigned char));
          }
        }
      }
    }
  }
  // transpose
  sdecomp_transpose_execute(plan, buf_bef, buf_aft);
  // check aft
  {
    for(int j = 0; j < mysizes_aft[1]; j++){
      for(int i = 0; i < mysizes_aft[0]; i++){
        for(int k = 0; k < mysizes_aft[2]; k++){
          int index = j * mysizes_aft[0] * mysizes_aft[2] + i * mysizes_aft[2] + k;
          VALUE(aft)
          for(size_t n = 0; n < size_elem; n++){
            assert(memcmp(&buf_aft[size_elem * (size_t)index + n], &value, sizeof(unsigned char)) == 0);
          }
        }
      }
    }
  }
  TAIL
  return 0;
}

/* 3d forward tests */

int sdecomp_test_transpose_3d_x1_to_y1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_X1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Y1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_x_to_y(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_x2_to_y2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_X2PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Y2PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_x_to_y(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_y1_to_z1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Y1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Z1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_y_to_z(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_y2_to_z2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Y2PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Z2PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_y_to_z(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_z1_to_x2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Z1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_X2PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_z_to_x(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_z2_to_x1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Z2PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_X1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_z_to_x(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

/* 3d backward tests */

int sdecomp_test_transpose_3d_x1_to_z2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_X1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Z2PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_x_to_z(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_x2_to_z1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_X2PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Z1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_x_to_z(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_z1_to_y1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Z1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Y1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_z_to_y(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_z2_to_y2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Z2PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_Y2PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_z_to_y(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_y1_to_x1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Y1PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_X1PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_y_to_x(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

int sdecomp_test_transpose_3d_y2_to_x2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan){
  const sdecomp_pencil_t pencil_bef = SDECOMP_Y2PENCIL;
  const sdecomp_pencil_t pencil_aft = SDECOMP_X2PENCIL;
  const int *glsizes = plan->gsizes;
  const size_t size_elem = plan->size_elem;
  assert(pencil_bef == plan->pencil_bef);
  assert(pencil_aft == plan->pencil_aft);
  sdecomp_test_transpose_3d_y_to_x(sdecomp, plan, glsizes, size_elem, pencil_bef, pencil_aft);
  MSG
  return 0;
}

#undef HEAD
#undef VALUE
#undef TAIL
#undef MSG

