#include <assert.h>
#include <mpi.h>
#include "domain.h"
#include "arrays/fluid.h"


#if NDIMS == 2

/**
 * @brief update halo cells of ux-like array
 * @param[in   ] domain  : information about domain decomposition and size
 * @param[inout] ux : array whose shape is ux-like
 * @return          : error code
 */
int domain_communicate_halo_ux_like(const domain_t *domain, double *ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  /* ! update y halo values ! 29 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int ymrank, yprank;
    MPI_Cart_shift(comm, 1, 1, &ymrank, &yprank);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_contiguous(isize+1, MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    {
      int size;
      MPI_Type_size(dtype, &size);
      assert((isize+1) * sizeof(double) == (size_t)size);
    }
    // send in positive direction
    MPI_Sendrecv(
      &UX(1,   jsize), 1, dtype, yprank, 0,
      &UX(1,       0), 1, dtype, ymrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &UX(1,       1), 1, dtype, ymrank, 0,
      &UX(1, jsize+1), 1, dtype, yprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief update halo cells of ux-like array
 * @param[in   ] domain  : information about domain decomposition and size
 * @param[inout] ux : array whose shape is ux-like
 * @return          : error code
 */
int domain_communicate_halo_ux_like(const domain_t *domain, double *ux){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! update y halo values ! 29 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int ymrank, yprank;
    MPI_Cart_shift(comm, 1, 1, &ymrank, &yprank);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(ksize, isize+1, (isize+1)*(jsize+2), MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    {
      int size;
      MPI_Type_size(dtype, &size);
      assert((isize+1) * ksize * sizeof(double) == (size_t)size);
    }
    // send in positive direction
    MPI_Sendrecv(
      &UX(1,   jsize, 1), 1, dtype, yprank, 0,
      &UX(1,       0, 1), 1, dtype, ymrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &UX(1,       1, 1), 1, dtype, ymrank, 0,
      &UX(1, jsize+1, 1), 1, dtype, yprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  /* ! update z halo values ! 29 ! */
  {
    const MPI_Comm comm = domain->sdecomp->comm_cart;
    // check neighbour ranks, negative and positive
    int zmrank, zprank;
    MPI_Cart_shift(comm, 2, 1, &zmrank, &zprank);
    // create datatype
    MPI_Datatype dtype;
    MPI_Type_vector(jsize, isize+1, isize+1, MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    {
      int size;
      MPI_Type_size(dtype, &size);
      assert((isize+1) * jsize * sizeof(double) == (size_t)size);
    }
    // send in positive direction
    MPI_Sendrecv(
      &UX(1, 1,   ksize), 1, dtype, zprank, 0,
      &UX(1, 1,       0), 1, dtype, zmrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &UX(1, 1,       1), 1, dtype, zmrank, 0,
      &UX(1, 1, ksize+1), 1, dtype, zprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  return 0;
}

#endif // NDIMS
