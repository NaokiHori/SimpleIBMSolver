#include <assert.h>
#include <mpi.h>
#include "domain.h"
#include "arrays/fluid.h"


#if NDIMS == 2

/**
 * @brief update halo cells of arrays whose shape is NOT ux-like
 * @param[in   ] domain  : information about domain decomposition and size
 * @param[inout] uy : array whose shape is NOT ux-like
 * @return          : error code
 */
static int communicate_halo_others(const domain_t *domain, double *uy){
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
    MPI_Type_contiguous(isize, MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    {
      int size;
      MPI_Type_size(dtype, &size);
      assert(isize * sizeof(double) == (size_t)size);
    }
    // send in positive direction
    MPI_Sendrecv(
      &UY(1,   jsize), 1, dtype, yprank, 0,
      &UY(1,       0), 1, dtype, ymrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &UY(1,       1), 1, dtype, ymrank, 0,
      &UY(1, jsize+1), 1, dtype, yprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  return 0;
}

#else // NDIMS == 3

/**
 * @brief update halo cells of arrays whose shape is NOT ux-like
 * @param[in   ] domain  : information about domain decomposition and size
 * @param[inout] uy : array whose shape is NOT ux-like
 * @return          : error code
 */
static int communicate_halo_others(const domain_t *domain, double *uy){
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
    MPI_Type_vector(ksize+2, isize, (isize+2)*(jsize+2), MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    {
      int size;
      MPI_Type_size(dtype, &size);
      assert(isize * (ksize+2) * sizeof(double) == (size_t)size);
    }
    // send in positive direction
    MPI_Sendrecv(
      &UY(1,   jsize, 0), 1, dtype, yprank, 0,
      &UY(1,       0, 0), 1, dtype, ymrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &UY(1,       1, 0), 1, dtype, ymrank, 0,
      &UY(1, jsize+1, 0), 1, dtype, yprank, 0,
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
    MPI_Type_vector(jsize+2, isize, isize+2, MPI_DOUBLE, &dtype);
    MPI_Type_commit(&dtype);
    {
      int size;
      MPI_Type_size(dtype, &size);
      assert(isize * (jsize+2) * sizeof(double) == (size_t)size);
    }
    // send in positive direction
    MPI_Sendrecv(
      &UY(1, 0,   ksize), 1, dtype, zprank, 0,
      &UY(1, 0,       0), 1, dtype, zmrank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // send in negative direction
    MPI_Sendrecv(
      &UY(1, 0,       1), 1, dtype, zmrank, 0,
      &UY(1, 0, ksize+1), 1, dtype, zprank, 0,
      comm, MPI_STATUS_IGNORE
    );
    // clean-up used datatype
    MPI_Type_free(&dtype);
  }
  return 0;
}

#endif // NDIMS

/**
 * @brief update halo cells of uy-like array
 * @param[in   ] domain  : information about domain decomposition and size
 * @param[inout] uy_like : array whose shape is uy-like
 * @return               : error code
 */
int domain_communicate_halo_uy_like(const domain_t *domain, double *uy_like){
  // just pass arguments to a common function defined above,
  // since uy, uz, p have same the shape
  communicate_halo_others(domain, uy_like);
  return 0;
}

#if NDIMS == 3

/**
 * @brief update halo cells of uz-like array
 * @param[in   ] domain  : information about domain decomposition and size
 * @param[inout] uz_like : array whose shape is uz-like
 * @return               : error code
 */
int domain_communicate_halo_uz_like(const domain_t *domain, double *uz_like){
  // just pass arguments to a common function defined above,
  // since uy, uz, p have same the shape
  communicate_halo_others(domain, uz_like);
  return 0;
}

#endif

/**
 * @brief update halo cells of p-like array
 * @param[in   ] domain  : information about domain decomposition and size
 * @param[inout] p_like : array whose shape is p-like
 * @return               : error code
 */
int domain_communicate_halo_p_like(const domain_t *domain, double *p_like){
  // just pass arguments to a common function defined above,
  // since uy, uz, p have same the shape
  communicate_halo_others(domain, p_like);
  return 0;
}

