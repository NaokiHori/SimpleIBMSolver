#if !defined(SDECOMP_H)
#define SDECOMP_H

#include <stddef.h> // size_t
#include <mpi.h>

/* ! definition of sdecomp_t ! 10 ! */
/** @struct sdecomp_t
 *  @brief struct storing information about domain decomposition
 *  @var ndims     : number of spatial dimensions
 *  @var comm_cart : default Cartesian communicator (x1 pencil)
 */
typedef struct {
  MPI_Comm comm_cart;
  int ndims;
  int padding; // to suppress warning -Wpadded
} sdecomp_t;

// direction
typedef enum {
  SDECOMP_XDIR = 0,
  SDECOMP_YDIR = 1,
  SDECOMP_ZDIR = 2
} sdecomp_dir_t;

typedef enum {
  SDECOMP_X1PENCIL = 0,
  SDECOMP_Y1PENCIL = 1,
  SDECOMP_Z1PENCIL = 2,
  SDECOMP_X2PENCIL = 3,
  SDECOMP_Y2PENCIL = 4,
  SDECOMP_Z2PENCIL = 5
} sdecomp_pencil_t;

/* constructor and destructor */
extern sdecomp_t *sdecomp_init(const MPI_Comm comm_default, const int ndims, const int *dims, const int *periods);
extern int sdecomp_finalise(sdecomp_t *sdecomp);

/* get number of processes in one dimension */
extern int sdecomp_get_nprocs(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir);
/* get position of my process */
extern int sdecomp_get_myrank(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir);

/* get local size and offset of my pencil */
extern int sdecomp_get_pencil_mysize(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const int num);
extern int sdecomp_get_pencil_offset(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const sdecomp_dir_t dir, const int num);

/* kernel functions to decide local size / offset of pencils */
extern int sdecomp_kernel_get_mysize(const int num_total, const int nprocs, const int myrank);
extern int sdecomp_kernel_get_offset(const int num_total, const int nprocs, const int myrank);

/* for parallel matrix transpose */
// declaration of struct
typedef struct sdecomp_transpose_t_ sdecomp_transpose_t;
// constructors
extern sdecomp_transpose_t *sdecomp_transpose_fwrd_init(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const int *gsizes, const size_t size_elem, const MPI_Datatype mpi_datatype);
extern sdecomp_transpose_t *sdecomp_transpose_bwrd_init(const sdecomp_t *sdecomp, const sdecomp_pencil_t pencil, const int *gsizes, const size_t size_elem, const MPI_Datatype mpi_datatype);
// executor
extern int sdecomp_transpose_execute(sdecomp_transpose_t *plan, const void * restrict sendbuf, void * restrict recvbuf);
// finaliser
extern int sdecomp_transpose_finalise(sdecomp_transpose_t *plan);

/* test functions */
// 2d parallel transpose
extern int sdecomp_test_transpose_2d_x1_to_y1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_2d_y1_to_x1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
// 3d paralsdecomp_lel transpose, forward
extern int sdecomp_test_transpose_3d_x1_to_y1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_y1_to_z1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_z1_to_x2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_x2_to_y2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_y2_to_z2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_z2_to_x1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
// 3d paralsdecomp_lel transpose, backward
extern int sdecomp_test_transpose_3d_x1_to_z2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_z2_to_y2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_y2_to_x2(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_x2_to_z1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_z1_to_y1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);
extern int sdecomp_test_transpose_3d_y1_to_x1(const sdecomp_t *sdecomp, sdecomp_transpose_t *plan);

/* memory management */
extern void *sdecomp_calloc(const size_t count, const size_t size);
extern void sdecomp_free(void *ptr);

#endif // SDECOMP_H
