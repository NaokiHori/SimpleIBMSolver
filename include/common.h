#if !defined(COMMON_H)
#define COMMON_H

#if !defined(NDIMS)
#error "define NDIMS: -DNDIMS=2 or -DNDIMS=3"
#endif

#include <stdlib.h> // size_t
#include <math.h>
// C99 does not specify M_PI in math.h
#if !defined(M_PI)
#define M_PI 3.1415926535897932
#endif

/* Runge-Kutta configurations */
#define RKSTEPMAX 3
typedef struct {
  double alpha;
  double beta;
  double gamma;
} rkcoef_t;
extern const rkcoef_t RKCOEFS[RKSTEPMAX];

/* directory names */
#define NDIGITS_STEP 10
#define DIRNAME_SAVE "output/save"
#define DIRNAME_LOG  "output/log"
#define DIRNAME_STAT "output/stat"

extern void *common_calloc(const size_t count, const size_t size);
extern void common_free(void *ptr);
extern double common_get_wtime(void);
extern int common_get_ndigits(int num);

#endif // COMMON_H
