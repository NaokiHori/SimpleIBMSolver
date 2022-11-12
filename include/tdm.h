#if !defined(TDM_H)
#define TDM_H

#include <stdbool.h>
#include <fftw3.h>

extern int tdm_solve_double(const int n, const int m, const bool is_periodic, const double * restrict l, const double * restrict c, const double * restrict u, double * restrict q);
extern int tdm_solve_fftw_complex(const int n, const int m, const bool is_periodic, const double * restrict l, const double * restrict c, const double * restrict u, fftw_complex * restrict q);
extern int tdm_cleanup(void);

#endif // TDM_H
