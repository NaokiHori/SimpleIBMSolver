#include "common.h"
#include "domain.h"
#include "fileio.h"
#include "internal.h"



/**
 * @brief save members in domain_t
 * @param[in] dirname : name of directory to which *.npy files will be written
 * @param[in] domain  : domain information
 * @return            : error code
 */
int domain_save(const char dirname[], const domain_t *domain){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const double lx = domain->lengths[0];
  const double ly = domain->lengths[1];
  const double *xf = domain->xf;
  const double *xc = domain->xc;
  const double dy = domain->dy;
  // generate y coordinate, which might be useful in post-processings
  double *yf = common_calloc(gljsize+1, sizeof(double));
  double *yc = common_calloc(gljsize,   sizeof(double));
  {
    for(int j = 1; j <= gljsize+1; j++){
      yf[j-1] = 1.*(j-1)*dy;
    }
    for(int j = 1; j <= gljsize; j++){
      yc[j-1] = 0.5*(yf[j-1]+yf[j  ]);
    }
  }
  fileio_w_0d_serial(dirname, "glisize", NPYIO_INT,    sizeof(int),    &glisize);
  fileio_w_0d_serial(dirname, "gljsize", NPYIO_INT,    sizeof(int),    &gljsize);
  fileio_w_0d_serial(dirname, "lx",      NPYIO_DOUBLE, sizeof(double), &lx     );
  fileio_w_0d_serial(dirname, "ly",      NPYIO_DOUBLE, sizeof(double), &ly     );
  fileio_w_0d_serial(dirname, "dy",      NPYIO_DOUBLE, sizeof(double), &dy     );
  fileio_w_1d_serial(dirname, "xf",      NPYIO_DOUBLE, sizeof(double), glisize+1, xf);
  fileio_w_1d_serial(dirname, "xc",      NPYIO_DOUBLE, sizeof(double), glisize+2, xc);
  fileio_w_1d_serial(dirname, "yf",      NPYIO_DOUBLE, sizeof(double), gljsize+1, yf);
  fileio_w_1d_serial(dirname, "yc",      NPYIO_DOUBLE, sizeof(double), gljsize,   yc);
  common_free(yf);
  common_free(yc);
  return 0;
}

