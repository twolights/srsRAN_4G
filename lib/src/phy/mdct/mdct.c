//
// Created by ykchen on 9/8/24.
//

#include <stdlib.h>
#include "srsran/phy/mdct/mdct.h"
#include "srsran/phy/mdct/differential_prod.h"


// XXX x_tilde with different psi's
static cf_t calculate_D(const cf_t* y, const cf_t* x_tilde, cf_t* temp, cf_t* output, uint32_t n, uint32_t Q, uint32_t psi)
{
  int32_t d = (int32_t)(1 + (Q * psi));
  differential_product(y, temp, d, n);
  srsran_vec_prod_conj_ccc(temp, x_tilde, output, n);
  return srsran_vec_acc_cc(output, n);
}

static cf_t calculate_C(const cf_t* y, const cf_t* x, cf_t* temp, cf_t* output, uint32_t n, uint32_t Q, uint32_t PSI, uint32_t tau)
{
  cf_t result = 0.0;
  for (int i = 0; i < PSI; i++) {
//    cf_t* x_tilde = x + i * n;
    result += calculate_D(y, x, temp, output, n, Q, i);
  }
  return result / (float)PSI;
}

int srsran_prepare_pss_mdct(srsran_pss_mdct_t* mdct, uint32_t n, uint32_t Q, uint32_t PSI)
{
  // TODO allocate x_tilde
  mdct->output = (cf_t*)malloc(n * sizeof(cf_t));
  if (mdct->output == NULL) {
    return SRSRAN_ERROR;
  }
  mdct->temp = (cf_t*)malloc(n * sizeof(cf_t));
  if (mdct->temp == NULL) {
    free(mdct->output);
    return SRSRAN_ERROR;
  }
  mdct->n = n;
  mdct->Q = Q;
  mdct->PSI = PSI;
  return SRSRAN_SUCCESS;
}

int srsran_destroy_pss_mdct(srsran_pss_mdct_t* mdct)
{
  if (mdct == NULL) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  // TODO free x_tilde
  free(mdct->temp);
  mdct->temp = NULL;
  free(mdct->output);
  mdct->output = NULL;
  return SRSRAN_SUCCESS;
}

SRSRAN_API int mdct_detect_pss(const srsran_pss_mdct_t* mdct,
                               const cf_t* in, uint32_t nof_samples,
                               srsran_pss_mdct_detect_res_t* result)
{
  return SRSRAN_SUCCESS;
}
