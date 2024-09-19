//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/mdct.h"
#include "srsran/phy/dft/dft.h"
#include "srsran/phy/mdct/differential_prod.h"
#include <stdlib.h>
#include <string.h>

// XXX copied from pss_nr.c, should be in a common header
#define PSS_NR_SUBC_BEGIN 56

static cf_t calculate_D(const cf_t* y, const cf_t* x_tilde, cf_t* temp, cf_t* output, uint32_t n, uint32_t Q, uint32_t psi)
{
  int32_t d = (int32_t)(1 + (Q * psi));
  differential_product(y, temp, d, n);
  srsran_vec_prod_conj_ccc(temp, x_tilde, output, n);
  return srsran_vec_acc_cc(output, n);
}

static cf_t calculate_C(const srsran_pss_mdct_t* mdct, const cf_t* y, uint32_t N_id_2, int32_t tau)
{
  cf_t result = 0.0;
  for (int psi = 0; psi < mdct->PSI; psi++) {
    result += calculate_D(y, mdct->x_tilde[N_id_2][psi], mdct->temp, mdct->output, mdct->symbol_sz, mdct->Q, psi);
  }
  return result / (float)mdct->PSI;
}

// x_tilde_{n,d}^{r}: d: delay, r: NID2
//static void fill_x_tilde(cf_t* out, uint32_t symbol_sz, uint32_t n, int32_t d, uint32_t r)
static void fill_x_tilde(srsran_pss_mdct_t* mdct, uint32_t r, uint32_t psi)
{
  uint32_t d = (uint32_t)(1 + (mdct->Q * psi));
  // TODO: Should be optimized by multiplying the constant phase instead of computing the diff. product for each N_id_2
  differential_product(mdct->pss_x[r], mdct->x_tilde[r][psi], (int)d, mdct->symbol_sz);
}

static void prepare_pss_x(srsran_pss_mdct_t* mdct, int32_t f_offset)
{
  cf_t ssb_grid[SRSRAN_SSB_NOF_RE];
  cf_t* pss_in_ssb = &ssb_grid[SRSRAN_PSS_NR_SYMBOL_IDX * SRSRAN_SSB_BW_SUBC];
  int N_id_2;
  srsran_dft_plan_t ifft_plan;
  srsran_dft_plan_c(&ifft_plan, mdct->symbol_sz, SRSRAN_DFT_BACKWARD);
  // TODO: Could be optimized by multiplying the constant phase instead of computing the IFFT for each N_id_2
  for (N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    mdct->pss_x[N_id_2] = (cf_t*)malloc(mdct->symbol_sz * sizeof(cf_t));
    srsran_vec_cf_zero(mdct->temp, mdct->symbol_sz);
    srsran_pss_nr_put(ssb_grid, N_id_2, 1.0f);
    srsran_vec_cf_copy(&mdct->temp[0],
                       &pss_in_ssb[SRSRAN_SSB_BW_SUBC / 2 - f_offset],
                       SRSRAN_SSB_BW_SUBC / 2 + f_offset);
    srsran_vec_cf_copy(&mdct->temp[mdct->symbol_sz - SRSRAN_SSB_BW_SUBC / 2 + f_offset],
                       &pss_in_ssb[0],
                       SRSRAN_SSB_BW_SUBC / 2 - f_offset);
    srsran_dft_run_c(&ifft_plan, mdct->temp, mdct->pss_x[N_id_2]);
  }
  srsran_dft_plan_free(&ifft_plan);
}

int srsran_prepare_pss_mdct(srsran_pss_mdct_t* mdct,
                            uint32_t symbol_sz, int32_t f_offset,
                            uint32_t Q, uint32_t PSI)
{
  int i, j;
  mdct->symbol_sz = symbol_sz;
  mdct->f_offset = f_offset;
  mdct->Q = Q;
  mdct->PSI = PSI;
  mdct->output = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
  mdct->debug = false;
  if (mdct->output == NULL) {
    return SRSRAN_ERROR;
  }
  mdct->temp = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
  if (mdct->temp == NULL) {
    free(mdct->output);
    return SRSRAN_ERROR;
  }
  prepare_pss_x(mdct, f_offset);
  for (i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    mdct->x_tilde[i] = (cf_t**)malloc(PSI * sizeof(cf_t*));
    if (mdct->x_tilde[i] == NULL) {
      for (j = 0; j < i; j++) {
        free(mdct->x_tilde[j]);
      }
      free(mdct->temp);
      free(mdct->output);
      return SRSRAN_ERROR;
    }
    for(j = 0; j < PSI; j++) {
      mdct->x_tilde[i][j] = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
      // TODO error handling
      fill_x_tilde(mdct, i, j);
    }
  }
  return SRSRAN_SUCCESS;
}

int srsran_destroy_pss_mdct(srsran_pss_mdct_t* mdct)
{
  if (mdct == NULL) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  for (int i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    for (int j = 0; j < mdct->PSI; j++) {
      free(mdct->x_tilde[i][j]);
    }
    free(mdct->x_tilde[i]);
    free(mdct->pss_x[i]);
  }
  free(mdct->temp);
  mdct->temp = NULL;
  free(mdct->output);
  mdct->output = NULL;
  return SRSRAN_SUCCESS;
}

SRSRAN_API int correlation_detect_pss(const srsran_pss_mdct_t* mdct,
                                      const cf_t* in, uint32_t nof_samples,
                                      uint32_t window_sz,
                                      srsran_pss_detect_res_t* result)
{
  float peak = -1 * INFINITY;
  for (uint32_t N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    for (int32_t tau = 0; tau < nof_samples - mdct->symbol_sz; tau += (int)window_sz) {
      srsran_vec_prod_conj_ccc(&in[tau], mdct->pss_x[N_id_2], mdct->temp, mdct->symbol_sz);
      float corr_mag = srsran_vec_acc_cc(mdct->temp, mdct->symbol_sz);
      if (corr_mag > peak) {
        peak = corr_mag;
        result->tau = tau;
        result->N_id_2 = N_id_2;
        result->peak_value = peak;
      }
    }
  }
  return SRSRAN_SUCCESS;
}

SRSRAN_API int mdct_detect_pss(const srsran_pss_mdct_t* mdct,
                               const cf_t* in, uint32_t nof_samples,
                               uint32_t window_sz,
                               srsran_pss_detect_res_t* result)
{
  float peak = -1 * INFINITY;
  // TODO: Should be optimized by multiplying the constant phase instead of computing MDCT for each N_id_2
  for (uint32_t N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    for (int32_t tau = 0; tau < nof_samples - mdct->symbol_sz; tau += (int)window_sz) {
      cf_t corr = calculate_C(mdct, in + tau, N_id_2, tau);
      float corr_mag = cabsf(corr);
      if(mdct->debug) {
        printf("N_id_2=%d, tau=%d, corr_mag=%f\n", N_id_2, tau, corr_mag);
      }
      if (corr_mag > peak) {
        peak = corr_mag;
        result->tau = tau;
        result->N_id_2 = N_id_2;
        result->peak_value = peak;
      }
    }
  }
  return SRSRAN_SUCCESS;
}
