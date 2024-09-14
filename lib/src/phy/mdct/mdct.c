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

static cf_t time_domain_pss[SRSRAN_NOF_NID_2_NR][SRSRAN_MDCT_PSS_FFT_SIZE];

// Prepare all time domain PSS sequences for all NID2's
static void prepare_all_nr_pss_x() {
  static bool initialized = false;
  if (initialized) {
    return;
  }

  uint32_t N_id_2;
  cf_t ssb_grid[SRSRAN_SSB_NOF_RE];
  uint32_t pss_loc_in_ssb = SRSRAN_PSS_NR_SYMBOL_IDX * SRSRAN_SSB_BW_SUBC + PSS_NR_SUBC_BEGIN;
  srsran_dft_plan_t ifft_plan;
  cf_t temp[SRSRAN_MDCT_PSS_FFT_SIZE];
  int factor = floor(SRSRAN_MDCT_PSS_FFT_SIZE / SRSRAN_PSS_NR_LEN);

  // TODO: Could be optimized by multiplying the constant phase instead of computing the IFFT for each N_id_2
  for (N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    srsran_pss_nr_put(ssb_grid, N_id_2, 1.0f);
    memset(temp, 0, SRSRAN_MDCT_PSS_FFT_SIZE * sizeof(cf_t));
    for (int i = 0; i < SRSRAN_PSS_NR_LEN; i++) {
      temp[i * factor] = ssb_grid[pss_loc_in_ssb + i];
    }
//    memcpy(temp, &ssb_grid[pss_loc_in_ssb], SRSRAN_PSS_NR_LEN * sizeof(cf_t));
    srsran_dft_plan_c(&ifft_plan, SRSRAN_MDCT_PSS_FFT_SIZE, SRSRAN_DFT_BACKWARD);
    srsran_dft_run_c(&ifft_plan, temp, time_domain_pss[N_id_2]);
    srsran_dft_plan_free(&ifft_plan);
  }
  initialized = true;
}

// XXX for testing
void get_time_domain_pss(cf_t* out, uint32_t N_id_2)
{
  prepare_all_nr_pss_x();
  memcpy(out, time_domain_pss[N_id_2], SRSRAN_PSS_NR_LEN * sizeof(cf_t));
}

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
    result += calculate_D(y, mdct->x_tilde[N_id_2][psi], mdct->temp, mdct->output, mdct->n, mdct->Q, psi);
  }
  return result / (float)mdct->PSI;
}

// x_tilde_{n,d}^{r}: d: delay, r: NID2
static void get_x_tilde(cf_t* out, uint32_t n, int32_t d, uint32_t r)
{
  prepare_all_nr_pss_x();
  // TODO: Should be optimized by multiplying the constant phase instead of computing the diff. product for each N_id_2
  differential_product(time_domain_pss[r], out, d, n);
}

int srsran_prepare_pss_mdct(srsran_pss_mdct_t* mdct, uint32_t n, uint32_t Q, uint32_t PSI)
{
  int i, j;
  int32_t d;
  mdct->output = (cf_t*)malloc(n * sizeof(cf_t));
  if (mdct->output == NULL) {
    return SRSRAN_ERROR;
  }
  mdct->temp = (cf_t*)malloc(n * sizeof(cf_t));
  if (mdct->temp == NULL) {
    free(mdct->output);
    return SRSRAN_ERROR;
  }
  for (i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    mdct->x_tilde[i] = (cf_t**)malloc(PSI * sizeof(cf_t*));
    if (mdct->x_tilde[i] == NULL) {
      for (j = 0; j < i; j++) {
        free(mdct->x_tilde[j]);
        free(mdct->temp);
        free(mdct->output);
      }
      return SRSRAN_ERROR;
    }
    for(j = 0; j < PSI; j++) {
      mdct->x_tilde[i][j] = (cf_t*)malloc(n * sizeof(cf_t));
      // TODO error handling
      d = (int32_t)(1 + Q * j);
      get_x_tilde(mdct->x_tilde[i][j], n, d, i);
    }
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
  for (int i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    for (int j = 0; j < mdct->PSI; j++) {
      free(mdct->x_tilde[i][j]);
    }
    free(mdct->x_tilde[i]);
  }
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
  float peak = -1 * INFINITY;
  // TODO: Should be optimized by multiplying the constant phase instead of computing MDCT for each N_id_2
  for (uint32_t N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    for (int32_t tau = 0; tau < nof_samples - mdct->n; tau++) {
      cf_t corr = calculate_C(mdct, in + tau, N_id_2, tau);
      float corr_mag = cabsf(corr);
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
