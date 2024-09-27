//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/mdct.h"
#include "srsran/phy/dft/dft.h"
#include "srsran/phy/mdct/differential_prod.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


// TODO move this to a common place
void unwrap_phase(const float* phase, float* target, size_t length) {
    if (length < 2) {
        return;
    }

    const float TWO_PI = 2.0f * M_PI;
    float offset = 0.0f;

    target[0] = phase[0];
    for (size_t i = 1; i < length; i++) {
        float delta = phase[i] - phase[i - 1];

        if (delta > M_PI) {
            offset -= TWO_PI;
        } else if (delta < -M_PI) {
            offset += TWO_PI;
        }

        target[i] = phase[i] + offset;
    }
}

// TODO make this return phase
static inline cf_t
calculate_D(const srsran_pss_mdct_t* mdct, uint32_t N_id_2, uint32_t psi)
{
  srsran_vec_prod_conj_ccc(mdct->y_tilde[psi], mdct->x_tilde[N_id_2][psi], mdct->temp, mdct->symbol_sz);
  cf_t result = srsran_vec_acc_cc(mdct->temp, mdct->symbol_sz);
//  printf("d=%d (Q=%u, psi=%u), atan2f(%f, %f) = %f\n",
//         d, Q, psi,
//         crealf(result), cimagf(result),
//         atan2f(cimagf(result), crealf(result))
//  );
  return result;
}

static inline cf_t
calculate_C(const srsran_pss_mdct_t* mdct, const cf_t* y, uint32_t N_id_2)
{
  cf_t result = 0.0;
  for (int psi = 0; psi < mdct->PSI; psi++) {
    result += calculate_D(mdct, N_id_2, psi);
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
  mdct->debug = false;
  mdct->temp = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
  if (mdct->temp == NULL) {
    return SRSRAN_ERROR;
  }
  mdct->y_tilde = (cf_t**)malloc(mdct->PSI * sizeof(cf_t*));
  if (mdct->y_tilde == NULL) {
    free(mdct->temp);
    return SRSRAN_ERROR;
  }
  for (i = 0; i < mdct->PSI; i++) {
      mdct->y_tilde[i] = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
      if (mdct->y_tilde[i] == NULL) {
        for (j = 0; j < i; j++) {
          free(mdct->y_tilde[j]);
        }
        free(mdct->y_tilde);
        free(mdct->temp);
        return SRSRAN_ERROR;
      }
  }

  prepare_pss_x(mdct, f_offset);
  for (i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    mdct->x_tilde[i] = (cf_t**)malloc(PSI * sizeof(cf_t*));
    if (mdct->x_tilde[i] == NULL) {
      for (j = 0; j < i; j++) {
        free(mdct->x_tilde[j]);
      }
      free(mdct->temp);
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
  int i, j;
  if (mdct == NULL) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  for (i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    for (j = 0; j < mdct->PSI; j++) {
      free(mdct->x_tilde[i][j]);
    }
    free(mdct->x_tilde[i]);
    free(mdct->pss_x[i]);
  }
  free(mdct->temp);
  mdct->temp = NULL;
  for (i = 0; i < mdct->PSI; i++) {
    free(mdct->y_tilde[i]);
    mdct->y_tilde[i] = NULL;
  }
  free(mdct->y_tilde);
  return SRSRAN_SUCCESS;
}

SRSRAN_API int srsran_detect_pss_correlation(const srsran_pss_mdct_t* mdct,
                                             const cf_t* in, uint32_t nof_samples,
                                             uint32_t window_sz,
                                             srsran_pss_detect_res_t* result)
{
  float peak = -1 * INFINITY;
  for (uint32_t N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    for (int32_t tau = 0; tau < nof_samples - mdct->symbol_sz; tau += (int)window_sz) {
      srsran_vec_prod_conj_ccc(&in[tau], mdct->pss_x[N_id_2], mdct->temp, mdct->symbol_sz);
      float corr_mag = cabsf(srsran_vec_acc_cc(mdct->temp, mdct->symbol_sz));
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

static void prepare_y_tilde(const srsran_pss_mdct_t* mdct, const cf_t* in, uint32_t tau)
{
  int32_t d;
  for (int psi = 0; psi < mdct->PSI; psi++) {
    d = (int32_t)(1 + (mdct->Q * psi));
    differential_product(&in[tau], mdct->y_tilde[psi], d, mdct->symbol_sz);
  }
}

static int mdct_detect_pss_with_nid2_set(const srsran_pss_mdct_t* mdct,
                                                uint32_t min_N_id_2, uint32_t max_N_id_2,
                                                const cf_t* in, uint32_t nof_samples,
                                                uint32_t window_sz,
                                                srsran_pss_detect_res_t* result)
{
  float peak = -1 * INFINITY;
  cf_t** y_tilde = (cf_t**)malloc(mdct->PSI * sizeof(cf_t*));
  memset(y_tilde, 0, mdct->PSI * sizeof(cf_t*));

  // TODO: Should be optimized by multiplying the constant phase instead of computing MDCT for each N_id_2
  for (int32_t tau = 0; tau < nof_samples - mdct->symbol_sz; tau += (int)window_sz) {
    prepare_y_tilde(mdct, in, tau);
    for (uint32_t N_id_2 = min_N_id_2; N_id_2 <= max_N_id_2; N_id_2++) {
      cf_t corr = calculate_C(mdct, in + tau, N_id_2);
      float corr_mag = cabsf(corr);
      if (mdct->debug) {
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
  for (int i = 0; i < mdct->PSI; i++) {
    if(y_tilde[i] == NULL) {
      continue;
    }
    free(y_tilde[i]);
  }
  free(y_tilde);
  return SRSRAN_SUCCESS;
}

SRSRAN_API int srsran_detect_pss_mdct(const srsran_pss_mdct_t* mdct,
                                      const cf_t* in, uint32_t nof_samples,
                                      uint32_t window_sz,
                                      srsran_pss_detect_res_t* result)
{
  return mdct_detect_pss_with_nid2_set(mdct,
                                       0, SRSRAN_NOF_NID_2_NR - 1,
                                       in, nof_samples,
                                       window_sz,
                                       result);
}

SRSRAN_API int srsran_find_pss_mdct(const srsran_pss_mdct_t* mdct,
                                    uint32_t N_id_2,
                                    const cf_t* in, uint32_t nof_samples,
                                    uint32_t window_sz,
                                    srsran_pss_detect_res_t* result)
{
  return mdct_detect_pss_with_nid2_set(mdct,
                                       N_id_2, N_id_2,
                                       in, nof_samples,
                                       window_sz,
                                       result);
}
