//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/mdct.h"
#include "srsran/phy/dft/dft.h"
#include "srsran/phy/utils/vector.h"
//#include "srsran/phy/mdct/differential_prod.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


inline static void srsran_vec_cf_copy_reversed(const cf_t* a, cf_t* b, uint32_t n)
{
  for (uint32_t i = 0; i < n; i++) {
    b[i] = a[n - i - 1];
  }
}

inline static void differential_product(const cf_t* a, cf_t* c, uint32_t d, uint32_t n)
{
  //  ssb_vec_prod_conj_circ_shift(a, a, c, n, d);
  uint32_t half = n / 2, half_d = half - d;
  if (d == 0) {
    srsran_vec_prod_conj_ccc(a, a, c, n);
    return;
  }
  srsran_vec_prod_conj_ccc(&a[d], &a[0], &c[0], half_d);
  srsran_vec_prod_conj_ccc(&a[0], &a[half_d], &c[half_d], d);
  srsran_vec_cf_copy_reversed(&c[0], &c[n - half], half);
  //  // If n is odd, the middle element is multiplied by its conjugate
  //  c[half + 1] = a[half + 1 + d]  * conj(a[half + 1]);
  //  srsran_vec_prod_conj_ccc(&a[d], &a[0], &c[0], n - d);
  //  srsran_vec_prod_conj_ccc(&a[0], &a[n - d], &c[n - d], d);
}

// TODO move this to a common place
void unwrap_phase(const float* phase, float* target, size_t length) {
  if (phase == NULL || target == NULL || length < 2) {
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

static inline uint32_t get_d(const srsran_pss_mdct_t* mdct, uint32_t psi)
{
  return (1 + (mdct->Q * psi));
}

static inline cf_t
calculate_D(const srsran_pss_mdct_t* mdct, uint32_t N_id_2, uint32_t psi)
{
  srsran_vec_prod_conj_ccc(mdct->y_tilde[psi], mdct->x_tilde[N_id_2][psi], mdct->temp, mdct->symbol_sz);
  cf_t result = srsran_vec_acc_cc(mdct->temp, mdct->symbol_sz);
  return result;
}

static inline cf_t
calculate_C(const srsran_pss_mdct_t* mdct, uint32_t N_id_2)
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
  int32_t d = get_d(mdct, psi);
  differential_product(mdct->pss_x[r], mdct->x_tilde[r][psi], (int)d, mdct->symbol_sz);
}

static inline void preserve_best_y_tilde(srsran_pss_mdct_t* mdct)
{
//  for (int psi = 0; psi < mdct->PSI; psi++) {
//    memcpy(mdct->y_tilde_best[psi], mdct->y_tilde[psi], mdct->symbol_sz * sizeof(cf_t));
//  }
  mdct->y_tilde_best = mdct->y_tilde;
}

static inline void switch_y_tilde(srsran_pss_mdct_t* mdct)
{
  mdct->y_tilde_current_index = !mdct->y_tilde_current_index;
  mdct->y_tilde = mdct->y_tilde_buffers[mdct->y_tilde_current_index];
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

SRSRAN_API int srsran_prepare_pss_mdct(srsran_pss_mdct_t* mdct,
                                       uint32_t srate_hz, uint32_t symbol_sz, int32_t f_offset,
                                       uint32_t Q, uint32_t PSI)
{
  int i, j;
  mdct->srate_hz = srate_hz;
  mdct->symbol_sz = symbol_sz;
  mdct->f_offset = f_offset;
  mdct->Q = Q;
  mdct->PSI = PSI;
  mdct->debug = false;
  mdct->temp = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
  if (mdct->temp == NULL) {
    return SRSRAN_ERROR;
  }
  mdct->phase = (float*)malloc(symbol_sz * sizeof(float));
  for (i = 0; i < 2; i++) {
    mdct->y_tilde_buffers[i] = (cf_t**)malloc(PSI * sizeof(cf_t*));
    if (mdct->y_tilde_buffers[0] == NULL) {
      free(mdct->y_tilde_buffers[0]);
      free(mdct->temp);
      free(mdct->phase);
      return SRSRAN_ERROR;
    }
  }
  mdct->y_tilde_current_index = 0;
  mdct->y_tilde = mdct->y_tilde_buffers[0];
  mdct->y_tilde_best = mdct->y_tilde_buffers[1];
  for (i = 0; i < mdct->PSI; i++) {
      mdct->y_tilde[i] = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
      mdct->y_tilde_best[i] = (cf_t*)malloc(symbol_sz * sizeof(cf_t));
      if (mdct->y_tilde[i] == NULL || mdct->y_tilde_best[i] == NULL) {
        for (j = 0; j < i; j++) {
          free(mdct->y_tilde[j]);
          free(mdct->y_tilde_best[j]);
        }
        free(mdct->y_tilde);
        free(mdct->y_tilde_best);
        free(mdct->temp);
        free(mdct->phase);
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
      free(mdct->phase);
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

SRSRAN_API int srsran_destroy_pss_mdct(srsran_pss_mdct_t* mdct)
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
    free(mdct->y_tilde_best[i]);
    mdct->y_tilde[i] = NULL;
    mdct->y_tilde_best[i] = NULL;
  }
  free(mdct->y_tilde);
  free(mdct->y_tilde_best);
  return SRSRAN_SUCCESS;
}

int estimate_coarse_cfo_with_mdct(const srsran_pss_mdct_t* mdct, srsran_pss_detect_res_t* res)
{
  if (mdct == NULL || res == NULL || res->peak_value < 0 || res->tau < 0) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  int i, psi;
  uint32_t nof_samples_to_process;
  uint32_t d;
  float theta_D, f_D;
  float* unwrapped = &mdct->phase[mdct->PSI];
  for (psi = 0; psi < mdct->PSI; psi++) {
    srsran_vec_div_ccc(mdct->y_tilde_best[psi], mdct->x_tilde[res->N_id_2][psi], mdct->temp, mdct->symbol_sz);
    d = get_d(mdct, psi);
    theta_D = 0;
    nof_samples_to_process = mdct->symbol_sz - d;
    // Get average f_D over all but the tail d the samples
    for (i = 0; i < nof_samples_to_process; i++) {
      theta_D += cargf(mdct->temp[i]);
    }
    mdct->phase[psi] = theta_D / (float)nof_samples_to_process;
  }
  unwrap_phase(mdct->phase, unwrapped, mdct->PSI);
  for (psi = 0; psi < mdct->PSI; psi++) {
    d = get_d(mdct, psi);
    unwrapped[psi] /= (float)(2 * M_PI * d);
  }
  // Average over all the symbols & psi's
  f_D = srsran_vec_acc_ff(unwrapped, mdct->PSI) / (float)mdct->PSI;
  f_D *= (float)mdct->srate_hz;
  res->coarse_cfo = f_D;
  return SRSRAN_SUCCESS;
}

int estimate_coarse_cfo(const srsran_pss_mdct_t* mdct,
                        const cf_t* in, uint32_t nof_samples,
                        srsran_pss_detect_res_t* res)
{
  if (mdct == NULL || in == NULL || res == NULL || res->peak_value < 0 || res->tau < 0) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  if (nof_samples < mdct->symbol_sz) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  srsran_vec_prod_conj_ccc(&in[res->tau], mdct->pss_x[res->N_id_2], mdct->temp, mdct->symbol_sz);
  float* target = (float*)malloc(mdct->symbol_sz * sizeof(float));
  for (int i = 0; i < mdct->symbol_sz; i++) {
    mdct->phase[i] = cargf(mdct->temp[i]);
  }
  unwrap_phase(mdct->phase, target, mdct->symbol_sz);
  float phase_diff = mdct->phase[mdct->symbol_sz - 1] - mdct->phase[0];
//  res->coarse_cfo = 75000;
  res->coarse_cfo = (float)(phase_diff / (2 * M_PI * mdct->symbol_sz / mdct->srate_hz)) * (float)mdct->symbol_sz;
  free(target);
  return SRSRAN_SUCCESS;
}

int estimate_cfo_by_half_pss(const srsran_pss_mdct_t* mdct,
                             const cf_t* in, uint32_t nof_samples,
                             srsran_pss_detect_res_t* res)
{
  if (mdct == NULL || in == NULL || res == NULL || res->peak_value < 0 || res->tau < 0) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  if (nof_samples < mdct->symbol_sz) {
    return SRSRAN_ERROR_INVALID_INPUTS;
  }
  cf_t c1 = srsran_vec_dot_prod_ccc(&in[res->tau], mdct->pss_x[res->N_id_2], mdct->symbol_sz / 2);
  cf_t c2 = srsran_vec_dot_prod_ccc(&in[res->tau + mdct->symbol_sz / 2],
                                    &mdct->pss_x[res->N_id_2][mdct->symbol_sz / 2],
                                    mdct->symbol_sz / 2);
  res->coarse_cfo = cargf(conjf(c1) * c2) / M_PI;
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
  estimate_coarse_cfo(mdct, in, nof_samples, result);
  return SRSRAN_SUCCESS;
}

static inline void prepare_y_tilde(const srsran_pss_mdct_t* mdct, const cf_t* in, uint32_t tau)
{
  uint32_t d;
  for (int psi = 0; psi < mdct->PSI; psi++) {
    d = get_d(mdct, psi);
    differential_product(&in[tau], mdct->y_tilde[psi], d, mdct->symbol_sz);
  }
}

static inline int mdct_detect_pss_with_nid2_set(srsran_pss_mdct_t* mdct,
                                                uint32_t min_N_id_2, uint32_t max_N_id_2,
                                                const cf_t* in, uint32_t nof_samples,
                                                uint32_t window_sz,
                                                bool estimate_cfo,
                                                srsran_pss_detect_res_t* result)
{
  float peak = -1 * INFINITY;
  bool y_tilde_switched;

  // TODO: Should be optimized by multiplying the constant phase instead of computing MDCT for each N_id_2
  for (int32_t tau = 0; tau < nof_samples - mdct->symbol_sz; tau += (int)window_sz) {
    y_tilde_switched = false;
    prepare_y_tilde(mdct, in, tau);
    for (uint32_t N_id_2 = min_N_id_2; N_id_2 <= max_N_id_2; N_id_2++) {
      cf_t corr = calculate_C(mdct, N_id_2);
      float corr_mag = cabsf(corr);
#if 0
      if (mdct->debug) {
        printf("N_id_2=%d, tau=%d, corr_mag=%f\n", N_id_2, tau, corr_mag);
      }
#endif
      if (corr_mag > peak) {
        peak = corr_mag;
        result->tau = tau;
        result->N_id_2 = N_id_2;
        result->peak_value = peak;
        preserve_best_y_tilde(mdct);
        y_tilde_switched = true;
      }
    }
    if(y_tilde_switched) {
      switch_y_tilde(mdct);
    }
  }
  if (estimate_cfo) {
    result->coarse_cfo = 0;
//    estimate_coarse_cfo_with_mdct(mdct, result);
  }
  return SRSRAN_SUCCESS;
}

SRSRAN_API int srsran_detect_pss_mdct(srsran_pss_mdct_t* mdct,
                                      const cf_t* in, uint32_t nof_samples,
                                      uint32_t window_sz,
                                      bool estimate_cfo,
                                      srsran_pss_detect_res_t* result)
{
  return mdct_detect_pss_with_nid2_set(mdct,
                                       0, SRSRAN_NOF_NID_2_NR - 1,
                                       in, nof_samples,
                                       window_sz,
                                       estimate_cfo,
                                       result);
}

SRSRAN_API int srsran_find_pss_mdct(srsran_pss_mdct_t* mdct,
                                    uint32_t N_id_2,
                                    const cf_t* in, uint32_t nof_samples,
                                    uint32_t window_sz,
                                    bool estimate_cfo,
                                    srsran_pss_detect_res_t* result)
{
  return mdct_detect_pss_with_nid2_set(mdct,
                                       N_id_2, N_id_2,
                                       in, nof_samples,
                                       window_sz,
                                       estimate_cfo,
                                       result);
}
