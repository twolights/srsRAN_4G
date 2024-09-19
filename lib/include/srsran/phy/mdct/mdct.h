//
// Created by ykchen on 9/8/24.
//

#ifndef SRSRAN_MDCT_H
#define SRSRAN_MDCT_H

#include <complex.h>
#include "srsran/phy/sync/pss_nr.h"

#define SRSRAN_MDCT_RECOMMENDED_Q 1
#define SRSRAN_MDCT_RECOMMENDED_PSI 6

typedef struct SRSRAN_API {
  cf_t* pss_x[SRSRAN_NOF_NID_2_NR];
  cf_t** x_tilde[SRSRAN_NOF_NID_2_NR];
  cf_t* temp;
  cf_t* output;
  uint32_t symbol_sz;
  uint32_t Q;
  uint32_t PSI;
  bool debug;
} srsran_pss_mdct_t;

typedef struct SRSRAN_API {
  int32_t   tau;        // time offset
  uint32_t  N_id_2;
  float     peak_value;
} srsran_pss_detect_res_t;

SRSRAN_API int srsran_prepare_pss_mdct(srsran_pss_mdct_t* mdct, uint32_t symbol_sz, uint32_t Q, uint32_t PSI);
SRSRAN_API int srsran_destroy_pss_mdct(srsran_pss_mdct_t* mdct);
SRSRAN_API int correlation_detect_pss(const srsran_pss_mdct_t* mdct,
                                      const cf_t* in, uint32_t nof_samples,
                                      uint32_t window_sz,
                                      srsran_pss_detect_res_t* result);
// TODO rename
SRSRAN_API int mdct_detect_pss(const srsran_pss_mdct_t* mdct,
                               const cf_t* in, uint32_t nof_samples,
                               uint32_t window_sz,
                               srsran_pss_detect_res_t* result);

#endif // SRSRAN_MDCT_H
