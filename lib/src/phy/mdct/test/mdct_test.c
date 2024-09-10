//
// Created by ykchen on 9/12/24.
//
#include "srsran/phy/mdct/mdct.h"
#include <string.h>

static cf_t time_domain_pss[SRSRAN_NOF_NID_2_NR][SRSRAN_PSS_NR_LEN];

#define NOF_SAMPLES 150

void get_time_domain_pss(cf_t* out, uint32_t N_id_2);

int main() {
  int i;
  for (i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    get_time_domain_pss(time_domain_pss[i], i);
  }
  cf_t buffer[NOF_SAMPLES];
  srsran_pss_mdct_t mdct;
  srsran_pss_mdct_detect_res_t res;
  memset(buffer, 0, NOF_SAMPLES * sizeof(cf_t));
  memcpy(buffer + 11, time_domain_pss[1], SRSRAN_PSS_NR_LEN * sizeof(cf_t));
  srsran_prepare_pss_mdct(&mdct, SRSRAN_MDCT_PSS_FFT_SIZE, 1, 6);
  mdct_detect_pss(&mdct, buffer, NOF_SAMPLES, &res);
  printf("MDCT: detected tau=%d, N_id_2=%d, peak_value=%f\n", res.tau, res.N_id_2, res.peak_value);
  srsran_destroy_pss_mdct(&mdct);
  return 0;
}
