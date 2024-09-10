//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/differential_prod.h"
#include "srsran/phy/sync/pss_nr.h"
#include "srsran/phy/dft/dft.h"
#include <complex.h>
#include <stdio.h>
#include <string.h>

#define N 10
#define FFT_SIZE (SRSRAN_PSS_NR_LEN + 1)
#define PSS_NR_SUBC_BEGIN 56

bool test_on_sequence()
{
  cf_t a[N] = {
      1.0 + 1.0*I,
      2.0 + 2.0*I,
      3.0 + 3.0*I,
      4.0 + 4.0*I,
      5.0 + 5.0*I,
      6.0 + 6.0*I,
      7.0 + 7.0*I,
      8.0 + 8.0*I,
      9.0 + 9.0*I,
      10.0 + 10.0*I
    };
  cf_t c[N];
  int32_t d = 1;
  differential_product(a, c, d, N);
  return true;
}

void prepare_all_pss_nr(cf_t all_pss[SRSRAN_NOF_NID_2_NR][FFT_SIZE])
{
  uint32_t N_id_2;
  cf_t ssb_grid[SRSRAN_SSB_NOF_RE];
  uint32_t pss_loc_in_ssb = SRSRAN_PSS_NR_SYMBOL_IDX * SRSRAN_SSB_BW_SUBC + PSS_NR_SUBC_BEGIN;

  for (N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    srsran_pss_nr_put(ssb_grid, N_id_2, 1.0f);
    memcpy(all_pss[N_id_2], &ssb_grid[pss_loc_in_ssb], SRSRAN_PSS_NR_LEN * sizeof(cf_t));
    all_pss[N_id_2][SRSRAN_PSS_NR_LEN] = 0.0f;
  }
}

bool test_on_pss_nr()
{
  cf_t all_pss[SRSRAN_NOF_NID_2_NR][FFT_SIZE];
  cf_t all_pss_time_domain[SRSRAN_NOF_NID_2_NR][FFT_SIZE];
  srsran_dft_plan_t ifft_plan;

  prepare_all_pss_nr(all_pss);

  for (uint32_t N_id_2 = 0; N_id_2 < SRSRAN_NOF_NID_2_NR; N_id_2++) {
    srsran_dft_plan_c(&ifft_plan, FFT_SIZE, SRSRAN_DFT_BACKWARD);
    srsran_dft_run_c(&ifft_plan, all_pss[N_id_2], all_pss_time_domain[N_id_2]);
    srsran_dft_plan_free(&ifft_plan);
  }

//  for (int i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
//    differential_product(all_pss[i], all_pss_diff_prod[i], 1, SRSRAN_PSS_NR_LEN);
//  }

  return true;
}

int main() {
  test_on_sequence();
  test_on_pss_nr();
  return 0;
}