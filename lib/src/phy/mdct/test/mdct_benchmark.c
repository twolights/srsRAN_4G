//
// Created by ykchen on 11/23/24.
//

#include "mdct_test_common.h"
#include "srsran/phy/mdct/mdct.h"
#include "srsran/phy/utils/vector.h"
#include <time.h>

#define NUM_MDCT_TESTS 10000000

int main() {
  srsran_pss_mdct_t mdct;
  srsran_prepare_pss_mdct(&mdct,
                          SAMPLING_FREQUENCY,
                          SYMBOL_SIZE, -30,
                          SRSRAN_MDCT_RECOMMENDED_Q * 12,
                          SRSRAN_MDCT_RECOMMENDED_PSI);
  cf_t buffer[SYMBOL_SIZE];
  srsran_pss_detect_res_t res;
  prepare_mocked_received_samples(&mdct, buffer, SYMBOL_SIZE, 0, 0, false);

  clock_t start = clock(), end;

  for (int i = 0; i < NUM_MDCT_TESTS; i++) {
    srsran_detect_pss_mdct(&mdct, buffer, SYMBOL_SIZE, 1, false, &res);
  }

  end = clock();
  printf("Time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

  srsran_destroy_pss_mdct(&mdct);
  return 0;
}
