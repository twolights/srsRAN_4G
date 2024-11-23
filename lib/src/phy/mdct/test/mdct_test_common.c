//
// Created by ykchen on 11/23/24.
//

#include "mdct_test_common.h"
#include <string.h>

void append_pss(srsran_pss_mdct_t* mdct, cf_t* buffer, uint32_t N_id_2, int32_t tau, int beta)
{
  for (int i = 0; i < mdct->symbol_sz; i++) {
    buffer[tau + i] += mdct->pss_x[N_id_2][i] * beta / 100.0;
  }
}

void prepare_mocked_received_samples(srsran_pss_mdct_t* mdct, cf_t* buffer, size_t size, uint32_t N_id_2, int32_t tau, bool add_noise)
{
  if (add_noise) {
    memset(buffer, 0, size * sizeof(cf_t));
    // TODO fill buffer with noise
  } else {
    memset(buffer, 0, size * sizeof(cf_t));
  }
  append_pss(mdct, buffer, N_id_2, tau, 100);
}

