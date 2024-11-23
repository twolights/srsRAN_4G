//
// Created by ykchen on 11/23/24.
//

#ifndef SRSRAN_MDCT_TEST_COMMON_H
#define SRSRAN_MDCT_TEST_COMMON_H

#include "srsran/phy/mdct/mdct.h"

#define DETECTION_METHOD_CORRELATION 0
#define DETECTION_METHOD_MDCT 1

#define SAMPLING_FREQUENCY 23.04e6
#define SYMBOL_SIZE 1536
#define NOF_SAMPLES (SYMBOL_SIZE * 2)

void append_pss(srsran_pss_mdct_t* mdct, cf_t* buffer, uint32_t N_id_2, int32_t tau, int beta);
void prepare_mocked_received_samples(srsran_pss_mdct_t* mdct, cf_t* buffer, size_t size, uint32_t N_id_2, int32_t tau, bool add_noise);

#endif // SRSRAN_MDCT_TEST_COMMON_H
