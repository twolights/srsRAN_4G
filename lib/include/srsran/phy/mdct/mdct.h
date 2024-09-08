//
// Created by ykchen on 9/8/24.
//

#ifndef SRSRAN_MDCT_H
#define SRSRAN_MDCT_H

#include <complex.h>
#include "srsran/phy/utils/vector.h"

float mdct(const cf_t* y, const cf_t* x, uint32_t n, uint32_t Q, uint32_t PSI, int32_t tau);

#endif // SRSRAN_MDCT_H
