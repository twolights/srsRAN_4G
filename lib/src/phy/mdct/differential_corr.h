//
// Created by ykchen on 9/8/24.
//

#ifndef SRSRAN_DIFFERENTIAL_CORR_H
#define SRSRAN_DIFFERENTIAL_CORR_H

#include <complex.h>
#include "srsran/phy/utils/vector.h"

static void differential_product(const cf_t* a, const cf_t* b, cf_t* c, uint32_t d, uint32_t n);

#endif // SRSRAN_DIFFERENTIAL_CORR_H
