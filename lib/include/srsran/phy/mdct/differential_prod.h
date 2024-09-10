//
// Created by ykchen on 9/8/24.
//

#ifndef SRSRAN_DIFFERENTIAL_PROD_H
#define SRSRAN_DIFFERENTIAL_PROD_H

#include <complex.h>
#include "srsran/phy/utils/vector.h"

void srsran_vec_prod_conj_cyclic_ccc(const cf_t* a, const cf_t* b, cf_t* c, int32_t shifts, uint32_t n);
void differential_product(const cf_t* a, cf_t* c, int32_t d, uint32_t n);

#endif // SRSRAN_DIFFERENTIAL_PROD_H
