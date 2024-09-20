//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/differential_prod.h"
#include "srsran/phy/utils/vector.h"
#include <stdlib.h>

// TODO move this to vector.h
void ssb_vec_prod_conj_circ_shift(const cf_t* a, const cf_t* b, cf_t* c, uint32_t n, int shift);

// TODO make this inline?
void differential_product(const cf_t* a, cf_t* c, int32_t d, uint32_t n)
{
  ssb_vec_prod_conj_circ_shift(a, a, c, n, d);
}