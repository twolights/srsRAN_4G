//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/differential_prod.h"
#include "srsran/phy/utils/vector.h"
#include <stdlib.h>

void srsran_vec_prod_conj_cyclic_ccc(const cf_t* a, const cf_t* b, cf_t* c, int32_t shifts, uint32_t n)
{
  uint32_t offset = (uint32_t)abs(shifts);

  // Avoid negative of samples
  if (offset > n) {
    srsran_vec_cf_zero(c, n);
    return;
  }

  if (shifts > 0) {
    srsran_vec_prod_conj_ccc(&a[0], &b[offset], &c[0], n - offset);
    srsran_vec_prod_conj_ccc(&a[n - offset], &b[0], &c[n - offset], offset);
    return;
  }

  if (shifts < 0) {
    srsran_vec_prod_conj_ccc(&a[offset], &b[0], &c[0], n - offset);
    srsran_vec_prod_conj_ccc(&a[0], &b[n - offset], &c[n - offset], offset);
    return;
  }

  srsran_vec_prod_conj_ccc(a, b, c, n);
}

void differential_product(const cf_t* a, cf_t* c, int32_t d, uint32_t n)
{
  srsran_vec_prod_conj_cyclic_ccc(a, a, c, d, n);
}