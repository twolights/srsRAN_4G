//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/differential_prod.h"
#include "srsran/phy/utils/vector.h"
#include <stdlib.h>

void differential_product(const cf_t* a, cf_t* c, int32_t d, uint32_t n)
{
  uint32_t offset = (uint32_t)abs(d);

  // Avoid negative of samples
  if (d > n) {
    srsran_vec_cf_zero(c, n);
    return;
  }

  // Do inner product on cyclic shifted version of a
  if (d > 0) {
    srsran_vec_prod_conj_ccc(&a[0], &a[offset], &c[0], n - offset);
    srsran_vec_prod_conj_ccc(&a[n - offset], &a[0], &c[n - offset], offset);
    return;
  }

  if (d < 0) {
    srsran_vec_prod_conj_ccc(&a[offset], &a[0], &c[0], n - offset);
    srsran_vec_prod_conj_ccc(&a[0], &a[n - offset], &c[n - offset], offset);
    return;
  }

  srsran_vec_prod_conj_ccc(a, a, c, n);
}