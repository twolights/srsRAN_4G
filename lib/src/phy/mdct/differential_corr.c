//
// Created by ykchen on 9/8/24.
//

#include "differential_corr.h"
#include "srsran/phy/utils/vector.h"
#include "srsran/phy/utils/vector_simd.h"

static void differential_product(const cf_t* a, const cf_t* b, cf_t* c, uint32_t d, uint32_t n) {
  // Avoid negative of samples
  if (d > n) {
    srsran_vec_cf_zero(c, n);
    return;
  }

  if (d == 0) {
    srsran_vec_prod_conj_ccc(a, b, c, n);
    return;
  }

  // Do inner product on cyclic shifted vectors
  srsran_vec_prod_conj_ccc(&a[d], &b[0], &c[0], n - d);
  srsran_vec_prod_conj_ccc(&a[0], &b[d], &c[n - d], d);
}