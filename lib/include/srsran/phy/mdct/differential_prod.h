//
// Created by ykchen on 9/8/24.
//

#ifndef SRSRAN_DIFFERENTIAL_PROD_H
#define SRSRAN_DIFFERENTIAL_PROD_H

#include <complex.h>
#include "srsran/phy/utils/vector.h"

inline void differential_product(const cf_t* a, cf_t* c, uint32_t d, uint32_t n)
{
  //  ssb_vec_prod_conj_circ_shift(a, a, c, n, d);
  //  uint32_t half = n / 2, half_d = half - d;
  if (d == 0) {
    srsran_vec_prod_conj_ccc(a, a, c, n);
    return;
  }
  //  srsran_vec_prod_conj_ccc(&a[d], &a[0], &c[0], half_d);
  //  srsran_vec_prod_conj_ccc(&a[0], &a[half_d], &c[half_d], d);
  //  srsran_vec_cf_copy_reversed(&c[0], &c[n - half], half);
  //  // If n is odd, the middle element is multiplied by its conjugate
  //  c[half + 1] = a[half + 1 + d]  * conj(a[half + 1]);
  srsran_vec_prod_conj_ccc(&a[d], &a[0], &c[0], n - d);
  srsran_vec_prod_conj_ccc(&a[0], &a[n - d], &c[n - d], d);
}

#endif // SRSRAN_DIFFERENTIAL_PROD_H
