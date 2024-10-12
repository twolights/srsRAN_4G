//
// Created by ykchen on 9/8/24.
//

#include "srsran/phy/mdct/differential_prod.h"
#include "srsran/phy/utils/simd.h"
#include "srsran/phy/utils/vector.h"
#include "srsran/phy/utils/vector_simd.h"

// TODO move this to vector.h
//void ssb_vec_prod_conj_circ_shift(const cf_t* a, const cf_t* b, cf_t* c, uint32_t n, int shift);


void srsran_vec_cf_copy_reversed(const cf_t* a, cf_t* b, uint32_t n)
{
  for (uint32_t i = 0; i < n; i++) {
    b[i] = a[n - i - 1];
  }
}

void srsran_differential_product_simd(const cf_t* x, cf_t* z, uint32_t d, uint32_t n)
{
  int i = 0;
  int32_t half = (int32_t)n / 2;

#if SRSRAN_SIMD_CF_SIZE
  if (SRSRAN_IS_ALIGNED(x) && SRSRAN_IS_ALIGNED(z)) {
    for (; i < half - SRSRAN_SIMD_CF_SIZE + 1; i += SRSRAN_SIMD_CF_SIZE) {
      simd_cf_t a = srsran_simd_cfi_load(&x[(i + d) % half]);
      simd_cf_t b = srsran_simd_cfi_load(&x[i]);

      simd_cf_t r = srsran_simd_cf_conjprod(a, b);

      srsran_simd_cfi_store(&z[i], r);
      srsran_simd_cfi_store(&z[n - i - 1], r);
    }
  } else {
    for (; i < half - SRSRAN_SIMD_CF_SIZE + 1; i += SRSRAN_SIMD_CF_SIZE) {
      simd_cf_t a = srsran_simd_cfi_loadu(&x[(i + d) % half]);
      simd_cf_t b = srsran_simd_cfi_loadu(&x[i]);

      simd_cf_t r = srsran_simd_cf_conjprod(a, b);

      srsran_simd_cfi_storeu(&z[i], r);
      srsran_simd_cfi_storeu(&z[n - i - 1], r);
    }
  }
#endif

  for (; i < half; i++) {
    z[i] = z[n - i - 1] = x[(i + d) % half] * conjf(x[i]);
  }
}

//#define DP_USE_SIMD
//#define DP_USE_SRSRAN

// TODO make this inline?
void differential_product(const cf_t* a, cf_t* c, uint32_t d, uint32_t n)
{
//  ssb_vec_prod_conj_circ_shift(a, a, c, n, d);
  if (d == 0) {
    srsran_vec_prod_conj_ccc(a, a, c, n);
    return;
  }
#ifdef DP_USE_SIMD
  srsran_differential_product_simd(a, c, d, n);
#elif defined(DP_USE_SRSRAN)
  uint32_t half = n / 2;
  uint32_t half_d = half - d;
  srsran_vec_prod_conj_ccc(&a[d], &a[0], &c[0], half_d);
  srsran_vec_prod_conj_ccc(&a[0], &a[half_d], &c[half_d], d);
  srsran_vec_cf_copy_reversed(&c[0], &c[n - half], half);
#else
  uint32_t half = n / 2;
  cf_t x_d, x;
  for (uint32_t i = 0; i < half; i++) {
    x_d = a[(i + d) % half];
    x = a[i];
    c[i] = c[n - i - 1] = crealf(x_d) * crealf(x) + cimagf(x_d) * cimagf(x) + I * (cimag(x_d) * creal(x) - creal(x_d) * cimag(x));
  }
#endif
  // If n is odd, the middle element is multiplied by its conjugate
//  c[half + 1] = a[half + 1 + d]  * conj(a[half + 1]);
}