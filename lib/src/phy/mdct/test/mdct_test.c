//
// Created by ykchen on 9/12/24.
//
#include "srsran/phy/mdct/mdct.h"
#include <stdlib.h>
#include <string.h>

static srsran_pss_mdct_t mdct;

#define DETECTION_METHOD_CORRELATION 0
#define DETECTION_METHOD_MDCT 1

#define SAMPLING_FREQUENCY 23.04e6
#define SYMBOL_SIZE 1536
#define NOF_SAMPLES (SYMBOL_SIZE * 2)
#define NUM_SINGLE_CELL_TESTS 9
#define NUM_MULTIPLE_CELL_TESTS 9
#define ADJACENT_CELL_TEST_START 2
#define ADJACENT_CELL_TEST_SIZE 3
#define NUM_ADJACENT_CELLS 2
#define CFO_TO_TEST (15000 * 3)

// Test data for single cell tests, each row contains {N_id_2, tau}
static int TEST_DATA[NUM_SINGLE_CELL_TESTS][2] = {
    {0, 0},
    {1, 0},
    {2, 0},
    {0, 200},
    {1, 100},
    {2, 150},
    {0, 400},
    {1, 700},
    {2, 800},
};

// Test data for multiple cell tests, each row contains:
// {Correct N_id_2, tau, N_id_2, tau, beta, N_id_2, tau, beta}
// Beta represents signal strength, ranges from 0 to 99
static int TEST_DATA_MULTIPLE_CELLS[NUM_MULTIPLE_CELL_TESTS][8] = {
  {0, 0, 1, 100, 10, 2, 0, 10},
  {0, 150, 1, 100, 20, 2, 100, 30},
  {0, 300, 1, 100, 10, 2, 500, 10},
  {1, 0, 0, 0, 30, 2, 0, 30},
  {1, 100, 0, 100, 10, 2, 0, 0},
  {1, 140, 0, 0, 0, 2, 400, 0},
  {2, 0, 0, 100, 10, 1, 0, 10},
  {2, 300, 0, 100, 10, 1, 500, 30},
  {2, 600, 0, 100, 30, 1, 0, 10},
};

static void append_pss(srsran_pss_mdct_t* mdct, cf_t* buffer, uint32_t N_id_2, int32_t tau, int beta)
{
  for (int i = 0; i < mdct->symbol_sz; i++) {
    buffer[tau + i] += mdct->pss_x[N_id_2][i] * beta / 100.0;
  }
}

static void prepare_mocked_received_samples(srsran_pss_mdct_t* mdct, cf_t* buffer, uint32_t N_id_2, int32_t tau, bool add_noise)
{
  if (add_noise) {
    memset(buffer, 0, NOF_SAMPLES * sizeof(cf_t));
    // TODO fill buffer with noise
  } else {
    memset(buffer, 0, NOF_SAMPLES * sizeof(cf_t));
  }
  append_pss(mdct, buffer, N_id_2, tau, 100);
}

static void apply_frequency_offset(cf_t* buffer, uint32_t nof_samples, int offset_in_hz, int sampling_frequency_in_hz)
{
  for (int i = 0; i < nof_samples; i++) {
    buffer[i] *= cexpf(I * 2 * M_PI * offset_in_hz * i / sampling_frequency_in_hz);
  }
}

static bool print_test_result(int N_id_2, int tau, srsran_pss_detect_res_t* res)
{
  if (res->tau == tau && res->N_id_2 == N_id_2) {
    printf("passed for N_id_2=%d, tau=%d\n", N_id_2, tau);
    return true;
  }
  printf("failed:");
  if (res->tau != tau) {
    printf(" expected tau=%d, got tau=%d", tau, res->tau);
  }
  if (res->N_id_2 != N_id_2) {
    printf(" expected N_id_2=%d, got N_id_2=%d", N_id_2, res->N_id_2);
  }
  printf("\n");
  return false;
}

static void detect_pss(const cf_t* in, uint32_t nof_samples,
                       uint32_t window_sz, srsran_pss_detect_res_t* result,
                       int method)
{
  if(method == DETECTION_METHOD_MDCT) {
    mdct_detect_pss(&mdct, in, nof_samples, window_sz, result);
  } else {
    correlation_detect_pss(&mdct, in, nof_samples, window_sz, result);
  }
}

static bool test_single_cell(int cfo, int method)
{
  cf_t buffer[NOF_SAMPLES];
  int i, N_id_2, tau;
  srsran_pss_detect_res_t res;
  bool result = true;

  if(method == DETECTION_METHOD_MDCT) {
    printf("MDCT: Running single cell tests, without noise, CFO=%dHz\n", cfo);
  } else {
    printf("Correlation: Running single cell tests, without noise, CFO=%dHz\n", cfo);
  }
  for (i = 0; i < NUM_SINGLE_CELL_TESTS; i++) {
    N_id_2 = TEST_DATA[i][0];
    tau    = TEST_DATA[i][1];
    printf("Test %d/%d ", i + 1, NUM_SINGLE_CELL_TESTS);
    prepare_mocked_received_samples(&mdct, buffer, N_id_2, tau, false);
    if (cfo != 0) {
      apply_frequency_offset(buffer, NOF_SAMPLES, cfo, SAMPLING_FREQUENCY);
    }
    detect_pss(buffer, NOF_SAMPLES, 1, &res, method);
    if (!print_test_result(N_id_2, tau, &res)) {
      result = false;
    }
  }
  return result;
}

static bool test_multiple_cells(int cfo, int method)
{
  cf_t buffer[NOF_SAMPLES];
  int i, j, N_id_2, tau, adjacent_N_id_2, adjacent_tau, beta;
  srsran_pss_detect_res_t res;
  bool result = true;
  int* row;
  int* adjacent_cells[2];

  if(method == DETECTION_METHOD_MDCT) {
    printf("MDCT: Running multiple cells tests, without noise, CFO=%dHz\n", cfo);
  } else {
    printf("Correlation: Running multiple cells tests, without noise, CFO=%dHz\n", cfo);
  }
  for(i = 0; i < NUM_MULTIPLE_CELL_TESTS; i++) {
    row = TEST_DATA_MULTIPLE_CELLS[i];
    N_id_2 = row[0];
    tau = row[1];
    printf("Test %d/%d ", i + 1, NUM_MULTIPLE_CELL_TESTS);
    prepare_mocked_received_samples(&mdct, buffer, N_id_2, tau, true);
    adjacent_cells[0] = row + ADJACENT_CELL_TEST_START;
    adjacent_cells[1] = row + ADJACENT_CELL_TEST_START + ADJACENT_CELL_TEST_SIZE;
    for (j = 0; j < NUM_ADJACENT_CELLS; j++) {
      adjacent_N_id_2 = adjacent_cells[j][0];
      adjacent_tau = adjacent_cells[j][1];
      beta = adjacent_cells[j][2];
      if (beta <= 0) {
        continue;
      }
      append_pss(&mdct, buffer, adjacent_N_id_2, adjacent_tau, beta);
    }
    if (cfo != 0) {
      apply_frequency_offset(buffer, NOF_SAMPLES, cfo, SAMPLING_FREQUENCY);
    }
    detect_pss(buffer, NOF_SAMPLES, 1, &res, method);
    if (!print_test_result(N_id_2, tau, &res)) {
      result = false;
    }
  }
  return result;
}

static void test_mdct_on_samples_file(const char* filename, int cfo)
{
  size_t l = 29953;
  cf_t* buffer = (cf_t*)malloc(l * sizeof(cf_t));
  FILE* fp = fopen(filename, "rb");
  fread(buffer, sizeof(cf_t), l, fp);
  fclose(fp);

  if(cfo != 0) {
    apply_frequency_offset(buffer, l, cfo, SAMPLING_FREQUENCY);
  }
  srsran_pss_detect_res_t res;
  mdct_detect_pss(&mdct, buffer, l, 1, &res);
  printf("MDCT[CFO=%d] Detected N_id_2=%d, tau=%d, peak=%lf\n", cfo, res.N_id_2, res.tau, res.peak_value);
  correlation_detect_pss(&mdct, buffer, l, 1, &res);
  printf("Corr[CFO=%d]: Detected N_id_2=%d, tau=%d, peak=%lf\n", cfo, res.N_id_2, res.tau, res.peak_value);
  free(buffer);
}

static int test_cells(bool perform) {
  if(!perform) {
    return 0;
  }
  int result = 0;
  if(!test_single_cell(0, DETECTION_METHOD_MDCT)) {
    result = -1;
  };
  if(!test_single_cell(CFO_TO_TEST, DETECTION_METHOD_MDCT)) {
    result = -1;
  }
  if(!test_multiple_cells(0, DETECTION_METHOD_MDCT)) {
    result = -1;
  }
  if(!test_multiple_cells(CFO_TO_TEST, DETECTION_METHOD_MDCT)) {
    result = -1;
  }
  if(!test_single_cell(0, DETECTION_METHOD_CORRELATION)) {
    result = -1;
  };
  if(true == test_single_cell(CFO_TO_TEST, DETECTION_METHOD_CORRELATION)) {
    result = -1;
  }
  if(!test_multiple_cells(0, DETECTION_METHOD_CORRELATION)) {
    result = -1;
  }
  if(true == test_multiple_cells(CFO_TO_TEST, DETECTION_METHOD_CORRELATION)) {
    result = -1;
  }
  return result;
}

static int test_mdct_on_samples(bool perform) {
  if (!perform) {
    return 0;
  }
  test_mdct_on_samples_file("ssb_1726762194-NID2-1-offset-3600.dat", 0);
  test_mdct_on_samples_file("ssb_1726762194-NID2-1-offset-3600.dat", 15000 * 10);
  test_mdct_on_samples_file("ssb_1727003781-NID2-1-offset-3608.dat", 0);
  test_mdct_on_samples_file("ssb_1727003781-NID2-1-offset-3608.dat", 15000 * 10);
  return 0;  // TODO
}

int main() {
  int result = 0;
  srsran_prepare_pss_mdct(&mdct,
                          SYMBOL_SIZE, -30,
                          SRSRAN_MDCT_RECOMMENDED_Q * 12,
                          SRSRAN_MDCT_RECOMMENDED_PSI);
  result = test_cells(true) == 0;
  result = result == 0 && test_mdct_on_samples(true) == 0 ? 0 : -1;
  srsran_destroy_pss_mdct(&mdct);
  return result;
}
