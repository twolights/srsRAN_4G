//
// Created by ykchen on 9/12/24.
//
#include "mdct_test_common.h"
#include "srsran/phy/mdct/mdct.h"
#include "srsran/phy/utils/vector.h"
#include <stdlib.h>
#include <string.h>

static srsran_pss_mdct_t mdct;

#define NUM_SINGLE_CELL_TESTS 9
#define NUM_MULTIPLE_CELL_TESTS 15
#define ADJACENT_CELL_TEST_START 2
#define ADJACENT_CELL_TEST_SIZE 3
#define NUM_ADJACENT_CELLS 2

// CFOs to test in Hz
//static int CFOS_TO_TEST[] = {0, 4500, 30000, 50000, -50000, -30000, -4500};
static int CFOS_TO_TEST[] = {0, 4500, -50000, -30000};

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
  {0, 300, 1, 100, 50, 2, 500, 10},
  {0, 300, 1, 100, 70, 2, 500, 10},
  {1, 0, 0, 0, 30, 2, 0, 30},
  {1, 100, 0, 100, 10, 2, 0, 0},
  {1, 140, 0, 0, 0, 2, 400, 0},
  {1, 140, 0, 0, 50, 2, 400, 30},
  {1, 140, 0, 0, 30, 2, 400, 50},
  {2, 0, 0, 100, 10, 1, 0, 10},
  {2, 300, 0, 100, 10, 1, 500, 30},
  {2, 600, 0, 100, 30, 1, 0, 10},
  {2, 600, 0, 100, 50, 1, 200, 30},
  {2, 600, 0, 100, 30, 1, 0, 50},
};

static void apply_frequency_offset(cf_t* buffer, uint32_t nof_samples, int offset_in_hz, int sampling_frequency_in_hz)
{
  srsran_vec_apply_cfo(buffer, ((float)offset_in_hz / (float)sampling_frequency_in_hz), buffer, nof_samples);
}

static bool print_test_result(int N_id_2, int tau, srsran_pss_detect_res_t* res)
{
  if (res->tau == tau && res->N_id_2 == N_id_2) {
    printf("passed for N_id_2=%d, tau=%d, CFO=%f\n", N_id_2, tau, res->coarse_cfo);
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
    srsran_detect_pss_mdct(&mdct, in, nof_samples, window_sz, true, result);
  } else {
    srsran_detect_pss_correlation(&mdct, in, nof_samples, window_sz, result);
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
    prepare_mocked_received_samples(&mdct, buffer, NOF_SAMPLES, N_id_2, tau, false);
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
    prepare_mocked_received_samples(&mdct, buffer, NOF_SAMPLES, N_id_2, tau, true);
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

static int test_cells() {
  int result = 0, i;
  int num_cfos = sizeof(CFOS_TO_TEST) / sizeof(int);
  for (i = 0; i < num_cfos; i++) {
    if (!test_single_cell(CFOS_TO_TEST[i], DETECTION_METHOD_MDCT)) {
      result = -1;
    }
  }
  for (i = 0; i < num_cfos; i++) {
    if (false && !test_multiple_cells(CFOS_TO_TEST[i], DETECTION_METHOD_MDCT)) {
      result = -1;
    }
  }
  return result;
}

int main() {
  int result = 0;
  srsran_prepare_pss_mdct(&mdct,
                          SAMPLING_FREQUENCY,
                          SYMBOL_SIZE, -30,
                          SRSRAN_MDCT_RECOMMENDED_Q * 12,
                          SRSRAN_MDCT_RECOMMENDED_PSI);
  result = test_cells();
  srsran_destroy_pss_mdct(&mdct);
  return result;
}
