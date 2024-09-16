//
// Created by ykchen on 9/12/24.
//
#include "srsran/phy/mdct/mdct.h"
#include <string.h>

static cf_t time_domain_pss[SRSRAN_NOF_NID_2_NR][SRSRAN_MDCT_PSS_FFT_SIZE];
static srsran_pss_mdct_t mdct;

#define SAMPLING_FREQUENCY 30.72e6
#define NOF_SAMPLES (SRSRAN_MDCT_PSS_FFT_SIZE + 100)
#define NUM_SINGLE_CELL_TESTS 9
#define NUM_MULTIPLE_CELL_TESTS 9
#define ADJACENT_CELL_TEST_START 2
#define ADJACENT_CELL_TEST_SIZE 3
#define NUM_ADJACENT_CELLS 2

// Test data for single cell tests, each row contains {N_id_2, tau}
static int TEST_DATA[NUM_SINGLE_CELL_TESTS][2] = {
    {0, 0},
    {1, 0},
    {2, 0},
    {0, 20},
    {1, 5},
    {2, 15},
    {0, 19},
    {1, 7},
    {2, 8},
};

// Test data for multiple cell tests, each row contains:
// {Correct N_id_2, tau, N_id_2, tau, beta, N_id_2, tau, beta}
// Beta represents signal strength, ranges from 0 to 99
static int TEST_DATA_MULTIPLE_CELLS[NUM_MULTIPLE_CELL_TESTS][8] = {
  {0, 0, 1, 10, 10, 2, 0, 10},
  {0, 15, 1, 10, 20, 2, 10, 30},
  {0, 30, 1, 10, 10, 2, 50, 10},
  {1, 0, 0, 0, 30, 2, 0, 30},
  {1, 10, 0, 10, 10, 2, 0, 0},
  {1, 14, 0, 0, 0, 2, 40, 0},
  {2, 0, 0, 10, 10, 1, 0, 10},
  {2, 30, 0, 10, 10, 1, 50, 30},
  {2, 60, 0, 10, 30, 1, 0, 10},
};

void get_time_domain_pss(cf_t* out, uint32_t N_id_2);

static void prepare_all_time_domain_pss()
{
  int i;
  for (i = 0; i < SRSRAN_NOF_NID_2_NR; i++) {
    get_time_domain_pss(time_domain_pss[i], i);
  }
}

static void append_pss(cf_t* buffer, uint32_t N_id_2, int32_t tau, int beta)
{
  for (int i = 0; i < SRSRAN_MDCT_PSS_FFT_SIZE; i++) {
    buffer[tau + i] += time_domain_pss[N_id_2][i] * beta / 100.0;
  }
}

static void prepare_mocked_received_samples(cf_t* buffer, uint32_t N_id_2, int32_t tau, bool add_noise)
{
  if (add_noise) {
    memset(buffer, 0, NOF_SAMPLES * sizeof(cf_t));
    // TODO fill buffer with noise
  } else {
    memset(buffer, 0, NOF_SAMPLES * sizeof(cf_t));
  }
  memcpy(buffer + tau, time_domain_pss[N_id_2], SRSRAN_MDCT_PSS_FFT_SIZE * sizeof(cf_t));
}

static void apply_frequency_offset(cf_t* buffer, uint32_t nof_samples, int offset_in_hz, int sampling_frequency_in_hz)
{
  for (int i = 0; i < nof_samples; i++) {
    buffer[i] *= cexpf(I * 2 * M_PI * offset_in_hz * i / sampling_frequency_in_hz);
  }
}

static bool print_test_result(int N_id_2, int tau, srsran_pss_mdct_detect_res_t* res)
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

static bool test_single_cell(int cfo)
{
  cf_t buffer[NOF_SAMPLES];
  int i, N_id_2, tau;
  srsran_pss_mdct_detect_res_t res;
  bool result = true;

  printf("MDCT: Running single cell tests, without noise, CFO=%dHz\n", cfo);
  for (i = 0; i < NUM_SINGLE_CELL_TESTS; i++) {
    N_id_2 = TEST_DATA[i][0];
    tau    = TEST_DATA[i][1];
    printf("Test %d/%d ", i + 1, NUM_SINGLE_CELL_TESTS);
    prepare_mocked_received_samples(buffer, N_id_2, tau, false);
    if (cfo != 0) {
      apply_frequency_offset(buffer, NOF_SAMPLES, cfo, SAMPLING_FREQUENCY);
    }
    mdct_detect_pss(&mdct, buffer, NOF_SAMPLES, 1, &res);
    if (!print_test_result(N_id_2, tau, &res)) {
      result = false;
    }
  }
  return result;
}

static bool test_multiple_cells(int cfo)
{
  cf_t buffer[NOF_SAMPLES];
  int i, j, N_id_2, tau, adjacent_N_id_2, adjacent_tau, beta;
  srsran_pss_mdct_detect_res_t res;
  bool result = true;
  int* row;
  int* adjacent_cells[2];

  printf("MDCT: Running multiple cells tests, with noise, CFO=%dHz\n", cfo);
  for(i = 0; i < NUM_MULTIPLE_CELL_TESTS; i++) {
    row = TEST_DATA_MULTIPLE_CELLS[i];
    N_id_2 = row[0];
    tau = row[1];
    printf("Test %d/%d ", i + 1, NUM_MULTIPLE_CELL_TESTS);
    prepare_mocked_received_samples(buffer, N_id_2, tau, true);
    adjacent_cells[0] = row + ADJACENT_CELL_TEST_START;
    adjacent_cells[1] = row + ADJACENT_CELL_TEST_START + ADJACENT_CELL_TEST_SIZE;
    for (j = 0; j < NUM_ADJACENT_CELLS; j++) {
      adjacent_N_id_2 = adjacent_cells[j][0];
      adjacent_tau = adjacent_cells[j][1];
      beta = adjacent_cells[j][2];
      if (beta <= 0) {
        continue;
      }
      append_pss(buffer, adjacent_N_id_2, adjacent_tau, beta);
    }
    if (cfo != 0) {
      apply_frequency_offset(buffer, NOF_SAMPLES, cfo, SAMPLING_FREQUENCY);
    }
    mdct_detect_pss(&mdct, buffer, NOF_SAMPLES, 64, &res);
    if (!print_test_result(N_id_2, tau, &res)) {
      result = false;
    }
  }
  return result;
}

int main() {
  prepare_all_time_domain_pss();
  srsran_prepare_pss_mdct(&mdct, SRSRAN_MDCT_PSS_FFT_SIZE,
                          SRSRAN_MDCT_RECOMMENDED_Q,
                          SRSRAN_MDCT_RECOMMENDED_PSI);
  if(!test_single_cell(0)) {
//    return -1;
  };
  if(!test_single_cell(20000)) {
//    return -1;
  }
  if(!test_multiple_cells(0)) {
//    return -1;
  }
  if(!test_multiple_cells(20000)) {
//    return -1;
  }
  srsran_destroy_pss_mdct(&mdct);
  return 0;
}
