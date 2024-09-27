//
// Created by ykchen on 9/27/24.
//

#include <assert.h>
#include <math.h>
#include <stdio.h>

#define NUMBER_OF_SAMPLES 200

void unwrap_phase(const float* phase, float* target, size_t length);

static const float wrapped_phase[NUMBER_OF_SAMPLES] = {
    0.00000f, 0.02045f, 0.04091f, 0.06136f, 0.08182f, 0.10227f, 0.12272f, 0.14318f,
    0.16363f, 0.18409f, 0.20454f, 0.22500f, 0.24545f, 0.26590f, 0.28636f, 0.30681f,
    0.32727f, 0.34772f, 0.36817f, 0.38863f, 0.40908f, 0.42954f, 0.44999f, 0.47044f,
    0.49090f, 0.51135f, 0.53181f, 0.55226f, 0.57271f, 0.59317f, 0.61362f, 0.63408f,
    0.65453f, 0.67499f, 0.69544f, 0.71589f, 0.73635f, 0.75680f, 0.77726f, 0.79771f,
    0.81816f, 0.83862f, 0.85907f, 0.87953f, 0.89998f, 0.92043f, 0.94089f, 0.96134f,
    0.98180f, 1.00225f, 1.02270f, 1.04316f, 1.06361f, 1.08407f, 1.10452f, 1.12498f,
    1.14543f, 1.16588f, 1.18634f, 1.20679f, 1.22725f, 1.24770f, 1.26815f, 1.28861f,
    1.30906f, 1.32952f, 1.34997f, 1.37042f, 1.39088f, 1.41133f, 1.43179f, 1.45224f,
    1.47270f, 1.49315f, 1.51360f, 1.53406f, 1.55451f, 1.57497f, 1.59542f, 1.61587f,
    1.63633f, 1.65678f, 1.67724f, 1.69769f, 1.71814f, 1.73860f, 1.75905f, 1.77951f,
    1.79996f, 1.82041f, 1.84087f, 1.86132f, 1.88178f, 1.90223f, 1.92269f, 1.94314f,
    1.96359f, 1.98405f, 2.00450f, 2.02496f, 2.04541f, 2.06586f, 2.08632f, 2.10677f,
    2.12723f, 2.14768f, 2.16813f, 2.18859f, 2.20904f, 2.22950f, 2.24995f, 2.27041f,
    2.29086f, 2.31131f, 2.33177f, 2.35222f, 2.37268f, 2.39313f, 2.41358f, 2.43404f,
    2.45449f, 2.47495f, 2.49540f, 2.51585f, 2.53631f, 2.55676f, 2.57722f, 2.59767f,
    2.61812f, 2.63858f, 2.65903f, 2.67949f, 2.69994f, 2.72040f, 2.74085f, 2.76130f,
    2.78176f, 2.80221f, 2.82267f, 2.84312f, 2.86357f, 2.88403f, 2.90448f, 2.92494f,
    2.94539f, 2.96584f, 2.98630f, 3.00675f, 3.02721f, 3.04766f, 3.06811f, 3.08857f,
    3.10902f, 3.12948f, -3.13325f, -3.11280f, -3.09235f, -3.07189f, -3.05144f, -3.03098f, -3.01053f, -2.99008f,
    -2.96962f, -2.94917f, -2.92871f, -2.90826f, -2.88780f, -2.86735f, -2.84690f, -2.82644f,
    -2.80599f, -2.78553f, -2.76508f, -2.74463f, -2.72417f, -2.70372f, -2.68326f, -2.66281f,
    -2.64236f, -2.62190f, -2.60145f, -2.58099f, -2.56054f, -2.54009f, -2.51963f, -2.49918f,
    -2.47872f, -2.45827f, -2.43781f, -2.41736f, -2.39691f, -2.37645f, -2.35600f, -2.33554f,
    -2.31509f, -2.29464f, -2.27418f, -2.25373f, -2.23327f, -2.21282f
};

static const float unwrapped_phase[NUMBER_OF_SAMPLES] = {
    0.00000f, 0.02045f, 0.04091f, 0.06136f, 0.08182f, 0.10227f, 0.12272f, 0.14318f,
    0.16363f, 0.18409f, 0.20454f, 0.22500f, 0.24545f, 0.26590f, 0.28636f, 0.30681f,
    0.32727f, 0.34772f, 0.36817f, 0.38863f, 0.40908f, 0.42954f, 0.44999f, 0.47044f,
    0.49090f, 0.51135f, 0.53181f, 0.55226f, 0.57271f, 0.59317f, 0.61362f, 0.63408f,
    0.65453f, 0.67499f, 0.69544f, 0.71589f, 0.73635f, 0.75680f, 0.77726f, 0.79771f,
    0.81816f, 0.83862f, 0.85907f, 0.87953f, 0.89998f, 0.92043f, 0.94089f, 0.96134f,
    0.98180f, 1.00225f, 1.02270f, 1.04316f, 1.06361f, 1.08407f, 1.10452f, 1.12498f,
    1.14543f, 1.16588f, 1.18634f, 1.20679f, 1.22725f, 1.24770f, 1.26815f, 1.28861f,
    1.30906f, 1.32952f, 1.34997f, 1.37042f, 1.39088f, 1.41133f, 1.43179f, 1.45224f,
    1.47270f, 1.49315f, 1.51360f, 1.53406f, 1.55451f, 1.57497f, 1.59542f, 1.61587f,
    1.63633f, 1.65678f, 1.67724f, 1.69769f, 1.71814f, 1.73860f, 1.75905f, 1.77951f,
    1.79996f, 1.82041f, 1.84087f, 1.86132f, 1.88178f, 1.90223f, 1.92269f, 1.94314f,
    1.96359f, 1.98405f, 2.00450f, 2.02496f, 2.04541f, 2.06586f, 2.08632f, 2.10677f,
    2.12723f, 2.14768f, 2.16813f, 2.18859f, 2.20904f, 2.22950f, 2.24995f, 2.27041f,
    2.29086f, 2.31131f, 2.33177f, 2.35222f, 2.37268f, 2.39313f, 2.41358f, 2.43404f,
    2.45449f, 2.47495f, 2.49540f, 2.51585f, 2.53631f, 2.55676f, 2.57722f, 2.59767f,
    2.61812f, 2.63858f, 2.65903f, 2.67949f, 2.69994f, 2.72040f, 2.74085f, 2.76130f,
    2.78176f, 2.80221f, 2.82267f, 2.84312f, 2.86357f, 2.88403f, 2.90448f, 2.92494f,
    2.94539f, 2.96584f, 2.98630f, 3.00675f, 3.02721f, 3.04766f, 3.06811f, 3.08857f,
    3.10902f, 3.12948f, 3.14993f, 3.17039f, 3.19084f, 3.21129f, 3.23175f, 3.25220f,
    3.27266f, 3.29311f, 3.31356f, 3.33402f, 3.35447f, 3.37493f, 3.39538f, 3.41583f,
    3.43629f, 3.45674f, 3.47720f, 3.49765f, 3.51811f, 3.53856f, 3.55901f, 3.57947f,
    3.59992f, 3.62038f, 3.64083f, 3.66128f, 3.68174f, 3.70219f, 3.72265f, 3.74310f,
    3.76355f, 3.78401f, 3.80446f, 3.82492f, 3.84537f, 3.86582f, 3.88628f, 3.90673f,
    3.92719f, 3.94764f, 3.96810f, 3.98855f, 4.00900f, 4.02946f, 4.04991f, 4.07037
};

int main() {
    float target[NUMBER_OF_SAMPLES];
    unwrap_phase(wrapped_phase, target, NUMBER_OF_SAMPLES);
    for (int i = 0; i < NUMBER_OF_SAMPLES; i++) {
      assert(fabsf(target[i] - unwrapped_phase[i]) < 1e-5);
    }
    return 0;
}
