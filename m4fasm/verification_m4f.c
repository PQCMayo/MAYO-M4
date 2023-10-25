#include <stdint.h>

#include "mayo.h"
#include "m4f_asm.h"
#include "bitsliced_arithmetic.h"


void bitsliced_m_calculate_SPS_m4f_asm(uint32_t *, const uint32_t *, const unsigned char *);
void bitsliced_multiply_bins_asm(uint32_t *, uint32_t *, const int);

#if N_MAX <= 78
void bitsliced_m_calculate_PS_m4f_asm(uint32_t *, const uint32_t *, const unsigned char *, const int);
#else
void bitsliced_m_calculate_PS_m4f_stack_asm(uint32_t *, const uint32_t *, const unsigned char *, const int);
void bitsliced_multiply_bins_stack_asm(uint32_t *, uint32_t *, const int);
#endif


// compute P * S^t = [ P1  P2 ] * [S1] = [P1*S1 + P2*S2] in bitsliced form
//                   [  0  P3 ]   [S2]   [        P3*S2]
void bitsliced_m_calculate_PS_m4f(const uint32_t *P1, const uint32_t *P2, const uint32_t *P3, const unsigned char *S, uint32_t *PS) {
    const int m = M_MAX;
    const int v = V_MAX;
    const int o = O_MAX;
    const int k = K_MAX;
    const int n = o + v;
    const int m_legs = m / 32;

    #if N_MAX > 78
    uint32_t accumulator[16 * M_MAX * N_MAX / 8] = {0};
    int P1_used;
    int P3_used;
    for (int col = 0; col < k; col++) {
        memset(accumulator, 0, sizeof accumulator);
        P1_used = 0;
        for (int row = 0; row < v; row++) {
                bitsliced_m_calculate_PS_m4f_stack_asm(accumulator + ( row * 16 )*m_legs * 4 ,  P1 + P1_used * m_legs * 4, S + col*n + row, v-row);
                P1_used += v-row;
                bitsliced_m_calculate_PS_m4f_stack_asm(accumulator + ( row * 16  )*m_legs * 4 , P2 + (row * o)*m_legs * 4, S + (col * n) + v, o);
        }

        P3_used = 0;
        for (int row = v; row < n; row++) {
            bitsliced_m_calculate_PS_m4f_stack_asm(accumulator + ( row * 16 )*m_legs * 4 , P3 + P3_used * m_legs * 4, S + col * n + row, n-row);
            P3_used += (n-row);

        }
        bitsliced_multiply_bins_stack_asm(accumulator, PS + col * m_legs * 4, n);
    }
    #else
    uint32_t accumulator[16 * M_MAX * K_MAX * N_MAX / 8] = {0};
    int P1_used = 0;

    //bitsliced_m_calculate_PS_m4f_P1_asm(accumulator, P1, S);
    for (int row = 0; row < v; row++) {

        bitsliced_m_calculate_PS_m4f_asm(accumulator + ( (row * k) * 16)*m_legs * 4,  P1 + P1_used * m_legs * 4, S+row, v-row);
        P1_used += (v-row);
        bitsliced_m_calculate_PS_m4f_asm(accumulator + ( (row * k) * 16)*m_legs * 4, P2 + (row * o)*m_legs * 4, S+v, o);
    }

    int P3_used = 0;
    for (int row = v; row < n; row++) {
        bitsliced_m_calculate_PS_m4f_asm(accumulator + ( (row * k) * 16)*m_legs * 4,  P3 + P3_used * m_legs * 4, S+row, n-row);
        P3_used += (n-row);
    }

    // multiply stuff according to the bins of the accumulator and add to PS.
    bitsliced_multiply_bins_asm(accumulator, PS, n*k);
    #endif
}


void bitsliced_m_calculate_SPS_m4f(const uint32_t *PS, const unsigned char *S, uint32_t *SPS){
    const int k = K_MAX;
    const int m = M_MAX;
    uint32_t accumulator[16*M_MAX*K_MAX*K_MAX/8] = {0};
    const int m_legs = m/32;

    for (int col = 0; col < k; col += 1) {
        bitsliced_m_calculate_SPS_m4f_asm(accumulator + ( (col) * 16 )*m_legs * 4, PS + (col) * m_legs * 4, S);
    }

    bitsliced_multiply_bins_asm(accumulator, SPS, k*k);
}