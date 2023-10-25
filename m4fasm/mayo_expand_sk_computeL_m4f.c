#include "mayo.h"

#include "bitsliced_arithmetic.h"

void mayo_expand_sk_computeL_inner1_m4f(uint32_t *p2, const uint32_t *p1, const uint8_t *O, const uint8_t *O_max);
void mayo_expand_sk_computeL_inner2_m4f(uint32_t *p2, const uint32_t *p1, const uint8_t *O, size_t c);

void mayo_expand_sk_computeL_m4f(const uint32_t *restrict bitsliced_P1, const uint8_t *restrict O, uint32_t *restrict bitsliced_P2)
{
    const int param_m = M_MAX;
    const int m_legs = param_m / 32;
    const int param_o = O_MAX;
    const int param_v = V_MAX;
    int bs_mat_entries_used;

    for (int i = 0; i < m_legs; i++)
    {
        bs_mat_entries_used = 0;
        for (int r = 0; r < param_v - 1; r++)
        {
            bs_mat_entries_used++;
            int c = r + 1;
            mayo_expand_sk_computeL_inner1_m4f(bitsliced_P2 + m_legs * 4 * (r * param_o) + i, bitsliced_P1 + m_legs * 4 * bs_mat_entries_used + i, O + c * param_o, O + param_v * param_o);
            bs_mat_entries_used += param_v - (r + 1);
        }
    }

    for (int i = 0; i < m_legs; i++)
    {
        bs_mat_entries_used = 0;
        for (int c = 1; c < param_v; c++)
        {
            mayo_expand_sk_computeL_inner2_m4f(bitsliced_P2 + m_legs * 4 * (c * param_o) + i, bitsliced_P1 + m_legs * 4 * (bs_mat_entries_used + 1) + i , O, c);
            bs_mat_entries_used += 1;
        }
    }

}
