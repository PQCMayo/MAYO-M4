#ifndef ARITHMETIC_M4F_H
#define ARITHMETIC_M4F_H

#include <stdint.h>
#include "mayo.h"

// M4R
// TODO: remove all the _notbitsliced_ from the names
void multiply_P1_right_m4f_O_asm2(uint32_t *table, const unsigned char *O, int col);
void multiply_P1t_right_notbitsliced_m4f_V_V_O_asm(uint32_t *acc, const uint32_t *table, const uint64_t *P1, int rows);
void multiply_P1t_right_notbitsliced_m4f_first_V_V_O_asm(uint32_t *acc, const uint32_t *table, const uint64_t *P1);
void multiply_P1_right_notbitsliced_m4f_V_V_O_asm(uint32_t *acc, const uint32_t *table, const uint64_t *P1, int rows);

void multiply_P1_right_m4f_K_asm2_transposed(uint32_t *table, const unsigned char *O, int col);
void multiply_P1_right_notbitsliced_m4f_V_V_K_asm(uint32_t *acc, const uint32_t *table, const uint64_t *P1, int rows);


// EF
void ef_bitslice_asm(uint32_t *a_bs, unsigned char *a);
void ef_unbitslice_asm(unsigned char *a, uint32_t *a_bs);
uint8_t ef_inner1_asm(uint32_t *pivot_row, uint32_t *a_bs_row, int pivot_col, int pivot_row_lower_bound, int pivot_row_ctr, uint8_t *pivot_is_zero, int row_upper_bound);
void    ef_inner2_asm(uint32_t *pivot_row2, uint32_t *pivot_row, uint8_t pivot);
void    ef_inner3_asm(uint32_t *a_bs_row, uint32_t *pivot_row, int pivot_row_ctr, int pivot_row_lower_bound, int pivot_row_upper_bound, int pivot_is_zero);
int     ef_inner4_asm(uint32_t *a_bs, uint32_t *pivot_row, int pivot_row_lower_bound, int pivot_row_ctr, int pivot_col, int pivot_is_nonzero);
void    backsub_inner_asm(uint8_t *a1, uint8_t *a2, uint8_t *a3, int correct_column, int r);

// matmul
void mul_add_mat_x_m_mat_m4f_K_V_O_triangular_asm(uint64_t *acc, const unsigned char *mat, const uint64_t *bs_mat);
void mul_add_mat_x_m_mat_m4f_K_V_K_triangular_asm(uint64_t *acc, const unsigned char *mat, const uint64_t *bs_mat);
void mul_add_mat_trans_x_m_mat_m4f_V_O_O_asm(uint64_t *acc, const unsigned char *mat, const uint64_t *bs_mat);


// verification
// TODO: remove or implment
//#if N_MAX <= 78
//void calculate_PS_m4f_asm(uint32_t *acc, const uint32_t *P, const unsigned char *S, const int cols);
//#else
void calculate_PS_m4f_stack_asm(uint64_t *acc, const uint64_t *P, const unsigned char *S, const int cols);
void multiply_bins_stack_asm(uint64_t *out, uint64_t *bins, const int cols);
//#endif

void multiply_bins_asm(uint64_t *out, uint64_t *bins, const int cols);
void calculate_SPS_m4f_asm(uint64_t *acc, const uint64_t *PS, const unsigned char * S);

#endif