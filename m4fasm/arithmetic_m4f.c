#include "arithmetic_m4f.h"
#include "arithmetic.h"
#include "simple_arithmetic.h"
#include <stdalign.h>
#include <string.h>
#include "mayo.h"

#define MAYO_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAYO_MIN(x, y) (((x) < (y)) ? (x) : (y))

// -----------------------------------------------------------------------------------------
// NON-OPTIMIZED ARITHMETIC
// -----------------------------------------------------------------------------------------


void m_upper(const mayo_params_t* p, const uint64_t *in, uint64_t *out, int size) {
    #ifndef ENABLE_PARAMS_DYNAMIC
    (void) p;
    #endif
    // Look into AVX2'ing this 
    const int m_vec_limbs = PARAM_m_vec_limbs(p);
    int m_vecs_stored = 0;
    for (int r = 0; r < size; r++) {
        for (int c = r; c < size; c++) {
            m_vec_copy(in + m_vec_limbs * (r * size + c), out + m_vec_limbs * m_vecs_stored );
            if (r != c) {
                m_vec_add(in + m_vec_limbs * (c * size + r), out + m_vec_limbs * m_vecs_stored );
            }
            m_vecs_stored ++;
        }
    }
}

// TODO: optimize
/*
// multiplies the transpose of a single matrix with m matrices and adds result to acc
static void mul_add_mat_trans_x_m_mat(int m_legs, const unsigned char *mat, const uint64_t *bs_mat, uint64_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {
    (void) m_legs;
    (void) mat_rows;
    (void) mat_cols;
    (void) bs_mat_cols;
    mul_add_mat_trans_x_m_mat_m4f_V_O_O_asm(acc, mat, bs_mat);

}
*/


// -----------------------------------------------------------------------------------------
// OPTIMIZED ARITHMETIC
// -----------------------------------------------------------------------------------------

// multiplies a single matrix with m matrices and adds result to acc
// TODO: optimize
/*
void mul_add_mat_x_m_mat(int m_legs, const unsigned char *mat, const uint32_t *bs_mat, uint32_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {
  (void) m_legs;
  if(mat_rows == K_MAX && mat_cols == V_MAX && bs_mat_cols == O_MAX){
    mul_add_mat_x_m_mat_m4f_K_V_O_triangular_asm(acc, mat, bs_mat);
  } else if(mat_rows == K_MAX && mat_cols == V_MAX && bs_mat_cols == K_MAX){
    mul_add_mat_x_m_mat_m4f_K_V_K_triangular_asm(acc, mat, bs_mat);
  }
}
*/

// TODO: optimize
/*
static void repack_add(uint32_t *out, const uint32_t *P2, const int dim0, const int dim1){
  const int v = dim0;
  const int o = dim1;
  const int m = M_MAX;
  const int o_size = (o+7)/8;

  for (int row = 0; row < v; row++)
  {
    for(int col = 0; col< o; col+=2){
      for(int k=0; k < m ; k += 8){
        unsigned char byte0 = P2[(row * m + k)*o_size + (col/8)] >> (4*(col % 8));
        unsigned char byte1 = P2[(row * m + k + 1)*o_size + (col/8)] >> (4*(col % 8));
        unsigned char byte2 = P2[(row * m + k + 2)*o_size + (col/8)] >> (4*(col % 8));
        unsigned char byte3 = P2[(row * m + k + 3)*o_size + (col/8)] >> (4*(col % 8));
        unsigned char byte4 = P2[(row * m + k + 4)*o_size + (col/8)] >> (4*(col % 8));
        unsigned char byte5 = P2[(row * m + k + 5)*o_size + (col/8)] >> (4*(col % 8));
        unsigned char byte6 = P2[(row * m + k + 6)*o_size + (col/8)] >> (4*(col % 8));
        unsigned char byte7 = P2[(row * m + k + 7)*o_size + (col/8)] >> (4*(col % 8));

        uint32_t out0 = (byte0 & 0xF) | 
                        ((byte1 & 0xF) << (1*4) ) |
                        ((byte2 & 0xF) << (2*4) ) |
                        ((byte3 & 0xF) << (3*4) ) |
                        ((byte4 & 0xF) << (4*4) ) |
                        ((byte5 & 0xF) << (5*4) ) |
                        ((byte6 & 0xF) << (6*4) ) |
                        ((byte7 & 0xF) << (7*4) );

        uint32_t out1 = (byte0 >> 4) | 
                        (byte1 & 0xF0) |
                        (byte2 & 0xF0) << (1*4) |
                        (byte3 & 0xF0) << (2*4) |
                        (byte4 & 0xF0) << (3*4) |
                        (byte5 & 0xF0) << (4*4) |
                        (byte6 & 0xF0) << (5*4) |
                        (byte7 & 0xF0) << (6*4);

        out[0]   ^= out0;
        if(col+1 < o)
          out[m/8] ^= out1;
        out += 1;
      }
      if(col+1 < o)
        out += m/8;
    }
  }
}




static void multiply_P1P1t_right_notbitsliced_m4f(uint32_t *P1_O, const uint32_t *P1, const unsigned char *O, const int v, const int o, const int m){
    const int o_size = (o+7)/8;
    const int m_legs = m/32;
    const int m_bytes = m/2;

    uint32_t table[o_size*256];

    const uint32_t *P1t_ptr = P1;
    for (int col = -(v%2); col < v; col += 2 ){

        for (int i = 0; i < o_size; i++)
        {
            table[i] = 0;
        }

        // build table.
        multiply_P1_right_m4f_O_asm2(table, O + (col*o/2), col);
        


        // do pairs of field elements (P1t)
        if(col >= 0){
            multiply_P1t_right_notbitsliced_m4f_V_V_O_asm(&P1_O[((col)*m_legs*32) * o_size], table, P1t_ptr + (m_bytes/4) * (col), v-col-1);
        } else {
            multiply_P1t_right_notbitsliced_m4f_first_V_V_O_asm(P1_O, table, P1);
        }
        
        
        if(col >= 0)
            P1t_ptr += (m_bytes/4)*(v-col-1);
        P1t_ptr += (m_bytes/4)*(v-(col+1)-1);
        // do pairs of field elements (P1)
        multiply_P1_right_notbitsliced_m4f_V_V_O_asm(P1_O, table, P1 + (m_bytes/4) * col, col+1);
    }
}

void P1P1t_times_O(const mayo_params_t* p, const uint64_t* P1, const unsigned char* O, uint64_t* acc) {
    (void)p;
    uint32_t P1_O[(O_MAX + 7)/8 * M_MAX * V_MAX] = {0};

    // TODO: do packing outside of this function (should be able to just change the representation)
    unsigned char O_packed[V_MAX*O_MAX/2];
    for(unsigned int i=0; i < sizeof O_packed; i++){
        O_packed[i] = O[2*i] ^ (O[2*i+1] << 4);
    }

    multiply_P1P1t_right_notbitsliced_m4f(P1_O, (uint32_t *)P1, O_packed, V_MAX, O_MAX, M_MAX);
    repack_add((uint32_t *)acc, P1_O, V_MAX, O_MAX);
}
*/

// TODO: optimize
/*
static void multiply_P1_right_transposed_notbitsliced_m4f(uint32_t *P1_O, const uint32_t *P1, const unsigned char *O, const int v, const int o, const int m){
    const int o_size = (o+7)/8;
    const int m_bytes = m / 2;

    uint32_t table[o_size*256];
    for (int col = -(v%2); col < v; col += 2 ){

        for (int i = 0; i < o_size; i++)
        {
            table[i] = 0;
        }

        // build table.
        multiply_P1_right_m4f_K_asm2_transposed(table, O + col, col);
        // do pairs of field elements (P1)
        multiply_P1_right_notbitsliced_m4f_V_V_K_asm(P1_O, table, P1 + (m_bytes/4) * col, col+1);
    }

}


void P1_times_Vt(const mayo_params_t* p, const uint32_t* P1, const unsigned char* V, uint32_t* acc){
    (void)p;

    uint32_t P1_Vt[(K_MAX + 7)/8 * M_MAX * V_MAX] = {0};

    multiply_P1_right_transposed_notbitsliced_m4f(P1_Vt, P1, V, V_MAX, K_MAX, M_MAX);

    repack_add(acc, P1_Vt, V_MAX, K_MAX);
}
*/

// TODO: optimize
/*
static void multiply_P1_right_notbitsliced_m4f(uint32_t *P2, const uint32_t *P1, const unsigned char *O, const int v, const int o, const int m){
    const int o_size = (o+7)/8;
    const int m_bytes = m / 2;

    uint32_t table[o_size*256];

    for (int col = -(v%2); col < v; col += 2 ){

        for (int i = 0; i < o_size; i++)
        {
            table[i] = 0;
        }
        // build table.
        multiply_P1_right_m4f_O_asm2(table, O + (col*o/2), col);

        // do pairs of field elements (P1)
        multiply_P1_right_notbitsliced_m4f_V_V_O_asm(P2, table, P1 + (m_bytes/4) * col, col+1);

    }
}

void P1_times_O(const mayo_params_t* p, const uint32_t* P1, const unsigned char* O, uint32_t* acc){
    (void)p;
    uint32_t P1_O[(O_MAX + 7)/8 * M_MAX * V_MAX] = {0};

    // TODO: do packing outside of this function (should be able to just change the representation)
    unsigned char O_packed[V_MAX*O_MAX/2];
    for(unsigned int i=0; i < sizeof O_packed; i++){
        O_packed[i] = O[2*i] ^ (O[2*i+1] << 4);
    }
    
    multiply_P1_right_notbitsliced_m4f(P1_O, P1, O_packed, V_MAX, O_MAX, M_MAX);
    repack_add(acc, P1_O, V_MAX, O_MAX);

}
*/




// put matrix in row echelon form with ones on first nonzero entries *in
// constant time*
// TODO: optimize
/*
static void EF_m4f(unsigned char *A, int nrows, int ncols) {

    uint32_t _pivot_row[(K_MAX * O_MAX + 1 + 31) / 32 * 4];
    uint32_t _pivot_row2[(K_MAX * O_MAX + 1 + 31) / 32 * 4];
    uint32_t bitsliced_A[((K_MAX * O_MAX + 1 + 31) / 32) * 4 * M_MAX];

    int legs = (ncols + 31) / 32;


    ef_bitslice_asm(bitsliced_A, A);

    // pivot row is secret, pivot col is not

    //unsigned char inverse;
    int pivot_row = 0;
    for (int pivot_col = 0; pivot_col < ncols; pivot_col++) {

        int pivot_row_lower_bound = MAYO_MAX(0, pivot_col + nrows - ncols);
        int pivot_row_upper_bound = MAYO_MIN(nrows - 1, pivot_col);
        // the pivot row is guaranteed to be between these lower and upper bounds if
        // A has full rank

        // zero out pivot row
        for (int i = 0; i < legs * 4; i++) {
            _pivot_row[i] = 0;
        }

        // try to get a pivot row in constant time
        unsigned char pivot = 0;
        unsigned char pivot_is_zero = 1;


        // make sure pivot is non-zero
        pivot = ef_inner1_asm(_pivot_row, bitsliced_A + pivot_row_lower_bound*legs*4, pivot_col, pivot_row_lower_bound, pivot_row, &pivot_is_zero, MAYO_MIN(nrows - 1, pivot_row_upper_bound + 16));

        ef_inner2_asm(_pivot_row2, _pivot_row, pivot);

        // conditionally write pivot row to the correct row, if there is a nonzero
        // pivot
        ef_inner3_asm(bitsliced_A + pivot_row_lower_bound * legs * 4, _pivot_row2, pivot_row, pivot_row_lower_bound, pivot_row_upper_bound, pivot_is_zero ? 1 : 0);

        // eliminate entries below pivot
        pivot_row = ef_inner4_asm(bitsliced_A, _pivot_row2, pivot_row_lower_bound, pivot_row, pivot_col, pivot_is_zero ? 0 : 1);
    }

    // unbitslice the matrix A
    for (int i = 0; i < nrows; i++) {
        ef_unbitslice_asm(A + i*ncols, bitsliced_A + i * legs * 4);
    }
} */

// sample a solution x to Ax = y, with r used as randomness
// require:
// - A is a matrix with m rows and k*o+1 collumns (values in the last collum are
// not important, they will be overwritten by y) in row major order
// - y is a vector with m elements
// - r and x are k*o bytes long
// return: 1 on success, 0 on failure
// TODO: optimize
/*
int sample_solution(const mayo_params_t *p, unsigned char *A,
                           const unsigned char *y, const unsigned char *r,
                           unsigned char *x, int _k, int _o, int _m, int _A_cols ) {
   int finished;
    int col_upper_bound;
    int correct_column;

    (void) p;
    (void) _k;
    (void) _o;
    (void) _m;
    (void) _A_cols;


    const int k = K_MAX;
    const int o = O_MAX;
    const int m = M_MAX;
    const int A_cols = (K_MAX * O_MAX + 1);




    // x <- r
    for (int i = 0; i < k * o; i++) {
        x[i] = r[i];
    }

    // compute Ar;
    unsigned char Ar[M_MAX];
    for (int i = 0; i < m; i++) {
        A[k * o + i * (k * o + 1)] = 0; // clear last col of A
    }
    mat_mul(A, r, Ar, k * o + 1, m, 1);

    // move y - Ar to last column of matrix A
    for (int i = 0; i < m; i++) {
        A[k * o + i * (k * o + 1)] = sub_f(y[i], Ar[i]);
    }

    EF_m4f(A, m, k * o + 1);

    // check if last row of A (excluding the last entry of y) is zero
    unsigned char full_rank = 0;
    for (int i = 0; i < A_cols - 1; i++) {
        full_rank |= A[(m - 1) * A_cols + i];
    }


    if (full_rank == 0) {
        return 0;
    }

    // back substitution in constant time
    // the index of the first nonzero entry in each row is secret, which makes
    // things less efficient
    for (int u = m - 1; u >= 0; u--) {
        finished = 0;
        col_upper_bound = MAYO_MIN(u + (32/(m-u)), k*o);
        // the first nonzero entry in row r is between r and col_upper_bound with probability at least ~1-q^{-32}

        for (int col = u; col <= col_upper_bound; col++) {
            correct_column = (A[u * A_cols + col] != 0) && !finished;

            x[col] ^= correct_column * A[u * A_cols + A_cols - 1];

            if(u>0) {
                backsub_inner_asm(A + u * A_cols + A_cols - 1, A + col, A + A_cols - 1, correct_column, u);
            }

            finished = finished || correct_column;
        }
    }
    return 1;
}

*/


// compute P * S^t = [ P1  P2 ] * [S1] = [P1*S1 + P2*S2]
//                   [  0  P3 ]   [S2]   [        P3*S2]
// TODO: optimize
/*
void m_calculate_PS(const uint32_t *P1, const uint32_t *P2, const uint32_t *P3, const unsigned char *S,
                              const int _m, const int _v, const int _o, const int _k, uint32_t *PS) {
    (void) _m;
    (void) _v;
    (void) _o;
    (void) _k;

    const int m = M_MAX;
    const int v = V_MAX;
    const int o = O_MAX;
    const int k = K_MAX;

    const int m_legs = m / 32;                                
    const int n = o + v;

      // use more stack efficient version for MAYO_3 and MAYO_5
    #if defined(PQM4) && N_MAX > 78
    uint32_t accumulator[16 * M_MAX * N_MAX / 8];
    int P1_used;
    int P3_used;
    for (int col = 0; col < k; col++) {
        #if 0
        for(unsigned int i = 0; i < sizeof(accumulator)/4; i++) {
            accumulator[i] = 0;
        }
        #else
        memset(accumulator, 0, sizeof accumulator);
        #endif
        P1_used = 0;
        for (int row = 0; row < v; row++) {
                calculate_PS_m4f_stack_asm(accumulator + ( row * 16 )*m_legs * 4 ,  P1 + P1_used * m_legs * 4, S + col*n + row, v-row);
                P1_used += v-row;
                calculate_PS_m4f_stack_asm(accumulator + ( row * 16  )*m_legs * 4 , P2 + (row * o)*m_legs * 4, S + (col * n) + v, o);
        }

        P3_used = 0;
        for (int row = v; row < n; row++) {
            calculate_PS_m4f_stack_asm(accumulator + ( row * 16 )*m_legs * 4 , P3 + P3_used * m_legs * 4, S + col * n + row, n-row);
            P3_used += (n-row);

        }
        
        //bitsliced_multiply_bins_stack_asm(accumulator, PS + col * m_legs * 4, n);
        multiply_bins_stack_asm(PS + col * m_legs * 4, accumulator, n);
    }


    #else

    alignas (32) uint32_t accumulator[16 * M_MAX * K_MAX * N_MAX / 8] = {0};
        int P1_used = 0;

    //calculate_PS_m4f_P1_asm(accumulator, P1, S);
    for (int row = 0; row < v; row++) {

        calculate_PS_m4f_asm(accumulator + ( (row * k) * 16)*m_legs * 4,  P1 + P1_used * m_legs * 4, S+row, v-row);
        P1_used += (v-row);
        calculate_PS_m4f_asm(accumulator + ( (row * k) * 16)*m_legs * 4, P2 + (row * o)*m_legs * 4, S+v, o);
    }

    int P3_used = 0;
    for (int row = v; row < n; row++) {
        calculate_PS_m4f_asm(accumulator + ( (row * k) * 16)*m_legs * 4,  P3 + P3_used * m_legs * 4, S+row, n-row);
        P3_used += (n-row);
    }


    // multiply stuff according to the bins of the accumulator and add to PS.
    multiply_bins_asm(PS, accumulator, n*k);

    #endif
} 
*/

// TODO: optimize
/*
void m_calculate_SPS(const uint32_t *PS, const unsigned char *S, int _m, int _k, int  _n, uint32_t *SPS){
    (void) _n;
    (void) _m;
    (void) _k;

    const int m = M_MAX;
    const int k = K_MAX;

    uint32_t accumulator[16*M_MAX*K_MAX*K_MAX/8] = {0};
    const int m_legs = m/32;

    for (int col = 0; col < k; col += 1) {
        calculate_SPS_m4f_asm(accumulator + ( (col) * 16 )*m_legs * 4, PS + (col) * m_legs * 4, S);
    }

    multiply_bins_asm(SPS, accumulator, k*k);
}
*/


// TODO: optimize
/*
void m_calculate_PS_SPS(const uint64_t *P1, const uint64_t *P2, const uint64_t *P3, const unsigned char *S,
                              const int m, const int v, const int o, const int k, uint64_t *SPS) {

    (void) m;
    (void) v;
    (void) o;
    (void) k;
    // compute P * S^t = {(P1, P2), (0, P3)} * S^t = {(P1*S1 + P2*S2), (P3 * S2)}
    alignas (32) uint32_t PS[N_MAX * K_MAX * M_MAX / 8];

    m_calculate_PS((uint32_t *)P1, (uint32_t *)P2, (uint32_t *)P3, S, m, v, o, k, PS);
    m_calculate_SPS(PS, S, m, k, v+o, (uint32_t *)SPS);
}
*/


// TODO: optimize
/*
void V_times_L__V_times_P1_times_Vt(const mayo_params_t* p, const uint64_t* L, const unsigned char* V, uint64_t* M, const uint64_t* P1, uint64_t* Y) {
        alignas (32) uint32_t Pv[N_MINUS_O_MAX * K_MAX * M_MAX / 8] = {0};
        const int param_k = K_MAX;
        const int param_m = M_MAX;
        const int param_n = N_MAX;
        const int param_o = O_MAX;
        

        mul_add_mat_x_m_mat(param_m / 32, V, (uint32_t *)L, (uint32_t *)M,
                                      param_k, param_n - param_o, param_o);
        // compute all the v_i^t * P^(1) * v_j

        P1_times_Vt(p, (uint32_t *)P1, V, Pv);
        mul_add_mat_x_m_mat(param_m / 32, V, Pv,
                                      (uint32_t *) Y, param_k, param_n - param_o,
                                      param_k);
}
*/


// TODO: optimize
/*
void Ot_times_P1O_P2(const mayo_params_t* p, const uint64_t* P1, const unsigned char* O, uint64_t* P1O_P2, uint64_t* P3){
    const int param_m = M_MAX;
    const int m_legs = param_m / 32;
    const int param_v = V_MAX;
    const int param_o = O_MAX;
    P1_times_O(p, (uint32_t *) P1, O, (uint32_t *) P1O_P2);

    // compute P3 = O^t * (P1*O + P2)
    mul_add_mat_trans_x_m_mat(m_legs, O, P1O_P2, P3, param_v, param_o, param_o);
}
*/

// TODO: remove those generic functions

void m_calculate_SPS(const uint64_t *PS, const unsigned char *S, int m, int k, int  n, uint64_t *SPS) {
    mayo_generic_m_calculate_SPS(PS, S, m, k, n, SPS);
}

void m_calculate_PS(const uint64_t *P1, const uint64_t *P2, const uint64_t *P3, const unsigned char *S,
                    const int m, const int v, const int o, const int k, uint64_t *PS) {
    mayo_generic_m_calculate_PS(P1, P2, P3, S, m, v, o, k, PS);
}



// a > b -> b - a is negative
// returns 0xFFFFFFFF if true, 0x00000000 if false
static inline uint64_t ct_64_is_greater_than(int a, int b) {
    int64_t diff = ((int64_t) b) - ((int64_t) a);
    return (uint64_t) (diff >> (8*sizeof(uint64_t)-1));
}

// if a == b -> 0x00, else 0xFF
static inline unsigned char ct_compare_8(unsigned char a, unsigned char b) {
    return (int8_t)((-(int32_t)(a ^ b)) >> (8*sizeof(uint32_t)-1));
}

// if a == b -> 0x0000000000000000, else 0xFFFFFFFFFFFFFFFF
static inline uint64_t ct_compare_64(int a, int b) {
    return (uint64_t)((-(int64_t)(a ^ b)) >> (8*sizeof(uint64_t)-1));
}

static inline unsigned char
m_extract_element(const uint64_t *in, int index) {
    const int leg = index / 16;
    const int offset = index % 16;

    return (in[leg] >> (offset*4)) & 0xF;
}

static inline void
ef_pack_m_vec(const unsigned char *in, uint64_t *out, int ncols) {
    int i;
    unsigned char *out8 = (unsigned char *)out;
    for(i = 0; i+1 < ncols; i += 2){
        out8[i/2]  = (in[i+0] << 0) | (in[i+1] << 4);
    }
    if (ncols % 2 == 1){
        out8[i/2]  = (in[i+0] << 0);
    }
}

static inline void
ef_unpack_m_vec(int legs, const uint64_t *in, unsigned char *out) {
    const unsigned char *in8 = (const unsigned char *)in;
    for(int i = 0; i < legs * 16; i += 2){
        out[i]   = (in8[i/2]) & 0xF;
        out[i+1] = (in8[i/2] >> 4);
    }
}



// put matrix in row echelon form with ones on first nonzero entries *in
// constant time*
static inline void EF(unsigned char *A, int nrows, int ncols) {

    uint64_t _pivot_row[(K_MAX * O_MAX + 1 + 15) / 16];
    uint64_t _pivot_row2[(K_MAX * O_MAX + 1 + 15) / 16];
    uint64_t packed_A[((K_MAX * O_MAX + 1 + 15) / 16) * M_MAX] = {0};

    int row_len = (ncols + 15) / 16;

    // nibbleslice the matrix A
    for (int i = 0; i < nrows; i++) {
        ef_pack_m_vec(A + i * ncols, packed_A + i * row_len, ncols);
    }

    // pivot row is secret, pivot col is not

    unsigned char inverse;
    int pivot_row = 0;
    for (int pivot_col = 0; pivot_col < ncols; pivot_col++) {

        int pivot_row_lower_bound = MAYO_MAX(0, pivot_col + nrows - ncols);
        int pivot_row_upper_bound = MAYO_MIN(nrows - 1, pivot_col);
        // the pivot row is guaranteed to be between these lower and upper bounds if
        // A has full rank

        // zero out pivot row
        for (int i = 0; i < row_len; i++) {
            _pivot_row[i] = 0;
            _pivot_row2[i] = 0;
        }

        // try to get a pivot row in constant time
        unsigned char pivot = 0;
        uint64_t pivot_is_zero = -1;
        for (int row = pivot_row_lower_bound;
                row <= MAYO_MIN(nrows - 1, pivot_row_upper_bound + 32); row++) {

            uint64_t is_pivot_row = ~ct_compare_64(row, pivot_row);
            uint64_t below_pivot_row = ct_64_is_greater_than(row, pivot_row);

            for (int j = 0; j < row_len; j++) {
                _pivot_row[j] ^= (is_pivot_row | (below_pivot_row & pivot_is_zero)) &
                                 packed_A[row * row_len + j];
            }
            pivot = m_extract_element(_pivot_row, pivot_col);
            pivot_is_zero = ~ct_compare_64((int) pivot, 0);
        }

        // multiply pivot row by inverse of pivot
        inverse = inverse_f(pivot);
        vec_mul_add_u64(row_len, _pivot_row, inverse, _pivot_row2);

        // conditionally write pivot row to the correct row, if there is a nonzero
        // pivot
        for (int row = pivot_row_lower_bound; row <= pivot_row_upper_bound; row++) {
            uint64_t do_copy = ~ct_compare_64(row, pivot_row) & ~pivot_is_zero;
            uint64_t do_not_copy = ~do_copy;
            for (int col = 0; col < row_len; col++) {
                packed_A[row * row_len + col] =
                    (do_not_copy & packed_A[row * row_len + col]) +
                    (do_copy & _pivot_row2[col]);
            }
        }

        // eliminate entries below pivot
        for (int row = pivot_row_lower_bound; row < nrows; row++) {
            unsigned char below_pivot = (row > pivot_row);
            unsigned char elt_to_elim = m_extract_element(packed_A + row * row_len, pivot_col);

            vec_mul_add_u64(row_len, _pivot_row2, below_pivot * elt_to_elim,
                                    packed_A + row * row_len);                            
        }

        pivot_row += (-(int32_t)(~pivot_is_zero));
    }

    unsigned char temp[(O_MAX * K_MAX + 1 + 15)];

    // unbitslice the matrix A
    for (int i = 0; i < nrows; i++) {
        ef_unpack_m_vec(row_len, packed_A + i * row_len, temp);
        for (int j = 0; j < ncols; j++) {
            A[i * ncols + j] = temp[j];
        }
    }
}
inline void vec_mul_add_u64(const int legs, const uint64_t *in, unsigned char a, uint64_t *acc) {
    uint32_t tab = mul_table(a);

    uint64_t lsb_ask = 0x1111111111111111ULL;

    for(int i=0; i < legs; i++){
        acc[i] ^= ( in[i]       & lsb_ask) * (tab & 0xff)
                ^ ((in[i] >> 1) & lsb_ask) * ((tab >> 8)  & 0xf)
                ^ ((in[i] >> 2) & lsb_ask) * ((tab >> 16) & 0xf)
                ^ ((in[i] >> 3) & lsb_ask) * ((tab >> 24) & 0xf);
    }
}


// sample a solution x to Ax = y, with r used as randomness
// require:
// - A is a matrix with m rows and k*o+1 collumns (values in the last collum are
// not important, they will be overwritten by y) in row major order
// - y is a vector with m elements
// - r and x are k*o bytes long
// return: 1 on success, 0 on failure
int sample_solution(const mayo_params_t *p, unsigned char *A,
                           const unsigned char *y, const unsigned char *r,
                           unsigned char *x, int k, int o, int m, int A_cols) {
    #ifdef MAYO_VARIANT
    (void) p;
    #endif

    unsigned char finished;
    int col_upper_bound;
    unsigned char correct_column;

    // x <- r
    for (int i = 0; i < k * o; i++) {
        x[i] = r[i];
    }

    // compute Ar;
    unsigned char Ar[M_MAX];
    for (int i = 0; i < m; i++) {
        A[k * o + i * (k * o + 1)] = 0; // clear last col of A
    }
    mat_mul(A, r, Ar, k * o + 1, m, 1);

    // move y - Ar to last column of matrix A
    for (int i = 0; i < m; i++) {
        A[k * o + i * (k * o + 1)] = sub_f(y[i], Ar[i]);
    }

    EF(A, m, k * o + 1);

    // check if last row of A (excluding the last entry of y) is zero
    unsigned char full_rank = 0;
    for (int i = 0; i < A_cols - 1; i++) {
        full_rank |= A[(m - 1) * A_cols + i];
    }

// It is okay to leak if we need to restart or not
#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(&full_rank, 1);
#endif

    if (full_rank == 0) {
        return 0;
    }

    // back substitution in constant time
    // the index of the first nonzero entry in each row is secret, which makes
    // things less efficient

    for (int row = m - 1; row >= 0; row--) {
        finished = 0;
        col_upper_bound = MAYO_MIN(row + (32/(m-row)), k*o);
        // the first nonzero entry in row r is between r and col_upper_bound with probability at least ~1-q^{-32}

        for (int col = row; col <= col_upper_bound; col++) {

            // Compare two chars in constant time.
            // Returns 0x00 if the byte arrays are equal, 0xff otherwise.
            correct_column = ct_compare_8((A[row * A_cols + col]), 0) & ~finished;

            unsigned char u = correct_column & A[row * A_cols + A_cols - 1];
            x[col] ^= u;

            for (int i = 0; i < row; i += 8) {
                uint64_t tmp = ( (uint64_t) A[ i    * A_cols + col] <<  0) ^ ( (uint64_t) A[(i+1) * A_cols + col] <<  8)
                             ^ ( (uint64_t) A[(i+2) * A_cols + col] << 16) ^ ( (uint64_t) A[(i+3) * A_cols + col] << 24)
                             ^ ( (uint64_t) A[(i+4) * A_cols + col] << 32) ^ ( (uint64_t) A[(i+5) * A_cols + col] << 40)
                             ^ ( (uint64_t) A[(i+6) * A_cols + col] << 48) ^ ( (uint64_t) A[(i+7) * A_cols + col] << 56);

                tmp = mul_fx8(u, tmp);

                A[ i    * A_cols + A_cols - 1] ^= (tmp      ) & 0xf;
                A[(i+1) * A_cols + A_cols - 1] ^= (tmp >> 8 ) & 0xf;
                A[(i+2) * A_cols + A_cols - 1] ^= (tmp >> 16) & 0xf;
                A[(i+3) * A_cols + A_cols - 1] ^= (tmp >> 24) & 0xf;
                A[(i+4) * A_cols + A_cols - 1] ^= (tmp >> 32) & 0xf;
                A[(i+5) * A_cols + A_cols - 1] ^= (tmp >> 40) & 0xf;
                A[(i+6) * A_cols + A_cols - 1] ^= (tmp >> 48) & 0xf;
                A[(i+7) * A_cols + A_cols - 1] ^= (tmp >> 56) & 0xf;
            }

            finished = finished | correct_column;
        }
    }
    return 1;
}