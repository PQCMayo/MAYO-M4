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


void m_upper(int m_legs, const uint64_t *in, uint64_t *out, int size) {
    int m_vecs_stored = 0;
    for (int r = 0; r < size; r++) {
        for (int c = r; c < size; c++) {
            m_vec_copy(m_legs, in + m_legs * 2 * (r * size + c), out + m_legs * 2 * m_vecs_stored );
            if (r != c) {
                m_vec_add(m_legs, in + m_legs * 2 * (c * size + r), out + m_legs * 2 * m_vecs_stored );
            }
            m_vecs_stored ++;
        }
    }
}



// multiplies the transpose of a single matrix with m matrices and adds result to acc
static void mul_add_mat_trans_x_m_mat(int m_legs, const unsigned char *mat, const uint64_t *bs_mat, uint64_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {
    (void) m_legs;
    (void) mat_rows;
    (void) mat_cols;
    (void) bs_mat_cols;
    mul_add_mat_trans_x_m_mat_m4f_V_O_O_asm(acc, mat, bs_mat);

}


// -----------------------------------------------------------------------------------------
// OPTIMIZED ARITHMETIC
// -----------------------------------------------------------------------------------------

// multiplies a single matrix with m matrices and adds result to acc
void mul_add_mat_x_m_mat(int m_legs, const unsigned char *mat, const uint32_t *bs_mat, uint32_t *acc, int mat_rows, int mat_cols, int bs_mat_cols) {
  (void) m_legs;
  if(mat_rows == K_MAX && mat_cols == V_MAX && bs_mat_cols == O_MAX){
    mul_add_mat_x_m_mat_m4f_K_V_O_triangular_asm(acc, mat, bs_mat);
  } else if(mat_rows == K_MAX && mat_cols == V_MAX && bs_mat_cols == K_MAX){
    mul_add_mat_x_m_mat_m4f_K_V_K_triangular_asm(acc, mat, bs_mat);
  }
}


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





// put matrix in row echelon form with ones on first nonzero entries *in
// constant time*
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




// compute P * S^t = [ P1  P2 ] * [S1] = [P1*S1 + P2*S2]
//                   [  0  P3 ]   [S2]   [        P3*S2]
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

    /* Old approach which is constant time but doesn't have to be
    unsigned char S1[V_MAX*K_MAX];
    unsigned char S2[O_MAX*K_MAX];
    unsigned char *s1_write = S1;
    unsigned char *s2_write = S2;
    for (int r=0; r < k; r++)
    {
        for (int c = 0; c < n; c++)
        {
            if(c < v){
                *(s1_write++) = S[r*n + c];
            } else {
                *(s2_write++) = S[r*n + c];
            }
        }
    }

    mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(m_legs, P1, S1, PS, v, v, k, 1); // P1 * S1
    mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(m_legs, P2, S2, PS, v, o, k, 0); // P2 * S2
    mul_add_bitsliced_m_upper_triangular_mat_x_mat_trans(m_legs, P3, S2, PS + v*k*m_legs*4, o, o, k, 1); // P3 * S2
    */

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


void Ot_times_P1O_P2(const mayo_params_t* p, const uint64_t* P1, const unsigned char* O, uint64_t* P1O_P2, uint64_t* P3){
    const int param_m = M_MAX;
    const int m_legs = param_m / 32;
    const int param_v = V_MAX;
    const int param_o = O_MAX;
    P1_times_O(p, (uint32_t *) P1, O, (uint32_t *) P1O_P2);

    // compute P3 = O^t * (P1*O + P2)
    mul_add_mat_trans_x_m_mat(m_legs, O, P1O_P2, P3, param_v, param_o, param_o);
}

