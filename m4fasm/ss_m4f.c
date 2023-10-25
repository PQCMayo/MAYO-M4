
#include <stdint.h>
#include "mayo.h"

#define MAYO_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MAYO_MIN(x, y) (((x) < (y)) ? (x) : (y))

#define PARAM_m(p) PARAM_NAME(m)
#define PARAM_n(p) PARAM_NAME(n)
#define PARAM_o(p) PARAM_NAME(o)
#define PARAM_v(p) PARAM_NAME(v)
#define PARAM_A_cols(p) PARAM_NAME(A_cols)
#define PARAM_k(p) PARAM_NAME(k)
#define PARAM_q(p) PARAM_NAME(q)
#define PARAM_m_bytes(p) PARAM_NAME(m_bytes)
#define PARAM_O_bytes(p) PARAM_NAME(O_bytes)
#define PARAM_v_bytes(p) PARAM_NAME(v_bytes)
#define PARAM_r_bytes(p) PARAM_NAME(r_bytes)
#define PARAM_P1_bytes(p) PARAM_NAME(P1_bytes)
#define PARAM_P2_bytes(p) PARAM_NAME(P2_bytes)
#define PARAM_P3_bytes(p) PARAM_NAME(P3_bytes)
#define PARAM_csk_bytes(p) PARAM_NAME(csk_bytes)
#define PARAM_esk_bytes(p) PARAM_NAME(esk_bytes)
#define PARAM_cpk_bytes(p) PARAM_NAME(cpk_bytes)
#define PARAM_epk_bytes(p) PARAM_NAME(epk_bytes)
#define PARAM_sig_bytes(p) PARAM_NAME(sig_bytes)
// static const unsigned char f_tail[] = PARAM_NAME(f_tail);
#define PARAM_salt_bytes(p) PARAM_NAME(salt_bytes)
#define PARAM_sk_seed_bytes(p) PARAM_NAME(sk_seed_bytes)
#define PARAM_digest_bytes(p) PARAM_NAME(digest_bytes)
#define PARAM_pk_seed_bytes(p) PARAM_NAME(pk_seed_bytes)
#define PARAM_f_tail(p) f_tail


// GF(16) multiplication mod x^4 + x + 1
static inline unsigned char mul_f(unsigned char a, unsigned char b) {
    unsigned char p = 0;
    int mask;
    for (int i = 0; i <= 3; ++i) {
        p ^= -(b & 1) & a;
        mask = -((a >> 3) & 1);
        a = ((a << 1) ^ (0x13 & mask));
        b >>= 1;
    }
    return p;
}

// GF(16) addition
static inline unsigned char add_f(unsigned char a, unsigned char b) {
    return a ^ b;
}

// GF(16) subtraction
static inline unsigned char sub_f(unsigned char a, unsigned char b) {
    return a ^ b;
}

static inline unsigned char lincomb(const unsigned char *a,
                                    const unsigned char *b, int n, int m) {
    unsigned char ret = 0;
    for (int i = 0; i < n; ++i, b += m) {
        ret = add_f(mul_f(a[i], *b), ret);
    }
    return ret;
}

static void mat_mul(const unsigned char *a, const unsigned char *b,
                    unsigned char *c, int colrow_ab, int row_a, int col_b) {
    for (int i = 0; i < row_a; ++i, a += colrow_ab) {
        for (int j = 0; j < col_b; ++j, ++c) {
            *c = lincomb(a, b + j, colrow_ab, col_b);
        }
    }
}


// assembly
void ef_bitslice_asm(uint32_t *a_bs, unsigned char *a);
void ef_unbitslice_asm(unsigned char *a, uint32_t *a_bs);
uint8_t ef_inner1_asm(uint32_t *pivot_row, uint32_t *a_bs_row, int pivot_col, int pivot_row_lower_bound, int pivot_row_ctr, uint8_t *pivot_is_zero, int row_upper_bound);
void    ef_inner2_asm(uint32_t *pivot_row2, uint32_t *pivot_row, uint8_t pivot);
void    ef_inner3_asm(uint32_t *a_bs_row, uint32_t *pivot_row, int pivot_row_ctr, int pivot_row_lower_bound, int pivot_row_upper_bound, int pivot_is_zero);
int     ef_inner4_asm(uint32_t *a_bs, uint32_t *pivot_row, int pivot_row_lower_bound, int pivot_row_ctr, int pivot_col, int pivot_is_nonzero);
void    backsub_inner_asm(uint8_t *a1, uint8_t *a2, uint8_t *a3, int correct_column, int r);


// put matrix in row echelon form with ones on first nonzero entries *in
// constant time*
static void EF(unsigned char *A, int nrows, int ncols) {

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
int sample_solution_m4f(const mayo_params_t *p, unsigned char *A,
                           const unsigned char *y, const unsigned char *r,
                           unsigned char *x) {
    int finished;
    int col_upper_bound;
    int correct_column;
    const int k = PARAM_k(p);
    const int o = PARAM_o(p);
    const int m = PARAM_m(p);
    const int A_cols = PARAM_A_cols(p);

    (void) p;

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


            #if 1
            if(u>0) {
                backsub_inner_asm(A + u * A_cols + A_cols - 1, A + col, A + A_cols - 1, correct_column, u);
            }
            #else
            for (int i = 0; i < u; i++) {
                A[i * A_cols + A_cols - 1] ^=
                    correct_column *
                    mul_f(A[u * A_cols + A_cols - 1], A[i * A_cols + col]);
            }
            #endif

            // for (int i = 0; i < r; i++) {
                // A[i * A_cols + A_cols - 1] ^=
                //     correct_column *
                //     mul_f(A[r * A_cols + A_cols - 1], A[i * A_cols + col]);

            // }


            finished = finished || correct_column;
        }
    }
    return 1;
}


