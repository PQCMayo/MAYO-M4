// SPDX-License-Identifier: Apache-2.0

#include <stdint.h>
#include <stddef.h>
#include <stdalign.h>
#include <string.h>
#include <stdbool.h>

#include "arithmetic.h"
#include "aes_ctr.h"

void compute_P3(uint8_t P3[P3_BYTES],
                const uint8_t seed_pk[PK_SEED_BYTES],
                const uint8_t O[O_BYTES]) {
    // aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    aes128ctr_ctx ctx;
    int row, col;

    aes128ctr_init(&ctx, seed_pk);  // initialize AES-128-CTR with seed_pk
    
    // Here we use the notation from the specification document.
    // So the parameters are in lowercase letters instead of uppercase macros.
    // e.g. n = PARAM_N, o = PARAM_O, m = PARAM_M, etc.

    // Compute P3_h = Upper(-O^T * P1_h * O - O^T * P2_h) minimizing temporary storage
    // The entries of the matrices P1_h (resp. P2_h) are obtained following the order (i, j, h) where
    // - i is the row index
    // - j is the column index
    // - h is the matrix index
    // So the first m entries are the (0, 0) coefficients of all m matrices,
    // then the next m entries are the (0, 1) coefficients of all m matrices, and so on.
    
    // Each P1_h is a (n - o) * (n - o) upper-triangular matrix
    // Each P2_h is a (n - o) * o full matrix
    // Each P3_h is a o * o upper-triangular matrix (transformed to upper-triangular in place)
    // O is a (n - o) * o full matrix
    
    // Compute -O^T * P1_h * O and accumulate into P3_h
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = row; col < PARAM_N - PARAM_O; col++) {
            // compute entry P1_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, &ctx);
            // accumulate into P3_h the value O^T[:, row] * P1_h[row, col] * O[col, :]
            compute_P3_from_P1(P3, entry, O, row, col);
        }
    }
    // Compute -O^T * P2_h and accumulate into P3_h
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = 0; col < PARAM_O; col++) {
            // compute entry P2_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, &ctx);
            // accumulate into P3_h the value O^T[:, row] * P2_h[row, col]
            compute_P3_from_P2(P3, entry, O, row, col);
        }
    }
}

void compute_L(uint8_t L[L_BYTES],
               const uint8_t O[O_BYTES],
               const uint8_t P1[P1_BYTES],
               aes128ctr_ctx *ctx) {
    // aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    int row, col;
    uint8_t *ptr;
    
    memset(L, 0, L_BYTES);  // initialize L to zero

    // Use P1_h to accumulate into L
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        // skip computation for row == col since the diagonal of P1_h + P1_h^T is zero
        for (col = row + 1; col < PARAM_N - PARAM_O; col++) {
            // retrieve entry P1_h[row, col] for all m matrices
            ptr = m_utmatrix_entry_ptr((uint8_t *)P1, row, col, PARAM_N - PARAM_O);
            memcpy(entry, ptr, M_ENTRIES_BYTES);  // no leak since P1 is public
            compute_L_from_P1(L, entry, O, row, col);
        }
    }
    // Use P2_h to accumulate into L
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = 0; col < PARAM_O; col++) {
            // compute entry P2_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, ctx);
            // accumulate into L the value P2_h[row, col]
            ptr = m_matrix_entry_ptr(L, row, col, PARAM_O);
            vec_add_u8(ptr, entry, M_ENTRIES_BYTES);
        }
    }
}

void build_linear_system(uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                         uint8_t y[M_ENTRIES_BYTES],
                         const uint8_t seed_pk[PK_SEED_BYTES],
                         const uint8_t O[O_BYTES],
                         const uint8_t V[PARAM_K * V_BYTES]) {
    // aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    // buf contains either all M matrices (stored transposed) or all u vectors at a given time
    // buf size is the larger of the two
#if PARAM_K * M_ENTRIES_BYTES * PARAM_O > BINOMIAL2(PARAM_K + 1) * M_ENTRIES_BYTES
    uint8_t buf[PARAM_K * M_ENTRIES_BYTES * PARAM_O] = {0};
#else
    uint8_t buf[BINOMIAL2(PARAM_K + 1) * M_ENTRIES_BYTES] = {0};
#endif
    aes128ctr_ctx ctx;
    size_t i;
    int row, col;

    // Here we use the notation from the specification document.
    // So the parameters are in lowercase letters instead of uppercase macros.
    // e.g. n = PARAM_N, o = PARAM_O, m = PARAM_M, etc.

    // Build linear system A and target vector y minimizing temporary storage
    // The entries of the matrices P1_j (resp. P2_j) are obtained as described in compute_P3.
    
    // For all i = 0, ..., k - 1, M_i is a m * o matrix defined as:
    // M_i[j, :] = v_i^T * (P1_j + P1_j^T) * O + v_i^T * P2_j
    // Similarly, for all i = 0, ..., k - 1 and j = k - 1, ..., i, u_i,j is a m-dimensional vector defined as:
    // [ v_i^T * P1_a * v_j ]_{a=0}^{m-1}  if i == j
    // [ v_i^T * P1_a * v_j + v_j^T * P1_a * v_i ]_{a=0}^{m-1}  if i != j
    
    // We now switch the indexing of P1_j (resp. P2_j) to P1_h (resp. P2_h) to avoid confusion.

    aes128ctr_init(&ctx, seed_pk);  // initialize AES-128-CTR with seed_pk
    // Compute v_i^T * (P1_h + P1_h^T) * O + v_i^T * P2_h and accumulate into P3_h
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = row; col < PARAM_N - PARAM_O; col++) {
            // compute entry P1_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, &ctx);
            // compute M and u from P1_h[row, col]
            compute_M_from_P1(
                // explicit cast to avoid compiler warning about pointer alignment
                (uint8_t (*)[M_ENTRIES_BYTES * PARAM_O])buf, 
                entry, O, V, row, col);
        }
    }
    // Compute v_i^T * P2_h and accumulate into M_i
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = 0; col < PARAM_O; col++) {
            // compute entry P2_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, &ctx);
            // compute M from P2_h[row, col]
            compute_M_from_P2(
                // explicit cast to avoid compiler warning about pointer alignment
                (uint8_t (*)[M_ENTRIES_BYTES * PARAM_O])buf,
                entry, V, row, col);
        }
    }
    // Build A from M (stored in buf)
    build_A(A, (const uint8_t (*)[M_ENTRIES_BYTES * PARAM_O])buf);

    for (i = 0; i < sizeof(buf); i++) {  // clear buf
        buf[i] = 0;
    }

    aes128ctr_init(&ctx, seed_pk);  // initialize AES-128-CTR with seed_pk
    // Compute v_i^T * P1_h * v_j and accumulate into u_i
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = row; col < PARAM_N - PARAM_O; col++) {
            // compute entry P1_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, &ctx);
            // compute u from P1_h[row, col]
            compute_u_from_Pi(
                // explicit cast to avoid compiler warning about pointer alignment
                (uint8_t (*)[M_ENTRIES_BYTES])buf, 
                entry, V, PARAM_N - PARAM_O, row, col);
        }
    }
    // Build y from u (stored in buf)
    build_y(y, (const uint8_t (*)[M_ENTRIES_BYTES])buf);

    for (i = 0; i < sizeof(buf); i++) {  // clear buf
        buf[i] = 0;
    }
}

bool sample_solution(uint8_t x[KO_ENTRIES_BYTES],
                     uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                     uint8_t y[M_ENTRIES_BYTES]) {
    int i, j;
    uint8_t full_rank;
    
    // At the start, x is initialized with r random vector
    // y = y - A * r
    for (i = 0; i < PARAM_M; i++) {
        for (j = 0; j < PARAM_K * PARAM_O; j++) {
            // y[i] -= A[i][j] * x[j]
            add_nibble(y, i, mul_f(nibble(A[i], j), nibble(x, j)));
        }
    }
    
    echelon_form(A, y);

    // check if last row of A is zero
    full_rank = 0;
    for (i = 0; i < KO_ENTRIES_BYTES; i++) {
        full_rank |= A[PARAM_M - 1][i];
    }
    if (full_rank == 0) {
        return false;  // A is not full rank
    }

    backward_substitution(A, y, x);

    return true;  // A is full rank and solution has been found
}

void compute_s(uint8_t s[SIG_BYTES],
               const uint8_t x[KO_ENTRIES_BYTES],
               const uint8_t O[O_BYTES],
               const uint8_t V[PARAM_K * V_BYTES]) {
    int i, row, col;
    uint8_t e0, e1;

    for (i = 0; i < SIG_BYTES; i++) {  // clear s for accumulation
        s[i] = 0;
    }

    for (i = 0; i < PARAM_K; i++) {
        // handle vinegar part
        for (row = 0; row < PARAM_N - PARAM_O; row++) {
            e0 = matrix_entry_val(V, i, row, PARAM_N - PARAM_O);  // v_i[row]
            add_nibble(s, i * PARAM_N + row, e0);
        }
        // handle oil part and x concatenation
        for (col = 0; col < PARAM_O; col++) {
            e0 = nibble(x, i * PARAM_O + col); // x[i * PARAM_O + col]
            add_nibble(s, i * PARAM_N + PARAM_N - PARAM_O + col, e0);
            for (row = 0; row < PARAM_N - PARAM_O; row++) {
                e1 = matrix_entry_val(O, row, col, PARAM_O);  // O[row, col]
                add_nibble(s, i * PARAM_N + row, mul_f(e0, e1));
            }
        }
    }
}

void quadratic_map(uint8_t y[M_ENTRIES_BYTES],
                   const uint8_t s[SIG_BYTES],
                   const uint8_t seed_pk[PK_SEED_BYTES],
                   const uint8_t P3[P3_BYTES]) {
    alignas(8) uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    uint8_t u[BINOMIAL2(PARAM_K + 1)][M_ENTRIES_BYTES] = {0};
    aes128ctr_ctx ctx;
    int row, col;
    uint8_t *ptr;

    aes128ctr_init(&ctx, seed_pk);  // initialize AES-128-CTR with seed_pk

    // Compute P1_h and accumulate into u
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = row; col < PARAM_N - PARAM_O; col++) {
            // compute entry P1_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, &ctx);
            // accumulate into u the value s_:[row] * P1_h[row, col] * s_:[col]
            compute_u_from_Pi(u, entry, s, PARAM_N, row, col);
        }
    }
    // Compute P2_h and accumulate into u
    for (row = 0; row < PARAM_N - PARAM_O; row++) {
        for (col = 0; col < PARAM_O; col++) {
            // compute entry P2_h[row, col] for all m matrices
            aes128ctr_stream(entry, M_ENTRIES_BYTES, &ctx);
            // accumulate into u the value s_:[row] * P2_h[row, col] * s_:[PARAM_N - PARAM_O + col]
            compute_u_from_Pi(u, entry, s, PARAM_N, row, PARAM_N - PARAM_O + col);
        }
    }
    // Use P3_h to accumulate into u
    for (row = 0; row < PARAM_O; row++) {
        for (col = row; col < PARAM_O; col++) {
            // entry P3[row, col]
            ptr = m_utmatrix_entry_ptr((uint8_t *)P3, row, col, PARAM_O);
            memcpy(entry, ptr, M_ENTRIES_BYTES);
            // accumulate into u the value s_:[PARAM_N - PARAM_O + row] * P3_h[row, col] * s_:[PARAM_N - PARAM_O + col]
            compute_u_from_Pi(u, entry, s, PARAM_N, PARAM_N - PARAM_O + row, PARAM_N - PARAM_O + col);
        }
    }

    // Build y from u
    // clear y since build_y expects it to be initialized to t, and we want t = 0 for verification
    memset(y, 0, M_ENTRIES_BYTES);
    build_y(y, (const uint8_t (*)[M_ENTRIES_BYTES])u);
}
