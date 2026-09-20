// SPDX-License-Identifier: Apache-2.0

/*
Copyright 2023-2025 the MAYO team. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*
This file contains some functions for arithmetic operations over GF(16)
taken from the MAYO reference implementation. 
These are: 
- mul_f;
- mul_table;
- vec_scalarmul_u64, modified from vec_mul_add_u64;
- echelon_form, modified from EF, but the constant-time idea is still the same;
- back_substitution, modified from the last part of sample_solution, similarly to echelon_form;
*/

#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include <stdint.h>
#include <stddef.h>
#include <stdalign.h>
#include <string.h>
#include <stdbool.h>

#include "params.h"
#include "aes_ctr.h"
#include "blocker.h"
#include "comparison.h"

/**
 * Computes P3 part of the public key from seed_pk and oil_subspace.
 * 
 * @param[out] P3 Pointer to the output buffer for P3 of size P3_BYTES
 * @param[in] seed_pk Pointer to the seed_pk of size PK_SEED_BYTES
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 */
#define compute_P3 MAYO_NAMESPACE(compute_P3)
void compute_P3(uint8_t P3[P3_BYTES],
                const uint8_t seed_pk[PK_SEED_BYTES],
                const uint8_t O[O_BYTES]);

/**
 * Computes L from O, P1 and P2 (generated on-the-fly with AES-128-CTR).
 * 
 * @param[out] L Output buffer to store L of size PARAM_M * (PARAM_N - PARAM_O) * PARAM_O / 2 bytes
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 * @param[in] P1 Pointer to the P1 part of the public key of size P1_BYTES
 * @param[in] ctx Pointer to an initialized aes128ctr_ctx for generating P2
 */
#define compute_L MAYO_NAMESPACE(compute_L)
void compute_L(uint8_t L[(PARAM_N - PARAM_O) * PARAM_O * M_ENTRIES_BYTES],
               const uint8_t O[O_BYTES],
               const uint8_t P1[P1_BYTES],
               aes128ctr_ctx *ctx);

/**
 * Builds the linear system A and target vector y for signing.
 * 
 * @param[out] A Pointer to the linear system matrix (assumed to be zero-initialized)
 * @param[out] y Pointer to the target vector (assumed to be initialized to t)
 * @param[in] seed_pk Pointer to the seed_pk of size PK_SEED_BYTES
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 * @param[in] V Pointer to the vinegar variables of size PARAM_K * V_BYTES
 */
#define build_linear_system MAYO_NAMESPACE(build_linear_system)
void build_linear_system(uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                         uint8_t y[M_ENTRIES_BYTES],
                         const uint8_t seed_pk[PK_SEED_BYTES],
                         const uint8_t O[O_BYTES],
                         const uint8_t V[PARAM_K * V_BYTES]);

/**
 * Samples a solution x to the linear system Ax = y.
 * 
 * @param[out] x Pointer to the output solution vector of size KO_ENTRIES_BYTES
 * initially set to the random vector r
 * @param[in] A Pointer to the linear system matrix
 * @param[in] y Pointer to the target vector
 * @return true if a valid solution was found, false otherwise
 */
#define sample_solution MAYO_NAMESPACE(sample_solution)
bool sample_solution(uint8_t x[KO_ENTRIES_BYTES],
                     uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                     uint8_t y[M_ENTRIES_BYTES]);

/**
 * Computes the signature from the solution x, oil subspace O and vinegar variables V.
 * 
 * @param[out] s Pointer to the output signature of size SIG_BYTES
 * @param[in] x Pointer to the solution vector of size KO_ENTRIES_BYTES
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 * @param[in] V Pointer to the vinegar variables of size PARAM_K * V_BYTES
 */
#define compute_s MAYO_NAMESPACE(compute_s)
void compute_s(uint8_t s[SIG_BYTES],
               const uint8_t x[KO_ENTRIES_BYTES],
               const uint8_t O[O_BYTES],
               const uint8_t V[PARAM_K * V_BYTES]);

/**
 * Evaluates the quadratic map to compute y from s, seed_pk, and P3,
 * i.e. computes y = P(s) where P is the quadratic map defined by P1, P2 and P3.
 * P1 and P2 are generated on-the-fly with AES-128-CTR using seed_pk,
 * while P3 is given as input since it is part of the public key.
 * 
 * @param[out] y Pointer to the output vector of size M_ENTRIES_BYTES to store the result
 * @param[in] s Pointer to the input signature of size SIG_BYTES
 * @param[in] seed_pk Pointer to the seed_pk of size PK_SEED_BYTES
 * @param[in] P3 Pointer to the P3 part of the public key of size P3_BYTES
 */
#define quadratic_map MAYO_NAMESPACE(quadratic_map)
void quadratic_map(uint8_t y[M_ENTRIES_BYTES],
                   const uint8_t s[SIG_BYTES],
                   const uint8_t seed_pk[PK_SEED_BYTES],
                   const uint8_t P3[P3_BYTES]);


/**************************************************************************************
 *   Implementation of arithmetic operations in GF(16) as static inline functions     *
 **************************************************************************************/

/**
 * Carryless multiplication in GF(2^4) with modulus x^4 + x + 1
 * 
 * @param[in] a First operand (expected to be in 0..15)
 * @param[in] b Second operand (expected to be in 0..15)
 * @return uint8_t Result of the multiplication
 */
static inline uint8_t mul_f(uint8_t a, uint8_t b) {
    // carryless multiply
    uint8_t p;

#if !(((defined(__clang__) && __clang_major__ < 15) || (!defined(__clang__) && defined(__GNUC__) && __GNUC__ <= 12)) && (defined(__x86_64__) || defined(_M_X64)))
    a ^= uint8_t_blocker;  // prevent GCC/Clang from optimizing this code, see https://github.com/PQCMayo/MAYO-C/pull/8/changes/3dace0e22d376451fac99899cbb3252d05712995#diff-c84e4f0b977868b489bda2319658059b92f056ac5793aceba06712d2dcd13c7f
#endif

    p  = (a & 1)*b;
    p ^= (a & 2)*b;
    p ^= (a & 4)*b;
    p ^= (a & 8)*b;

    // reduce mod x^4 + x + 1
    uint8_t top_p = p & 0xf0;
    uint8_t out = (p ^ (top_p >> 4) ^ (top_p >> 3)) & 0x0f;
    return out;
}

/**
 * Computes the multiplicative inverse in GF(16).
 * 
 * @param[in] a Input element (expected to be in 0..15)
 * @return uint8_t Multiplicative inverse of a in GF(16)
 */
static inline uint8_t inverse_f(uint8_t a) {
    // multiplicative inverse in GF(16) with modulus x^4 + x + 1
    // using exponentiation by squaring: a^(2^4 - 2) = a^14
    uint8_t a2 = mul_f(a, a);         // a^2
    uint8_t a4 = mul_f(a2, a2);       // a^4
    uint8_t a8 = mul_f(a4, a4);       // a^8
    uint8_t a12 = mul_f(a8, a4);      // a^12
    uint8_t a14 = mul_f(a12, a2);     // a^14
    return a14;
}

/**
 * Returns the nibble at position pos in vec
 * 
 * @param[in] vec Pointer to a vector of nibbles packed in uint8_t
 * @param[in] pos Position of the nibble
 * @return uint8_t The nibble at the specified position
 */
static inline uint8_t nibble(const uint8_t *vec, const int pos) {
    int byte_pos = pos >> 1;
    int nibble_pos = pos & 1;
    return (vec[byte_pos] >> (4 * nibble_pos)) & 0x0f;
}

/**
 * Sets the nibble at position pos in vec to val
 * 
 * @param[out] vec Pointer to a vector of nibbles packed in uint8_t
 * @param[in] pos Position of the nibble to set
 * @param[in] val Value to set the nibble to (expected to be in 0..15)
 */
static inline void set_nibble(uint8_t *vec, const int pos, const uint8_t val) {
    int byte_pos = pos >> 1;
    int nibble_pos = pos & 1;
    vec[byte_pos] = (vec[byte_pos] & ~(0x0f << (4 * nibble_pos))) | ((val & 0x0f) << (4 * nibble_pos));
}

/**
 * Adds the nibble at position pos in vec with val
 * 
 * @param[out] vec Pointer to a vector of nibbles packed in uint8_t
 * @param[in] pos Position of the nibble to add
 * @param[in] val Value to add to the nibble (expected to be in 0..15)
 */
static inline void add_nibble(uint8_t *vec, const int pos, const uint8_t val) {
    int byte_pos = pos >> 1;
    int nibble_pos = pos & 1;
    vec[byte_pos] ^= (val & 0x0f) << (4 * nibble_pos);
}

/**
 * Retrieves an entry from a GF(16) matrix stored in a packed format.
 * 
 * @param[in] matrix Pointer to the matrix data (each entry is 4 bits, two entries per byte)
 * @param[in] row Row index of the entry to retrieve
 * @param[in] col Column index of the entry to retrieve
 * @param[in] ncols Number of columns in the matrix
 * @return uint8_t The retrieved matrix entry (in 0..15)
 */
static inline uint8_t matrix_entry_val(const uint8_t *matrix,
                                       const int row, const int col,
                                       const int ncols) {
    return nibble(matrix, row * ncols + col);
}

/**
 * Retrieves a pointer to an entry of PARAM_M GF(16) elements of an full matrix stored in a packed format.
 * 
 * @param[in] matrix Pointer to the matrix data (each entry is M_ENTRIES_BYTES bytes)
 * @param[in] row Row index of the entry
 * @param[in] col Column index of the entry
 * @param[in] ncols Number of columns in the full matrix
 * @return uint8_t* Pointer to the byte containing the requested matrix entry
 */
static inline uint8_t* m_matrix_entry_ptr(uint8_t *matrix,
                                          const int row, const int col,
                                          const int ncols) {
    return matrix + M_ENTRIES_BYTES * (row * ncols + col);
}

/**
 * Retrieves a pointer to an entry of PARAM_M GF(16) elements of an upper-triangular matrix stored in a packed format.
 * 
 * @param[in] matrix Pointer to the matrix data (each entry is M_ENTRIES_BYTES bytes)
 * @param[in] row Row index of the entry
 * @param[in] col Column index of the entry
 * @param[in] dim Number of rows/columns in the upper-triangular matrix
 * @return uint8_t* Pointer to the byte containing the requested matrix entry
 */
static inline uint8_t* m_utmatrix_entry_ptr(uint8_t *matrix,
                                            const int row, const int col,
                                            const int dim) {
    return matrix + M_ENTRIES_BYTES * (row * dim + col - row * (row + 1) / 2);
}

/**
 * Precomputed multiplication table for GF(16) elements.
 * Each byte contains four 4-bit results for multiplying by 1, 2, 4, and 8,
 * i.e. by x^0, x^1, x^2, and x^3 in GF(16).
 * 
 * @param[in] b GF(16) element to multiply with (expected to be in 0..15)
 * @return uint32_t Packed multiplication results
 */
static inline uint32_t mul_table(uint8_t b) {
    uint32_t x = ((uint32_t) b) * 0x08040201;
    uint32_t high_half = x & 0xf0f0f0f0;
    return (x ^ (high_half >> 4) ^ (high_half >> 3));
}

/**
 * Adds two vectors of GF(16) elements packed in uint8_t.
 * 
 * @param[out] out Output vector to store the result
 * @param[in] in Input vector to add
 * @param[in] len Length of the vectors in bytes
 */
static inline void vec_add_u8(uint8_t *out, const uint8_t *in, const size_t len) {
    for (size_t i = 0; i < len; i++) {
        out[i] ^= in[i];
    }
}

/**
 * Scalarly multiplies a vector of GF(16) elements packed in uint64_t with a single GF(16) element.
 * 
 * @param[out] out Output vector to store the result (assumed to be aligned to 8 bytes)
 * @param[in] outlen Length of the output vector in uint64_t elements
 * @param[in] in Input vector to multiply of length at least outlen (assumed to be aligned to 8 bytes)
 * @param[in] a GF(16) scalar to multiply with (expected to be in 0..15)
 */
static inline void vec_scalarmul_u64(uint8_t *out, const int outlen, 
                                     const uint8_t *in, 
                                     const uint8_t a) {
    uint32_t tab = mul_table(a);
    uint64_t lsb_mask = 0x1111111111111111ULL;
    uint64_t *in_u64 = (uint64_t *)in;
    uint64_t *out_u64 = (uint64_t *)out;

    for(int i = 0; i < outlen; i++) {
        out_u64[i] = ( in_u64[i]       & lsb_mask) * (tab & 0xf)
                   ^ ((in_u64[i] >> 1) & lsb_mask) * ((tab >> 8)  & 0xf)
                   ^ ((in_u64[i] >> 2) & lsb_mask) * ((tab >> 16) & 0xf)
                   ^ ((in_u64[i] >> 3) & lsb_mask) * ((tab >> 24) & 0xf);
    }
}

/**************************************************************************************
 *                      Helper functions for computing P3                             *
 **************************************************************************************/

/**
 * Computes the contribution to P3 from an entry of P1, i.e. O^T[:, row] * P1[row, col] * O[col, :].
 * 
 * @param[out] P3 Pointer to the output buffer for P3 of size P3_BYTES
 * @param[in] entry Pointer to the input entry of P1, which is a vector of PARAM_M GF(16) elements 
 * packed in uint64_t (assumed to be aligned to 8 bytes)
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 * @param[in] row Row index of the P1 entry
 * @param[in] col Column index of the P1 entry
 */
static inline void compute_P3_from_P1(uint8_t P3[P3_BYTES],
                                      const uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8],
                                      const uint8_t O[O_BYTES],
                                      const uint8_t row, const uint8_t col) {
    // temporary buffers aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t tmp0[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    alignas(8) uint8_t tmp1[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    int i, j;
    uint8_t *ptr;
    uint8_t coeff;

    for (i = 0; i < PARAM_O; i++) {
        coeff = matrix_entry_val(O, row, i, PARAM_O);  // O^T[i, row]
        // compute O^T[i, row] * P1[row, col]
        vec_scalarmul_u64(tmp0, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);

        for (j = 0; j < PARAM_O; j++) {
            coeff = matrix_entry_val(O, col, j, PARAM_O);  // O[col, j]
            // compute O^T[i, row] * P1[row, col] * O[col, j]
            vec_scalarmul_u64(tmp1, TO_QWORDS(M_ENTRIES_BYTES), tmp0, coeff);
            // accumulate into P3[i, j] the value O^T[i, row] * P1[row, col] * O[col, j]
            ptr = m_utmatrix_entry_ptr(P3, (i <= j) ? i : j, (i <= j) ? j : i, PARAM_O);
            vec_add_u8(ptr, tmp1, M_ENTRIES_BYTES);
        }
    }
}

/**
 * Computes the contribution to P3 from an entry of P2, i.e. O^T[:, row] * P2[row, col].
 * 
 * @param[out] P3 Pointer to the output buffer for P3 of size P3_BYTES
 * @param[in] entry Pointer to the input entry of P2, which is a vector of PARAM_M GF(16) elements 
 * packed in uint64_t (assumed to be aligned to 8 bytes)
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 * @param[in] row Row index of the P2 entry
 * @param[in] col Column index of the P2 entry
 */
static inline void compute_P3_from_P2(uint8_t P3[P3_BYTES],
                                      const uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8],
                                      const uint8_t O[O_BYTES],
                                      const uint8_t row, const uint8_t col) {
    // temporary buffer aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t tmp[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    int i;
    uint8_t *ptr;
    uint8_t coeff;

    for (i = 0; i < PARAM_O; i++) {
        coeff = matrix_entry_val(O, row, i, PARAM_O);  // O^T[i, row]
        // compute O^T[i, row] * P2[row, col]
        vec_scalarmul_u64(tmp, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);
        // accumulate into P3[i, j] the value O^T[i, row] * P2[row, col]
        ptr = m_utmatrix_entry_ptr(P3, (i <= col) ? i : col, (i <= col) ? col : i, PARAM_O);
        vec_add_u8(ptr, tmp, M_ENTRIES_BYTES);
    }
}

/**************************************************************************************
 *                 Helper functions for computing L, M and u                          *
 **************************************************************************************/

/**
 * Computes the contribution to L from an entry of P1, 
 * i.e. P1[row, col] * O[col, :] and P1^T[col, row] * O[row, :].
 * 
 * @param[out] L Output buffer to store L of size L_BYTES
 * @param[in] entry Pointer to the input entry of P1, which is a vector of PARAM_M GF(16) elements
 * packed in uint64_t (assumed to be aligned to 8 bytes)
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 * @param[in] row Row index of the P1 entry
 * @param[in] col Column index of the P1 entry
 */
static inline void compute_L_from_P1(uint8_t L[L_BYTES],
                                     const uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8],
                                     const uint8_t O[O_BYTES],
                                     const uint8_t row, const uint8_t col) {
    // temporary buffer aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t tmp[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    int i;
    uint8_t *ptr;
    uint8_t coeff;

    if (row == col) {
        return;  // no contribution to L from the diagonal entries of P1
    }
    for (i = 0; i < PARAM_O; i++) {
        coeff = matrix_entry_val(O, col, i, PARAM_O);  // O[col, i]
        // compute P1[row, col] * O[col, i]
        vec_scalarmul_u64(tmp, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);
        // accumulate into L[row, i] the value P1[row, col] * O[col, i]
        ptr = m_matrix_entry_ptr(L, row, i, PARAM_O);
        vec_add_u8(ptr, tmp, M_ENTRIES_BYTES);

        coeff = matrix_entry_val(O, row, i, PARAM_O);  // O[row, i]
        // compute P1^T[col, row] * O[row, i]
        vec_scalarmul_u64(tmp, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);
        // accumulate into L[col, i] the value P1^T[col, row] * O[row, i]
        ptr = m_matrix_entry_ptr(L, col, i, PARAM_O);
        vec_add_u8(ptr, tmp, M_ENTRIES_BYTES);
    }
}

#define U_INDEX_AUX(i, j) ((i) * PARAM_K - (i) * ((i) - 1) / 2 + PARAM_K - (j) - 1)
// Determine the ell value in the signature/verification loop from i and j
#define U_INDEX(i, j) (U_INDEX_AUX(((i) <= (j) ? (i) : (j)), ((i) <= (j) ? (j) : (i))))

/**
 * Computes the contribution to M from a given entry of the matrix P1, 
 * i.e. v_i[row] * P1[row, col] * O[col, :] and v_i[col] * P1^T[col, row] * O[row, :].
 * 
 * @param[out] M Output matrix M to store the result
 * @param[in] entry Input entry of the matrix P1 (PARAM_M GF(16) elements packed in uint64_t)
 * @param[in] O Pointer to the oil subspace of size O_BYTES
 * @param[in] V Pointer to the vinegar variables of size M * V_BYTES
 * @param[in] row Row index of the entry
 * @param[in] col Column index of the entry
 */
static inline void compute_M_from_P1(uint8_t M[PARAM_K][M_ENTRIES_BYTES * PARAM_O],
                                     const uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8],
                                     const uint8_t O[O_BYTES], const uint8_t V[PARAM_K * V_BYTES],
                                     const uint8_t row, const uint8_t col) {
    // temporary buffers aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t tmp0[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    alignas(8) uint8_t tmp1[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    int i, j;
    uint8_t coeff;
    uint8_t *ptr;

    if (row == col) {
        return;  // skip computation since the diagonal of P1 + P1^T is zero
    }
    for (i = 0; i < PARAM_K; i++) {
        // Handle v_i^T * P1 * O
        coeff = matrix_entry_val(V, i, row, PARAM_N - PARAM_O);  // v_i[row]
        // compute v_i[row] * P1[row, col]
        vec_scalarmul_u64(tmp0, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);
        for (j = 0; j < PARAM_O; j++) {
            coeff = matrix_entry_val(O, col, j, PARAM_O);  // O[col, j]
            vec_scalarmul_u64(tmp1, TO_QWORDS(M_ENTRIES_BYTES), tmp0, coeff);
            // accumulate into M_i the value v_i[row] * P1[row, col] * O[col, j]
            ptr = m_matrix_entry_ptr((uint8_t *)M, i, j, PARAM_O);
            vec_add_u8(ptr, tmp1, M_ENTRIES_BYTES);
        }
        
        // Handle v_i^T * P1^T * O
        coeff = matrix_entry_val(V, i, col, PARAM_N - PARAM_O);  // v_i[col]
        // compute v_i[col] * P1^T[col, row]
        vec_scalarmul_u64(tmp0, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);
        for (j = 0; j < PARAM_O; j++) {
            coeff = matrix_entry_val(O, row, j, PARAM_O);  // O[row, j]
            vec_scalarmul_u64(tmp1, TO_QWORDS(M_ENTRIES_BYTES), tmp0, coeff);
            // accumulate into M_i the value v_i[col] * P1^T[col, row] * O[row, j]
            ptr = m_matrix_entry_ptr((uint8_t *)M, i, j, PARAM_O);
            vec_add_u8(ptr, tmp1, M_ENTRIES_BYTES);
        }
    }
}

/**
 * Computes M from a given entry of the matrix P2.
 * 
 * @param[out] M Output matrix M to store the result
 * @param[in] entry Input entry of the matrix P2 (PARAM_M GF(16) elements packed in uint64_t)
 * @param[in] V Pointer to the vinegar variables of size M * V_BYTES
 * @param[in] row Row index of the entry
 * @param[in] col Column index of the entry
 */
static inline void compute_M_from_P2(uint8_t M[PARAM_K][M_ENTRIES_BYTES * PARAM_O],
                                     const uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8],
                                     const uint8_t V[PARAM_K * V_BYTES],
                                     const uint8_t row, const uint8_t col) {
    // temporary buffer aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t tmp[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    int i;
    uint8_t coeff;
    uint8_t *ptr;

    for (i = 0; i < PARAM_K; i++) {
        coeff = matrix_entry_val(V, i, row, PARAM_N - PARAM_O);  // v_i[row]
        // compute v_i[row] * P2[row, col]
        vec_scalarmul_u64(tmp, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);
        // see M as a matrix with entries of M_ENTRIES_BYTES bytes
        ptr = m_matrix_entry_ptr((uint8_t *)M, i, col, PARAM_O);
        // accumulate tmp into M_i[:, col]
        vec_add_u8(ptr, tmp, M_ENTRIES_BYTES);
    }
}

/**
 * Computes u from a given entry of the matrix Pk, i.e. P1, P2 or P3.
 * 
 * @param[out] u Output vector u to store the result
 * @param[in] entry Input entry of the matrix Pk (PARAM_M GF(16) elements packed in uint64_t)
 * @param[in] W Pointer to a matrix with PARAM_K columns, and usually 
 * PARAM_N - PARAM_O rows if they are vinegar variables, or
 * PARAM_N rows if they are signature variables
 * @param[in] W_ncols Number of columns in W
 * @param[in] row Row index of the entry
 * @param[in] col Column index of the entry
 */
static inline void compute_u_from_Pi(uint8_t u[BINOMIAL2(PARAM_K + 1)][M_ENTRIES_BYTES],
                                     const uint8_t entry[TO_QWORDS(M_ENTRIES_BYTES) * 8],
                                     const uint8_t *W, const int W_ncols,
                                     const uint8_t row, const uint8_t col) {
    // temporary buffers aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t tmp0[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    alignas(8) uint8_t tmp1[TO_QWORDS(M_ENTRIES_BYTES) * 8];
    int i, j, ell;
    uint8_t coeff;

    for (i = 0; i < PARAM_K; i++) {
        coeff = matrix_entry_val(W, i, row, W_ncols);  // v_i[row]
        // compute v_i[row] * Pk[row, col]
        vec_scalarmul_u64(tmp0, TO_QWORDS(M_ENTRIES_BYTES), entry, coeff);
        // accumulate into u[ell] the value v_i^T * Pk * v_j
        for (j = 0; j < PARAM_K; j++) {
            coeff = matrix_entry_val(W, j, col, W_ncols);  // v_j[col]
            // compute v_i[row] * Pk[row, col] * v_j[col]
            vec_scalarmul_u64(tmp1, TO_QWORDS(M_ENTRIES_BYTES), tmp0, coeff);
            // accumulate into u[ell]
            ell = U_INDEX(i, j);
            vec_add_u8(u[ell], tmp1, M_ENTRIES_BYTES);
        }
    }
}

/**************************************************************************************
 *           Helper functions for building the linear system from M and u             *
 **************************************************************************************/

/**
 * Multiplies A by z, where the columns of A are seen as polynomials mod f(z).
 * 
 * @param[in,out] A Input/output matrix A
 */
static inline void A_mul_z(uint8_t A[PARAM_M][KO_ENTRIES_BYTES]) {
    int row, col;
    // temporary buffers aligned to 8 bytes as required by vec_scalarmul_u64
    alignas(8) uint8_t tmp0[TO_QWORDS(KO_ENTRIES_BYTES) * 8];
    alignas(8) uint8_t tmp1[TO_QWORDS(KO_ENTRIES_BYTES) * 8];
    // shift A up by one
    for (col = 0; col < KO_ENTRIES_BYTES; col++) {  // copy last row of A to tmp0
        tmp0[col] = A[PARAM_M - 1][col];
    }
    for (row = PARAM_M - 1; row > 0; row--) {  // shift A up by one
        for (col = 0; col < KO_ENTRIES_BYTES; col++) {
            A[row][col] = A[row - 1][col];
        }
    }
    for (col = 0; col < KO_ENTRIES_BYTES; col++) {  // clear first row of A
        A[0][col] = 0;
    }
    // add the contribution of the modulus
    vec_scalarmul_u64(tmp1, TO_QWORDS(KO_ENTRIES_BYTES), tmp0, F_TAIL0);
    vec_add_u8(A[0], tmp1, KO_ENTRIES_BYTES);
    vec_scalarmul_u64(tmp1, TO_QWORDS(KO_ENTRIES_BYTES), tmp0, F_TAIL1);
    vec_add_u8(A[1], tmp1, KO_ENTRIES_BYTES);
    vec_scalarmul_u64(tmp1, TO_QWORDS(KO_ENTRIES_BYTES), tmp0, F_TAIL2);
    vec_add_u8(A[2], tmp1, KO_ENTRIES_BYTES);
    vec_scalarmul_u64(tmp1, TO_QWORDS(KO_ENTRIES_BYTES), tmp0, F_TAIL3);
    vec_add_u8(A[3], tmp1, KO_ENTRIES_BYTES);
}

/**
 * Accumulates the contributions of M_i and M_j into the matrix A.
 * 
 * @param[out] A Output matrix A
 * @param[in] M Input matrices M
 * @param[in] i Index i of M
 * @param[in] j Index j of M
 */
static inline void A_add_M(uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                           const uint8_t M[PARAM_K][M_ENTRIES_BYTES * PARAM_O],
                           const int i, const int j) {
    int row, col;
    uint8_t coeff;

    // accumulate M_j into A
    for (row = 0; row < PARAM_M; row++) {
        for (col = 0; col < PARAM_O; col++) {
            // M_j[row, col] (M is stored transposed)
            coeff = matrix_entry_val(M[j], col, row, PARAM_M);
            add_nibble(A[row], i * PARAM_O + col, coeff);
        }
    }
    if (i == j) {
        return;  // skip accumulation of M_i since it would be added twice otherwise
    }
    // accumulate M_i into A
    for (row = 0; row < PARAM_M; row++) {
        for (col = 0; col < PARAM_O; col++) {
            // M_i[row, col] (M is stored transposed)
            coeff = matrix_entry_val(M[i], col, row, PARAM_M);
            add_nibble(A[row], j * PARAM_O + col, coeff);
        }
    }
}

/**
 * Builds the matrix A from the matrices M.
 * 
 * @param[out] A Output matrix A (assumed to be zero-initialized)
 * @param[in] M Input matrices M
 */
static inline void build_A(uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                           const uint8_t M[PARAM_K][M_ENTRIES_BYTES * PARAM_O]) {
    int i, j;
    
    // Build A from maximum ell to minimum ell
    A_add_M(A, M, PARAM_K - 1, PARAM_K - 1);  // maximum ell
    // for all ell from second maximum to minimum
    for (i = PARAM_K - 2; i >= 0; i--) {
        for (j = i; j < PARAM_K; j++) {
            // multiply by z mod f(z)
            A_mul_z(A);
            // accumulate M into A
            A_add_M(A, M, i, j);
        }
    }
}

/**
 * Multiplies y by z mod f(z).
 * 
 * @param[in,out] y Input/output vector y
 */
static inline void y_mul_z(uint8_t y[M_ENTRIES_BYTES]) {
    int i;
    uint8_t last;

    last = y[M_ENTRIES_BYTES - 1] >> 4;  // get last GF(16) element of y
    // shift y down by one GF(16) element
    for (i = M_ENTRIES_BYTES - 1; i > 0; i--) {
        y[i] = (y[i] << 4) | (y[i - 1] >> 4);
    }
    y[0] = (y[0] << 4);
    // add the contribution of the modulus
    y[0] ^= (mul_f(last, F_TAIL0) | (mul_f(last, F_TAIL1) << 4));
    y[1] ^= (mul_f(last, F_TAIL2) | (mul_f(last, F_TAIL3) << 4));
}

/**
 * Builds the target vector y from the vectors u.
 * 
 * @param[out] y Output target vector y (assumed to be initialized to t)
 * @param[in] u Input vectors u
 */
static inline void build_y(uint8_t y[M_ENTRIES_BYTES],
                           const uint8_t u[BINOMIAL2(PARAM_K + 1)][M_ENTRIES_BYTES]) {
    uint8_t t[M_ENTRIES_BYTES];  // holds the target vector y initial value
    int i, j, ell;

    // No leak since the value of y is public
    memcpy(t, y, M_ENTRIES_BYTES);  // copy initial value of target vector y
    memset(y, 0, M_ENTRIES_BYTES);  // clear y to accumulate into it

    // Build y from maximum ell to minimum ell
    ell = U_INDEX(PARAM_K - 1, PARAM_K - 1);
    vec_add_u8(y, u[ell], M_ENTRIES_BYTES);  // accumulate u[ell_max] into y
    // for all ell from second maximum to minimum
    for (i = PARAM_K - 2; i >= 0; i--) {
        for (j = i; j < PARAM_K; j++) {
            // multiply by z mod f(z)
            y_mul_z(y);
            // accumulate u[ell] into y
            ell = U_INDEX(i, j);
            vec_add_u8(y, u[ell], M_ENTRIES_BYTES);
        }
    }
    // finally add the initial value t
    vec_add_u8(y, t, M_ENTRIES_BYTES);
}

/**************************************************************************************
 *           Implementation of Gaussian elimination in constant-time                  *
 **************************************************************************************/

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

/**
 * Transforms the augmented matrix [A | y] into echelon form using Gaussian elimination
 * in constant-time.
 * 
 * @param[in,out] A Input/output matrix A
 * @param[in,out] y Input/output target vector y
 */
static inline void echelon_form(uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                                uint8_t y[M_ENTRIES_BYTES]) {
    alignas(8) uint8_t tmp0[TO_QWORDS(KO_ENTRIES_BYTES + 1) * 8];
    alignas(8) uint8_t tmp1[TO_QWORDS(KO_ENTRIES_BYTES + 1) * 8];
    int row, byte;  // row and byte indices
    int pivot_row, pivot_col;
    // lower and upper bounds for pivot row at each step
    int pivot_row_lb, pivot_row_ub;
    uint8_t pivot, pivot_inv;
    
    pivot_row = 0;
    for (pivot_col = 0; pivot_col < PARAM_K * PARAM_O; pivot_col++) {        
        // the pivot row is between these bounds if A has full rank
        pivot_row_lb = MAX(0, pivot_col + PARAM_M - PARAM_K * PARAM_O - 1);
        pivot_row_ub = MIN(PARAM_M - 1, pivot_col);

        // clear tmp0 and tmp1
        for (byte = 0; byte < KO_ENTRIES_BYTES + 1; byte++) {
            tmp0[byte] = 0;
            tmp1[byte] = 0;
        }

        // Iterate to find a pivot row considering at most 32 rows after the upper bound
        // since the probability of not finding a pivot after 32 rows is negligible (2^-128 for q = 16).
        // The pivot row is contained in tmp0 at the end of this loop
        uint8_t pivot_is_zero = 0xFF;
        pivot = 0;  // avoid compiler warning
        for (row = pivot_row_lb; row <= MIN(PARAM_M - 1, pivot_row_ub + 32); row++) {
            // The idea is to conditionally add the candidate pivot row into tmp0:
            // - if it is the current pivot_row, or 
            // - if we have not found a pivot yet and the index is below row,
            //   i.e. the left entries, w.r.t. pivot_col, of the row are all zero.
            // We continue this process until we find a nonzero pivot, then we
            // do not add any more rows, even if we don't stop the loop.
            // Then tmp0 will contain the sum of all these rows, with a nonzero pivot.

            uint8_t is_pivot_row = ~ct_compare_32(row, pivot_row);
            uint8_t below_pivot_row = ct_is_greater_than(row, pivot_row);
            
            // conditionally add row into tmp0
            for (byte = 0; byte < KO_ENTRIES_BYTES; byte++) {
                tmp0[byte] ^= A[row][byte] &
                    (is_pivot_row | (below_pivot_row & pivot_is_zero));
            }
            tmp0[KO_ENTRIES_BYTES] ^= nibble(y, row) &
                (is_pivot_row | (below_pivot_row & pivot_is_zero));

            // update pivot_is_zero
            pivot = nibble(tmp0, pivot_col);
            pivot_is_zero = ~ct_compare_8(pivot, 0);
        }

        // multiply pivot row by inverse of pivot and store in tmp1
        pivot_inv = inverse_f(pivot);
        vec_scalarmul_u64(tmp1, TO_QWORDS(KO_ENTRIES_BYTES + 1), tmp0, pivot_inv);

        // conditionally write pivot row to the correct row, i.e. pivot_row,
        // if the pivot is non-zero
        for (row = pivot_row_lb; row <= pivot_row_ub; row++) {
            uint8_t do_copy = ~ct_compare_32(row, pivot_row) & ~pivot_is_zero;
            uint8_t do_not_copy = ~do_copy;
            for (byte = 0; byte < KO_ENTRIES_BYTES; byte++) {
                A[row][byte] = (do_not_copy & A[row][byte]) + (do_copy & tmp1[byte]);
            }
            uint8_t new_nibble = (do_not_copy & nibble(y, row)) + 
                (do_copy & nibble(tmp1, PARAM_K * PARAM_O));
            set_nibble(y, row, new_nibble);
        }

        // eliminate entries below pivot
        for (row = pivot_row_lb; row < PARAM_M; row++) {
            // negate comparison result to get 1 if row < pivot_row, 0 otherwise
            uint8_t below_pivot = -ct_is_greater_than(row, pivot_row);
            uint8_t elt_to_elim = nibble(A[row], pivot_col);
            vec_scalarmul_u64(tmp0, TO_QWORDS(KO_ENTRIES_BYTES + 1), tmp1, below_pivot * elt_to_elim);
            vec_add_u8(A[row], tmp0, KO_ENTRIES_BYTES);
            add_nibble(y, row, nibble(tmp0, PARAM_K * PARAM_O));
        }

        // update pivot_row, i.e. increment iff pivot is non-zero
        pivot_row += (-(int8_t)(~pivot_is_zero));
    }
}

/**
 * Solves the upper-triangular system A * x = y using backward substitution
 * in constant-time.
 * 
 * @param[in] A Input upper-triangular matrix A
 * @param[in,out] y Input target vector y (its value will be modified during the computation)
 * @param[out] x Output solution vector x initialised with a randomized initial value
 */
static inline void backward_substitution(uint8_t A[PARAM_M][KO_ENTRIES_BYTES],
                                         uint8_t y[M_ENTRIES_BYTES],
                                         uint8_t x[KO_ENTRIES_BYTES]) {
    int i;
    int row, col;
    int column_ub;
    uint8_t finished;
    uint8_t correct_column;
    uint8_t first_nonzero;

    for (row = PARAM_M - 1; row >= 0; row--) {
        finished = 0;  // indicates whether we have found the first non-zero column in this row

        // The first nonzero entry in row r is between r and k*o.
        // We can limit the search to r + 32, and we will still find it
        // except with probability q^-32 (for q = 16 we have 2^-128 probability of failure).
        column_ub = MIN(PARAM_K * PARAM_O, row + 32);

        // find first non-zero column in row in constant-time
        for (col = row; col < column_ub; col++) {
            // Compare two bytes in constant time.
            correct_column = ct_compare_8(nibble(A[row], col), 0) & ~finished;
            first_nonzero = correct_column & nibble(y, row);

            // x[col] += first_nonzero
            add_nibble(x, col, first_nonzero);
            // y -= first_nonzero * A[:, col]
            for (i = 0; i < row; i++) {
                add_nibble(y, i, mul_f(first_nonzero, nibble(A[i], col)));
            }
            
            // mark as finished if we found the first non-zero column
            finished = finished | correct_column;
        }
    }
}

#endif  // ARITHMETIC_H
