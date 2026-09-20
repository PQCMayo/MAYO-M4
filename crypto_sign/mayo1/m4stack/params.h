// SPDX-License-Identifier: Apache-2.0

/*
The structure of this file is inspired by https://github.com/pq-crystals/dilithium.git
*/

#ifndef PARAMS_H
#define PARAMS_H

#include "config.h"

#define PARAM_Q 16  // cardinality of the field
#define PK_SEED_BYTES 16  // number of bytes in the seed used to generate P1 and P2

#if MAYO_MODE == 1
#define PARAM_N 86
#define PARAM_M 78
#define PARAM_O 8
#define PARAM_K 10
#define SALT_BYTES 24  // number of bytes in the salt used for signing
#define DIGEST_BYTES 32  // number of bytes in the digest of the message to be signed
#define F_TAIL0 8  // constant coefficient in f(z) = z^78 + z^2 + z + x^3
#define F_TAIL1 1  // coefficient of z in f(z) = z^78 + z^2 + z + x^3
#define F_TAIL2 1  // coefficient of z^2 in f(z) = z^78 + z^2 + z + x^3
#define F_TAIL3 0  // coefficient of z^3 in f(z) = z^78 + z^2 + z + x^3

#elif MAYO_MODE == 2
#define PARAM_N 81
#define PARAM_M 64
#define PARAM_O 17
#define PARAM_K 4
#define SALT_BYTES 24  // number of bytes in the salt used for signing
#define DIGEST_BYTES 32  // number of bytes in the digest of the message to be signed
#define F_TAIL0 8  // constant coefficient in f(z) = z^64 + x^3*z^3 + x*z^2 + x^3
#define F_TAIL1 0  // coefficient of z in f(z) = z^64 + x^3*z^3 + x*z^2 + x^3
#define F_TAIL2 2  // coefficient of z^2 in f(z) = z^64 + x^3*z^3 + x*z^2 + x^3
#define F_TAIL3 8  // coefficient of z^3 in f(z) = z^64 + x^3*z^3 + x*z^2 + x^3

#elif MAYO_MODE == 3
#define PARAM_N 118
#define PARAM_M 108
#define PARAM_O 10
#define PARAM_K 11
#define SALT_BYTES 32  // number of bytes in the salt used for signing
#define DIGEST_BYTES 48  // number of bytes in the digest of the message to be signed
#define F_TAIL0 8  // constant coefficient in f(z) = z^108 + (x^2 + x + 1)*z^3 + z^2 + x^3
#define F_TAIL1 0  // coefficient of z in f(z) = z^108 + (x^2 + x + 1)*z^3 + z^2 + x^3
#define F_TAIL2 1  // coefficient of z^2 in f(z) = z^108 + (x^2 + x + 1)*z^3 + z^2 + x^3
#define F_TAIL3 7  // coefficient of z^3 in f(z) = z^108 + (x^2 + x + 1)*z^3 + z^2 + x^3

#elif MAYO_MODE == 5
#define PARAM_N 154
#define PARAM_M 142
#define PARAM_O 12
#define PARAM_K 12
#define SALT_BYTES 40  // number of bytes in the salt used for signing
#define DIGEST_BYTES 64  // number of bytes in the digest of the message to be signed
#define F_TAIL0 4  // constant coefficient in f(z) = z^142 + z^3 + x^3*z^2 + x^2
#define F_TAIL1 0  // coefficient of z in f(z) = z^142 + z^3 + x^3*z^2 + x^2
#define F_TAIL2 8  // coefficient of z^2 in f(z) = z^142 + z^3 + x^3*z^2 + x^2
#define F_TAIL3 1  // coefficient of z^3 in f(z) = z^142 + z^3 + x^3*z^2 + x^2

#else
#error "Unsupported MAYO_MODE"

#endif

#if PARAM_N != PARAM_M + PARAM_O
// sanity check since the implementation assumes N = M + O
#error "PARAM_N must be equal to PARAM_M + PARAM_O"
#endif

#if PARAM_M % 2 != 0
// sanity check since the implementation assumes PARAM_M is even
#error "PARAM_M must be even"
#endif

#if PARAM_O % 2 != 0 && PARAM_K % 2 != 0
// sanity check since the implementation assumes at least one of PARAM_O or PARAM_K is even
#error "PARAM_O and PARAM_K cannot both be odd"
#endif

#define BINOMIAL2(n) (((n) * ((n) - 1)) / 2)

#define SK_SEED_BYTES SALT_BYTES  // number of bytes in the seed used to generate the secret key
#define R_BYTES SALT_BYTES  // number of bytes in the random coins used for signing
#define O_BYTES (((PARAM_N - PARAM_O)*PARAM_O + 1)/2)  // number of bytes to represent the O matrix
#define V_BYTES ((PARAM_N - PARAM_O + 1)/2)  // number of bytes to store the vinegar variables
#define P1_BYTES (PARAM_M * BINOMIAL2(PARAM_N - PARAM_O + 1) / 2)  // number of bytes to represent the P1 matrices
#define P2_BYTES (PARAM_M * ((PARAM_N - PARAM_O)*PARAM_O) / 2)  // number of bytes to represent the P2 matrices
#define P3_BYTES (PARAM_M * BINOMIAL2(PARAM_O + 1) / 2)  // number of bytes to represent the P3 matrices
#define L_BYTES (PARAM_M * (PARAM_N - PARAM_O)*PARAM_O / 2)  // number of bytes to represent the L matrices

#define CSK_BYTES SK_SEED_BYTES  // number of bytes in the compact representation of a secret key
#define ESK_BYTES (SK_SEED_BYTES + O_BYTES + P1_BYTES + L_BYTES)  // number of bytes in the expanded representation of a secret key
#define CPK_BYTES (PK_SEED_BYTES + P3_BYTES)  // number of bytes in the compact representation of a public key
#define EPK_BYTES (P1_BYTES + P2_BYTES + P3_BYTES)  // number of bytes in the expanded representation of a public key
#define SIG_BYTES ((PARAM_N * PARAM_K + 1)/2 + SALT_BYTES)  // number of bytes in a signature

#define M_ENTRIES_BYTES (PARAM_M / 2)  // number of bytes to store M GF(16) elements
#define KO_ENTRIES_BYTES ((PARAM_K * PARAM_O) / 2)  // number of bytes to store k*o GF(16) elements

#define TO_BYTES(len) (((len) + 1) / 2)
#define TO_QWORDS(bytes) (((bytes) + 7) / 8)

#endif  // PARAMS_H