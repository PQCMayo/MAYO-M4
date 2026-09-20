// SPDX-License-Identifier: Apache-2.0

// standard headers
#include <stddef.h>
#include <stdint.h>
#include <string.h>
// project headers
#include <mayo.h>
#include <randombytes.h>
// local headers
#include "params.h"
#include "fips202.h"
#include "aes_ctr.h"
#include "arithmetic.h"

int mayo_keypair_internal(const mayo_params_t *p, uint8_t *pk, uint8_t *seed_sk) {
    (void)p;  // unused parameter

    // contiguous buffer for seed_pk and oil_subspace
    uint8_t buf[PK_SEED_BYTES + O_BYTES];

    // derive seed_pk and oil_subspace from seed_sk
    shake256(buf, PK_SEED_BYTES + O_BYTES, seed_sk, SK_SEED_BYTES);

    // copy seed_pk to pk
    memcpy(pk, buf, PK_SEED_BYTES);
    // zero initialize P3 part of the public key
    memset(pk + PK_SEED_BYTES, 0, P3_BYTES);

    // compute P3 part of the public key
    compute_P3(pk + PK_SEED_BYTES, buf, buf + PK_SEED_BYTES);

    // clear sensitive data
    memset(buf, 0, sizeof(buf));

    return MAYO_OK;
}

int mayo_keypair(const mayo_params_t *p, uint8_t *pk, uint8_t *sk) {
    (void)p;  // unused parameter
    uint8_t seed_sk[SK_SEED_BYTES];

    // pick seed_sk securely at random
    if (randombytes(seed_sk, SK_SEED_BYTES) != 0) {
        return MAYO_ERR;  // error in randomness source
    }

    // generate keypair from seed_sk
    int ret = mayo_keypair_internal(p, pk, seed_sk);
    if (ret != MAYO_OK) {
        memset(pk, 0, CPK_BYTES);
        memset(seed_sk, 0, SK_SEED_BYTES);
        return ret;  // error in keypair generation
    }

    // copy seed_sk to sk
    memcpy(sk, seed_sk, SK_SEED_BYTES);

    // clear sensitive data
    memset(seed_sk, 0, SK_SEED_BYTES);

    return MAYO_OK;
}

int mayo_expand_pk(const mayo_params_t *p, const uint8_t *cpk, uint64_t *epk) {
    (void)p;  // unused parameter
    (void)cpk;  // unused parameter
    (void)epk;  // unused parameter
    return MAYO_OK;

    /*
    aes128ctr_ctx ctx;
    aes128ctr_init(&ctx, cpk);  // initialize AES-128-CTR with seed_pk
    aes128ctr_stream(epk, P1_BYTES + P2_BYTES, &ctx);  // generate P1 and P2
    memcpy(epk + P1_BYTES + P2_BYTES, cpk + PK_SEED_BYTES, P3_BYTES);  // copy P3
    return MAYO_OK;
    */
}

int mayo_expand_sk(const mayo_params_t *p, const uint8_t *csk, sk_t *esk) {
    (void)p;  // unused parameter
    (void)csk;  // unused parameter
    (void)esk;  // unused parameter
    return MAYO_OK;

    /*
    uint8_t buf[PK_SEED_BYTES + O_BYTES];  // contiguous buffer for seed pk and oil subspace
    aes128ctr_ctx ctx;
    
    // append seed_sk to esk
    memcpy(esk, csk, SK_SEED_BYTES);

    // derive seed_pk and oil_subspace from seed_sk
    shake256(buf, PK_SEED_BYTES + O_BYTES, csk, SK_SEED_BYTES);
    // append O to esk
    memcpy(esk + SK_SEED_BYTES, buf + PK_SEED_BYTES, O_BYTES);

    // generate P1 using AES-128-CTR with seed_pk
    aes128ctr_init(&ctx, buf);  // initialize AES-128-CTR with seed_pk
    aes128ctr_stream(esk + SK_SEED_BYTES + O_BYTES, P1_BYTES, &ctx);  // generate P1

    // compute L using O, P1 and P2 (generated on-the-fly with AES-128-CTR)
    compute_L(esk + SK_SEED_BYTES + O_BYTES + P1_BYTES,
              buf + PK_SEED_BYTES,
              esk + SK_SEED_BYTES + O_BYTES,
              &ctx);
    
    // clear sensitive data
    memset(buf, 0, sizeof(buf));
    memset(&ctx, 0, sizeof(ctx));

    return MAYO_OK;
    */
}

int mayo_sign_signature_internal(const mayo_params_t *p, 
                                 uint8_t *sig, size_t *siglen,
                                 const uint8_t *m, size_t mlen,
                                 const uint8_t rnd[R_BYTES],
                                 const uint8_t *sk) {
    (void)p;  // unused parameter
    int ctr;
    uint8_t buf[PK_SEED_BYTES + O_BYTES];  // contiguous buffer for seed pk and oil subspace
    // contiguous buffer for message digest, salt (or randomness), seed_sk and counter
    uint8_t tmp[DIGEST_BYTES + SALT_BYTES + SK_SEED_BYTES + 1];
    
    // vinegar variables and r (after sample_solution the solution x is in place of r)
    uint8_t V[PARAM_K * V_BYTES + KO_ENTRIES_BYTES];
    uint8_t A[PARAM_M][KO_ENTRIES_BYTES] = {0};  // linear system matrix
    uint8_t y[M_ENTRIES_BYTES];  // linear system target vector (initialized to t)
    
    // derive seed pk and oil subspace from seed_sk
    shake256(buf, PK_SEED_BYTES + O_BYTES, sk, SK_SEED_BYTES);

    // compute message digest
    shake256(tmp, DIGEST_BYTES, m, mlen);
    // append randomness and seed_sk
    memcpy(tmp + DIGEST_BYTES, rnd, R_BYTES);
    memcpy(tmp + DIGEST_BYTES + R_BYTES, sk, SK_SEED_BYTES);
    // compute salt (overlapping buffer, works for this implementation of shake256)
    shake256(tmp + DIGEST_BYTES, SALT_BYTES, tmp, DIGEST_BYTES + R_BYTES + SK_SEED_BYTES);
    
    // fiat-shamir with aborts
    for (ctr = 0; ctr < 256; ctr++) {
        // derive vinegar variables and r
        tmp[DIGEST_BYTES + SALT_BYTES + SK_SEED_BYTES] = (uint8_t)ctr;
        shake256(V, 
                 PARAM_K * V_BYTES + KO_ENTRIES_BYTES, 
                 tmp, 
                 DIGEST_BYTES + SALT_BYTES + SK_SEED_BYTES + 1);
        
        // initialize linear system A and target vector y
        shake256(y, M_ENTRIES_BYTES, tmp, DIGEST_BYTES + SALT_BYTES);
        build_linear_system(A, y, buf, buf + PK_SEED_BYTES, V);

        // solve linear system Ax = y
        if (sample_solution(V + PARAM_K * V_BYTES, A, y)) {
            break;
        } else {
            memset(A, 0, sizeof(A));  // reset A for next iteration
            memset(y, 0, sizeof(y));  // reset y for next iteration
        }
    }
    // ctr == 256 with negligible probability
    
    // compute s from solution x, oil subspace O and vinegar variables V
    compute_s(sig, V + PARAM_K * V_BYTES, buf + PK_SEED_BYTES, V);
    // append salt to signature
    memcpy(sig + SIG_BYTES - SALT_BYTES, tmp + DIGEST_BYTES, SALT_BYTES);
    // set signature length
    *siglen = SIG_BYTES;

    // clear sensitive data
    memset(buf, 0, sizeof(buf));
    memset(tmp, 0, sizeof(tmp));
    memset(V, 0, sizeof(V));
    memset(A, 0, sizeof(A));
    memset(y, 0, sizeof(y));

    return MAYO_OK;
}

int mayo_sign_signature(const mayo_params_t *p, 
                        uint8_t *sig, size_t *siglen,
                        const uint8_t *m, size_t mlen,
                        const uint8_t *sk) {
    (void)p;  // unused parameter
    uint8_t rnd[R_BYTES];
#if MAYO_RANDOMIZED_SIGNING
    // pick randomness for signing
    if (randombytes(rnd, R_BYTES) != 0) {
        return MAYO_ERR;  // error in randomness source
    }
#else
    // set randomness for signing to zero
    memset(rnd, 0, R_BYTES);
#endif
    return mayo_sign_signature_internal(p, sig, siglen, m, mlen, rnd, sk);
}

int mayo_sign(const mayo_params_t *p, 
              uint8_t *sm, size_t *smlen, 
              const uint8_t *m, size_t mlen, 
              const uint8_t *sk) {
    (void)p;  // unused parameter
    size_t siglen = 0;
    int ret = mayo_sign_signature(p, sm, &siglen, m, mlen, sk);
    if (ret != MAYO_OK || siglen != SIG_BYTES) {
        memset(sm, 0, siglen);
        *smlen = 0;
        return ret;
    }
    memcpy(sm + SIG_BYTES, m, mlen);
    *smlen = SIG_BYTES + mlen;
    return ret;
}

int mayo_verify(const mayo_params_t *p, 
                const uint8_t *m, size_t mlen,
                const uint8_t *sig,
                const uint8_t *pk) {
    (void)p;  // unused parameter
    // contiguous buffer for message digest and salt, or
    // result vector y of the quadratic map
    // size is the larger of the two
#if DIGEST_BYTES + SALT_BYTES > M_ENTRIES_BYTES
    uint8_t tmp[DIGEST_BYTES + SALT_BYTES];
#else
    uint8_t tmp[M_ENTRIES_BYTES];
#endif
    uint8_t t[M_ENTRIES_BYTES];  // target vector for verification
    
    // compute message digest
    shake256(tmp, DIGEST_BYTES, m, mlen);
    // append salt from signature
    memcpy(tmp + DIGEST_BYTES, sig + SIG_BYTES - SALT_BYTES, SALT_BYTES);
    // compute target vector t for linear system
    shake256(t, M_ENTRIES_BYTES, tmp, DIGEST_BYTES + SALT_BYTES);

    // evaluate quadratic map at signature and store result y in tmp
    quadratic_map(tmp, sig, pk, pk + PK_SEED_BYTES);

    if (memcmp(tmp, t, M_ENTRIES_BYTES) == 0) {
        return MAYO_OK;  // signature is valid
    } else {
        return MAYO_ERR;  // signature is invalid
    }
}

int mayo_open(const mayo_params_t *p,
              uint8_t *m, size_t *mlen,
              const uint8_t *sm, size_t smlen,
              const uint8_t *pk) {
    (void)p;  // unused parameter
    if (smlen < SIG_BYTES) {
        *mlen = 0;
        return MAYO_ERR;  // signature too short to be valid
    }
    int ret = mayo_verify(p, sm + SIG_BYTES, smlen - SIG_BYTES, sm, pk);
    if (ret != MAYO_OK) {
        *mlen = 0;
        return ret;  // signature verification failed
    }
    // signature verification succeeded, copy message to output
    memcpy(m, sm + SIG_BYTES, smlen - SIG_BYTES);
    *mlen = smlen - SIG_BYTES;
    return MAYO_OK;
}
