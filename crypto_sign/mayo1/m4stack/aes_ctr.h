// SPDX-License-Identifier: Apache-2.0

#ifndef AES128CTR_H
#define AES128CTR_H

#include <stdint.h>
#include <stddef.h>

// AES-128-CTR context structure
typedef struct aes128ctr_ctx {
    uint8_t round_keys[11 * 16];  // expanded AES key
    uint8_t state[16];  // current output block
    uint64_t pos;  // position in the current keystream in bytes (lower 4 bits are offset in block, upper bits are block counter)
} aes128ctr_ctx;

/**
 * Initializes the AES-128-CTR context with the given key.
 * 
 * @param[out] ctx Pointer to the AES-128-CTR context to initialize
 * @param[in] key 16-byte AES key
 */
void aes128ctr_init(aes128ctr_ctx *ctx, const uint8_t key[16]);

/**
 * Generates AES-128-CTR keystream and writes it to the output buffer.
 * 
 * @param[out] out Output buffer to write the keystream
 * @param[in] outlen Length of the output buffer in bytes
 * @param[in,out] ctx Pointer to the AES-128-CTR context
 */
void aes128ctr_stream(uint8_t *out, int outlen, aes128ctr_ctx *ctx);

/**
 * AES ECB encryption function. Used by randombytes_ctrdrbg.c for KAT testing.
 */
void aes256ecb(const uint8_t *key, const uint8_t *input, uint8_t *output);
#define AES_ECB_encrypt(input, key, output) aes256ecb(key, input, output);

#endif  // AES128CTR_H
