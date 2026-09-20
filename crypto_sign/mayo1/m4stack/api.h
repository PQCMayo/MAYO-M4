#ifndef API_H
#define API_H

#include <stdint.h>
#include <stddef.h>

#include "params.h"

// NIST/SUPERCOP/PQClean API for signature schemes

#define CRYPTO_SECRETKEYBYTES CSK_BYTES
#define CRYPTO_PUBLICKEYBYTES CPK_BYTES
#define CRYPTO_BYTES SIG_BYTES

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);
int crypto_sign(unsigned char *sm, size_t *smlen, 
                const unsigned char *m, size_t mlen, 
                const unsigned char *sk);
int crypto_sign_open(unsigned char *m, size_t *mlen, 
                     const unsigned char *sm, size_t smlen, 
                     const unsigned char *pk);

#endif  // API_H
