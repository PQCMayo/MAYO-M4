// SPDX-License-Identifier: Apache-2.0

/*
The structure of this file is inspired by https://github.com/pq-crystals/dilithium.git
*/

#ifndef CONFIG_H
#define CONFIG_H

#ifndef MAYO_RANDOMIZED_SIGNING
#define MAYO_RANDOMIZED_SIGNING 1
#endif

#if MAYO_RANDOMIZED_SIGNING != 0 && MAYO_RANDOMIZED_SIGNING != 1
#error "MAYO_RANDOMIZED_SIGNING must be either 0 or 1"
#endif

#ifndef MAYO_MODE
#define MAYO_MODE 2
#endif

#if MAYO_MODE == 1
#define CRYPTO_ALGNAME "MAYO_1"
#define MAYO_NAMESPACETOP mayo1
#define MAYO_NAMESPACE(s) mayo1_##s
#elif MAYO_MODE == 2
#define CRYPTO_ALGNAME "MAYO_2"
#define MAYO_NAMESPACETOP mayo2
#define MAYO_NAMESPACE(s) mayo2_##s
#elif MAYO_MODE == 3
#define CRYPTO_ALGNAME "MAYO_3"
#define MAYO_NAMESPACETOP mayo3
#define MAYO_NAMESPACE(s) mayo3_##s
#elif MAYO_MODE == 5
#define CRYPTO_ALGNAME "MAYO_5"
#define MAYO_NAMESPACETOP mayo5
#define MAYO_NAMESPACE(s) mayo5_##s
#else
#error "Unsupported MAYO_MODE"
#endif // MAYO_MODE

#ifndef USE_BLOCKERS
#define USE_BLOCKERS 1
#endif

#if USE_BLOCKERS != 0 && USE_BLOCKERS != 1
#error "USE_BLOCKERS must be either 0 or 1"
#endif

#endif  // CONFIG_H
