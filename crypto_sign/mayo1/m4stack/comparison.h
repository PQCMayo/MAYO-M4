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
The following constant-time comparisons are directly taken from the MAYO
reference implementation, and adapted to be used in this project. The original code can be found
in arithmetic.h.
*/

#ifndef COMPARISON_H
#define COMPARISON_H

#include <stdint.h>
#include <stddef.h>

#include "blocker.h"

#if !(((!defined(__clang__) && defined(__GNUC__) && __GNUC__ <= 12)) && (defined(__x86_64__) || defined(_M_X64)))
// a > b -> b - a is negative
// returns 0xFFFFFFFF if true, 0x00000000 if false
static inline uint32_t ct_is_greater_than(int a, int b) {
    int32_t diff = b - a;
    return ((uint32_t) (diff >> (8*sizeof(uint32_t)-1)) ^ uint32_t_blocker);
}

// a > b -> b - a is negative
// returns 0xFFFFFFFFFFFFFFFF if true, 0x0000000000000000 if false
static inline uint64_t ct_64_is_greater_than(int a, int b) {
    int64_t diff = ((int64_t) b) - ((int64_t) a);
    return ((uint64_t) (diff >> (8*sizeof(uint64_t)-1)) ^ uint64_t_blocker);
}

// if a == b -> 0x00000000, else 0xFFFFFFFF
static inline uint32_t ct_compare_32(int a, int b) {
    return ((uint32_t)((-(int32_t)(a ^ b)) >> (8*sizeof(uint32_t)-1)) ^ uint32_t_blocker);
}

// if a == b -> 0x0000000000000000, else 0xFFFFFFFFFFFFFFFF
static inline uint64_t ct_compare_64(int a, int b) {
    return ((uint64_t)((-(int64_t)(a ^ b)) >> (8*sizeof(uint64_t)-1)) ^ uint64_t_blocker);
}

// if a == b -> 0x00, else 0xFF
static inline uint8_t ct_compare_8(uint8_t a, uint8_t b) {
    return ((int8_t)((-(int32_t)(a ^ b)) >> (8*sizeof(uint32_t)-1)) ^ uint8_t_blocker);
}
#else
// a > b -> b - a is negative
// returns 0xFFFFFFFF if true, 0x00000000 if false
static inline uint32_t ct_is_greater_than(int a, int b) {
    int32_t diff = b - a;
    return ((uint32_t) (diff >> (8*sizeof(uint32_t)-1)));
}

// a > b -> b - a is negative
// returns 0xFFFFFFFF if true, 0x00000000 if false
static inline uint64_t ct_64_is_greater_than(int a, int b) {
    int64_t diff = ((int64_t) b) - ((int64_t) a);
    return ((uint64_t) (diff >> (8*sizeof(uint64_t)-1)));
}

// if a == b -> 0x00000000, else 0xFFFFFFFF
static inline uint32_t ct_compare_32(int a, int b) {
    return ((uint32_t)((-(int32_t)(a ^ b)) >> (8*sizeof(uint32_t)-1)));
}

// if a == b -> 0x0000000000000000, else 0xFFFFFFFFFFFFFFFF
static inline uint64_t ct_compare_64(int a, int b) {
    return ((uint64_t)((-(int64_t)(a ^ b)) >> (8*sizeof(uint64_t)-1)));
}

// if a == b -> 0x00, else 0xFF
static inline uint8_t ct_compare_8(uint8_t a, uint8_t b) {
    return ((int8_t)((-(int32_t)(a ^ b)) >> (8*sizeof(uint32_t)-1)));
}
#endif

#endif  // COMPARISON_H