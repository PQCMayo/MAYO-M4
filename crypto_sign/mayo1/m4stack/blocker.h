// SPDX-License-Identifier: Apache-2.0

#ifndef BLOCKER_H
#define BLOCKER_H

#include <stdint.h>
#include <stddef.h>

#include "config.h"

// Variables to prevent certain compiler optimizations
// see https://github.com/PQCMayo/MAYO-C/pull/8/changes/3dace0e22d376451fac99899cbb3252d05712995#diff-c84e4f0b977868b489bda2319658059b92f056ac5793aceba06712d2dcd13c7f
extern volatile uint8_t uint8_t_blocker_variable;
extern volatile uint32_t uint32_t_blocker_variable;
extern volatile uint64_t uint64_t_blocker_variable;

#if USE_BLOCKERS == 1
#define uint8_t_blocker uint8_t_blocker_variable
#define uint32_t_blocker uint32_t_blocker_variable
#define uint64_t_blocker uint64_t_blocker_variable
#elif USE_BLOCKERS == 0
#define uint8_t_blocker 0
#define uint32_t_blocker 0
#define uint64_t_blocker 0
#else
#error "Invalid value for USE_BLOCKERS, expected 0 or 1"
#endif

#endif  // BLOCKER_H