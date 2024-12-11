.syntax unified
.cpu cortex-m4
.thumb

#include "verification.i"
#include "asm_params.h"


// For MAYO3 and MAYO5 we need to compute PS incrementally
// TODO: evaluate if those are still worthwhile, alternatively remove
//#if N <= 78
//.global calculate_PS_m4f_asm
//.type calculate_PS_m4f_asm, %function
//.align 2
//calculate_PS_m4f_asm:
//    macro_calculate_PS_m4f_asm (M/32), N, K
//#else
.global calculate_PS_m4f_stack_asm
.type calculate_PS_m4f_stack_asm, %function
.align 2
calculate_PS_m4f_stack_asm:
    macro_calculate_PS_m4f_stack_asm ((M + 15)/16), N, K

.global multiply_bins_stack_asm
.type multiply_bins_stack_asm, %function
.align 2
multiply_bins_stack_asm:
    macro_multiply_bins_stack_asm M, K

//#endif

.global calculate_SPS_m4f_asm
.type calculate_SPS_m4f_asm, %function
.align 2
calculate_SPS_m4f_asm:
    macro_calculate_SPS_m4f_asm ((M + 15)/16), N, K


.global multiply_bins_asm
.type multiply_bins_asm, %function
.align 2
multiply_bins_asm:
    macro_multiply_bins_asm M




