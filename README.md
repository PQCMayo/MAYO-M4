# MAYO-M4
Arm Cortex-M4 implementation of [MAYO](https://pqmayo.org/).
The code in this repository implements the Round-3 version of the MAYO specification.
It is based on the nibble-sliced implementation as described in the paper **Nibbling MAYO: Optimized Implementations for AVX2 and Cortex-M4** available [here](https://eprint.iacr.org/2023/1683.pdf), but with adapted parameters.

This repository is based on [pqm4](https://github.com/mupq/pqm4) and you will find the usual `test.py`, `testvectors.py`, and `benchmarks.py` scripts. 
Please follow the installation steps in pqm4. 
We target the [NUCLEO-L4R5ZI board](https://www.st.com/en/evaluation-tools/nucleo-l4r5zi.html), but tests can also be performed using qemu.

```
git clone --recurse-submodules https://github.com/PQCMayo/MAYO-M4.git
cd MAYO-M4
```

## Running tests and benchmarks
```
# run tests using qemu
./test.py -p mps2-an386 mayo1 mayo2 mayo3
# run testvectors using qemu
./testvectors.py -p mps2-an386 mayo1 mayo2 mayo3

# run tests on the board
./test.py -p nucleo-l4r5zi -u /dev/ttyACM0 mayo1 mayo2 mayo3
# run testvectors on the board
./testvectors.py -p nucleo-l4r5zi -u /dev/ttyACM0 mayo1 mayo2 mayo3
# run benchmarks on the board
./benchmarks.py -p nucleo-l4r5zi -u /dev/ttyACM0 mayo1 mayo2 mayo3 -i 1000
# print benchmarks
./convert_benchmarks.py md
```

## Benchmarks

Cycle counts on the NUCLEO-L4R5ZI, averaged over 1000 runs.
Built with `arm-none-eabi-gcc` 13.2.1 at `-O3`.

| Scheme | KeyGen | ExpandSK | ExpandPK | Sign | Verify | ExpandSK + Sign | ExpandPK + Verify |
| ------ | ------ | -------- | -------- | ---- | ------ | --------------- | ----------------- |
| MAYO-1 | 9,631,297 | 9,786,056 | 6,911,431 | 8,322,944 | 4,627,416 | 19,223,231 | 11,539,985 |
| MAYO-2 | 8,898,780 | 8,444,298 | 5,207,345 | 3,620,090 | 1,681,028 | 10,312,833 | 6,888,289 |
| MAYO-3 | 27,169,341 | 29,677,229 | 18,077,967 | 17,535,468 | 10,397,669 | 43,708,718 | 28,475,402 |


