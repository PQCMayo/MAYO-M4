# MAYO-M4
Arm Cortex-M4 implementation of [MAYO](https://pqmayo.org/).

This repository includes the Arm Cortex-M4 implementations as described in the paper **Nibbling MAYO: Optimized Implementations for AVX2 and Cortex-M4** available at [here](https://eprint.iacr.org/2023/1683.pdf).

It consists of two variants:
1. A slower version compatible with the round-1 specification of MAYO using bitsliced representation. This version is contained in the [main branch](https://github.com/PQCMayo/MAYO-M4/tree/main) of this repository.
2. A faster version that changes representation of keys and PRNG output to nibble-sliced representation compatible with the [nibbling-mayo branch](https://github.com/PQCMayo/MAYO-C/tree/nibbling-mayo) of the reference implementation. This version is contained in the [nibbling-mayo branch](https://github.com/PQCMayo/MAYO-M4/tree/nibbling-mayo) of this repository.

This repository is based on [pqm4](https://github.com/mupq/pqm4) and you will find the usual `test.py`, `testvectors.py`, and `benchmarks.py` scripts. 
Please follow the installation steps in pqm4. 
We target the [NUCLEO-L4R5ZI board](https://www.st.com/en/evaluation-tools/nucleo-l476rg.html), but tests can also be performed using qemu.

```
git clone --recurse-submodules https://github.com/PQCMayo/MAYO-M4.git
cd MAYO-M4
```

## Running tests and benchmarks for round-1 MAYO
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
./benchmarks.py -p nucleo-l4r5zi -u /dev/ttyACM0 mayo1 mayo2 mayo3 -i 10
# print benchmarks
./convert_benchmarks.py md
```

## Running tests and benchmarks for nibble-sliced MAYO
```
# switch to nibbling-mayo branch including dependencies
git checkout nibbling-mayo --recurse-submodules

# run tests using qemu
./test.py -p mps2-an386 mayo1 mayo2 mayo3
# run testvectors using qemu
./testvectors.py -p mps2-an386 mayo1 mayo2 mayo3

# run tests on the board
./test.py -p nucleo-l4r5zi -u /dev/ttyACM0 mayo1 mayo2 mayo3
# run testvectors on the board
./testvectors.py -p nucleo-l4r5zi -u /dev/ttyACM0 mayo1 mayo2 mayo3
# run benchmarks on the board
./benchmarks.py -p nucleo-l4r5zi -u /dev/ttyACM0 mayo1 mayo2 mayo3 -i 10
# print benchmarks
./convert_benchmarks.py md
```
