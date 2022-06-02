# DreyFVS

This is a submission for the Heuristic Track of PACE 2022. This is an iterative solver applying two alogithms to improve a solution generating with a Sinkhorn algorithm.

## Installation

To compile code on a Debian system you can use:
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

## Using the program

The program uses standard input to read the input graph and print the solution in the standard output just before finishing. Some debug information is printed in the standard error. The program can be used like that:
```bash
./build/DFVS < in.txt > out.txt
```