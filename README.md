# DreyFVS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6638217.svg)](https://doi.org/10.5281/zenodo.6638217)

DreyFVS is a solver for the Directed Feedback Vertex Set problem. It was submitted to the Heuristic Track of PACE Challenge 2022. This is an iterative solver applying two local search algorithms to improve a solution generated with a Sinkhorn-Knopp algorithm. More information can be found in the [solver description](solver_description.pdf)

## Requirements

* A 64-bit Linux operating system.
* A ![C++17](https://camo.githubusercontent.com/f6214f517f9a411ceb0a232746fc26dcefeaf7161239126ef7db4135b19fcd4d/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f432b2b2d31372d626c75652e7376673f7374796c653d666c6174) compiler, preferably GCC
* The [cmake](http://www.cmake.org/) build system (>= 3.5.1).

## Installation

To compile the code on a linux system you can use:
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

## Using the program

The program uses standard input to read the input graph and print the solution in the standard output just before finishing. Some debug information is printed in the standard error. The input and output formats are described on the [PACE 2022 - Heuristic Track Details](https://pacechallenge.org/2022/tracks/#heuristic-track). 

The program can be used like that:
```bash
./build/DFVS < in.txt > out.txt
```
Then the solver runs for 10 minutes and stops when it receives a signal SIGTERM.
