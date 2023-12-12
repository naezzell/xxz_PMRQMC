#!/bin/bash
#
# This program is introduced in the paper:
# Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
#
# This program is licensed under a Creative Commons Attribution 4.0 International License:
# http://creativecommons.org/licenses/by/4.0/
#
#

# clean previous .bin file away
./clean.sh

# prepare Hamiltonian and observable helper files
g++ -O3 -std=c++11 -o xxz_prepare.bin xxz_prepare.cpp
./xxz_prepare.bin 10 0.1 0.5

# run QMC simulation
g++ -O3 -std=c++11 -o xxz_PMRQMC.bin xxz_PMRQMC.cpp
./xxz_PMRQMC.bin 10 0.1 0.5


