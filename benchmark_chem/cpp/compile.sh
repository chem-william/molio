#!/bin/bash

g++ -std=c++17 \
    -O3 \
    read_chem.cpp \
    -I/PATH/TO/CHEMFILES/INSTALLATION/chemfiles_install/include \
    -L/PATH/TO/CHEMFILES/INSTALLATION/chemfiles_install/lib \
    -lchemfiles \
    -o read_chem
