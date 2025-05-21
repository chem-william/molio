#include <iostream>
#include <string>
#include "chemfiles.hpp"
using namespace chemfiles;

int main(int argc, char* argv[]) {
    std::string filename = "./water.pdb";
    
    if (argc > 1) {
        filename = argv[1];
    }
    
    auto file = Trajectory(filename);

    Frame frame;
    volatile int total_atoms = 0;
    while (!file.done()) {
        frame = file.read();
        total_atoms += frame.size();
    }
    
    volatile int throwaway = file.size();
    return 0;
} 
