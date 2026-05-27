#!/bin/bash

echo "Compiling waveform processing tools..."

# Set boost path
export ANACONDA=/home/robertkralik/anaconda3/include

# Compile the core library
echo "  - Compiling waveform_core.cpp..."
g++ -O2 -Wall -fPIC -pthread -std=c++17 -m64 -I$ANACONDA -I$ROOTSYS/include -c waveform_core.cpp

# Compile the processing tool
echo "  - Compiling process_waveforms..."
g++ -O2 -Wall -fPIC -pthread -std=c++17 -m64 -I$ANACONDA -I$ROOTSYS/include -c process_waveforms.cpp
g++ -O2 -m64 -std=c++17 process_waveforms.o waveform_core.o -lm -L$ROOTSYS/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -o process_waveforms

# Compile the drawing tool
echo "  - Compiling draw_waveforms..."
g++ -O2 -Wall -fPIC -pthread -std=c++17 -m64 -I$ANACONDA -I$ROOTSYS/include -c draw_waveforms.cpp
g++ -O2 -m64 -std=c++17 draw_waveforms.o waveform_core.o -lm -L$ROOTSYS/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -o draw_waveforms

# Clean up object files
rm -f *.o

echo ""
echo "Compilation complete!"
echo "Usage:"
echo "  ./process_waveforms [--gate N] [--pregate N] input_file output_file"
echo "  ./draw_waveforms [--gate N] [--pregate N] [--num N] [--channel N] input_file output_file"
