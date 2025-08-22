export BOOST=$HOME/anaconda3/enva/my-root-environment/include/boost

g++  -O2 -Wall -fPIC -pthread -std=c++17 -m64 -I$BOOST -I$ROOTSYS/include -c process_waveforms.cpp
g++  -O2 -m64 -std=c++17 process_waveforms.o -lm -L$ROOTSYS/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -o process_waveforms
rm process_waveforms.o