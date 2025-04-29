export BOOST=/usr/local/opt/boost  # Boost path
export ROOTSYS=/usr/local/opt/root  # ROOT path
export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH

g++ -O2 -Wall -fPIC -pthread -std=c++20 -m64 \
    -I$BOOST/include \
    -I$(root-config --incdir) \
    -c main.cpp 

# Link the object file with ROOT and Boost libraries
g++ -O2 -m64 -std=c++20 main.o \
    -L$ROOTSYS/lib/root \
    -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
    -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc \
    -pthread -lm -ldl -rdynamic \
    -L$BOOST/lib \
    -lboost_system -lboost_filesystem \
    -o waveconvert

# Clean up the object file
rm main.o

