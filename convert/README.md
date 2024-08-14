This code was originally developed by [@mahditaani](https://github.com/mahditaani) and later updated by [@rkralik5](https://github.com/rkralik5)

The waveconvert program uses the CAEN xml readout file to create a root file for faster processing. 

The program requires ROOT and BOOST

To compile the program make sure you have the ROOTSYS and BOOST environment variables set.
The ROOTSYS variable is usually set by sourcing the `thisroot.sh` file from your ROOT installation
The BOOST environment variable should be set to the directory containing the `boost` directory.

Then source `compile.sh` and the program should compile. 

To use the program you should type in the following:
    waveconvert inputfile.xml outputfile.root 

This will create two trees with the following structure

device
    frequency   - frequency of sampling the wave
    NSamples  - total number of waves 
    resolution  - digitiser ADC resolution
    voltLow     - minimum voltage i.e. value from a reading of 0
    voltHigh    - maximum voltage i.e. value from a readion of `resolution`
data
    wavex - vector of x axis values for each wave representing time in ns 
    wavey - vector of y axis values for each wave representing recorded voltage in mV
    clocktime - time of recording of the waveform in units of Unix time in seconds
