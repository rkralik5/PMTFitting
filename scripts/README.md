# Useful scripts to analyse waveforms

Inputs are ROOT files created by scripts in the [convert](../convert) folder

All the scripts are based off of ROOT and require a ROOT version installed

For convenience some of the scripts can be run using a bash script.
    For example:
        ./RunFitSingle.sh inputFile.root - will run the Fit.C script on a single file.
        ./RunFit.sh /path/to/files - will look through every file in a directory and run the Fit.C script on the root files it finds. 

Originally created by [@mahditaani](https://github.com/mahditaani) and updated by [@rkralik5](https://github.com/rkralik5)

## [DrawWaveform](DrawWaveform.C)

Draw a random selection of N waveforms from an input ROOT file using TGraph, save them into an output ROOT file

Run as:
    root 'DrawWaveform.C("inFileName","outFileName",N)'

## [Fit](Fit.C)

This script reads in the waveforms, finds a minimum peak and then integrates around that peak. The results are plotted on a histogram and then a fit is run to characterise the  PMT. 

The parameters of the fit can be found in: 

The result of the fit is written to a text file `pmt.txt`:
    filename    q0  s0  q1  s1  mu  gain    photoelectronResolution peakToValleyRatio

## [DarkRate](DarkRate.C)

This script reads in the waveforms and then looks looks for the number of times the charge is higher than a threshold value in a sliding time window. 
The threshold is given by the single pe integrated voltage (q1) multiplied by some factor. The takes several measurements then calculates a mean and standard deviation.
e.g.
    DarkRate.C("test.root",5,0.3)
        This means that the single pe has a mean integrated voltage of 5mV and it will seat a threshold of 0.3pe
The result of the fit is written to a text file `darkrates.txt`:
    filename    meanDR  stdDevDR