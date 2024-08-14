This program contains the code necessary to process waveform data from the CAEN desktop digitiser.

This is a copy of a public repository with the same name from [@mahditaani](https://github.com/mahditaani), updated by [@rkralik5](https://github.com/rkralik5) to read waveforms from CAEN CoMPASS and improved fitting for multi-PE data


waveconvert program is in the convert directory
    This takes the CAEN xml readout file as an input and outputs a root file.

The scripts used to analyse the waveforms are in the scripts directory
    These scripts use the converted file (.root) as an input.
