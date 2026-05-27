This repository contains ROOT-based code designed for processing and analysing
PMT waveform data from the CAEN desktop digitiser.

The code currently works with outputs from the CAEN Scope software. In the
future, the plan is to update it to also work with CoMPASS.

This repository was copied from an original public repository with the same
name from [@mahditaani][mahdiTaaniLink] and upgraded by
[@rkralik5][robertkralik]

The structure of the repository is as follows:
# [Waveform Processing](WaveformProcessing)
Converts CAEN digitizer readout files (`.xml` from CAENScope, `.dat` from QA DAQ)
into ROOT files. Two tools are provided:

- **`process_waveforms`** — processes all waveforms in bulk, extracting charge,
  baseline, and peak voltage into a ROOT TTree
- **`draw_waveforms`** — draws annotated sample waveforms for visual inspection

Both tools share a common physics library (`waveform_core`) and support
configurable integration gate and pre-gate widths via command-line options.
See the [WaveformProcessing README](WaveformProcessing/README.md) for full usage.

# [scripts](scripts)
Analysis scripts for waveform integration, dark rate analysis, or plotting.
These scripts use the ROOT files from [convert](convert) as inputs

# PMT Filename Parser

The PMT measurements are saved in a structured format containing the necessary
information. The filename parser extracts the PMT manufacturer, model, input
voltage, channel, and measurement information from filenames such as:
```
Ch1_Hamamatsu_R12860_1300V_Ch2_Hamamatsu_R12860_1300V_CrosstalkStudy.root
```


[robertkralik]: https://github.com/rkralik5
[mahdiTaaniLink]: https://github.com/mahditaani