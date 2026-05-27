This code was originally developed by [@mahditaani](https://github.com/mahditaani) and later updated by [@rkralik5](https://github.com/rkralik5).

Converts CAEN digitizer readout files (`.xml` or `.dat`) into ROOT files with processed waveform data. Two tools are provided: `process_waveforms` for bulk processing, and `draw_waveforms` for visual inspection of individual waveforms.

## Source files

- `waveform_core.cpp/hpp` — shared physics library (baseline calculation, charge integration)
- `process_waveforms_new.cpp` — bulk waveform processing tool
- `draw_waveforms.cpp` — waveform drawing/inspection tool

## Dependencies

- [ROOT](https://root.cern/) — `ROOTSYS` must be set (typically via `source thisroot.sh`)
- [Boost](https://www.boost.org/) — the include path must be set in `compile.sh` (the `ANACONDA` variable should point to the directory containing the `boost/` subdirectory)

## Compilation

```bash
source compile.sh
```

This produces two executables: `process_waveforms` and `draw_waveforms`.

## Usage

### process_waveforms

Processes all waveforms in an input file and writes extracted features to a ROOT TTree.

```
./process_waveforms [OPTIONS] <input_file> <output_file>
```

| Option | Description | Default |
|---|---|---|
| `--gate N` | Integration gate width in time bins | 50 |
| `--pregate N` | Number of bins before peak to start integration | 5 |
| `--help` | Show usage information | |

Example:
```bash
./process_waveforms --gate 100 --pregate 10 data.xml output.root
```

### draw_waveforms

Draws a sample of waveforms with annotations (baseline, minimum, integration range) and saves them as ROOT canvases.

```
./draw_waveforms [OPTIONS] <input_file> <output_file>
```

| Option | Description | Default |
|---|---|---|
| `--gate N` | Integration gate width in time bins | 50 |
| `--pregate N` | Number of bins before peak to start integration | 5 |
| `--num N` | Number of waveforms to draw | 5 |
| `--channel N` | Restrict to channel N (can be repeated for multiple channels) | all channels |
| `--help` | Show usage information | |

Example:
```bash
./draw_waveforms --gate 100 --pregate 10 --num 3 --channel 8 data.xml waveforms.root
```

## Output structure (`process_waveforms`)

The output ROOT file contains a single `TTree` named `Data` with the following branches:

| Branch | Type | Description |
|---|---|---|
| `Channel` | `UShort_t` | Digitizer channel number |
| `Timestamp` | `Long64_t` | Timestamp from start of run |
| `Clocktime` | `Long64_t` | Unix timestamp of the event |
| `Charge` | `Float_t` | Integrated charge in pC |
| `Baseline` | `Float_t` | Baseline voltage in mV |
| `MinVoltage` | `Float_t` | Minimum (peak) voltage in mV |
| `MinTime` | `Float_t` | Time of minimum in ns |
