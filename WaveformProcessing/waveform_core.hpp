#ifndef WAVEFORM_CORE_HPP
#define WAVEFORM_CORE_HPP

#include <vector>
#include <string>

// Processing parameters (digitizer hardware constants)
extern const int fImpedance;          ///< Impedance in Ohm (fixed for digitizer)
extern const float fTimeBinWidth;     ///< Time bin width in ns (500MHz frequency)
extern const float fADCTomV;          ///< ADC to mV conversion (2000mV / 14 bit resolution)

/// @brief Find minimum value and its position in the waveform
/// @param adcSamples Vector of waveform samples in ADC counts
/// @param minADC Reference to store minimum ADC value
/// @param minTime Reference to store time of minimum
/// @return Position of minimum in the waveform
int FindMinimum(const std::vector<int>& adcSamples, int& minADC, float& minTime);

/// @brief Calculate baseline from waveform samples outside integration gate
/// @param adcSamples Vector of waveform samples in ADC counts
/// @param minPos Position of minimum in the waveform
/// @param preGate Number of bins before minimum for integration start
/// @param gate Integration gate width in bins
/// @return Baseline in ADC counts (truncated mean of 50% of values outside gate)
float CalculateBaseline(const std::vector<int>& adcSamples, int minPos, int preGate, int gate);

/// @brief Integrate waveform to get deposited charge
/// @param adcSamples Vector of waveform samples in ADC counts
/// @param baselineADC Baseline in ADC counts
/// @param minPos Position of minimum in the waveform
/// @param preGate Number of bins before minimum for integration start
/// @param gate Integration gate width in bins
/// @return Integrated charge in pC
float IntegrateCharge(const std::vector<int>& adcSamples, float baselineADC,
                      int minPos, int preGate, int gate);

/// @brief Process a single waveform and extract features
/// @param adcSamples Vector of raw ADC samples
/// @param charge Reference to store integrated charge in pC
/// @param baseline Reference to store baseline voltage in mV
/// @param minVolt Reference to store minimum voltage in mV
/// @param minTime Reference to store time of minimum in ns
/// @param preGate Number of bins before minimum for integration start
/// @param gate Integration gate width in bins
void ProcessSingleWaveform(const std::vector<int>& adcSamples, float& charge,
                          float& baseline, float& minVolt, float& minTime,
                          int preGate, int gate);

/// @brief Check file extension to determine file type
/// @param filename Input filename
/// @return File type: "xml", "dat", or "unknown"
std::string GetFileType(const std::string& filename);

#endif // WAVEFORM_CORE_HPP
