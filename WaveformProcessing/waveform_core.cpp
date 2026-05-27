#include "waveform_core.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>

// Processing parameters (digitizer hardware constants)
const int fImpedance = 50;          ///< Impedance in Ohm (fixed for digitizer)
const float fTimeBinWidth = 2.0;    ///< Time bin width in ns (500MHz frequency)
const float fADCTomV = 2000.0 / 16384; ///< ADC to mV conversion (2000mV / 14 bit resolution)

int FindMinimum(const std::vector<int>& adcSamples, int& minADC, float& minTime) {
    // Don't look within a 5 bin-wide buffer in beginning and end
    auto minIt = std::min_element(adcSamples.begin() + 5, adcSamples.end() - 5);
    int minPos = std::distance(adcSamples.begin(), minIt);
    
    minADC = *minIt;
    minTime = minPos * fTimeBinWidth;
    
    return minPos;
}

float CalculateBaseline(const std::vector<int>& adcSamples, int minPos, int preGate, int gate) {
    // Get the integration range
    int lowInt = minPos > preGate ? minPos - preGate : 0;
    int highInt = lowInt + gate;
    if(static_cast<size_t>(highInt) >= adcSamples.size()) {
        highInt = adcSamples.size() - 1;
        lowInt = adcSamples.size() - 1 - gate;
    }

    // Make a vector for baseline calculation (values outside the gate)
    std::vector<int> vBaseline;
    vBaseline.reserve(adcSamples.size() - gate);
    for(size_t i = 0; i < adcSamples.size(); i++) {
        if(static_cast<int>(i) >= lowInt && static_cast<int>(i) < highInt) continue;
        vBaseline.push_back(adcSamples[i]);
    }

    // Calculate the baseline as a truncated mean of 50% of values outside the gate
    std::sort(vBaseline.begin(), vBaseline.end());
    float baseline = 0;
    int nBaseline = 0;
    for(size_t iBln = vBaseline.size()/4; iBln < (vBaseline.size() - vBaseline.size()/4); iBln++) {
        baseline += vBaseline[iBln];
        nBaseline++;
    }
    return nBaseline > 0 ? baseline / nBaseline : 0;
}

float IntegrateCharge(const std::vector<int>& adcSamples, float baselineADC,
                      int minPos, int preGate, int gate) {
    // Get the integration range
    int lowInt = minPos > preGate ? minPos - preGate : 0;
    int highInt = lowInt + gate;
    if(static_cast<size_t>(highInt) >= adcSamples.size()) {
        highInt = adcSamples.size() - 1;
        lowInt = adcSamples.size() - 1 - gate;
    }

    // Integrate the charge (in ADC counts)
    float chargeADC = 0;
    for(int iSample = lowInt; iSample < highInt; ++iSample) {
        chargeADC += baselineADC - adcSamples[iSample];
    }
    
    // Convert to pC: charge = (ADC * ADCTomV) * time / impedance
    float charge = chargeADC * fADCTomV * fTimeBinWidth / fImpedance;
    return charge;
}

void ProcessSingleWaveform(const std::vector<int>& adcSamples, float& charge,
                          float& baseline, float& minVolt, float& minTime,
                          int preGate, int gate) {
    // Find minimum position and ADC value
    int minADC;
    int minPos = FindMinimum(adcSamples, minADC, minTime);
    
    // Calculate baseline in ADC counts
    float baselineADC = CalculateBaseline(adcSamples, minPos, preGate, gate);
    
    // Integrate charge (function handles ADC→mV conversion internally)
    charge = IntegrateCharge(adcSamples, baselineADC, minPos, preGate, gate);

    // Convert ADC values to mV for output
    baseline = baselineADC * fADCTomV;
    minVolt = minADC * fADCTomV;
}

std::string GetFileType(const std::string& filename) {
    if(filename.find(".xml") != std::string::npos) {
        return "xml";
    } else if(filename.find(".dat") != std::string::npos) {
        return "dat";
    }
    return "unknown";
}
