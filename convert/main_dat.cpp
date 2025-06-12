#include <iostream>
#include <fstream>
#include <cstdint>
#include "TFile.h"
#include "TTree.h"
#include "TArrayS.h"

void ConvertDAT(const std::string& inName, const std::string& outName) {
    std::cout << "Converting " << inName << " to " << outName << std::endl;

    //need to set these values manually since using wavedump outputs
    // If using CoMPASS outputs need to set these values manually
//#define fResolution 500 ///< Number of ADC bins (=2^ADCResolution)
//#define fVoltLow 0.f ///< Voltage range low
//#define fVoltHigh 2.f ///< Voltage range high
//#define fFrequency 500e6 ///< Frequency of the digitiser in Hz
//#define fWindowSize 600 ///< Size of the waveform (number of time bins)

    // Hardcoded digitizer settings
    float frequency = 500e6;      // 500 MHz
    int NSamples = 1030; // Number of samples in the waveform
    int WSize = 1030;    // Size of the waveform (number of time bins)
    int resolution = 16384;        // 14-bit ADC → 2^14
    float voltLow = 0;
    float voltHigh = 2.0;



    // Open input .dat file
    std::ifstream infile(inName, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open file: " << inName << std::endl;
        return;
    }

    // Read header
    uint32_t header[6]; // [EventSize, BoardID, Pattern, Channel, EventCounter, TriggerTimeTag]
    infile.read(reinterpret_cast<char*>(header), sizeof(header));
    if (infile.gcount() != sizeof(header)) {
        std::cerr << "Header read failed" << std::endl;
        return;
    }

    // Read waveform data
    uint16_t rawSamples[1030];
    infile.read(reinterpret_cast<char*>(rawSamples), sizeof(rawSamples));
    if (infile.gcount() != sizeof(rawSamples)) {
        std::cerr << "Data read failed" << std::endl;
        return;
    }

    // Create output ROOT file
    TFile* outFile = new TFile(outName.c_str(), "RECREATE");

    // ----- Device TTree -----
    TTree* tDevice = new TTree("Device", "Device Settings");
    tDevice->Branch("frequency", &frequency, "Frequency/F");
    tDevice->Branch("NSamples", &NSamples, "NSamples/I");
    tDevice->Branch("WSize", &WSize, "WaveformSize/I");
    tDevice->Branch("resolution", &resolution, "TimeResolution/I");
    tDevice->Branch("voltLow", &voltLow, "VoltageLow/F");
    tDevice->Branch("voltHigh", &voltHigh, "VoltageHigh/F");
    tDevice->Fill();

    // ----- Data TTree -----
    UShort_t Channel = header[3];
    Long64_t Timestamp = static_cast<Long64_t>(header[4]);
    Long64_t Clocktime = static_cast<Long64_t>(header[5]);

    TArrayS* Samples = new TArrayS(WSize);
    for (int i = 0; i < WSize; ++i) {
        Samples->SetAt(rawSamples[i], i);
    }

    TTree* tData = new TTree("Data", "Wave Data");
    tData->Branch("Channel", &Channel, "Channel/S");
    tData->Branch("Timestamp", &Timestamp, "Timestamp/L");
    tData->Branch("Clocktime", &Clocktime, "Clocktime/L");
    tData->Branch("Samples", &Samples);

    tData->Fill();

    // ----- Write and Close -----
    outFile->Write();
    outFile->Close();
    std::cout << "Conversion complete.\n";
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./converter input.dat output.root" << std::endl;
        return 1;
    }
    ConvertDAT(argv[1], argv[2]);
    return 0;
}