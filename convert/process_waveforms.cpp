// C++ Includes
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdint>

// BOOST Includes
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

// ROOT Includes
#include "TTree.h"
#include "TFile.h"

// Processing parameters (digitizer hardware constants)
const int fImpedance = 50;          ///< Impedance in Ohm (fixed for digitizer)
const float fFrequency = 500e6;     ///< Digitizer frequency in Hz (500 MS/s)
const float fTimeBinWidth = 2.0;    ///< Time bin width in ns (1000/500MHz = 2ns)
const int fResolutionBits = 14;     ///< ADC resolution in bits (always 14-bit)
const int fResolution = 16384;      ///< ADC resolution (2^14 = 16384)
const float fVoltLow = 0.0;         ///< Voltage range low in V (always 0V)
const float fVoltHigh = 2.0;        ///< Voltage range high in V (always 2V)
const float fADCTomV = 2000.0 / 16384; ///< ADC to mV conversion (2000mV / 2^14)

// Integration parameters
const int fPreGate = 5;         ///< Number of time bins before peak position to start integration
const int fGate = 50;           ///< Integration range for the waveform integration in time bins

/// @brief Find minimum value and its position in the waveform
/// @param adcSamples Vector of waveform samples in ADC counts
/// @param minADC Reference to store minimum ADC value
/// @param minTime Reference to store time of minimum
/// @param timeBinWidth Width of time bin in ns
/// @return Position of minimum in the waveform
int FindMinimum(const std::vector<int>& adcSamples, int& minADC, float& minTime, float timeBinWidth) {
    // Don't look within a 5 bin-wide buffer in beginning and end
    auto minIt = std::min_element(adcSamples.begin() + 5, adcSamples.end() - 5);
    int minPos = std::distance(adcSamples.begin(), minIt);
    
    minADC = *minIt;
    minTime = minPos * timeBinWidth;
    
    return minPos;
}

/// @brief Calculate baseline from waveform samples outside integration gate
/// @param adcSamples Vector of waveform samples in ADC counts
/// @param minPos Position of minimum in the waveform
/// @param preGate Number of bins before minimum for integration start
/// @param gate Integration gate width in bins
/// @return Baseline in ADC counts (truncated mean of 50% of values outside gate)
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

/// @brief Integrate waveform to get deposited charge
/// @param adcSamples Vector of waveform samples in ADC counts
/// @param baselineADC Baseline in ADC counts
/// @param minPos Position of minimum in the waveform
/// @param timeBinWidth Width of time bin in ns
/// @param preGate Number of bins before minimum for integration start
/// @param gate Integration gate width in bins
/// @return Integrated charge in pC
float IntegrateCharge(const std::vector<int>& adcSamples, float baselineADC, int minPos, 
                     float timeBinWidth, int preGate, int gate) {
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
    float charge = chargeADC * fADCTomV * timeBinWidth / fImpedance;
    return charge;
}

/// @brief Process a single waveform and extract features
/// @param adcSamples Vector of raw ADC samples
/// @param charge Reference to store integrated charge in pC
/// @param baseline Reference to store baseline voltage in mV
/// @param minVolt Reference to store minimum voltage in mV
/// @param minTime Reference to store time of minimum in ns
void ProcessSingleWaveform(const std::vector<int>& adcSamples, float& charge,
													 float& baseline, float& minVolt, float& minTime) {
    // Find minimum position and ADC value
    int minADC;
    int minPos = FindMinimum(adcSamples, minADC, minTime, fTimeBinWidth);
    
    // Calculate baseline in ADC counts
    float baselineADC = CalculateBaseline(adcSamples, minPos, fPreGate, fGate);
    
    // Integrate charge (function handles ADC→mV conversion internally)
    charge = IntegrateCharge(adcSamples, baselineADC, minPos, fTimeBinWidth, fPreGate, fGate);
    
    // Convert ADC values to mV for output
    baseline = baselineADC * fADCTomV;
    minVolt = minADC * fADCTomV;
}

/// @brief Parse waveform values directly from XML string to vector
/// @param line String of space separated waveform values from xml file
/// @param samples Output vector of waveform values (ADC counts)
/// @param expectedSize Expected number of samples for validation
void ParseXMLWaveformData(const std::string& line, std::vector<int>& samples, int expectedSize) {
    samples.clear();
    samples.reserve(expectedSize);
    
    std::stringstream ss(line);
    std::string value;
    
    // Parse space-separated ADC values
    while(getline(ss, value, ' ')) {
        if(!value.empty()) {
            samples.push_back(std::stoi(value));
        }
    }
}

/// @brief Process XML file and convert to ROOT format with processed data
/// @param inName Input XML filename
/// @param outName Output ROOT filename
void ProcessXMLFile(std::string inName, std::string outName) {
	std::cout << "Processing XML file " << inName << " to " << outName << std::endl;
	
	// Populate tree structure pt
  using boost::property_tree::ptree;
  ptree pt;

	std::cout << "Opening XML file..." << std::endl;
	std::ifstream is(inName);
	std::cout << "Reading XML content..." << std::endl;
  read_xml(is, pt);
	std::cout << "XML content read successfully." << std::endl;
	TFile *fOut = new TFile(outName.c_str(),"RECREATE");

	// Only read waveform size from XML - other parameters are hardware constants
	int WSize = -1; ///< Size of the waveform (number of time bins)
	
	// Store hardware constants in Device tree for reference (need non-const copies for ROOT)
	float frequency = fFrequency;
	int resolution = fResolution;
	float voltLow = fVoltLow;
	float voltHigh = fVoltHigh;
	float timeBinWidth = fTimeBinWidth;
	float ADCTomV = fADCTomV;
	
 	TTree *tDevice = new TTree("Device","Device Settings");
  tDevice->Branch("frequency",    &frequency,     "Frequency/F");
	tDevice->Branch("WSize",        &WSize,         "WaveformSize/I");
	tDevice->Branch("resolution",   &resolution,    "TimeResolution/I");
	tDevice->Branch("voltLow",      &voltLow,       "VoltageLow/F");
	tDevice->Branch("voltHigh",     &voltHigh,      "VoltageHigh/F");
	tDevice->Branch("timeBinWidth", &timeBinWidth,  "TimeBinWidth/F");
	tDevice->Branch("ADCTomV",      &ADCTomV,       "ADCTomV/F");

	std::cout << "Reading digitizer settings..." << std::endl;

	// Read only the waveform size from XML settings
	BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer.settings")){
		if(v.first == "window"){
			WSize = v.second.get<int>("<xmlattr>.size", -1);
		}
	}

	std::cout << "Waveform size: " << WSize << std::endl;
	
	// Validate that we found the waveform size
	if(WSize <= 0) {
		std::cerr << "Error: Could not read waveform size from XML file!" << std::endl;
		fOut->Close();
		return;
	}
	
	tDevice->Fill(); // Only fill once per file

	// Create processed data tree
	UShort_t Channel = -1;       ///< Channel number
	Long64_t Timestamp = -1;     ///< Timestamp of the event from start of run
	Long64_t Clocktime = -1;     ///< Clocktime of the event (in Unix time)
	Float_t Charge = -1;         ///< Integrated charge in pC
	Float_t Baseline = -1;       ///< Baseline voltage in mV
	Float_t MinVoltage = -1;     ///< Minimum voltage in mV
	Float_t MinTime = -1;        ///< Time of minimum in ns
	
	TTree *tProcessed = new TTree("ProcessedData","Processed Wave Data");
	tProcessed->Branch("Channel",    &Channel,    "Channel/s");
	tProcessed->Branch("Timestamp",  &Timestamp,  "Timestamp/L");
	tProcessed->Branch("Clocktime",  &Clocktime,  "Clocktime/L");
	tProcessed->Branch("Charge",     &Charge,     "Charge/F");
	tProcessed->Branch("Baseline",   &Baseline,   "Baseline/F");
	tProcessed->Branch("MinVoltage", &MinVoltage, "MinVoltage/F");
	tProcessed->Branch("MinTime",    &MinTime,    "MinTime/F");

	// Temporary storage for waveform samples during processing
	std::vector<int> adcSamples;

	std::cout << "Processing events..." << std::endl;

	// Loop over the events
	BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer")){
		if(v.first == "event"){
			int id = v.second.get<int>("<xmlattr>.id", -1);
			if(id%100 == 0)
				std::cout << "Event:\t" << id/1000 << "k\r" << std::flush;
				
			Timestamp = v.second.get<long long int>("<xmlattr>.timestamp", -1);
			Clocktime = v.second.get<long long int>("<xmlattr>.clocktime", -1);
			
			// Loop over the channels containing a waveform
			for(auto& t : v.second){
				if(t.first == "trace"){
					Channel = t.second.get<int>("<xmlattr>.channel", -1);
					
					// Parse waveform data directly to vector
					ParseXMLWaveformData(t.second.data(), adcSamples, WSize);
					
					// Process the waveform to extract features
					ProcessSingleWaveform(adcSamples, Charge, Baseline, MinVoltage, MinTime);
					
					// Fill the processed data tree
					tProcessed->Fill();
				}
			}
		}
  }
	
	std::cout << std::endl; // New line after progress updates
	fOut->Write();
	fOut->Close();
}

/// @brief Process DAT file and convert to ROOT format with processed data
/// @param inName Input DAT filename  
/// @param outName Output ROOT filename
void ProcessDATFile(std::string inName, std::string outName) {
    std::cout << "Processing DAT file " << inName << " to " << outName << std::endl;
    
    // Open .dat input file
    std::ifstream infile(inName, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open file: " << inName << std::endl;
        return;
    }

    // Determine waveform size from the first event
    uint32_t eventSize;
    infile.read(reinterpret_cast<char*>(&eventSize), sizeof(uint32_t));
    int WSize = (eventSize - 24) / 2; // 24 bytes header, rest is waveform data
    std::cout << "Waveform size: " << WSize << std::endl;
    infile.seekg(0, std::ios::beg); // Reset to the beginning for reading

    // Create ROOT file
    TFile* fOut = new TFile(outName.c_str(), "RECREATE");

    // Store hardware constants in Device tree for reference (need non-const copies for ROOT)
    float frequency = fFrequency;
    int resolution = fResolution;
    float voltLow = fVoltLow;
    float voltHigh = fVoltHigh;
    float timeBinWidth = fTimeBinWidth;
    float ADCTomV = fADCTomV;
    
    TTree *tDevice = new TTree("Device","Device Settings");
    tDevice->Branch("frequency",    &frequency,     "Frequency/F");
    tDevice->Branch("WSize",        &WSize,         "WaveformSize/I");
    tDevice->Branch("resolution",   &resolution,    "TimeResolution/I");
    tDevice->Branch("voltLow",      &voltLow,       "VoltageLow/F");
    tDevice->Branch("voltHigh",     &voltHigh,      "VoltageHigh/F");
    tDevice->Branch("timeBinWidth", &timeBinWidth,  "TimeBinWidth/F");
    tDevice->Branch("ADCTomV",      &ADCTomV,       "ADCTomV/F");
    tDevice->Fill(); // Only fill once per file

    // Create processed data tree
    UShort_t Channel = -1;       ///< Channel number
    Long64_t Timestamp = -1;     ///< Timestamp of the event from start of run
    Long64_t Clocktime = -1;     ///< Clocktime of the event (in Unix time)
    Float_t Charge = -1;         ///< Integrated charge in pC
    Float_t Baseline = -1;       ///< Baseline voltage in mV
    Float_t MinVoltage = -1;     ///< Minimum voltage in mV
    Float_t MinTime = -1;        ///< Time of minimum in ns
    
    TTree *tProcessed = new TTree("ProcessedData","Processed Wave Data");
    tProcessed->Branch("Channel",    &Channel,    "Channel/s");
    tProcessed->Branch("Timestamp",  &Timestamp,  "Timestamp/L");
    tProcessed->Branch("Clocktime",  &Clocktime,  "Clocktime/L");
    tProcessed->Branch("Charge",     &Charge,     "Charge/F");
    tProcessed->Branch("Baseline",   &Baseline,   "Baseline/F");
    tProcessed->Branch("MinVoltage", &MinVoltage, "MinVoltage/F");
    tProcessed->Branch("MinTime",    &MinTime,    "MinTime/F");

    // Allocate space for reading
    uint32_t header[6];
    uint16_t* rawSamples = new uint16_t[WSize];
    std::vector<int> adcSamples;
    adcSamples.reserve(WSize);
    int NSamples = 0;

    // Read loop
    while (infile.read(reinterpret_cast<char*>(header), sizeof(header))) {
        if (!infile.read(reinterpret_cast<char*>(rawSamples), WSize * sizeof(uint16_t)))
            break;

        Channel = static_cast<UShort_t>(header[3]);
        Timestamp = static_cast<Long64_t>(header[4]);
        Clocktime = static_cast<Long64_t>(header[5]);

        // Convert uint16_t array to std::vector<int> for processing
        adcSamples.clear();
        for (int i = 0; i < WSize; ++i) {
            adcSamples.push_back(static_cast<int>(rawSamples[i]));
        }

        // Process the waveform to extract features
        ProcessSingleWaveform(adcSamples, Charge, Baseline, MinVoltage, MinTime);
        
        // Fill the processed data tree
        tProcessed->Fill();

        if (++NSamples % 10000 == 0)
            std::cout << "Event:\t" << NSamples / 1000 << "k\r" << std::flush;
    }

    std::cout << std::endl; // New line after progress updates
    std::cout << "Total events processed: " << NSamples << std::endl;

    // Finalize
    fOut->Write();
    fOut->Close();
    delete[] rawSamples;
    
    std::cout << "Conversion complete." << std::endl;
}

/// @brief Check file extension to determine file type
/// @param filename Input filename
/// @return File type: "xml", "dat", or "unknown"
std::string GetFileType(const std::string& filename) {
    if(filename.find(".xml") != std::string::npos) {
        return "xml";
    } else if(filename.find(".dat") != std::string::npos) {
        return "dat";
    }
    return "unknown";
}

/// @brief Main conversion function that detects file type and processes accordingly
/// @param inName Input filename (XML or DAT)
/// @param outName Output ROOT filename
void ProcessWaveforms(std::string inName, std::string outName) {
    std::string fileType = GetFileType(inName);
    
    if(fileType == "xml") {
        ProcessXMLFile(inName, outName);
    } else if(fileType == "dat") {
        ProcessDATFile(inName, outName);
    } else {
        std::cerr << "Unsupported file type. Only .xml and .dat files are supported." << std::endl;
        std::cerr << "Input file: " << inName << std::endl;
    }
}

int main(int argc, char *argv[]){	
	if(argc != 3) {
		std::cout << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
		std::cout << "Supported input formats: .xml, .dat" << std::endl;
		return 1;
	}
	ProcessWaveforms(argv[1], argv[2]);
	return 0;
}
