// C++ Includes
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdint>

// BOOST Includes
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

// ROOT Includes
#include "TTree.h"
#include "TFile.h"

// Local includes
#include "waveform_core.hpp"

// Default integration parameters (can be overridden by command line arguments)
int fPreGate = 5;         ///< Number of time bins before peak position to start integration
int fGate = 50;           ///< Integration range for the waveform integration in time bins

/// @brief Print usage information
void PrintUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [OPTIONS] <input_file> <output_file>" << std::endl;
    std::cout << std::endl;
    std::cout << "Process PMT waveform data from XML or DAT files." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --gate N        Set integration gate width to N bins (default: " << fGate << ")" << std::endl;
    std::cout << "  --pregate N     Set pre-gate width to N bins (default: " << fPreGate << ")" << std::endl;
    std::cout << "  --help          Show this help message" << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  input_file      Input file (.xml or .dat)" << std::endl;
    std::cout << "  output_file     Output ROOT file" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << programName << " data.xml output.root" << std::endl;
    std::cout << "  " << programName << " --gate 100 --pregate 10 data.xml output.root" << std::endl;
}

/// @brief Parse command line arguments
/// @param argc Number of arguments
/// @param argv Argument array
/// @param inputFile Reference to store input filename
/// @param outputFile Reference to store output filename
/// @return true if parsing successful, false if should exit
bool ParseArguments(int argc, char* argv[], std::string& inputFile, std::string& outputFile) {
    std::vector<std::string> positionalArgs;
    
    for(int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if(arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return false;
        }
        else if(arg == "--gate") {
            if(i + 1 >= argc) {
                std::cerr << "Error: --gate requires a value" << std::endl;
                PrintUsage(argv[0]);
                return false;
            }
            try {
                fGate = std::stoi(argv[++i]);
                if(fGate <= 0) {
                    std::cerr << "Error: gate must be positive" << std::endl;
                    return false;
                }
            } catch(const std::exception& e) {
                std::cerr << "Error: invalid gate value '" << argv[i] << "'" << std::endl;
                return false;
            }
        }
        else if(arg == "--pregate") {
            if(i + 1 >= argc) {
                std::cerr << "Error: --pregate requires a value" << std::endl;
                PrintUsage(argv[0]);
                return false;
            }
            try {
                fPreGate = std::stoi(argv[++i]);
                if(fPreGate < 0) {
                    std::cerr << "Error: pregate must be non-negative" << std::endl;
                    return false;
                }
            } catch(const std::exception& e) {
                std::cerr << "Error: invalid pregate value '" << argv[i] << "'" << std::endl;
                return false;
            }
        }
        else if(arg[0] == '-') {
            std::cerr << "Error: unknown option '" << arg << "'" << std::endl;
            PrintUsage(argv[0]);
            return false;
        }
        else {
            positionalArgs.push_back(arg);
        }
    }
    
    // Check for required positional arguments
    if(positionalArgs.size() != 2) {
        std::cerr << "Error: exactly 2 arguments required (input_file, output_file)" << std::endl;
        PrintUsage(argv[0]);
        return false;
    }
    
    inputFile = positionalArgs[0];
    outputFile = positionalArgs[1];
    return true;
}

/// @brief Process XML file and convert to ROOT format with processed data
/// @param inName Input XML filename
/// @param outName Output ROOT filename
void ProcessXMLFile(std::string inName, std::string outName) {
    std::cout << "Processing XML file " << inName << " to " << outName << std::endl;
    
    // Populate tree structure pt
    using boost::property_tree::ptree;
    ptree pt;

    // Read the XML file into the property tree
    std::ifstream is(inName);
    read_xml(is, pt);
    
    // Read only the waveform size from XML settings
    int WSize = -1; ///< Size of the waveform (number of time bins)
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer.settings")){
        if(v.first == "window"){
            WSize = v.second.get<int>("<xmlattr>.size", -1);
        }
    }

    // Create ROOT file
    TFile *fOut = new TFile(outName.c_str(),"RECREATE");

    // Create processed data tree
    UShort_t Channel = -1;       ///< Channel number
    Long64_t Timestamp = -1;     ///< Timestamp of the event from start of run
    Long64_t Clocktime = -1;     ///< Clocktime of the event (in Unix time)
    Float_t Charge = -1;         ///< Integrated charge in pC
    Float_t Baseline = -1;       ///< Baseline voltage in mV
    Float_t MinVoltage = -1;     ///< Minimum voltage in mV
    Float_t MinTime = -1;        ///< Time of minimum in ns
    
    TTree *tData = new TTree("Data","Processed Wave Data");
    tData->Branch("Channel",    &Channel,    "Channel/s");
    tData->Branch("Timestamp",  &Timestamp,  "Timestamp/L");
    tData->Branch("Clocktime",  &Clocktime,  "Clocktime/L");
    tData->Branch("Charge",     &Charge,     "Charge/F");
    tData->Branch("Baseline",   &Baseline,   "Baseline/F");
    tData->Branch("MinVoltage", &MinVoltage, "MinVoltage/F");
    tData->Branch("MinTime",    &MinTime,    "MinTime/F");

    // Temporary storage for waveform samples during processing
    std::vector<int> adcSamples;

    std::cout << "Processing events..." << std::endl;

    // Loop over events (triggers)
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer")){
        if(v.first == "event"){
            int id = v.second.get<int>("<xmlattr>.id", -1);
            if(id%100 == 0)
                std::cout << "Event:\t" << id/1000 << "k\r" << std::flush;
                
            Timestamp = v.second.get<long long int>("<xmlattr>.timestamp", -1);
            Clocktime = v.second.get<long long int>("<xmlattr>.clocktime", -1);
            
            // Loop over all channels containing a waveform
            for(auto& t : v.second){
                if(t.first == "trace"){
                    Channel = t.second.get<int>("<xmlattr>.channel", -1);
                    
                    // Parse waveform data to vector
                    adcSamples.clear();
                    adcSamples.reserve(WSize);

                    std::stringstream ss(t.second.data()); // Hold waveform data
                    std::string value; // Hold each separate ADC value
    
                    // Parse space-separated ADC values
                    while(getline(ss, value, ' ')) {
                        if(!value.empty()) {
                            adcSamples.push_back(std::stoi(value));
                        }
                    }
                    
                    // Process the waveform to extract features
                    ProcessSingleWaveform(adcSamples, Charge, Baseline, MinVoltage, MinTime, fPreGate, fGate);
                    
                    // Fill the processed data tree
                    tData->Fill();
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
    infile.seekg(0, std::ios::beg); // Reset to the beginning for reading

    // Create ROOT file
    TFile* fOut = new TFile(outName.c_str(), "RECREATE");

    // Create processed data tree
    UShort_t Channel = -1;       ///< Channel number
    Long64_t Timestamp = -1;     ///< Timestamp of the event from start of run
    Long64_t Clocktime = -1;     ///< Clocktime of the event (in Unix time)
    Float_t Charge = -1;         ///< Integrated charge in pC
    Float_t Baseline = -1;       ///< Baseline voltage in mV
    Float_t MinVoltage = -1;     ///< Minimum voltage in mV
    Float_t MinTime = -1;        ///< Time of minimum in ns
    
    TTree *tData = new TTree("Data","Processed Wave Data");
    tData->Branch("Channel",    &Channel,    "Channel/s");
    tData->Branch("Timestamp",  &Timestamp,  "Timestamp/L");
    tData->Branch("Clocktime",  &Clocktime,  "Clocktime/L");
    tData->Branch("Charge",     &Charge,     "Charge/F");
    tData->Branch("Baseline",   &Baseline,   "Baseline/F");
    tData->Branch("MinVoltage", &MinVoltage, "MinVoltage/F");
    tData->Branch("MinTime",    &MinTime,    "MinTime/F");

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
        ProcessSingleWaveform(adcSamples, Charge, Baseline, MinVoltage, MinTime, fPreGate, fGate);
        
        // Fill the processed data tree
        tData->Fill();

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
    std::string inputFile, outputFile;
    
    // Parse command line arguments
    if(!ParseArguments(argc, argv, inputFile, outputFile)) {
        return 1;
    }
    
    // Print processing parameters
    std::cout << "Processing parameters:" << std::endl;
    std::cout << "  Gate: " << fGate << " bins" << std::endl;
    std::cout << "  PreGate: " << fPreGate << " bins" << std::endl;
    std::cout << std::endl;
    
    ProcessWaveforms(inputFile, outputFile);
    return 0;
}
