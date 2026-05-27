// C++ Includes
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <set>

// BOOST Includes
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

// ROOT Includes
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TAxis.h"

// Local includes
#include "waveform_core.hpp"

// Default parameters
int fPreGate = 5;         ///< Number of time bins before peak position to start integration
int fGate = 50;           ///< Integration range for the waveform integration in time bins
int fNumWaveformsToDraw = 5;    ///< Number of waveforms to draw
std::set<int> fChannelsToProcess; ///< Channels to process (empty = all channels)

/// @brief Print usage information
void PrintUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [OPTIONS] <input_file> <output_file>" << std::endl;
    std::cout << std::endl;
    std::cout << "Draw sample PMT waveforms from XML or DAT files." << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --gate N        Set integration gate width to N bins (default: " << fGate << ")" << std::endl;
    std::cout << "  --pregate N     Set pre-gate width to N bins (default: " << fPreGate << ")" << std::endl;
    std::cout << "  --num N         Number of waveforms to draw (default: " << fNumWaveformsToDraw << ")" << std::endl;
    std::cout << "  --channel N     Channel to process (can be used multiple times, default: all)" << std::endl;
    std::cout << "  --help          Show this help message" << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments:" << std::endl;
    std::cout << "  input_file      Input file (.xml or .dat)" << std::endl;
    std::cout << "  output_file     Output ROOT file with canvases" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << programName << " data.xml waveforms.root" << std::endl;
    std::cout << "  " << programName << " --gate 100 --pregate 10 --num 3 data.xml waveforms.root" << std::endl;
    std::cout << "  " << programName << " --channel 0 --channel 2 data.xml waveforms.root" << std::endl;
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
        else if(arg == "--num") {
            if(i + 1 >= argc) {
                std::cerr << "Error: --num requires a value" << std::endl;
                PrintUsage(argv[0]);
                return false;
            }
            try {
                fNumWaveformsToDraw = std::stoi(argv[++i]);
                if(fNumWaveformsToDraw <= 0) {
                    std::cerr << "Error: number of waveforms must be positive" << std::endl;
                    return false;
                }
            } catch(const std::exception& e) {
                std::cerr << "Error: invalid number value '" << argv[i] << "'" << std::endl;
                return false;
            }
        }
        else if(arg == "--channel") {
            if(i + 1 >= argc) {
                std::cerr << "Error: --channel requires a value" << std::endl;
                PrintUsage(argv[0]);
                return false;
            }
            try {
                int channel = std::stoi(argv[++i]);
                if(channel < 0) {
                    std::cerr << "Error: channel must be non-negative" << std::endl;
                    return false;
                }
                fChannelsToProcess.insert(channel);
            } catch(const std::exception& e) {
                std::cerr << "Error: invalid channel value '" << argv[i] << "'" << std::endl;
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

/// @brief Draw a single waveform with annotations
/// @param adcSamples Vector of ADC samples
/// @param plotNum Plot number for naming
/// @param channel Channel number
/// @param eventId Event ID
void DrawSingleWaveform(const std::vector<int>& adcSamples, int plotNum, int channel, int eventId) {
    // Process the waveform to get parameters
    float charge, baseline, minVolt, minTime;
    ProcessSingleWaveform(adcSamples, charge, baseline, minVolt, minTime, fPreGate, fGate);
    
    std::cout << "Waveform " << plotNum << " (Event " << eventId << ", Channel " << channel 
              << "): charge = " << charge << " pC" << std::endl;
    
    // Create canvas
    TCanvas* c = new TCanvas(Form("Waveform_%i_Channel_%i", plotNum, channel), 
                            Form("Waveform %i (Event %i, Channel %i)", plotNum, eventId, channel),
                            800, 600);
    gStyle->SetOptStat(0);
    
    // Create graph
    TGraph* grWaveform = new TGraph(adcSamples.size());
    for(size_t iSample = 0; iSample < adcSamples.size(); ++iSample){
        grWaveform->SetPoint(iSample, iSample * fTimeBinWidth, adcSamples[iSample] * fADCTomV);
    }
    
    grWaveform->SetTitle(Form("Waveform %i (Event %i, Channel %i); Time [ns]; Voltage [mV]", 
                             plotNum, eventId, channel));
    grWaveform->GetXaxis()->CenterTitle();
    grWaveform->GetYaxis()->CenterTitle();
    grWaveform->SetLineColor(kBlue);
    grWaveform->SetLineWidth(2);
    grWaveform->Draw("AL");
    
    // Add annotations
    c->Update();
    
    // Minimum marker
    TMarker* mMinimum = new TMarker(minTime, minVolt, 29);
    mMinimum->SetMarkerColor(kRed);
    mMinimum->SetMarkerSize(2);
    mMinimum->Draw("same");
    
    // Baseline
    TLine* lBaseline = new TLine(gPad->GetUxmin(), baseline, gPad->GetUxmax(), baseline);
    lBaseline->SetLineColor(kGreen);
    lBaseline->SetLineWidth(2);
    lBaseline->Draw("same");
    
    // Integration range
    float pregate = minTime - fPreGate * fTimeBinWidth;
    float gate = pregate + fGate * fTimeBinWidth;
    TLine* lPreGate = new TLine(pregate, gPad->GetUymin(), pregate, gPad->GetUymax());
    lPreGate->SetLineStyle(2);
    lPreGate->SetLineWidth(2);
    lPreGate->SetLineColor(kMagenta);
    lPreGate->Draw("same");
    
    TLine* lGate = new TLine(gate, gPad->GetUymin(), gate, gPad->GetUymax());
    lGate->SetLineStyle(2);
    lGate->SetLineWidth(2);
    lGate->SetLineColor(kMagenta);
    lGate->Draw("same");
    
    // Legend
    TLegend* leg = new TLegend(0.6, 0.15, 0.85, 0.4);
    leg->AddEntry(mMinimum, Form("Minimum (%.1f mV)", minVolt), "p");
    leg->AddEntry(lBaseline, Form("Baseline (%.1f mV)", baseline), "l");
    leg->AddEntry(lPreGate, Form("Integration range"), "l");
    leg->AddEntry((TObject*)0, Form("Charge: %.2f pC", charge), "");
    leg->Draw("same");
    
    c->Write();
}

/// @brief Draw sample waveforms from XML file
/// @param inName Input XML filename
/// @param outName Output ROOT filename for canvases
void DrawWaveformsXML(std::string inName, std::string outName) {
    std::cout << "Drawing " << fNumWaveformsToDraw << " waveforms from XML file " << inName << std::endl;
    
    // Populate tree structure pt
    using boost::property_tree::ptree;
    ptree pt;

    // Read the XML file into the property tree
    std::ifstream is(inName);
    read_xml(is, pt);
    
    // Read only the waveform size from XML settings
    int WSize = -1;
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer.settings")){
        if(v.first == "window"){
            WSize = v.second.get<int>("<xmlattr>.size", -1);
        }
    }

    // Create output ROOT file for canvases
    TFile *fOut = new TFile(outName.c_str(),"RECREATE");
    fOut->cd();

    // Random number generator for selecting events
    TRandom3 rand;
    rand.SetSeed(17);
    
    // Count total events and valid traces first
    int totalEvents = 0;
    int totalValidTraces = 0;
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer")){
        if(v.first == "event") {
            totalEvents++;
            for(auto& t : v.second){
                if(t.first == "trace"){
                    int channel = t.second.get<int>("<xmlattr>.channel", -1);
                    if(fChannelsToProcess.empty() || fChannelsToProcess.count(channel)) {
                        totalValidTraces++;
                    }
                }
            }
        }
    }
    
    std::cout << "Total events in file: " << totalEvents << std::endl;
    std::cout << "Total valid traces (matching channel filter): " << totalValidTraces << std::endl;
    
    // Generate random trace indices to draw
    std::vector<int> traceIndices;
    for(int i = 0; i < fNumWaveformsToDraw; i++) {
        traceIndices.push_back(rand.Integer(totalValidTraces));
    }
    std::sort(traceIndices.begin(), traceIndices.end());
    
    std::vector<int> adcSamples;
    int currentTraceIndex = 0;
    int plotIndex = 0;
    size_t nextTargetIndex = 0;
    
    // Loop over events and draw selected ones
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer")){
        if(v.first == "event"){
            int id = v.second.get<int>("<xmlattr>.id", -1);
            
            // Loop over all channels containing a waveform
            for(auto& t : v.second){
                if(t.first == "trace"){
                    int channel = t.second.get<int>("<xmlattr>.channel", -1);
                    
                    // Check channel filter
                    if(!fChannelsToProcess.empty() && !fChannelsToProcess.count(channel)) {
                        continue; // Skip this channel
                    }
                    
                    // Check if this trace should be drawn
                    if(nextTargetIndex < traceIndices.size() && currentTraceIndex == traceIndices[nextTargetIndex]) {
                        // Parse waveform data to vector
                        adcSamples.clear();
                        adcSamples.reserve(WSize);

                        std::stringstream ss(t.second.data());
                        std::string value;
        
                        // Parse space-separated ADC values
                        while(getline(ss, value, ' ')) {
                            if(!value.empty()) {
                                adcSamples.push_back(std::stoi(value));
                            }
                        }
                        
                        // Process and draw this waveform
                        DrawSingleWaveform(adcSamples, plotIndex + 1, channel, id);
                        plotIndex++;
                        nextTargetIndex++;
                        
                        if(plotIndex >= fNumWaveformsToDraw) break;
                    }
                    currentTraceIndex++;
                }
            }
            if(plotIndex >= fNumWaveformsToDraw) break;
        }
    }
    
    std::cout << "Drew " << plotIndex << " waveforms" << std::endl;
    fOut->Write();
    fOut->Close();
}

/// @brief Draw sample waveforms from DAT file
/// @param inName Input DAT filename
/// @param outName Output ROOT filename for canvases
void DrawWaveformsDAT(std::string inName, std::string outName) {
    std::cout << "Drawing " << fNumWaveformsToDraw << " waveforms from DAT file " << inName << std::endl;
    
    // Open .dat input file
    std::ifstream infile(inName, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open file: " << inName << std::endl;
        return;
    }

    // Determine waveform size and count valid traces
    uint32_t eventSize;
    infile.read(reinterpret_cast<char*>(&eventSize), sizeof(uint32_t));
    int WSize = (eventSize - 24) / 2;
    infile.seekg(0, std::ios::beg);
    
    // Count total valid traces
    int totalValidTraces = 0;
    uint32_t header[6];
    while (infile.read(reinterpret_cast<char*>(header), sizeof(header))) {
        infile.seekg(WSize * sizeof(uint16_t), std::ios::cur); // Skip waveform data
        int channel = static_cast<int>(header[3]);
        if(fChannelsToProcess.empty() || fChannelsToProcess.count(channel)) {
            totalValidTraces++;
        }
    }
    infile.clear();
    infile.seekg(0, std::ios::beg);
    
    std::cout << "Total valid traces (matching channel filter): " << totalValidTraces << std::endl;

    // Create output ROOT file for canvases
    TFile* fOut = new TFile(outName.c_str(), "RECREATE");
    fOut->cd();

    // Random number generator for selecting events
    TRandom3 rand;
    rand.SetSeed(17);
    
    // Generate random trace indices to draw
    std::vector<int> traceIndices;
    for(int i = 0; i < fNumWaveformsToDraw; i++) {
        traceIndices.push_back(rand.Integer(totalValidTraces));
    }
    std::sort(traceIndices.begin(), traceIndices.end());

    // Allocate space for reading
    uint16_t* rawSamples = new uint16_t[WSize];
    std::vector<int> adcSamples;
    adcSamples.reserve(WSize);
    
    int currentTraceIndex = 0;
    int plotIndex = 0;
    size_t nextTargetIndex = 0;

    // Read events and draw selected ones
    while (infile.read(reinterpret_cast<char*>(header), sizeof(header)) && plotIndex < fNumWaveformsToDraw) {
        if (!infile.read(reinterpret_cast<char*>(rawSamples), WSize * sizeof(uint16_t)))
            break;

        int channel = static_cast<int>(header[3]);
        
        // Check channel filter
        if(!fChannelsToProcess.empty() && !fChannelsToProcess.count(channel)) {
            continue; // Skip this channel
        }

        // Check if this trace should be drawn
        if(nextTargetIndex < traceIndices.size() && currentTraceIndex == traceIndices[nextTargetIndex]) {
            Long64_t timestamp = static_cast<Long64_t>(header[4]);
            Long64_t clocktime = static_cast<Long64_t>(header[5]);

            // Convert uint16_t array to std::vector<int> for processing
            adcSamples.clear();
            for (int i = 0; i < WSize; ++i) {
                adcSamples.push_back(static_cast<int>(rawSamples[i]));
            }

            // Process and draw this waveform
            DrawSingleWaveform(adcSamples, plotIndex + 1, channel, currentTraceIndex);
            plotIndex++;
            nextTargetIndex++;
        }
        currentTraceIndex++;
    }

    std::cout << "Drew " << plotIndex << " waveforms" << std::endl;
    fOut->Write();
    fOut->Close();
    delete[] rawSamples;
}

/// @brief Main drawing function that detects file type and draws accordingly
/// @param inName Input filename (XML or DAT)
/// @param outName Output ROOT filename
void DrawWaveforms(std::string inName, std::string outName) {
    std::string fileType = GetFileType(inName);
    
    if(fileType == "xml") {
        DrawWaveformsXML(inName, outName);
    } else if(fileType == "dat") {
        DrawWaveformsDAT(inName, outName);
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
    
    // Print drawing parameters
    std::cout << "Drawing parameters:" << std::endl;
    std::cout << "  Gate: " << fGate << " bins" << std::endl;
    std::cout << "  PreGate: " << fPreGate << " bins" << std::endl;
    std::cout << "  Number of waveforms: " << fNumWaveformsToDraw << std::endl;
    if(!fChannelsToProcess.empty()) {
        std::cout << "  Channels: ";
        for(auto it = fChannelsToProcess.begin(); it != fChannelsToProcess.end(); ++it) {
            if(it != fChannelsToProcess.begin()) std::cout << ", ";
            std::cout << *it;
        }
        std::cout << std::endl;
    } else {
        std::cout << "  Channels: all" << std::endl;
    }
    std::cout << std::endl;
    
    DrawWaveforms(inputFile, outputFile);
    return 0;
}
