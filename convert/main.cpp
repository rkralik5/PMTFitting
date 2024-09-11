// C++ Includes
#include <iostream>

// BOOST Includes
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

// ROOT Includes
#include "TTree.h"
#include "TFile.h"
#include "TArrayS.h"

/// @brief Convert output of CAENScopt waveform values to vectors
/// @param line String of space separated waveform values
/// @param Samples Output TArrayS of waveform values
void PlotPoints(std::string line, TArrayS &Samples){
	/// Which time bin
	int iBin = 0;
  /// Store the rest of the line that's yet to be processed
  std::stringstream RestOfLine(line);
  /// Hold each value
  std::string value;
  // Loop through space separated waveform bin values
  while(getline(RestOfLine, value, ' ')) Samples[iBin++] = std::stoi(value);
}

void Convert(std::string inName, std::string outName){
	// Populate tree structure pt
  using boost::property_tree::ptree;
  ptree pt;
	std::ifstream is(inName);
  read_xml(is, pt);
	TFile *fOut = new TFile(outName.c_str(),"RECREATE");

	float frequency = -1; ///< Frequency of the digitiser
	int NSamples = -1; ///< Total number of samples
	int WSize = -1; ///< Size of the waveform (number of time bins)
	int resolution = -1; ///< Number of ADC bins (=2^ADCResolution)
	float voltLow = -1; ///< Voltage range low
	float voltHigh = -1; ///< Voltage range high
 	TTree *tDevice = new TTree("Device","Device Settings");
  tDevice->Branch("frequency", &frequency, "Frequency/F");
	tDevice->Branch("NSamples",  &NSamples,  "NSamples/I");
	tDevice->Branch("WSize",     &WSize,     "WaveformSize/I");
	tDevice->Branch("resolution",&resolution,"TimeResolution/I");
	tDevice->Branch("voltLow",   &voltLow,   "VoltageLow/F");
	tDevice->Branch("voltHigh",  &voltHigh,  "VoltageHigh/F");

	// Read in some basic digitiser information
	BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer.digitizer")){
		if(v.first == "frequency"){
			frequency	= v.second.get<double>("<xmlattr>.hz", -1);
    }
		if(v.first == "maxsamples"){
	    NSamples = v.second.get<int>("<xmlattr>.maxsamples", -1);
    }
		if(v.first == "resolution"){
	    int resBits	= v.second.get<int>("<xmlattr>.bits", -1);
			resolution = pow(2, resBits);
    }
		if(v.first == "voltagerange"){
	    voltLow = v.second.get<double>("<xmlattr>.low", -1);
	    voltHigh = v.second.get<double>("<xmlattr>.hi", -1);
    }
	}
	BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer.settings")){
		if(v.first == "window"){
			WSize = v.second.get<int>("<xmlattr>.size", -1);
		}
	}
	tDevice->Fill(); // Only fill once per file

	short Channel = -1; ///< Channel number
	long long int Timestamp = -1; ///< Timestamp of the event from start of run
	long long int Clocktime = -1; ///< Clocktime of the event (in Unix time)
	TArrayS Samples; ///< Array of samples (waveform values)
	TTree *tWaves = new TTree("Data","Wave Data");
	tWaves->Branch("Channel",  &Channel,  "Channel/S");
	tWaves->Branch("Timestamp",&Timestamp,"Timestamp/I");
	tWaves->Branch("Clocktime",&Clocktime,"Clocktime/I");
	tWaves->Branch("Samples",  &Samples);

	// Loop over the events
	BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer")){
		if(v.first == "event"){
			int id = v.second.get<int>("<xmlattr>.id", -1);
			if(id%10000 == 0)
				std::cout << "Event:\t" << id/1000 << "k\r" << std::flush;

			Timestamp = v.second.get<int>("<xmlattr>.timestamp", -1);
			Clocktime = v.second.get<int>("<xmlattr>.clocktime", -1);
			// Loop over the channels containing a waveform
			for(auto& t : v.second){
				if(t.first == "trace"){
					Samples.Reset();
					Samples.Set(WSize);
					Channel = t.second.get<int>("<xmlattr>.channel", -1);
					PlotPoints(t.second.data(), Samples);
					tWaves->Fill();
				}
			}
		}
  }

	fOut->Write();
	fOut->Close();
 }

int main(int argc, char *argv[]){	
	Convert(argv[1], argv[2]);
}
