// C++ Includes
#include <iostream>

// BOOST Includes
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

// ROOT Includes
#include "TTree.h"
#include "TFile.h"

/// @brief Convert output of CAENScopt waveform values to vectors
/// @param line String of space separated waveform value
/// @param timeRes Resolution of the time bins
/// @param voltRes Resolution of the voltage bins
/// @param px Output vector of x values (time in ns)
/// @param py Output vector of y values (voltage in mV)
void PlotPoints(std::string line, float timeRes, float voltRes, std::vector<float> &px, std::vector<float> &py){
	/// Which time bin
	int timeBin = 0;
  /// Store the rest of the line that's yet to be processed
  std::stringstream RestOfLine(line);
  /// Hold each value
  std::string value;

  // Loop through space separated waveform bin values
  while(getline(RestOfLine, value, ' ')){
		//float val = std::stof(value);
		float time = (float)(timeBin)*timeRes;
		//float voltage = (float)val;
		float voltage = std::stof(value)*voltRes;
		px.push_back(time);
		py.push_back(voltage);
		timeBin++;
  }
}

void Convert(std::string inName, std::string outName){
	// Populate tree structure pt
  using boost::property_tree::ptree;
  ptree pt;
	std::ifstream is(inName);
  read_xml(is, pt);
	TFile *fOut = new TFile(outName.c_str(),"RECREATE");

	float frequency = -1;
	int NSamples = -1;
	int resolution = -1;
	float voltLow = -1;
	float voltHigh = -1;
	std::vector<float> wavex;
	std::vector<float> wavey;
	int clocktime = -1;

 	TTree *tDevice = new TTree("device","Device Settings");
  tDevice->Branch("frequency", &frequency, "Frequency/F");
	tDevice->Branch("NSamples",  &NSamples,  "NSamples/I");
	tDevice->Branch("resolution",&resolution,"TimeResolution/I");
	tDevice->Branch("voltLow",   &voltLow,   "VoltageLow/F");
	tDevice->Branch("voltHigh",  &voltHigh,  "VoltageHigh/F");

	TTree *tWaves = new TTree("data","Wave Data");
	tWaves->Branch("wavex",    &wavex);
	tWaves->Branch("wavey",    &wavey);
	tWaves->Branch("clocktime",&clocktime,"Clocktime/I");

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
	tDevice->Fill(); // Only fill once per file

	/// Convert time bins into units of ns
	float timeResolution = 1E9/frequency;
	/// Convert voltage bins into units of mV
	float voltResolution = (voltHigh - voltLow)*1000/resolution;

	int iEvent = 0;
	// Now we loop over the individual waveforms
	BOOST_FOREACH(ptree::value_type const& v, pt.get_child("caendigitizer")){
		iEvent++;
		if(iEvent%10000 == 0)
			std::cout << "Event:\t" << iEvent/1000 << "k\r" << std::flush;
		
		if(v.first == "event"){
			clocktime = v.second.get<int>("<xmlattr>.clocktime", -1);

			wavex.clear();
			wavey.clear();
			PlotPoints(v.second.get<std::string>("trace"), timeResolution, 		   voltResolution, wavex, wavey);
			tWaves->Fill(); // Fill for each event
    }
  }

	fOut->Write();
	fOut->Close();
 }

int main(int argc, char *argv[]){	
	Convert(argv[1], argv[2]);
}
