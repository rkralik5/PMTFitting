// C++ Includes
#include <iostream>


// BOOST Includes
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>


// ROOT Includes
#include "TTree.h"
#include "TFile.h"





void PlotPoints(std::string line, float hz, std::vector<float> &px, std::vector<float> &py){
	int counter = 0;
		
    // Vector of string to save tokens 
    //std::vector <TVector2> points; 
	
      
    // stringstream class check1 
    std::stringstream check1(line); 
      
    //std::string intermediate; 
      std::string intermediates;
	  
    // Tokenizing w.r.t. space ' ' 
    while(getline(check1, intermediates, ' ')) 
    { 
//		TVector2 p;
		float intermediate = std::stof(intermediates);
		float x = (float)(counter)*1e9/hz;
		float y = (float)intermediate;
		//p.Set( x, y );
		px.push_back(x);
		py.push_back(y);
		
        //points.push_back(p); 
		counter++;
    } 

	//return points;
}


void Convert(std::string inName, std::string outName){
	// populate tree structure pt
    using boost::property_tree::ptree;
    ptree pt;
	std::ifstream is(inName);
    read_xml(is, pt);
	TFile *output = new TFile(outName.c_str(),"RECREATE");

	float frequency = -1;
	int maxSamples = -1;
	int resolution = -1;
	float voltLow = -1;
	float voltHigh = -1;
	std::vector<float> wavex;
	std::vector<float> wavey;
	int time = -1;
		

 	TTree *tDevice = new TTree("device","Device Settings");
    tDevice->Branch("frequency",&frequency);
	tDevice->Branch("maxSamples",&maxSamples);
	tDevice->Branch("resolution",&resolution);
	tDevice->Branch("voltLow",&voltLow);
	tDevice->Branch("voltHigh",&voltHigh);


	TTree *tWaves = new TTree("data","Wave Data");
	tWaves->Branch("wavex", &wavex);
	tWaves->Branch("wavey", &wavey);
	tWaves->Branch("time", &time);

	// Read in some basic digitiser information
	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("caendigitizer.digitizer") ) {
		
		//d.id = v.first.get("<xmlattr>.id", "DEFAULT_ID");
		
		if( v.first == "frequency" ) {
			frequency		= v.second.get<double>("<xmlattr>.hz", -1);
        }
		if( v.first == "maxsamples" ) {
	            maxSamples	= v.second.get<int>("<xmlattr>.maxsamples", -1);
        }
		if( v.first == "resolution" ) {
	            int resBits	= v.second.get<int>("<xmlattr>.bits", -1);
			resolution = pow(2, resBits);
        }
		if( v.first == "voltagerange" ) {
	            voltLow	= v.second.get<double>("<xmlattr>.low", -1);
	            voltHigh	= v.second.get<double>("<xmlattr>.hi", -1);
        }

    }
	tDevice->Fill();
	int c = 0;
	//std::vector<point> a;
	// Now we loop over the events
	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("caendigitizer") ) {
		c++;
		std::cout << "Run:\t" << c << std::endl;
		//if (c > 1003)break; 
		//d.id = v.first.get("<xmlattr>.id", "DEFAULT_ID");
		
		if( v.first == "event" ) {
			wavex.clear();
			wavey.clear();

			PlotPoints(v.second.get<std::string>("trace"), frequency, wavex, wavey);
			//ProcessWaveform5(a, totCharge, timePeak, frequency);
			tWaves->Fill();

        }
		
    }

	output->Write();
	output->Close();
 }

int main(int argc, char *argv[]){
	
	Convert(argv[1], argv[2]);

}
