///
/// @file MultiPEFit.C
/// @brief Fit a gaussian to the multi-PE integrated charge distribution
///

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric> // for std::accumulate

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "rootlogon.C"

#include "filename_parser.hpp"

#define PI 3.141592654
int fNChargeBins = 500; ///< Set the number of bins for the charge histogram
double e = 1.602e-19;
// TODO: #12 Change the gain to be loaded from somewhere or somethin else smart
float gain = 4973130; ///< Gain for this specific PMT

/// @brief Trivial function to check if a csv file is empty (and needs a header)
bool CheckEmptyFile(std::string outFileName);

////////////////////////////////////////////////////////////////////////////////
/// START OF MAIN FUNCTION
///
/// @brief Main fit function that loads waves, defines fitting functions, does the fits, plots the result and saves it to a csv file
/// @param inFileName (string) Name of the input ROOT file
/// @param outFileName (string) Name of the output csv file
void MultiPEFit(std::string inFileName, std::string outFileName="Output.csv",
	              int iChannel = -1){
	TFile* inFile = new TFile(inFileName.c_str(),"READ");
	if(inFile->IsZombie()){
		std::cerr << "Could not open file " << inFileName << std::endl;
		return;
	}

	// Get the measurement description from the filename
	ParsedFile parsedFile = parseFilename(inFileName);
	std::cout << "Analysing measurement: "
						<< parsedFile.description << std::endl;
	std::cout << "This measurement contains " << parsedFile.devices.size()
						<< " PMTs:" << std::endl;
	for(auto &device : parsedFile.devices){
		std::cout << device.manufacturer << " PMT " << device.model
							<< " at " << device.voltage <<"V";
		if(device.channel.has_value()){
			std::cout << " on channel " << device.channel.value();
		}
		std::cout << std::endl;
	}

	// If we specified channel, figure out which device it corresponds to
	FileInfo PMTInfo;
	if(iChannel != -1){
		for(auto &device : parsedFile.devices){
			if(device.channel.has_value() && device.channel.value() == iChannel){
				PMTInfo = device;
				break;
			}
		}
	} else { // If no channel is specified then use the first device
		PMTInfo = parsedFile.devices[0];
	}
	std::cout << "Using PMT " << PMTInfo.manufacturer << " "
						<< PMTInfo.model << " at " << PMTInfo.voltage
						<< "V on channel " << iChannel << std::endl;

	TTree *tWaves = (TTree*)inFile->Get("ProcessedData");
	Short_t Channel = -1; ///< Channel number
	Float_t Charge = 0; ///< Integrated charge in pC
	tWaves->SetBranchAddress("Channel", &Channel);
	tWaves->SetBranchAddress("Charge", &Charge);

	// Calculate the integrated charges
	int NEntries = tWaves->GetEntries();
	std::vector<float> vecCharge;
	vecCharge.reserve(NEntries);

	for (int i = 0; i < NEntries; i++){ // fill histogram
		//if(i%1000 == 0)
		//	std::cout << "Integrated:\t" << i/1000 << "k waveforms\r" << std::flush;
		tWaves->GetEntry(i);
		if(Channel != iChannel) continue; // Only get the channel we want
		vecCharge.push_back(Charge);
	}
	float vecAverage = accumulate(vecCharge.begin(),vecCharge.end(),0.0)/vecCharge.size();
	// Find the minimum and maximum to fill a histogram	
	auto chargeMinmax = minmax_element(vecCharge.begin(),vecCharge.end());
	float minimum = (float)*chargeMinmax.first;
	float maximum = (float)*chargeMinmax.second;
	std::cout << "minimum: " << minimum << " and maximum: " << maximum << std::endl;

	TH1F *hCharge = new TH1F(("charge_"+PMTInfo.model).c_str(),
													 ";Integrated Charge [pC];Area Normalized (arb. units)", fNChargeBins, minimum, maximum);
	hCharge->GetXaxis()->CenterTitle();
	hCharge->GetYaxis()->CenterTitle();
	for(auto iCharge : vecCharge) hCharge->Fill(iCharge);
	hCharge->Scale(1/hCharge->Integral());

	TCanvas c("c","c");
	c.cd();
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	// Fit with first Gaussian to get the range for final fit
	hCharge->Fit("gaus","EM");
	TF1 *pmt = hCharge->GetFunction("gaus");
	double mean = pmt->GetParameter(1);
	double sigma = pmt->GetParameter(2);
	double NPE_calc = pow(mean/sigma,2);

	// Now fit again with a Gaussian but only the peak to get the correct mean
	hCharge->Fit("gaus","EM","",mean-1.3*sigma,mean+sigma);
	hCharge->Draw("hist");
	//pmt->SetRange(minimum,maximum);
	pmt->Draw("same");

	// Number of Photo Electrons is simply mean/charge of 1 PE
	double NPE_gain = pmt->GetParameter(1)*1e-12/(gain*e);
	std::cout << "NPE: " << NPE_gain << std::endl;
	double NPE_calc2 = pow(vecAverage/pmt->GetParameter(2),2);

	c.Update();

	CornerLabel(PMTInfo.manufacturer+" "+PMTInfo.model);
	std::string plotLabel = Form("%s/%s_%s_%.2fV_%s.pdf",
															 parsedFile.path.c_str(),
															 PMTInfo.manufacturer.c_str(),
															 PMTInfo.model.c_str(), PMTInfo.voltage,
															 parsedFile.description.c_str());
	std::cout << "PlotLabel: " << plotLabel << std::endl;
	c.SaveAs(plotLabel.c_str());

	c.Clear();

	// Lets output this to some file
	std::ofstream outfile;
  outfile.open(outFileName.c_str(),
							 std::ios_base::app); // append instead of overwrite
	// Check if the file is empty and add the header if not
	if(CheckEmptyFile(outFileName)){
		outfile << "path,description,channel,pmt,model,voltage,NEntries,chisqr,charge_avg,charge_full,sigma_full,charge,sigma,npe_calc,npe_err,npe_calc2\n";
	}
	outfile << parsedFile.path << "," << parsedFile.description << ","
					<< PMTInfo.channel.value_or(-1) << "," << PMTInfo.manufacturer << ","
					<< PMTInfo.model << ","	<< PMTInfo.voltage << "," 
					<< "," << NEntries << "," << pmt->GetChisquare()/pmt->GetNDF() << ","
					<< vecAverage << "," << mean << "," << sigma << ","
					<< pmt->GetParameter(1) << "," << pmt->GetParameter(2) << ","
					<< NPE_calc << "," << 0.1*NPE_calc << "," << NPE_calc2 << "\n";
	outfile.close();
	vecCharge.clear();
}

///////////////////////////////////////////////////////////////////////////////
///  START OF FUNCTION DEFINITIONS
bool CheckEmptyFile(std::string outFileName){
	std::ifstream checkfile(outFileName.c_str());
  return checkfile.peek() == std::ifstream::traits_type::eof();
}