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
int fImpedance = 50; ///< Impedance in Ohm
int fNChargeBins = 500; ///< Set the number of bins for the charge histogram
double e = 1.602e-19;
int fPreGate = 6; ///< Number of time bins before peak position to start int
int fGate = 150; ///< Integration range for the waveform integration in time bins
// TODO: #12 Change the gain to be loaded from somewhere or somethin else smart
float gain = 4973130; ///< Gain for this specific PMT

// If using CoMPASS outputs need to set these values manually
#define fResolution 500 ///< Number of ADC bins (=2^ADCResolution)
#define fVoltLow 0.f ///< Voltage range low
#define fVoltHigh 2.f ///< Voltage range high
#define fFrequency 500e6 ///< Frequency of the digitiser in Hz
#define fWindowSize 600 ///< Size of the waveform (number of time bins)

/// @brief Get the conversion factor of ADC bins to mV and the size of the waveform in seconds
/// @param inFile Input ROOT file - if it has a Device TTree, it will use the values from there
/// @param ADCTomV Return conversion factor from ADC bins to mV
/// @param timeBinWidth Return width of the time bin in nanoseconds
/// @param resolution Number of ADC bins (=2^ADCResolution)
/// @param voltLow Voltage range low
/// @param voltHigh Voltage range high
/// @param frequency Frequency of the digitiser in Hz
/// @param WindowSize Size of the waveform (number of time bins)
void GetParams(TFile *inFile, float &ADCTomV, float &timeBinWidth,
							 int resolution = -1, float voltLow = -1, float voltHigh = -1, float frequency = -1);

/// @brief Integrate waveform to get a deposited charge
/// @param Samples Input TArrayS of waveform values (output from waveconvert)
/// @param ADCTomV Conversion factor from ADC bins to mV
/// @param timeBinWidth Width of the time bin in nanoseconds
/// @param preGate Length of time to integrate before the peak in bins of time. Default 5
/// @param gate Length of time to integrate from preGate in bins of time. Default 50
/// @return Inegrated charge in pC
float IntegrateCharge(TArrayS *&Samples, float ADCTomV, float timeBinWidth,
										  int preGate=5, int gate=50);

/// @brief Trivial funciton to check if a csv file is empty (and needs a header)
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

	// Get the parameters of the digitiser
	float ADCTomV;
	float timeBinWidth;
	GetParams(inFile, ADCTomV, timeBinWidth,
					  fResolution, fVoltLow, fVoltHigh, fFrequency);
	
	TTree *tWaves = (TTree*)inFile->Get("Data");
	Short_t Channel = -1; ///< Channel number
	TArrayS *Samples = new TArrayS; ///< Array of samples (waveform values)
	tWaves->SetBranchAddress("Channel", &Channel);
	tWaves->SetBranchAddress("Samples", &Samples);

	// Calculate the integrated charges
	int NEntries = tWaves->GetEntries();
	std::vector<float> vecCharge;
	vecCharge.reserve(NEntries);
	for (int i = 0; i < NEntries; i++){ // fill histogram
		if(i%1000 == 0)
			std::cout << "Integrated:\t" << i/1000 << "k waveforms\r" << std::flush;
		tWaves->GetEntry(i);
		if(Channel != iChannel) continue; // Only get the channel we want
		vecCharge.push_back(IntegrateCharge(Samples, ADCTomV, timeBinWidth,
																				fPreGate,	fGate));
	}
	std::cout << std::endl;
	float vecAverage = accumulate(vecCharge.begin(),vecCharge.end(),0.0)/vecCharge.size();
	std::cout << "Vector average is " << vecAverage << std::endl;

	// Find the minimum and maximum to fill a histogram	
	auto chargeMinmax = minmax_element(vecCharge.begin(),vecCharge.end());
	float minimum = (float)*chargeMinmax.first;
	float maximum = (float)*chargeMinmax.second;

	TH1F *hCharge = new TH1F("charge",
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

	//TODO: #26 Adapt this code to actual print the PMT type, label, voltage as columns
	std::string plotLabel = parsedFile.path + "/"
		+ PMTInfo.manufacturer + "_" + PMTInfo.model + "_"
		+ std::to_string(PMTInfo.voltage) + "V"
		+ parsedFile.description + ".pdf";
	c.SaveAs(plotLabel.c_str());

	// Lets output this to some file
	std::ofstream outfile;
  outfile.open(outFileName.c_str(),
							 std::ios_base::app); // append instead of overwrite
	// Check if the file is empty and add the header if not
	if(CheckEmptyFile(outFileName)){
		outfile << "path,description,channel,pmt,model,voltage,intrange,chisqr,charge_avg,charge_full,sigma_full,charge,sigma,npe_calc,npe_err,npe_calc2\n";
	}
	outfile << parsedFile.path << "," << parsedFile.description << ","
					<< PMTInfo.channel.value_or(-1) << "," << PMTInfo.manufacturer << ","
					<< PMTInfo.model << ","	<< PMTInfo.voltage << "," 
					<< fGate << "," << pmt->GetChisquare()/pmt->GetNDF() << ","
					<< vecAverage << "," << mean << "," << sigma << ","
					<< pmt->GetParameter(1) << "," << pmt->GetParameter(2) << ","
					<< NPE_calc << "," << 0.1*NPE_calc << "," << NPE_calc2 << "\n";
	outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
///  START OF FUNCTION DEFINITIONS
void GetParams(TFile *inFile, float &ADCTomV, float &timeBinWidth,
							 int resolution = -1, float voltLow = -1, float voltHigh = -1, float frequency = -1){
	TTree *tDevice = (TTree*)inFile->Get("Device");
	// Check if the TTree was loaded correctly
	if(tDevice == nullptr){
		if(resolution == -1 || voltLow == -1 || voltHigh == -1 ||
			 frequency == -1){
			std::cerr << "Could not load the Device TTree" << std::endl;
			std::cerr << "You have to input the values manually" << std::endl;
			return;
		}else{
			ADCTomV = (float)(voltHigh - voltLow)*1000/(float)resolution;
			timeBinWidth = 1./frequency;
			return;
		}
	}

	tDevice->SetBranchAddress("resolution", &resolution);
	tDevice->SetBranchAddress("voltLow", &voltLow);
	tDevice->SetBranchAddress("voltHigh", &voltHigh);
	tDevice->SetBranchAddress("frequency", &frequency);
	tDevice->GetEntry(0);
	ADCTomV = (voltHigh - voltLow)*1000/resolution;
	timeBinWidth = 1e9/frequency;
}

float IntegrateCharge(TArrayS *&Samples, float ADCTomV, float timeBinWidth,
										  int preGate, int gate){
	// Convert waveform to vector for easier manipulation
	std::vector<float> vecSamples;
	vecSamples.reserve(Samples->GetSize());
	for(int iSample = 0; iSample < Samples->GetSize(); ++iSample){
		vecSamples.push_back(ADCTomV*Samples->At(iSample));
	}

	// Find the minimum voltage
	// Don't look within a 5 bin-wide buffer in beginning and end, which would 
	// have incomplete signal
	int minPos = std::distance(vecSamples.begin(),
														 std::min_element(vecSamples.begin()+5,
														 									vecSamples.end()-5));
	float minCharge = vecSamples.at(minPos);

	// Get the integration range
	int lowInt = minPos>preGate ? minPos - preGate : 0;
	int highInt = lowInt+gate;
	if(highInt >= vecSamples.size()){
		highInt = vecSamples.size()-1;
		lowInt = vecSamples.size()-1-gate;
	}

	// Make a vector for baseline calculation
	std::vector<float> vBaseline;
	vBaseline.reserve(vecSamples.size()-gate);
	for(auto i=0; i<vecSamples.size(); i++){
		if(i>=lowInt && i<highInt) continue;
		vBaseline.push_back(vecSamples.at(i));
	}

	// Calculate the baseline as a truncated mean of 50% of values outside the gate
	std::sort(vBaseline.begin(),vBaseline.end());
	float baseline = 0;
	int nBaseline = 0;
	for(int iBln=vBaseline.size()/4; iBln<(vBaseline.size()-vBaseline.size()/4); iBln++){
		baseline += vBaseline.at(iBln);
		nBaseline++;
	}
	baseline /= (float)nBaseline;

	// Finally integrate the output voltages into a charge in pC
	float charge = 0;
	for(int iSample = lowInt; iSample < highInt; ++iSample){
		charge += baseline - vecSamples.at(iSample);
	}
	charge = charge*timeBinWidth/fImpedance; // Charge is voltage*time/impedance

	return charge;
}