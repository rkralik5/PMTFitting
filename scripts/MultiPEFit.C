///
/// @file MultiPEFit.C
/// @brief Fit a gaussian to the multi-PE integrated charge distribution
///

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

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

#define PI 3.141592654
int fImpedance = 50; ///< Impedance in Ohm
int fNChargeBins = 500; ///< Set the number of bins for the charge histogram
double e = 1.602e-19;
int fPreGate = 6; ///< Number of time bins before peak position to start int
int fGate = 50; ///< Integration range for the waveform integration in time bins
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

/// @brief Get PMT label from the inFileName
/// @param inFileName Name of the input file with path and extension
/// @return PMT label
std::string GetPMTLabel(std::string inFileName);

/// @brief Main fit function that loads waves, defines fitting functions, does the fits, plots the result and saves it to a csv file
/// @param inFileName (string) Name of the input ROOT file
/// @param outFileName (string) Name of the output csv file
void MultiPEFit(std::string inFileName, std::string outFileName="Output.csv"){
	TFile* inFile = new TFile(inFileName.c_str(),"READ");
	// Get the parameters of the digitiser
	float ADCTomV;
	float timeBinWidth;
	GetParams(inFile, ADCTomV, timeBinWidth,
					  fResolution, fVoltLow, fVoltHigh, fFrequency);
	
	TTree *tWaves = (TTree*)inFile->Get("Data");
	TArrayS *Samples = new TArrayS; ///< Array of samples (waveform values)
	tWaves->SetBranchAddress("Samples", &Samples);

	// Calculate the integrated charges
	int NEntries = tWaves->GetEntries();
	std::vector<float> vecCharge;
	vecCharge.reserve(NEntries);
	for (int i = 0; i < NEntries; i++){ // fill histogram
		if(i%1000 == 0)
			std::cout << "Integrated:\t" << i/1000 << "k waveforms\r" << std::flush;
		tWaves->GetEntry(i);
		//vecCharge.push_back(IntegrateCharge(wavex, wavey, fPreGate, fGate));
		vecCharge.push_back(IntegrateCharge(Samples, ADCTomV, timeBinWidth,
																				fPreGate,	fGate));
	}

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

	// Now fit again with a Gaussian but only the peak to get the correct mean
	hCharge->Fit("gaus","EM","",mean-sigma,mean+sigma);
	hCharge->Draw("hist");
	pmt->Draw("same");

	// Number of Photo Electrons is simply mean/charge of 1 PE
	double NPE = pmt->GetParameter(1)*1e-12/(gain*e);
	std::cout << "NPE: " << NPE << std::endl;

	// Print the voltage
	std::string VLabel = inFileName.substr(0,inFileName.find_last_of("."));
	VLabel = VLabel.substr(0,VLabel.find_last_of("_"));
	VLabel = VLabel.substr(VLabel.find_last_of("_")+1);
	TText text(.7,.5,("Input: "+VLabel).c_str());
	text.SetNDC();
	text.Draw("same");

	c.Update();

	std::string PMTLabel = GetPMTLabel(inFileName);
	CornerLabel(PMTLabel);

	std::string plotFileName = inFileName.substr(0,inFileName.find_last_of("."))+"_Fit";
	c.SaveAs((plotFileName+"_Gate"+std::to_string(fGate)+".pdf").c_str());

	//c.SetLogy();
	//c.SaveAs((plotFileName+"_Logy"+"_Gate"+std::to_string(fGate)+".pdf").c_str());

	// Lets output this to some file
	std::ofstream outfile;
  outfile.open(outFileName.c_str(),
							 std::ios_base::app); // append instead of overwrite
	outfile << fGate << "," << inFileName << ","
					<< pmt->GetChisquare()/pmt->GetNDF() << ","
					<< pmt->GetParameter(1) << ","
					<< pmt->GetParameter(2) << ","
					<< NPE << ",";
	outfile.close();
	
	delete hCharge;
	delete pmt;
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

std::string GetPMTLabel(std::string inFileName){
	// Remove the path from inFileName
	std::string inFile = inFileName.substr(inFileName.find_last_of("/")+1);

	// Get everything before first _ into PMTBrand
  std::stringstream RestOfLine(inFile);
  std::string PMTBrand;
	getline(RestOfLine, PMTBrand, '_');
	
	// Second value before _ is the PMTLabel
	std::string PMTLabel;
	getline(RestOfLine, PMTLabel, '_');

	return PMTBrand+" "+PMTLabel;
}