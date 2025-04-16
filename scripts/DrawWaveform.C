#include <vector>
#include <iostream>
#include <numeric> // for std::accumulate
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TArrayS.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TROOT.h"

#include "filename_parser.hpp"
#include "rootlogon.C"

int fImpedance = 50; ///< Impedance in Ohm
int fNChargeBins = 500; ///< Set the number of bins for the charge histogram
double e = 1.602e-19;
int fPreGate = 6; ///< Number of time bins before peak position to start int
int fGate = 150; ///< Integration range for the waveform integration in time bins

// Values to estimate the initial shape of the pedestal
float fThreshold = 0.5; ///< Threshold for pedestal in mV (not used in fit)
int fPeakTime = 390; ///< Approx. time of the peak in ns (not used in fit)

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
	int preGate, int gate, float& minTime, float& minVolt, float& baseline);

// TODO: #29 Make sure the DrawWaveform function works as expected
/// @brief Draw a random selection of waveforms from a ROOT file
/// @param inFileName Name of the input ROOT file
/// @param outFileName Name of the output ROOT file. Default SampleWaveforms.root
/// @param NPlots Int number of waveforms to be drawn. Default is 5
void DrawWaveform(std::string inFileName,
								  std::string outFileName="SampleWaveforms.root",
									int NPlots=5){
	TFile inFile(inFileName.c_str(),"READ");

	// Get the parameters of the digitiser
	float ADCTomV;
	float timeBinWidth;
	GetParams(&inFile, ADCTomV, timeBinWidth,
						fResolution, fVoltLow, fVoltHigh, fFrequency);
		
	TTree *tWaves = (TTree*)inFile.Get("Data");
	TArrayS *Samples = new TArrayS; ///< Array of samples (waveform values)
	tWaves->SetBranchAddress("Samples", &Samples);

	TFile outFile(outFileName.c_str(),"RECREATE");
	outFile.cd();

	// Get a random waveform and plot it with a TGraph
	int NEntries = tWaves->GetEntries();
	TRandom3 rand;
	for (int iPlot = 1; iPlot < NPlots; iPlot++){ // fill histogram
		int iEntry = rand.Integer(NEntries);
		tWaves->GetEntry(iEntry);
		float mintime = 0; // get time of peak voltage in ns
		float minvolt = 0; // get peak voltage in mV
		float baseline = 0; // get baseline voltage in mV
		float charge = IntegrateCharge(Samples, ADCTomV, timeBinWidth,
																	 fPreGate,	fGate, mintime, minvolt,
																	 baseline);
    
		new TCanvas;
		TGraph grWaveform(Samples->GetSize());
		for(int iSample = 0; iSample < Samples->GetSize(); ++iSample){
			grWaveform.AddPoint(iSample*timeBinWidth,ADCTomV*Samples->At(iSample));
		}

		grWaveform.SetTitle(Form("Waveform_%i",iPlot));
		grWaveform.GetXaxis()->SetTitle("Time [ns]");
		grWaveform.GetXaxis()->CenterTitle();
		grWaveform.GetYaxis()->SetTitle("Output voltage [mV]");
		grWaveform.GetYaxis()->CenterTitle();
		grWaveform.SetName(Form("Waveform_%i",iPlot));
		//grWaveform.GetXaxis()->SetRangeUser(wavex->front(),wavex->back());
		grWaveform.Draw();
		grWaveform.Write();
	}
}

void GetParams(TFile *inFile, float &ADCTomV, float &timeBinWidth,
							 int resolution = -1, float voltLow = -1, float voltHigh = -1,
							 float frequency = -1){
	TTree *tDevice = (TTree*)inFile->Get("Device");
	// Check if the TTree was loaded correctly
	if(tDevice == nullptr){
		std::cerr << "Could not load the Device TTree" << std::endl;
		std::cerr << "Using the manually inputed values!" << std::endl;
		if(resolution == -1 || voltLow == -1 || voltHigh == -1 ||	frequency == -1){
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
	int preGate, int gate, float& minTime, float& minVolt,
	float& baseline){
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

// Minimum time is just the minPos multiplied by the time resolution
minTime = minPos*timeBinWidth;
minVolt = vecSamples.at(minPos);

// Get the integration range
int lowInt = minPos>preGate ? minPos - preGate : 0;
int highInt = lowInt+gate;
if(highInt >= vecSamples.size()){
highInt = vecSamples.size()-1;
lowInt = vecSamples.size()-1-gate;
}
/*
// Make a vector for baseline calculation
std::vector<float> vBaseline;
vBaseline.reserve(vecSamples.size()-gate);
for(auto i=0; i<vecSamples.size(); i++){
if(i>=lowInt && i<highInt) continue;
vBaseline.push_back(vecSamples.at(i));
}

// Calculate the baseline as a truncated mean of 50% of values outside the gate
std::sort(vBaseline.begin(),vBaseline.end());
baseline = 0;
int nBaseline = 0;
for(int iBln=vBaseline.size()/4; iBln<(vBaseline.size()-vBaseline.size()/4); iBln++){
baseline += vBaseline.at(iBln);
nBaseline++;
}
baseline /= (float)nBaseline;
*/

// 1st baseline was with fPeakTime-40
baseline = accumulate(vecSamples.begin()+5,
				 vecSamples.begin()+fPeakTime-100,0.0)/
				(fPeakTime-105);

// Correct the minVolt by subtracting the baseline
minVolt = baseline-minVolt;

// Finally integrate the output voltages into a charge in pC
float charge = 0;
for(int iSample = lowInt; iSample < highInt; ++iSample){
charge += baseline - vecSamples.at(iSample);
}
charge = charge*timeBinWidth/fImpedance; // Charge is voltage*time/impedance
return charge;
}