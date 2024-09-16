/// @brief Script to calculate the dark rate of a channel using the ROOT files from waveconvert or from CoMPASS
/// @author: Robert Kralik kralikrobo@gmail.com (2024) @rkralik5

#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TArrayS.h>

#define NThresholds 20 ///< Number of thresholds to check
#define fNSecPerBin 60 ///< Number of seconds per bin

// If using CoMPASS outputs need to set these values manually
#define fResolution 500 ///< Number of ADC bins (=2^ADCResolution)
#define fVoltLow 0.f ///< Voltage range low
#define fVoltHigh 2.f ///< Voltage range high
#define fFrequency 500e6 ///< Frequency of the digitiser in Hz
#define fWindowSize 600 ///< Size of the waveform (number of time bins)

/// @brief Get the conversion factor of ADC bins to mV and the size of the waveform in seconds
/// @param inFile Input ROOT file - if it has a Device TTree, it will use the values from there
/// @param ADCTomV Return conversion factor from ADC bins to mV
/// @param waveSize Return size of the waveform in seconds
/// @param resolution Number of ADC bins (=2^ADCResolution)
/// @param voltLow Voltage range low
/// @param voltHigh Voltage range high
/// @param frequency Frequency of the digitiser in Hz
/// @param WindowSize Size of the waveform (number of time bins)
void GetParams(TFile *inFile, float &ADCTomV, float &waveSize,
							 int resolution = -1, float voltLow = -1, float voltHigh = -1, float frequency = -1, int WindowSize = -1);

/// @brief Subtract the baseline from a waveform and transform to vector
/// @param Samples TArrayS of waveform values
/// @param ADCTomV Conversion factor from ADC bins to mV
/// @return vector of floats of waveform values with the baseline subtracted
std::vector<float> BaselineCorrection(TArrayS *&Samples, float ADCTomV);

/// @brief Main script that calculates the dark rate of a channel
/// @param inFileName input ROOT file with waveforms
/// @param outFileName output file name with histograms and dark rates for each threshold
/// @param iChannel number of the channel to calculate the dark rate for
void DarkRate(std::string inFileName,
							std::string outROOTName="DarkRateOutput.root",
							std::string outCSVName="DarkRateOutput.csv",
							int iChannel = 0){
	std::cout << "Calculating the dark rate for channel " << iChannel << std::endl;
	
	TFile* inFile = new TFile(inFileName.c_str(),"READ");

	// Get the parameters of the digitiser
	float ADCTomV;
	float waveSize;
	GetParams(inFile, ADCTomV, waveSize,
					  fResolution, fVoltLow, fVoltHigh, fFrequency, fWindowSize);

	// Load the waveform TTree
	TTree *tWaves = (TTree*)inFile->Get("Data");
	Short_t Channel = -1; ///< Channel number
	int Clocktime = -1; ///< Clocktime of the event (in Unix time)
	TArrayS *Samples = new TArrayS; ///< Array of samples (waveform values)
	tWaves->SetBranchAddress("Channel", &Channel);
	tWaves->SetBranchAddress("Clocktime", &Clocktime);
	tWaves->SetBranchAddress("Samples", &Samples);

	// Initiate the thresholds
	std::vector<float> thresholds;
	for(int i=1; i<=NThresholds; i++) thresholds.push_back(-0.5*i);
	int darkVals[NThresholds] = {0};
	int NWaveforms = 0;

	// Get the time range of the data
	float beginning = 0;
	float end = 0;
	tWaves->GetEntry(0);
	beginning = Clocktime;	
	tWaves->GetEntry(tWaves->GetEntries()-1);
	end = Clocktime;
	int nBins = std::ceil((end - beginning)/fNSecPerBin); // 1 minute per bin
	std::cout << "Time range:\t" << beginning << " - " << end << std::endl;
	std::cout << "Number of bins:\t" << nBins << std::endl;
	// Create dark rate histograms for each threshold
	TH1F *hDarkRate[NThresholds];
	for(int iThreshold = 0; iThreshold < NThresholds; iThreshold++){
		hDarkRate[iThreshold] = new TH1F(Form("hDarkRate_%d",iThreshold),
																		 Form("Clock Times for Threshold %.1f mV;Unix Time [s];Dark Rate [Hz]",thresholds[iThreshold]),
																		 nBins, beginning, end);
	}
	// Create a histogram to count the number of waveforms per second
	TH1I hWaveforms("hWaveforms","Waveforms",nBins,beginning,end);
	
	// Now we loop over the events
	int NEvents = tWaves->GetEntries();
	//int NEvents = 1;
	std::cout << "Entries:\t" << NEvents/1000 << "k" << std::endl;
	for(int iWave = 0; iWave < NEvents; iWave++){
		tWaves->GetEntry(iWave);
		if(Channel != iChannel) continue;
		hWaveforms.Fill(Clocktime);
		NWaveforms++;

		std::vector<float> waveform = BaselineCorrection(Samples, ADCTomV);
		// Check if the waveform has a signal above the threshold for each threshold
		for(int iThreshold = 0; iThreshold < NThresholds; iThreshold++){
			bool IsAboveThreshold = false;
			for(auto point : waveform){
				if(point > thresholds.at(iThreshold)){
					IsAboveThreshold = true;
				}else if(IsAboveThreshold){
					IsAboveThreshold = false; //Only count once per peak below threshold
					darkVals[iThreshold]++;
					hDarkRate[iThreshold]->Fill(Clocktime);
				}
			}
		}
		if(iWave%10000 == 0)
			std::cout << "Processed:\t" << iWave/10000 
								<< "0k waveforms\r" << std::flush;
	}

	// Transform the histogram from counts to dark rate per second and save
	TFile *outFileROOT = new TFile(outROOTName.c_str(),"RECREATE");
	for(auto h : hDarkRate){
		h->Sumw2();
		h->Divide(&hWaveforms);
		h->Scale(1./waveSize);
		h->Write();
	}
	outFileROOT->Close();

	// Print the dark rates to a csv file
	std::ofstream outFileCSV;
  outFileCSV.open(outCSVName.c_str(),
							 		std::ios_base::app); // append instead of overwrite
	outFileCSV << outROOTName << ",";
	for(int iThreshold = 0; iThreshold < NThresholds; iThreshold++){
		double darkRate = (double)darkVals[iThreshold]/(NWaveforms*waveSize);
		double darkRateError = TMath::Sqrt(darkVals[iThreshold])/(NWaveforms*waveSize);
		outFileCSV << darkRate << "," << darkRateError << ",";
		std::cout << "Threshold " << thresholds[iThreshold]
							<< " with a Dark Rate = (" << darkRate
							<< " +/- " << darkRateError
							<< ") Hz, corresponding to " << darkVals[iThreshold]
							<< std::endl;
	}
	outFileCSV << std::endl;
	delete tWaves;
	inFile->Close();
}

void GetParams(TFile *inFile, float &ADCTomV, float &waveSize,
							 int resolution = -1, float voltLow = -1, float voltHigh = -1, float frequency = -1, int WindowSize = -1){
	TTree *tDevice = (TTree*)inFile->Get("Device");
	// Check if the TTree was loaded correctly
	if(tDevice == nullptr){
		if(resolution == -1 || voltLow == -1 || voltHigh == -1 ||
			 frequency == -1 || WindowSize == -1){
			std::cerr << "Could not load the Device TTree" << std::endl;
			std::cerr << "You have to input the values manually" << std::endl;
			return;
		}else{
			ADCTomV = (float)(voltHigh - voltLow)*1000/(float)resolution;
			waveSize = (float)WindowSize/frequency;
			return;
		}
	}

	tDevice->SetBranchAddress("resolution", &resolution);
	tDevice->SetBranchAddress("voltLow", &voltLow);
	tDevice->SetBranchAddress("voltHigh", &voltHigh);
	tDevice->SetBranchAddress("frequency", &frequency);
	tDevice->SetBranchAddress("WSize", &WindowSize);
	tDevice->GetEntry(0);
	ADCTomV = (voltHigh - voltLow)*1000/resolution;
	waveSize = WindowSize/frequency;
}

std::vector<float> BaselineCorrection(TArrayS *&Samples, float ADCTomV){
	// Convert waveform to vector
		std::vector<float> vecSamples;
		vecSamples.reserve(Samples->GetSize());
		for(int iSample = 0; iSample < Samples->GetSize(); ++iSample){
			vecSamples.push_back(ADCTomV*Samples->At(iSample));
		}
	// Sort waveform from smallest to largest values
	std::sort(vecSamples.begin(),vecSamples.end());
	// Calculate the baseline as a truncated mean of 50% of values outside the gate
	float baseline = 0;
	int nBaseline = 0;
	for(int iBln=vecSamples.size()/4; iBln<vecSamples.size()/4*3; iBln++){
		baseline += vecSamples.at(iBln);
		nBaseline++;
	}
	baseline /= (float)nBaseline;

	std::vector<float> returnVec;
	for(int iSample = 0; iSample < Samples->GetSize(); ++iSample){
		returnVec.push_back(ADCTomV*Samples->At(iSample) - baseline);
	}
	return returnVec;
}