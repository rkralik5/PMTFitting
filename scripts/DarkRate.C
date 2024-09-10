// TODO: #3 Add a description of the script
// TODO: #4 Make sure this works

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

#define NThresholds 20

/// @brief Check whether a waveform has a signal above a threshold
/// @param wy waveform y values in mV
/// @param threshold threshold value in mV
/// @return boolean whether the waveform has a signal above the threshold
bool IsAboveThreshold(std::vector<float> *&wy, int threshold);

void DarkRate(std::string inFileName, std::string outFileName="OutputDR.csv"){
	TFile inFile(inFileName.c_str(),"READ");
	TTree *tWaves = (TTree*)inFile.Get("data");
	std::vector<float> *wavey = 0;
	int clocktime = -1;
	tWaves->SetBranchAddress("wavey", &wavey);
	tWaves->SetBranchAddress("clocktime", &clocktime);

	std::vector<int> thresholds;
	for(int i=1; i<NThresholds+1; i++) thresholds.push_back(0.5*i);
	
	std::cout << "Entries:\t" << tWaves->GetEntries() << std::endl;
	int darkVals[NThresholds];
	std::vector<int> clocktimes[NThresholds];
	// Now we loop over the events
	for(int iWave = 0; iWave < tWaves->GetEntries(); iWave++){
		tWaves->GetEntry(iWave);
		for(int iThreshold = 0; iThreshold < NThresholds; iThreshold++){
			if(IsAboveThreshold(wavey, thresholds.at(iThreshold))){
				darkVals[iThreshold]++;
				clocktimes[iThreshold].push_back(clocktime);
			}
		}
		if(iWave%10000 == 0)
			std::cout << "Processed:\t" << iWave/10000 
								<< "0k waveforms\r" << std::flush;
	}

	std::ofstream outfile;
  outfile.open(outFileName.c_str(),
							 std::ios_base::app); // append instead of overwrite
	TH1F hClockTimes[NThresholds];
	for(int iThreshold = 0; iThreshold < NThresholds; iThreshold++){
		double darkRate = (double)darkVals[iThreshold]/(double)tWaves->GetEntries();
		std::cout << "Threshold" << thresholds[iThreshold]
							<< " with a Dark Rate = " << darkRate
							<< std::endl;
		outfile << thresholds[iThreshold] << "," << darkRate << "\n";
		
		// Create a TH1 histogram with 1s binning from clocktimes
		int minTime = *std::min_element(clocktimes[iThreshold].begin(), clocktimes[iThreshold].end());
		int maxTime = *std::max_element(clocktimes[iThreshold].begin(), clocktimes[iThreshold].end());
		int nBins = maxTime - minTime;
		std::cout << "minTime = " << minTime << " and maxTime = " << maxTime << std::endl;
		hClockTimes[iThreshold] = TH1F(Form("hClockTimes_%d",iThreshold),
																	 Form("Clock Times for Threshold %d",thresholds[iThreshold]),
																	 nBins, minTime, maxTime);
	}
	outfile.close();
	inFile.Close();
	delete tWaves;
}

bool IsAboveThreshold(std::vector<float> *&wy, int threshold){
	// Calculate the baseline as a truncated mean of 50% of values outside the gate
	std::sort(wy->begin(),wy->end());
	float baseline = 0;
	int nBaseline = 0;
	for(int iBln=wy->size()/4; iBln<(3/4)*wy->size(); iBln++){
		baseline += wy->at(iBln);
		nBaseline++;
	}
	baseline /= (float)nBaseline;

	float minimum = wy->back();
	if(baseline - minimum > threshold) return true;
	return false;
}
