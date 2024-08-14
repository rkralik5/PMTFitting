#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TRandom3.h"

/// @brief Draw a random selection of waveforms from a ROOT file
/// @param inFileName Name of the input ROOT file
/// @param outFileName Name of the output ROOT file. Default SampleWaveforms.root
/// @param NPlots Int number of waveforms to be drawn. Default is 5
void DrawWaveform(std::string inFileName,
								  std::string outFileName="SampleWaveforms.root",
									int NPlots=5){
	TFile inFile(inFileName.c_str(),"READ");
	TTree *tDevice = (TTree*)inFile.Get("device");
	TTree *tWaves = (TTree*)inFile.Get("data");

	float frequency;
	int maxSamples;
	int resolution;
	float voltLow;
	float voltHigh;
	std::vector<float> *wavex = 0;
	std::vector<float> *wavey = 0;
	int clocktime;

	tDevice->SetBranchAddress("frequency",&frequency);
	tDevice->SetBranchAddress("NSamples",&maxSamples);
	tDevice->SetBranchAddress("resolution",&resolution);
	tDevice->SetBranchAddress("voltLow",&voltLow);
	tDevice->SetBranchAddress("voltHigh",&voltHigh);
	tDevice->GetEntry(0); // Device information only set once

	tWaves->SetBranchAddress("wavex", &wavex);
	tWaves->SetBranchAddress("wavey", &wavey);
	tWaves->SetBranchAddress("clocktime", &clocktime);

	TFile outFile(outFileName.c_str(),"RECREATE");
	outFile.cd();

	/// Graph to hold the waveforms
	

	// Get a random waveform and plot it with a TGraph
	int NEntries = tWaves->GetEntries();
	TRandom3 rand;
	for (int iPlot = 1; iPlot < NPlots; iPlot++){ // fill histogram
		int iEntry = rand.Integer(NEntries);
		tWaves->GetEntry(iEntry);
    
		new TCanvas;
		// Need to conver std::vector to array by calling 1st element
		TGraph* grWaveform = new TGraph(wavex->size(), &(wavex->at(0)), &(wavey->at(0)));
		grWaveform->SetTitle(Form("Waveform_%i",iEntry));
		grWaveform->GetXaxis()->SetTitle("Time [ns]");
		grWaveform->GetXaxis()->CenterTitle();
		grWaveform->GetYaxis()->SetTitle("Output voltage [mV]");
		grWaveform->GetYaxis()->CenterTitle();
		grWaveform->SetName(Form("Waveform_%i",iEntry));
		grWaveform->GetXaxis()->SetRangeUser(wavex->front(),wavex->back());
		grWaveform->Draw();
		grWaveform->Write();
	}
}