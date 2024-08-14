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

#define PI 3.141592654
int fImpedance = 50; ///< Impedance in Ohm
int fNChargeBins = 500; ///< Set the number of bins for the charge histogram

Double_t factorial(int a){
	if(a > 1){return a*factorial(a-1);} else return 1;
}

Double_t Pois(double mu, double n){
	return pow(mu,n)*exp(-mu)/factorial(n);
}

Double_t Gaus(double x, double q, double s, double n){
	return 1/(s*sqrt(2*PI*n))*exp(-pow(x-n*q,2)/(2*n*pow(s,2)));
}

Double_t Ignxe( double x, double q0, double s0, double q1, double s1, double a, double n ){
	double qn = q0 + n*q1;
	double sn = sqrt(pow(s0,2) + n*pow(s1,2));
	
	return (a/2)*exp( -a*(x-qn-0.5*a*pow(sn,2)) ) * (ROOT::Math::erf( fabs(q0 - qn - pow(sn,2)*a)/(sn*sqrt(2)) ) + TMath::Sign(1,x - qn - pow(sn,2)*a )*ROOT::Math::erf( fabs(x - qn - pow(sn,2)*a)/(sn*sqrt(2))));
}

Double_t PMTF(Double_t *x_, Double_t *par){
	
	double x = x_[0];
	double q0 = par[0];
	double s0 = par[1];
	double q1 = par[2];
	double s1 = par[3];
	double w = par[4];
	double a = par[5];
	double mu = par[6];

	return Pois(mu, 0)*((1-w)*Gaus(x,q0,s0,1) + w*Ignxe(x, q0, s0, q1, s1, a, 0))
		   + Pois(mu, 1)*((1-w)*Gaus(x,q1,s1,1) + w*Ignxe(x, q0, s0, q1, s1, a, 1))
			 + Pois(mu, 2)*((1-w)*Gaus(x,q1,s1,2) + w*Ignxe(x, q0, s0, q1, s1, a, 2))
		   + Pois(mu, 3)*((1-w)*Gaus(x,q1,s1,3) + w*Ignxe(x, q0, s0, q1, s1, a, 3))
		   + Pois(mu, 4)*((1-w)*Gaus(x,q1,s1,4) + w*Ignxe(x, q0, s0, q1, s1, a, 4));
}

Double_t PMTF0(Double_t *x_, Double_t *par){
	double x = x_[0];
	double q0 = par[0];
	double s0 = par[1];
	double q1 = par[2];
	double s1 = par[3];
	double w = par[4];
	double a = par[5];
	double mu = par[6];

	return Pois(mu, 0)*((1-w)*Gaus(x,q0,s0,1) + w*Ignxe(x, q0, s0, q1, s1, a, 0));
}

Double_t PMTF1(Double_t *x_, Double_t *par){
	double x = x_[0];
	double q0 = par[0];
	double s0 = par[1];
	double q1 = par[2];
	double s1 = par[3];
	double w = par[4];
	double a = par[5];
	double mu = par[6];

	return Pois(mu, 1)*((1-w)*Gaus(x,q1,s1,1) + w*Ignxe(x, q0, s0, q1, s1, a, 1));
}

Double_t PMTF2(Double_t *x_, Double_t *par){
	double x = x_[0];
	double q0 = par[0];
	double s0 = par[1];
	double q1 = par[2];
	double s1 = par[3];
	double w = par[4];
	double a = par[5];
	double mu = par[6];

	return Pois(mu, 2)*((1-w)*Gaus(x,q1,s1,2) + w*Ignxe(x, q0, s0, q1, s1, a, 2));
}

void PrintVector(std::vector<float> v, std::string s = "Printing Vector"){
	std::cout << "===================" << s << "============" ;
	for (int i = 0; i < v.size(); i++){
		std::cout << v[i] << std::endl;
	}
}

void WriteVector(std::vector<float> v, std::string s){
	ofstream myfile;
	myfile.open (s);
	
	for (int i = 0; i < v.size(); i++){
		myfile << v[i] << "\n";
	}
	myfile.close();
}

/// @brief Integrate waveform to get a deposited charge
/// @param wx Input vector of waveform x values (time in ns)
/// @param wy Input vector of waveform y values (voltage in mV)
/// @param preGate Length of time to integrate before the peak in bins of time. Default 5
/// @param gate Length of time to integrate from preGate in bins of time. Default 50
/// @return Inegrated charge in pC
float IntegrateCharge(std::vector<float> *wx, std::vector<float> *wy,
										  int preGate=5, int gate=50);

void FixFit(TF1 *p, TF1 *a ){
	p->FixParameter(0, a->GetParameter(0));
	p->FixParameter(1, a->GetParameter(1));
	p->FixParameter(2, a->GetParameter(2));
	p->FixParameter(3, a->GetParameter(3));
	p->FixParameter(4, a->GetParameter(4));
	p->FixParameter(5, a->GetParameter(5));
	p->FixParameter(6, a->GetParameter(6));
}

void Fit(std::string inFileName, std::string outFileName=""){
	TFile inFile(inFileName.c_str(),"READ");
	TTree *tWaves = (TTree*)inFile.Get("data");

	// If no output filename simply append FitOutput to the input file name
	if(outFileName == "")
		outFileName = inFileName.substr(0,inFileName.find_last_of("."))+
									"_FitOutput.root";
	TFile outFile(outFileName.c_str());
/*
	TH1F *hCharge = new TH1F("ChargeHistogram", "");
	// Get the histogram of charges, either load or calculate
	if(outFile.GetListOfKeys()->Contains("ChargeHistogram"))
		hCharge = (TH1F*)outFile.Get("ChargeHistogram");
	else{}*/

	std::vector<float> *wavex = 0;
	std::vector<float> *wavey = 0;
	tWaves->SetBranchAddress("wavex", &wavex);
	tWaves->SetBranchAddress("wavey", &wavey);

	// Calculate the integrated charges
	int NEntries = tWaves->GetEntries();
	//int NEntries = 10000;
	std::vector<float> vecCharge;
	vecCharge.reserve(NEntries);
	for (int i = 0; i < NEntries; i++){ // fill histogram
		if(i%1000 == 0)
			std::cout << "Integrated:\t" << i/1000 << "k waveforms\r" << std::flush;
		tWaves->GetEntry(i);
		vecCharge.push_back(IntegrateCharge(wavex, wavey));
	}

	// Find the minimum and maximum to fill a histogram
	auto chargeMinmax = minmax_element(vecCharge.begin(),vecCharge.end());
	auto chargeMax = max_element(vecCharge.begin(),vecCharge.end());
	std::cout << "minimum and maximum charge is " << *chargeMinmax.first
						<< " and maximum is " << *chargeMinmax.second
						<< " and alternative maximum is " << *chargeMax
						<<std::endl;

	float min = *chargeMinmax.first;
	float max = 3; //*chargeMinmax.second;
	
	TH1F *h1 = new TH1F("charge", "charge", fNChargeBins,
											min, max);
	h1->GetXaxis()->SetTitle("Integrated Charge [pC]");
	for(auto iCharge : vecCharge) h1->Fill(iCharge);
  
	h1->Scale(1/h1->Integral());

	//Use TSpectrum to find the peak candidates
  TSpectrum *s = new TSpectrum(2);
  Int_t nfound = s->Search(h1,0.1,"",0.00001);
	std::cout << "Found " << nfound << " peaks" << std::endl;

	Double_t *xpeaks = s->GetPositionX();
	double peaksX[nfound];
	double peaksY[nfound];
	std::fill(peaksX, peaksX + nfound, -1); // fill array with -1
	std::fill(peaksY, peaksY + nfound, -1); // fill array with -1

  for (int p=0;p < nfound;p++) {
    Double_t xp = xpeaks[p];
	  std::cout << "peak " << p << "  " << xp << std::endl;
		peaksX[p] = xp;
		peaksY[p] = h1->GetBinContent(h1->GetXaxis()->FindBin(xp));
	}

	double q0 = (peaksX[0] > 0) ? peaksX[0] : 0; // if first peak not found then set it to zero
	double q1 = (peaksX[1] > 0) ? peaksX[1] : 0.4; // if second peak not found then set it to ten

	// Let's find the valley now
	double valley[2] = {-1, -1};
	if (nfound == 2){ // two peaks found lets look for valley
		int startBin = h1->GetXaxis()->FindBin(peaksX[0]);
		int endBin = h1->GetXaxis()->FindBin(peaksX[1]);

		double minValX = -1;
		double minValY = 99999999;
		std::cout << "BIN DEBUG" << std::endl;
		std::cout << startBin << "  " << endBin << std::endl;
		for (int i = startBin; i < endBin; i++){
			double val = h1->GetBinContent(i);
			if (val < minValY){
				minValY = val;
				minValX = h1->GetXaxis()->GetBinCenter(i);
			}
		}
		valley[0] = minValX;
		valley[1] = minValY;
	}

	std::cout << "VALLEY" << std::endl;
	std::cout << valley[0] << "   " << valley[1] << std::endl;

	TF1 *pmt = new TF1("pmt",PMTF,min,max,7); // x in [0;10], 2 parameters
	TH1F *testhi = new TH1F("testhi","Poisson distribution",fNChargeBins,min,max);
	TF1 *pmt0 = new TF1("pmt0",PMTF0,min,max,7); // x in [0;10], 2 parameters
	TF1 *pmt1 = new TF1("pmt1",PMTF1,min,max,7); // x in [0;10], 2 parameters
	TF1 *pmt2 = new TF1("pmt2",PMTF2,min,max,7); // x in [0;10], 2 parameters
    //give the parameters meaningful names
  pmt->SetParNames("Q_{0}","#sigma_{0}","Q_{1}","#sigma_{1}", "W", "a", "#mu");
	
	pmt->SetParameter(0,q0);
	pmt->SetParameter(1,0.2);
	pmt->SetParameter(2,q1);
	pmt->SetParameter(3,0.5);
	pmt->SetParameter(4,0.4999);
	pmt->SetParameter(5,0.00001379);
	pmt->FixParameter(6,0.0585);

	pmt->SetParLimits(0,min,max);
	pmt->SetParLimits(1,-0.5,0.5);
	pmt->SetParLimits(2,min,max);
	pmt->SetParLimits(3,0,10);
	pmt->SetParLimits(4,0.,1);
	pmt->SetParLimits(5,0.,1);
	pmt->SetParLimits(6,0.001,1.5000);
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	h1->Fit("pmt");

	FixFit(pmt0, pmt );
	FixFit(pmt1, pmt );
	FixFit(pmt2, pmt );

	pmt0->Draw("SAME");
	pmt1->Draw("SAME");
	pmt2->Draw("SAME");

	pmt0->SetLineColor(kBlue);
	pmt1->SetLineColor(kGray);
	pmt2->SetLineColor(kMagenta);

	TCanvas *c2 = new TCanvas("hist1","hist1");
	h1->Draw("HIST");

	double charge = pmt->GetParameter(2);
	double e = 1.6e-19;
	double gain = charge/e;

	double peRes = pmt->GetParameter(3)/pmt->GetParameter(2); // s1/q1
	double peak2valley = peaksY[1]/valley[1]; // height at q1/ height at valley

	std::cout << "Q1 = " << charge << "pC" << std::endl;
	std::cout << "Which is equal to a gain of " << gain << std::endl;

	std::cout << "Pedastal " << h1->Integral(h1->FindFixBin(-50),h1->FindFixBin(10)) << std::endl;
	std::cout << "Signal " << h1->Integral(h1->FindFixBin(11),h1->FindFixBin(100)) << std::endl;
	std::cout << "PE Res " << peRes << std::endl;
	std::cout << "P2V " << peak2valley << std::endl;
	std::cout << pmt->GetParameter(0) << "\t" << pmt->GetParameter(1) << "\t" << pmt->GetParameter(2) << "\t" << pmt->GetParameter(6) << "\t" << pmt->GetParameter(3) << std::endl;

	// Lets output this to some file
	std::ofstream outfile;
	std::ofstream outfileParts;

  outfile.open("pmt.txt", std::ios_base::app); // append instead of overwrite

	outfile << inFileName << "\t" << pmt->GetParameter(0) << "\t" << pmt->GetParameter(1) << "\t" << pmt->GetParameter(2) << "\t" << pmt->GetParameter(3) << "\t" << gain << "\t" << peRes << "\t" << peak2valley << "\n";

	outfile.close();
}

float IntegrateCharge(std::vector<float> *wx, std::vector<float> *wy,
										  int preGate=5, int gate=50){
	int timeBinWidth = wx->at(1) - wx->at(0);

	// Find the minimum voltage
	// Don't look withing a buffer of 5 bins which would have icomplete signal
	int minPos = std::distance(wy->begin(),
														 std::min_element(wy->begin()+5,wy->end()-5));
	float minCharge = wy->at(minPos);

	// Get the integration range
	int lowInt = minPos>preGate ? minPos - preGate : 0;
	int highInt = lowInt+gate;
	if(highInt >= wy->size()){
		highInt = wy->size()-1;
		lowInt = wy->size()-51;
	}

	// Make a vector for baseline calculation
	std::vector<float> vBaseline;
	vBaseline.reserve(wy->size()-gate);
	for(auto i=0; i<wy->size(); i++){
		if(i>=lowInt && i<highInt) continue;
		vBaseline.push_back(wy->at(i));
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
		charge += baseline - wy->at(iSample);
	}
	charge = charge*timeBinWidth/fImpedance; // Charge is voltage*time/impedance

	return charge;
}