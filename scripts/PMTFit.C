///
/// @file PMTFit
/// @brief Fit the PMT response function to the integrated charge distribution
/// @details This script reads in a ROOT file with a TTree of
/// waveforms and fits the integrated charge distribution to a PMT response 
/// function. The input is a CoMPASS-style ROOT file, which can also be
/// generated for CAEN Scope xml outputs using the waveconvert script.
/// The input ROOT file should have a TTree called "Data" with the waveform
/// data saves as a TArrayS called "Samples".
/// The PMT response function derived from a physics-motivated description and 
/// consists a sum of Gaussians for 0, 1, 2, 3, 4 PE and an exponential 
/// background.
/// The fit is done using the ROOT fitting functions.
/// The script outputs the fit parameters to a csv file.
/// The script also outputs a plot of the charge distribution with the fit overlaid.
/// The script can be run from the command line with the following arguments:
/// root 'PMTFit.C("<input ROOT file>","<output csv file>")'
/// Alternatively, the script can be run by Running
/// bash RunPMTFit.sh <input ROOT file> <output csv file>
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

// If using CoMPASS outputs need to set these values manually
#define fResolution 500 ///< Number of ADC bins (=2^ADCResolution)
#define fVoltLow 0.f ///< Voltage range low
#define fVoltHigh 2.f ///< Voltage range high
#define fFrequency 500e6 ///< Frequency of the digitiser in Hz
#define fWindowSize 600 ///< Size of the waveform (number of time bins)

Double_t factorial(int a){
	if(a > 1){return a*factorial(a-1);} else return 1;
}

Double_t Pois(double mu, double n){
	return pow(mu,n)*exp(-mu)/factorial(n);
}

/// @brief Gauss for n PE centered at (x-q0)
Double_t Gaus(double x, double q0, double q, double s, double n){
	return 1/(s*sqrt(2*PI*n))*exp(-pow(x-q0-n*q,2)/(2*n*pow(s,2)));
}

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

/// @brief Full PMT background function (pedestal and exponential background)
/// @param x Input
/// @param q0 Mean of the pedestal
/// @param s0 Sigma of the pedestal
/// @param q1 Mean of the 1 PE peak
/// @param s1 Sigma of the 1 PE peak
/// @param a alpha = coefficient of the exponential decrease
/// @param n number of PE
Double_t Ignxe(double x, double q0, double s0, double q1, double s1, double a, double n);

/// @brief Full PMT response function as described in https://doi.org/10.1016/0168-9002(94)90183-X
Double_t PMTF(Double_t *x_, Double_t *par);

/// @brief 0 PE (pedestal) component of the full PMT response function
Double_t PMTF0(Double_t *x_, Double_t *par);

/// @brief 1 PE component of the full PMT response function
Double_t PMTF1(Double_t *x_, Double_t *par);

/// @brief 2 PE component of the full PMT response function
Double_t PMTF2(Double_t *x_, Double_t *par);

/// @brief Take parameters from PMT response function a and fix them onto p
/// @param p Function to fix the paramters onto
/// @param a Function to take the paramters from
void FixFit(TF1 *&p, TF1 *&a);

void SetParametersFromFitResult(TF1*& fOut, TF1*& fIn);

/// @brief Integrate waveform to get a deposited charge using old input - kept for backwards compatibility
/// @param wx Input vector of waveform x values (time in ns)
/// @param wy Input vector of waveform y values (voltage in mV)
/// @param preGate Length of time to integrate before the peak in bins of time. Default 5
/// @param gate Length of time to integrate from preGate in bins of time. Default 50
/// @return Inegrated charge in pC
float IntegrateCharge(std::vector<float> *&wx, std::vector<float> *&wy,
										  int preGate=5, int gate=50);

/// @brief Integrate waveform to get a deposited charge
/// @param Samples Input TArrayS of waveform values (output from waveconvert)
/// @param ADCTomV Conversion factor from ADC bins to mV
/// @param timeBinWidth Width of the time bin in nanoseconds
/// @param preGate Length of time to integrate before the peak in bins of time. Default 5
/// @param gate Length of time to integrate from preGate in bins of time. Default 50
/// @return Inegrated charge in pC
float IntegrateCharge(TArrayS *&Samples, float ADCTomV, float timeBinWidth,
										  int preGate=5, int gate=50);

std::pair<float,float> FindChargeMinmax(std::vector<float> &vecCharge,
																			  int minBinContent=3);

/// @brief Get PMT label from the inFileName
/// @param inFileName Name of the input file with path and extension
/// @return PMT label
std::string GetPMTLabel(std::string inFileName);



/// @brief Main fit function that loads waves, defines fitting functions, does the fits, plots the result and saves it to a csv file
/// @param inFileName (string) Name of the input ROOT file
/// @param outFileName (string) Name of the output csv file
void PMTFit(std::string inFileName, std::string outFileName="Output.csv"){
	TFile* inFile = new TFile(inFileName.c_str(),"READ");
	// Get the parameters of the digitiser
	float ADCTomV;
	float timeBinWidth;
	GetParams(inFile, ADCTomV, timeBinWidth,
					  fResolution, fVoltLow, fVoltHigh, fFrequency);
	
	TTree *tWaves = (TTree*)inFile->Get("Data");
	//std::vector<float> *wavex = 0;
	//std::vector<float> *wavey = 0;
	//tWaves->SetBranchAddress("wavex", &wavex);
	//tWaves->SetBranchAddress("wavey", &wavey);

	//Short_t Channel = -1; ///< Channel number
	int Clocktime = -1; ///< Clocktime of the event (in Unix time)
	TArrayS *Samples = new TArrayS; ///< Array of samples (waveform values)
	//tWaves->SetBranchAddress("Channel", &Channel);
	tWaves->SetBranchAddress("Clocktime", &Clocktime);
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
	auto chargeMinmax = FindChargeMinmax(vecCharge, 5);
	float min = chargeMinmax.first; float max = chargeMinmax.second;

	TH1F *hCharge = new TH1F("charge",
													 ";Integrated Charge [pC];Area Normalized (arb. units)", fNChargeBins, min, max);
	hCharge->GetXaxis()->CenterTitle();
	hCharge->GetYaxis()->CenterTitle();
	for(auto iCharge : vecCharge) hCharge->Fill(iCharge);
	hCharge->Scale(1/hCharge->Integral());

	TCanvas c("c","c");
	c.cd();

	hCharge->Draw("axis");

	//TODO: #6 Add option to do multi-PE fit as well

	//Use TSpectrum to find the peak candidates
  TSpectrum *s = new TSpectrum(2);
  Int_t nfound = s->Search(hCharge,5,"same",0.005);
	std::cout << "Found " << nfound << " peaks" << std::endl;
	Double_t *xpeaks = s->GetPositionX();
	for(int i=0; i<nfound; ++i)
		std::cout << "Peak " << i << " is " << xpeaks[i] << std::endl;
	double q0 = (xpeaks[0] > 0) ? xpeaks[0] : 0; // if pedestal peak not found then set it to zero
	double q1 = nfound>1 ? xpeaks[1] : 0.5*max; // if second peak not found then set it to half the charge distribution width

	TF1 *pmt = new TF1("pmt",PMTF,min,max,8);
	pmt->SetNpx(1000);
	TF1 *pmt0 = new TF1("pmt0",PMTF0,min,max,8);
	TF1 *pmt1 = new TF1("pmt1",PMTF1,min,max,8);
	TF1 *pmt2 = new TF1("pmt2",PMTF2,min,max,8);

	// TODO: #5 Get then number of PE from fraction of pedestal
  pmt->SetParNames("Q_{0}","#sigma_{0}","Q_{1}","#sigma_{1}", "w", "a", "#mu",
									 "Scaling factor");
	pmt->SetParameter(0,q0); // Mean of the pedestal
	pmt->SetParameter(1,0.1*(q1-q0)); // Sigma of the pedestal
	pmt->SetParameter(2,q1); // Mean of the SPE peak
	pmt->SetParameter(3,0.3*q1); // Expected SPE resolution is 30%
	pmt->SetParameter(4,0.01); // Exponentional background contribution
	pmt->SetParameter(5,1); // Background decay constant
	pmt->SetParameter(6,0.2); // "True" number of PE (should be < 1 for SPE)
	pmt->SetParameter(7,1); // Scaling factor

	pmt->SetParLimits(0,min,0.8*q1);
	pmt->SetParLimits(1,0,0.5*(q1-q0));
	pmt->SetParLimits(2,q0,max);
	pmt->SetParLimits(3,0.05*q1,q1);
	pmt->SetParLimits(4,1e-5,1);
	pmt->SetParLimits(5,0,10);
	pmt->SetParLimits(6,0.01,1.);
	pmt->SetParLimits(7,0.,100);
	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);

	hCharge->Fit("pmt","EMR"); //EM
	std::cout << "Chi^2/NDF = " << pmt->GetChisquare()/pmt->GetNDF() << std::endl;

	// If the fit did badly then try to fix the exponential background
	if(pmt->GetChisquare()/pmt->GetNDF() > 5){
		std::cout << "Fit failed - trying to fix the exponential background" << std::endl;
		if(pmt->GetParameter(2) < pmt->GetParameter(0) || pmt->GetParameter(0) > 10*q0){
			std::cout << "Single PE mean < Pedestal mean - re-do the entire fit one at a time" << std::endl;
			pmt->SetParameter(0,q0); // Mean of the pedestal
			pmt->SetParameter(1,0.1*(q1-q0)); // Sigma of the pedestal
			pmt->SetParameter(2,q1); // Mean of the SPE peak
			pmt->SetParameter(3,0.3*q1); // Expected SPE resolution is 30%
			pmt->FixParameter(4,0); // NO Exponentional background contribution
			pmt->SetParameter(6,0.1); // "True" number of PE (should be < 1 for SPE)
			pmt->SetParameter(7,0.1); // Scaling factor - should not be needed but is...
			hCharge->Fit("pmt","EMR"); //EM
			std::cout << "Chi^2/NDF = " << pmt->GetChisquare()/pmt->GetNDF() << std::endl;
		}
		// Now get the estimate of the exponential background from the histogram
		double valley_pos = pmt->GetMinimumX(pmt->GetParameter(0),
																				 pmt->GetParameter(2));
		double valley_fit = pmt->Eval(valley_pos);
		double valley_hist = hCharge->GetBinContent(hCharge->FindBin(valley_pos));
		double a = valley_hist > valley_fit ? -log(valley_hist-valley_fit) : 1;
		if(a>500) a = 50;
		if(valley_pos > 0.9*pmt->GetParameter(2)) a = 1;
		std::cout << "Pre-set alpha is " << a << std::endl;
		pmt->SetParameter(4,0.01); // Turn on Exponentional background contribution
		pmt->SetParLimits(4, 1e-5, 1);
		pmt->SetParameter(5, a); // Background decay constant
		pmt->SetParLimits(5,1e-5, pmt->Eval(pmt->GetParameter(0))/pmt->GetParameter(7));
		hCharge->Fit("pmt","EMR");
	}

	double chisqr = pmt->GetChisquare()/pmt->GetNDF();
	std::cout << "Chi^2/NDF = " << chisqr << std::endl;

	double SPECharge = pmt->GetParameter(2);
	double gain = (SPECharge-pmt->GetParameter(0))*1e-12/e;
	double peRes = pmt->GetParameter(3)/pmt->GetParameter(2); // s1/q1

	FixFit(pmt0, pmt);
	FixFit(pmt1, pmt);
	FixFit(pmt2, pmt);

	pmt0->Draw("same");
	pmt1->Draw("same");
	if(2*SPECharge < max) // Only draw if second PE peak below maximum
		pmt2->Draw("same");

	pmt->SetLineColor(kRed+1);
	pmt0->SetLineColor(kGray+1);
	pmt0->SetLineStyle(2);
	pmt1->SetLineColor(kBlue+1);
	pmt1->SetLineStyle(2);
	pmt2->SetLineColor(kMagenta+1);
	pmt2->SetLineStyle(2);

	TLegend leg(0.62,0.3,0.88,0.55);
	leg.AddEntry(hCharge,"Data","lep");
	leg.AddEntry(pmt,"Full PMT response","l");
	leg.AddEntry(pmt0,"Pedestal","l");
	leg.AddEntry(pmt1,"Single PE","l");
	if(2*SPECharge<max) leg.AddEntry(pmt2,"Two PE","l");
	leg.Draw("same");

	// Print the voltage
	std::string VLabel = inFileName.substr(0,inFileName.find_last_of("."));
	VLabel = VLabel.substr(VLabel.find_last_of("_")+1);
	TText text(.4,.8,("Input: "+VLabel).c_str());
	text.SetNDC();
	text.Draw("same");

	c.Update();

	std::string PMTLabel = GetPMTLabel(inFileName);
	CornerLabel(PMTLabel);

	std::string plotFileName = inFileName.substr(0,inFileName.find_last_of("."))+"_Fit";
	//c.SaveAs((plotFileName+"_Gate"+std::to_string(fGate)+".pdf").c_str());

	c.SetLogy();
	c.SaveAs((plotFileName+"_Logy"+"_Gate"+std::to_string(fGate)+".pdf").c_str());

	double peak = pmt->Eval(pmt->GetParameter(2));
	double valley = pmt->GetMinimum(pmt->GetParameter(0), pmt->GetParameter(2));
	double valleyPos = pmt->GetMinimumX(pmt->GetParameter(0), pmt->GetParameter(2));
	double peak2valley = peak/valley;
	double peakHist = hCharge->GetBinContent(hCharge->FindBin(q1));
	double valleyHist = 2*hCharge->GetBinContent(hCharge->FindBin(q0));
	int valleyBin = 1;
	for(int iBin=hCharge->FindBin(q0); iBin<hCharge->FindBin(q1); ++iBin){
		if(hCharge->GetBinContent(iBin)<valleyHist){
			valleyHist = hCharge->GetBinContent(iBin);
			valleyBin = iBin;
		}
	}
	for(int iBin=valleyBin; iBin<hCharge->FindBin(q1)+10; ++iBin){
		if(hCharge->GetBinContent(iBin)>peakHist)
			peakHist = hCharge->GetBinContent(iBin);
	}

	std::cout << "Q1 = " << SPECharge << "pC" << std::endl;
	std::cout << "Which is equal to a gain of " << gain << std::endl;
	std::cout << "Pedastal " << hCharge->Integral(0,hCharge->FindFixBin(valleyPos)) << std::endl;
	std::cout << "Signal " << hCharge->Integral(hCharge->FindFixBin(valleyPos)+1,fNChargeBins) << std::endl;
	std::cout << "PE Res " << peRes << std::endl;
	std::cout << "Peak " << peak << " and valley " << valley << std::endl;
	std::cout << "P2V " << peak2valley << std::endl;
	std::cout << "PeakHist " << peakHist << " and valleyHist " << valleyHist
						<< " P2VHist " << peakHist/valleyHist << std::endl;

	// Lets output this to some file
	std::ofstream outfile;
  outfile.open(outFileName.c_str(),
							 std::ios_base::app); // append instead of overwrite
	outfile << fGate << "," << inFileName << ","
					<< pmt->GetChisquare()/pmt->GetNDF() << ","
					<< pmt->GetParameter(2) << ","
					<< pmt->GetParameter(3) << ","
					<< gain << ","
					<< peak2valley << ","
					<< peakHist/valleyHist << ","
					<< peRes << ","
					<< pmt->GetParameter(6) << "\n";
	outfile.close();
	
	delete hCharge;
	delete s;
	delete pmt;
	delete pmt0;
	delete pmt1;
	delete pmt2;
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

float IntegrateCharge(std::vector<float> *&wx, std::vector<float> *&wy,
										  int preGate, int gate){
	int timeBinWidth = wx->at(1) - wx->at(0);

	// Find the minimum voltage
	// Don't look within a 5 bin-wide buffer in beginning and end, which would 
	// have incomplete signal
	int minPos = std::distance(wy->begin(),
														 std::min_element(wy->begin()+5,wy->end()-5));
	float minCharge = wy->at(minPos);

	// Get the integration range
	int lowInt = minPos>preGate ? minPos - preGate : 0;
	int highInt = lowInt+gate;
	if(highInt >= wy->size()){
		highInt = wy->size()-1;
		lowInt = wy->size()-1-gate;
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

std::pair<float,float> FindChargeMinmax(std::vector<float> &vecCharge,
																				int minBinContent=3){
	auto chargeMinmax = minmax_element(vecCharge.begin(),vecCharge.end());
	float minimum = (float)*chargeMinmax.first;
	float maximum = (float)*chargeMinmax.second;
	TH1F hTemp("temp","temp",fNChargeBins,
						*chargeMinmax.first,*chargeMinmax.second);
	for(auto charge : vecCharge) hTemp.Fill(charge);

	for(int iBin=1; iBin<=fNChargeBins; ++iBin){
		if(hTemp.GetBinContent(iBin)<minBinContent)
			minimum = hTemp.GetBinLowEdge(iBin);
		else break;
	}
	for(int iBin=fNChargeBins; iBin>0; --iBin){
		if(hTemp.GetBinContent(iBin)<minBinContent)
			maximum = hTemp.GetBinLowEdge(iBin)+hTemp.GetBinWidth(iBin);
		else break;
	}
	return std::make_pair(minimum,maximum);
}

Double_t Ignxe( double x, double q0, double s0, double q1, double s1, double a, double n ){
	double qn = q0 + n*q1;
	double sn = sqrt(pow(s0,2) + n*pow(s1,2));
	
	return (a/2)*exp(-a*(x-qn-0.5*a*pow(sn,2))) * (ROOT::Math::erf( fabs(q0 - qn - pow(sn,2)*a)/(sn*sqrt(2)) ) + TMath::Sign(1,x - qn - pow(sn,2)*a )*ROOT::Math::erf( fabs(x - qn - pow(sn,2)*a)/(sn*sqrt(2))));
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
	double scale = par[7];
	double step = x > q0 ? 1 : 0;

	return scale*Pois(mu, 0)*((1-w)*Gaus(x, 0,q0,s0,1) + step*w*a*exp(-a*(x-q0))) +
				 scale*Pois(mu, 1)*((1-w)*Gaus(x,q0,q1,s1,1) + w*Ignxe(x,q0,s0,q1,s1,a,1)) + 
				 scale*Pois(mu, 2)*((1-w)*Gaus(x,q0,q1,s1,2) + w*Ignxe(x,q0,s0,q1,s1,a,2)) +
		     scale*Pois(mu, 3)*((1-w)*Gaus(x,q0,q1,s1,3) + w*Ignxe(x,q0,s0,q1,s1,a,3)) + 
		     scale*Pois(mu, 4)*((1-w)*Gaus(x,q0,q1,s1,4) + w*Ignxe(x,q0,s0,q1,s1,a,4));
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
	double scale = par[7];
	double step = x > q0 ? 1 : 0;

	return scale*Pois(mu,0)*((1-w)*Gaus(x,0,q0,s0,1) + step*w*a*exp(-a*(x-q0)));
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
	double scale = par[7];

	return scale*Pois(mu,1)*((1-w)*Gaus(x,q0,q1,s1,1) + w*Ignxe(x,q0,s0,q1,s1,a,1));
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
	double scale = par[7];

	return scale*Pois(mu,2)*((1-w)*Gaus(x,q0,q1,s1,2) + w*Ignxe(x,q0,s0,q1,s1,a,2));
}

void FixFit(TF1 *&p, TF1 *&a){
	p->FixParameter(0, a->GetParameter(0));
	p->FixParameter(1, a->GetParameter(1));
	p->FixParameter(2, a->GetParameter(2));
	p->FixParameter(3, a->GetParameter(3));
	p->FixParameter(4, a->GetParameter(4));
	p->FixParameter(5, a->GetParameter(5));
	p->FixParameter(6, a->GetParameter(6));
	p->FixParameter(7, a->GetParameter(7));
}

void SetParametersFromFitResult(TF1*& fOut, TF1*& fIn){
	fOut->SetNpx(fIn->GetNpx());
	
	for(int i=0; i<fIn->GetNpar(); i++){
		fOut->SetParameter(i,fIn->GetParameter(i));
		Double_t parmin = 0.8*fIn->GetParameter(i);
		Double_t parmax = 1.2*fIn->GetParameter(i);
		//fIn->GetParLimits(i,parmin,parmax);
		std::cout << "Parameter " << i
							<< " is " << fIn->GetParameter(i)
							<< " between " << parmin
							<< " and " << parmax
							<< std::endl;
		fOut->SetParLimits(i,parmin,parmax);
	}
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