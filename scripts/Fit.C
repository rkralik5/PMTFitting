#include <vector>
#include <iostream>
#include <fstream>

#define PI 3.141592654


Double_t factorial(int a){

	if (a > 1) {return a*factorial(a-1);} else return 1;

}

Double_t Pois(double mu, double n){
	
	return pow(mu,n)*exp(-mu)/factorial(n);
	
}

Double_t Gaus(double x, double q, double s, double n){
		
	return 1/(s*sqrt(2*PI*n))*exp(-pow(x-n*q,2)/(2*n*pow(s,2)));
	
}

Double_t Ignxe( double x, double q0, double s0, double q1, double s1, double a, double n ){
	
	double qn = q0 + n*q1;
	double sn = sqrt( pow(s0,2) + n*pow(s1,2) );
	
	
	return (a/2)*exp( -a*(x-qn-0.5*a*pow(sn,2) ) ) * ( ROOT::Math::erf( fabs(q0 - qn - pow(sn,2)*a)/(sn*sqrt(2)) ) + TMath::Sign(1,x - qn - pow(sn,2)*a )*ROOT::Math::erf( fabs(x - qn - pow(sn,2)*a)/(sn*sqrt(2)) )  ) ;
	
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

	return Pois(mu, 0)*( (1-w)*Gaus(x,q0,s0,1) + w*Ignxe(x, q0, s0, q1, s1, a, 0) ) 
		+ Pois(mu, 1)*( (1-w)*Gaus(x,q1,s1,1) + w*Ignxe(x, q0, s0, q1, s1, a, 1) ) 
		+ Pois(mu, 2)*( (1-w)*Gaus(x,q1,s1,2) + w*Ignxe(x, q0, s0, q1, s1, a, 2) )  
		+ Pois(mu, 3)*( (1-w)*Gaus(x,q1,s1,3) + w*Ignxe(x, q0, s0, q1, s1, a, 3) ) 
		+ Pois(mu, 4)*( (1-w)*Gaus(x,q1,s1,4) + w*Ignxe(x, q0, s0, q1, s1, a, 4) );
	
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

	return Pois(mu, 0)*( (1-w)*Gaus(x,q0,s0,1) + w*Ignxe(x, q0, s0, q1, s1, a, 0) );
	
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

	return Pois(mu, 1)*( (1-w)*Gaus(x,q1,s1,1) + w*Ignxe(x, q0, s0, q1, s1, a, 1) );
	
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

	return Pois(mu, 2)*( (1-w)*Gaus(x,q1,s1,2) + w*Ignxe(x, q0, s0, q1, s1, a, 2) );
	
}


void PrintVector(std::vector<float> v, std::string s = "Printing Vector"){

	std::cout << "===================" << s << "============" ;
	for (int i = 0; i < v.size(); i++){
		std::cout << v[i] << std::endl;
	}
}
void WriteVector(std::vector<float> v, std::string s ){

	ofstream myfile;
	myfile.open (s);
	
		for (int i = 0; i < v.size(); i++){
		myfile << v[i] << "\n";
	}
	myfile.close();
}

// This function makes sure that the i-th value away from peak position m is still in the vector wy. 
// This also used to be used to ignore negative values, this would ensure that we stop integrating when the charge fell below the basecharge x.
// However this made fitting less reliable so I have commented out that line. Feel free to revisit as the digitizer has pushed us to do this!  
float StillSignal(float x, std::vector<float> *wy, int m, int i){
	float sum = 0;
	if (m+i >= 0 && m+i < wy->size()) {
		sum += (x - wy->at(m+i));
	}
	//if (sum > 0 ){ return sum;} else return 0; 
	return sum;
}
void ProcessWaveform(std::vector<float> *wx, std::vector<float> *wy, float &totCharge, float &timePeak, float hz, float adcFactor){

	// Read in points
	// integration window
	// find peak
	float widthNS = 40;
	int width = (int)(widthNS*(hz/1e9));
	int buffer = 5; // used to ignore the first and last entries of the vector as the signal would be incomplete if it is in this region.
	int halfSigWidth = 4;

	int minPos = 0;
	float minCharge = 10000;
	float maxCharge = -10000;
	float minTime = -10000;
	float fullCharge = 0;

	float sumCharge = 0;
	int sumInt = 0;

	float currentWindow = 0;
	float minimumWindow = 0;
	int minimumStart = 0;


	// Lets calculate the charge in a sliding window
	for (int i = 0 + buffer; i < wy->size() - buffer; i++){ // find minimum
		if (wy->at(i) < minCharge){minCharge = wy->at(i); minPos = i; minTime = wx->at(i);}
	} 

	// make sure can integrate charge properly
	int lowInt = minPos - halfSigWidth; 
	int highInt = minPos + halfSigWidth;

	if (lowInt < 0) lowInt = 0;
	if (highInt > wy->size()) highInt = wy->size();

	for (int i = 0; i < wy->size(); i++){
	
		float totCharge = 0;
		if (i < lowInt || i > highInt){
			sumCharge += wy->at(i);
			sumInt++;
		}	
	}

	bool keepPoint = true;

	float baseCharge = sumCharge/sumInt;
	float countCharge = 0;
	bool stop = false;
	
	countCharge += StillSignal(baseCharge, wy, minPos, 0);
	countCharge += StillSignal(baseCharge, wy, minPos, -1);
	if ( StillSignal(baseCharge, wy, minPos, -1) != 0 ) {  countCharge += StillSignal(baseCharge, wy, minPos, -2); }

	countCharge += StillSignal(baseCharge, wy, minPos, +1);
	if ( StillSignal(baseCharge, wy, minPos, +1) != 0 ) { countCharge += StillSignal(baseCharge, wy, minPos, +2);  }


	totCharge = countCharge*adcFactor; // voltage in [mV]
	timePeak = 1e9*minTime/hz; // time in [ns]

}
void FixFit(TF1 *p, TF1 *a ){

	p->FixParameter(0, a->GetParameter(0) );                                            
	p->FixParameter(1, a->GetParameter(1) );                                              
	p->FixParameter(2, a->GetParameter(2) );                                              
	p->FixParameter(3, a->GetParameter(3) );                                              
	p->FixParameter(4, a->GetParameter(4) );                                              
	p->FixParameter(5, a->GetParameter(5) ); 
	p->FixParameter(6, a->GetParameter(6) );  

}
void Fit(std::string inName){
	
   
	TF1 *pmt = new TF1("pmt",PMTF,0,200,7); // x in [0;10], 2 parameters	
	TH1F *testhi = new TH1F("testhi","Poisson distribution",300,0,300); 
	TF1 *pmt0 = new TF1("pmt0",PMTF0,0,200,7); // x in [0;10], 2 parameters
	TF1 *pmt1 = new TF1("pmt1",PMTF1,0,200,7); // x in [0;10], 2 parameters
	TF1 *pmt2 = new TF1("pmt2",PMTF2,0,200,7); // x in [0;10], 2 parameters
    //give the parameters meaningful names
    pmt->SetParNames("Q_{0}","#sigma_{0}","Q_{1}","#sigma_{1}", "W", "a", "#mu");       
	
	TFile *input = new TFile(inName.c_str(),"READ");

	TTree *tDevice = (TTree*) input->Get("device");
	TTree *tWaves = (TTree*) input->Get("data");

	float frequency;
	int maxSamples;
	int resolution;
	float voltLow;
	float voltHigh;
	std::vector<float> *wavex = 0;
	std::vector<float> *wavey = 0;
	int time;

	tDevice->SetBranchAddress("frequency",&frequency);
	tDevice->SetBranchAddress("maxSamples",&maxSamples);
	tDevice->SetBranchAddress("resolution",&resolution);
	tDevice->SetBranchAddress("voltLow",&voltLow);
	tDevice->SetBranchAddress("voltHigh",&voltHigh);

	tWaves->SetBranchAddress("wavex", &wavex);
	tWaves->SetBranchAddress("wavey", &wavey);
	tWaves->SetBranchAddress("time", &time);
	
	tDevice->GetEntry(0); // frequency and maxSamples etc only set once

	float adcFactor = (voltHigh - voltLow)*1000/resolution; // Multiply by adc counts to get voltage [mV]
	std::cout << "VoltLow: " << voltLow << std::endl;
	std::cout << "VoltHigh: " << voltHigh << std::endl;
	std::cout << "resolution: " << resolution << std::endl;


	TH1D *h1 = new TH1D("charge", "charge", 500, -50, 200); // define histogram 
	h1->GetXaxis()->SetTitle("Integrated Voltage [mV]");

	for (int i = 0; i < tWaves->GetEntries(); i++){ // fill histogram
		tWaves->GetEntry(i);
		float totCharge = 0;
		float timePeak = 0;

		ProcessWaveform(wavex, wavey, totCharge, timePeak, frequency, adcFactor);
		if (totCharge != -999){  h1->Fill( totCharge ); }
		//h2->Fill( timePeak );

	}



	//double entries = h1->GetEntries();
	double entries = h1->Integral();
	h1->Scale(1/entries);
	std::cout << "Entries: " << entries << std::endl;


	//Use TSpectrum to find the peak candidates
   	TSpectrum *s = new TSpectrum(2);
   	Int_t nfound = s->Search(h1,2,"",0.00001);
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
	double q1 = (peaksX[1] > 0) ? peaksX[1] : 10; // if second peak not found then set it to ten

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



	pmt->SetParameter(0,q0);                                            
	pmt->SetParameter(1,1.5);                                              
	pmt->SetParameter(2,q1);                                              
	pmt->SetParameter(3,10);                                              
	pmt->SetParameter(4,0.4999);                                              
	pmt->SetParameter(5,0.00001379); 
	pmt->FixParameter(6,0.0585);  

	pmt->SetParLimits(0,-3,3);
	pmt->SetParLimits(1,0,30);
	pmt->SetParLimits(2,0,40);
	pmt->SetParLimits(3,0,10);
	pmt->SetParLimits(4,0.,1);
	pmt->SetParLimits(5,0.,1);
	pmt->SetParLimits(6,0.001,1.5000);
	
	gStyle->SetOptStat(0);                                                      
	gStyle->SetOptFit(1);


	//h1->Fit("pmt","","",-10,100);    
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


	double voltage = (pmt->GetParameter(2))*1e-3; 
	double timeb = 4e-9; // time per bin
	double charge = timeb*voltage/50;
	double e = 1.6e-19;
	double gain = charge/e;

	double peRes = pmt->GetParameter(3)/pmt->GetParameter(2); // s1/q1
	double peak2valley = peaksY[1]/valley[1]; // height at q1/ height at valley

	std::cout << "Q1 = " << voltage << std::endl;
	std::cout << "That is equal to a charge of " << charge << "C" << std::endl;
	std::cout << "Which is equal to a gain of " << gain << std::endl;

	std::cout << "Pedastal" << h1->Integral(h1->FindFixBin(-50),h1->FindFixBin(10)) << std::endl;
	std::cout << "Signal" << h1->Integral(h1->FindFixBin(11),h1->FindFixBin(100)) << std::endl;
	std::cout << "PE Res " << peRes << std::endl;
	std::cout << "P2V " << peak2valley << std::endl;
	std::cout << pmt->GetParameter(0) << "\t" << pmt->GetParameter(1) << "\t" << pmt->GetParameter(2) << "\t" << pmt->GetParameter(6) << "\t" << pmt->GetParameter(3) << std::endl;

	// Lets output this to some file
	std::ofstream outfile;
	std::ofstream outfileParts;

    outfile.open("pmt.txt", std::ios_base::app); // append instead of overwrite

	outfile << inName << "\t" << pmt->GetParameter(0) << "\t" << pmt->GetParameter(1) << "\t" << pmt->GetParameter(2) << "\t" << pmt->GetParameter(3) << "\t" << gain << "\t" << peRes << "\t" << peak2valley << "\n";

	outfile.close();

}


void Gain(double v){

	double voltage = v*1e-3;
	double time = 4e-9; // time per bin
	double charge = time*voltage/50;
	double e = 1.6e-19;
	double gain = charge/e;
	std::cout << "Q1 = " << voltage << std::endl;
	std::cout << "That is equal to a charge of " << charge << "C" << std::endl;
	std::cout << "Which is equal to a gain of " << gain << std::endl;
}
