/// This code is old and haven't been updated yet - DO NOT USE!

#include <fstream>


void GetMeanStd(std::vector<double> v, double &m, double &d){

	double mean = 0;
	for (int i = 0; i < v.size(); i++){
		mean += v[i];
	}
	mean /= (double)v.size();

	double dev = 0;
	for (int i = 0; i < v.size(); i++){
		dev+= (v[i] - mean)*(v[i] - mean);
	}

	m = mean;
	d = sqrt(dev/v.size());

}

 void DarkRate(std::string inFile, double pe = 10.0, double thresh = 1.0){
	std::ofstream outfile;

  outfile.open("darkrates.txt", std::ios_base::app); // append instead of overwrite

	TFile *input = new TFile (inFile.c_str(), "READ");

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
	
	tDevice->GetEntry(0); // frequency and maxSamples only set once

	float adcFactor = (voltHigh - voltLow)*1000/resolution; // Multiply by adc counts to get voltage [mV]



	int width = 5;
	int dark = 0;
	double timef = 0;
	
	std::cout << "Entries:\t" << tWaves->GetEntries() << std::endl;
	
	std::vector<double> darkVals;
	// Now we loop over the events
	for (int c = 0; c < tWaves->GetEntries(); c++){

		tWaves->GetEntry(c);
		int size = wavey->size();
		timef += (double)(4e-9)*size;
		
		if (c % 100000 == 0 || c == tWaves->GetEntries() -1){
			darkVals.push_back((double)dark/timef);
	
			dark = 0;
			timef = 0;
		}

		bool belowThreshold = true;

		
		

		if (c == 0) {std::cout << "Size:\t" << size << std::endl;}

		double baseline = 0;

		for (int i = 0; i < size; i++){
			baseline += wavey->at(i);
		}
		baseline /= (double)(size);


		for (int i = 0; i < size - width; i++){

			double sum = 0;

			for (int j = 0; j < width; j++){
				sum += wavey->at(i+j) - baseline;
			}

			sum*=adcFactor;

			if (belowThreshold){
				if ( sum >= pe*thresh ) { dark++; belowThreshold = false; }
			} else {
				if (sum < pe*thresh){belowThreshold = true;}
			}
			
		}
			    
	}


	double rateMean = 0;
	double rateSTD = 0;
	GetMeanStd(  darkVals, rateMean, rateSTD);
	std::cout << "DR = " << rateMean << "+-" << rateSTD << std::endl;
	outfile << inFile << "\t" << rateMean << "\t" << rateSTD << "\n";  

 }

