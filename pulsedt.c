#include <TROOT.h>
#include <math.h>

#include <stdio.h>
#include <vector>
#include "Kalman.h"
#include "sgfilter.h"
float q = 1;
float r = .005;
float P = 1;
float K = 0;
Kalman myFilter = Kalman(q,r,P,K);//(q,r,p,k) = (process noise variance,measurement noise variance,estimated error covariance,kalman gain)
// Histogram parameters
double hist_ylo[4] = { -0.5, -0.5, -0.5, -0.5 }; // Volts
double hist_yhi[4] = { 0.5, 0.5, 0.5, 0.5 }; // Volts
double pulsearea_max = 700E9;
double pulsearea_range = 300E9;

// Settings
double cfd_frac = 0.5;// Fraction of height to measure to
int skip = 32; // Skip start and end of signal
double repeat_thresh = 0.495; // When checking for clipping, signal must be > this threshold
unsigned int repeat_max = 2; // Discard if (# repeated values) > repeat_max
//double signal_thresh[4] = { 0.025, 0.02, 0.025, 0.025 }; // Discard if whole signal < signal_thresh
double signal_thresh[4] = { 0.01, 0.01, 0.01	, 0.01 }; // Discard if whole signal < signal_thresh
int interp_type = 0; // 0 is linear, 1 is cubic spline
int smoothingbool = 1; // true to smooth with a kalman filter
//SINC Filtering Settings and array declaration

//int FIRFilterBool = 0; //Steven Edit - testing new methods
//double upsampledpoints[1024][100] = {{0}};
//double upsampledtimes[1024][100] = {{0}};
// Floating-point equality check
const double epsilon = 0.001;
#define floateq(A,B) (fabs(A-B) < epsilon)

unsigned int active_channels[] = {1,2};
size_t nactive_channels = 2; // Length of the above array

void print_graph(double t[1024], double wave[1024],float amp,int center){
		//printf("t: %f,y: %f",t,wave);
		//for (int m = 0;m<1024;m++){printf("t: %f,y: %f",t[m],wave[m]);}}
		double delay = t[center]-4;
		for (int m = 0;m<1024;m++){
			t[m] = t[m]-delay;
			printf("t: %f\n",t[m]);
		}
		TGraph *printed_graph = new TGraph(1024,t,wave);
		TF1 *func = new TF1("f1","[0] * exp(-(x - [1]) / [2]) * (1 - exp(-(x - [1]) / [3]))");
		func->SetParameter(0,amp);
		func->SetParName(0,"Amplitude");
		func->SetParameter(1,0);
		func->SetParName(1,"Delay");
		func->SetParameter(2,2.3);
		func->SetParName(2,"Fall Time");
		func->SetParameter(3,1.2);
		func->SetParName(3,"Rise Time");
		printed_graph->Draw("AC*");
		printed_graph->Fit("f1","r");
		func->Draw("same");
		return;
}
// double double_exp_fit(double t,double *par){
// 		return par[0] * exp(-(t - par[1]) / par[2]) * (1 - exp(-(t - par[1]) / par[3]));
// }
void pulsedt(const char *filename){
	TFile *f = new TFile(filename);
	TTree *tree;
	double raw_waveform[4][1024];
	double time[4][1024];
	f->GetObject("tree", tree);
	tree->SetBranchAddress("ch0waveforms", &raw_waveform[0][0]);
	tree->SetBranchAddress("ch1waveforms", &raw_waveform[1][0]);
	tree->SetBranchAddress("ch2waveforms", &raw_waveform[2][0]);
	tree->SetBranchAddress("ch3waveforms", &raw_waveform[3][0]);
	tree->SetBranchAddress("ch0time", &time[0][0]);
	tree->SetBranchAddress("ch1time", &time[1][0]);
	tree->SetBranchAddress("ch2time", &time[2][0]);
	tree->SetBranchAddress("ch3time", &time[3][0]);
	//TH2D *dtAmp = new TH2D("dtAmplitude", "dt vs. Amplitude", 1000, -1, 200, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch1waveforms = new TH2D("ch1waveforms", "Channel 1 Waveforms", 1024, -0.1, 200, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch2waveforms = new TH2D("ch2waveforms", "Channel 2 Waveforms", 1024, -0.1, 200, 512, hist_ylo[1], hist_yhi[1]);
	TH2D *ch3waveforms = new TH2D("ch3waveforms", "Channel 3 Waveforms", 1024, -0.1, 200, 512, hist_ylo[2], hist_yhi[2]);
	TH2D *ch4waveforms = new TH2D("ch4waveforms", "Channel 4 Waveforms", 1024, -0.1, 200, 512, hist_ylo[3], hist_yhi[3]);
	TH2D *hwaveforms[] = { ch1waveforms, ch2waveforms, ch3waveforms, ch4waveforms };

	TH1D *ch1rise_time = new TH1D("ch1rise_time", "Channel 1 Rise Time", 1000, 0, 10);
	TH1D *ch2rise_time = new TH1D("ch2rise_time", "Channel 2 Rise Time", 1000, 0, 10);
	TH1D *ch3rise_time = new TH1D("ch3rise_time", "Channel 3 Rise Time", 1000, 0, 10);
	TH1D *ch4rise_time = new TH1D("ch4rise_time", "Channel 4 Rise Time", 1000, 0, 10);
	TH1D *hrise_time[] = { ch1rise_time, ch2rise_time, ch3rise_time, ch4rise_time };

	// /TH1I *hDPeakIndex = new TH1I("hDPeakIndex", "hDPeakIndex", 1025, -512, 512);

	TH1D *hPulseHeight[4];
	for (int i=0; i<4; i++) {
		std::stringstream ssname, sslabel;
		ssname << "hPulseHeight_ch" << i+1;
		sslabel << "Pulse Height (Channel " << i+1 << ")";
			hPulseHeight[i] = new TH1D(ssname.str().c_str(), sslabel.str().c_str(), 512, hist_ylo[i], hist_yhi[i]);
	}

	TH1D *hPulseArea[4];
	for (int i=0; i<4; i++) {
		std::stringstream ssname, sslabel;
		ssname << "hPulseArea_ch" << i+1;
		sslabel << "Pulse Area (Channel " << i+1 << ") in unit pC";
		hPulseArea[i] = new TH1D(ssname.str().c_str(), sslabel.str().c_str(), 1000, pulsearea_max-pulsearea_range, pulsearea_max+pulsearea_range);
	}

	// Print settings
	printf("CFD Fraction:      %f\n", cfd_frac);
	printf("Skip Ticks:        %d\n", skip);
	printf("Max Repeat Values: %d\n", repeat_max);
	printf("Minimum Signal Threshold (ch1): %f V\n", signal_thresh[0]);
	printf("Minimum Signal Threshold (ch2): %f V\n", signal_thresh[1]);
	printf("Minimum Signal Threshold (ch3): %f V\n", signal_thresh[2]);
	printf("Minimum Signal Threshold (ch4): %f V\n", signal_thresh[3]);
	printf("Voltage Range (ch1): [%f, %f] V\n", hist_ylo[0], hist_yhi[0]);
	printf("Voltage Range (ch2): [%f, %f] V\n", hist_ylo[1], hist_yhi[1]);
	printf("Voltage Range (ch3): [%f, %f] V\n", hist_ylo[2], hist_yhi[2]);
	printf("Voltage Range (ch4): [%f, %f] V\n", hist_ylo[3], hist_yhi[3]);
	printf("Max. Pulse Area: %f Vns\n", pulsearea_max);

	// dt between each channel
	int combinations = 1; // nactive_channels choose 2
	if (nactive_channels == 2) {
		combinations = 1;
	} else if (nactive_channels == 3) {
		combinations = 3;
	} else if (nactive_channels == 4) {
		combinations = 6;
	} else if (nactive_channels > 4) {
		combinations = 6;
		puts("WARNING: Too many channels: dt calculations capped at 4");
	}

	std::vector<TH1D*> hdt;
	// Delta-time between channels A and B is named "dtAB"
	for (int i=0; i<nactive_channels; i++) {
		for (int j=i+1; j<nactive_channels; j++) {
			unsigned int c1 = active_channels[i];
			unsigned int c2 = active_channels[j];
			std::stringstream ss1;
			ss1 << "dt" << c1 << c2;
			std::stringstream ss2;
			ss2 << "#Deltat_{" << c1 << c2 << "}";
			//TH1D *h = new TH1D(&name[0], &name[0], 1024, 0, 200);
			//TH1D *h = new TH1D(&name[0], &name[0], 2048, -60, 60);
			TH1D *h = new TH1D(ss1.str().c_str(), ss2.str().c_str(), 10000, -30, 30);
			hdt.push_back(h);

			//printf("%d %d (%d %d) %s\n", i, j, c1, c2, name);
		}
	}

	unsigned int ndiscarded = 0;

	long long nentries = tree->GetEntries();
	printf("Processing %lld entries\n", nentries);

	// For interpolating to find the constant-fraction time
	int interp_pts_up = 9; // # points above the principle point
	int interp_pts_down = 9; // # points below the principle point
	int ninterp_pts = interp_pts_up + 1 + interp_pts_down;
	int interp_pts_peak = 6;
	TGraph *interp_graph = new TGraph(ninterp_pts);
	TGraph *interp_graphpeak = new TGraph(2*interp_pts_peak+1);

	// filtered
	double waveform[4][1024];

	for (int i=0; i<nentries; i++) {
		if ((i % (nentries/10)) == 0) {
			printf("%d%%\n", int(100.0*float(i)/float(nentries)));
		}

		tree->GetEntry(i);

		// filter
		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			if (smoothingbool) {
				myFilter.setParameters(q,r,P);
				//SGSmoothing::Smooth(1024, &raw_waveform[c][0], &waveform[c][0], 5, 3);
				for (int j = 0;j<1024;j++){
				 	waveform[c][j] = myFilter.getFilteredValue(raw_waveform[c][j]);
				}

			} else {
					for (int j=0; j<1024; j++) { waveform[c][j] = raw_waveform[c][j];}
			}
		}



		int signal_good[4] = { 0, 0, 0, 0 };
		// discard is logically boolean. Its value is the repeat channel + 1
		int discard = 0;

		double prevval[4] = { -10000, -10000, -10000, -10000 };
		unsigned int repeat[4] = { 0, 0, 0, 0 };

		// Loop over each data point
		for (int j=0+skip; j<1024-skip; j++) {
			// Check threshold
			for (int k=0; k<nactive_channels; k++) {
				unsigned int c = active_channels[k] - 1;
				if (fabs(waveform[c][j]) > signal_thresh[c]) {
					signal_good[c] = 1;
				}
				if (floateq(waveform[c][j], prevval[c])) {
					++repeat[c];
				} else {
					repeat[c] = 0;
					prevval[c] = waveform[c][j];
				}
			}

			// Check repeat
			for (int k=0; k<nactive_channels; k++) {
				unsigned int c = active_channels[k] - 1;
				if ((repeat[c] >= repeat_max) && (fabs(prevval[c]) > repeat_thresh)) {
					discard = c+1;
					break;
				}
			}

			if (discard) {
				break;
			}
		}

		int signal_discard = 0;
		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			if (signal_good[c] == 0) {
				signal_discard = c+1;
				break;
			}
		}
		if (signal_discard) {
			++ndiscarded;
			//printf("Discard entry %d: channel %d below threshold\n", i, signal_discard);
			continue;
		}

		if (discard) {
			++ndiscarded;
			int ch = discard-1;
			//printf("Discard entry %d: %f repeated %d times in channel %d\n", i, prevval[ch], repeat[ch], discard);
			continue;
		}

		// 	g->SetPoint(i, time, 1.0 + TMath::Sin(TMath::TwoPi() * i / 99.0));
		//g->Draw("AF");
		//std::cout << "g->Integral() = " << g->Integral() << std::endl;

		// Find peak

		double peak_val[4] = { -100, -100, -100, -100 };
		double pulse_area[4] = { 0, 0, 0, 0 };
		unsigned int peak_idx[4] = { 0, 0, 0, 0 };
		double rise_time[4] = { 0, 0, 0, 0 };


		double baseline[4];
		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;

			baseline[c] = 0;
			for (int j = skip;j<=4*skip;j++){
				baseline[c]+=waveform[c][j];
			}
			baseline[c] = baseline[c]/(3*skip);
		}
		for (int j=0+skip; j<1024-skip; j++) {
			for (int k=0; k<nactive_channels; k++) {
				unsigned int c = active_channels[k] - 1;
				hwaveforms[c]->Fill(time[c][j], waveform[c][j]);
				double w = fabs(waveform[c][j]-baseline[c]);
				if (w > peak_val[c]) {
					peak_val[c] = w;
					peak_idx[c] = j;
				}
			}
		}


		//for (int k=0; k<nactive_channels; k++) {
		//unsigned int c = active_channels[k] - 1;
		// 	upsampledata( waveform,c,peak_idx[c]-2,N, L, T);
		// 	upsampledata( waveform,c,peak_idx[c],N, L, T);
		// 	upsampledata( waveform,c,peak_idx[c]-1,N, L, T);
		// 	upsampledata( waveform,c,peak_idx[c]+1,N, L, T);
		// 	upsampledata( waveform,c,peak_idx[c]+2,N, L, T);
		// 	for (int j=peak_idx[c]-2;j<=peak_idx[c]+2;j++){
		// 		for(int m=0;m<int(N);m++){
		// 			if(fabs(upsampledpoints[j][m]-baseline[c])>fabs(peak_val[c])){
		// 				peak_val[c] = (fabs(upsampledpoints[j][m]-baseline[c]));
		// 			}
		// 		}
		// 	}

		//find Pulse Height
		int area_range = 50;
		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			for (int m=-interp_pts_peak;m<=interp_pts_peak;m++){
				interp_graphpeak->SetPoint(m+interp_pts_peak,time[c][peak_idx[c]+m],fabs(waveform[c][peak_idx[c]+m]));
			}
			TFitResultPtr r =interp_graphpeak->Fit("gaus","SQIF","ROB=.9");
			if (peak_val[c]<0){peak_val[c] = -1*r->Value(0);}
			else{peak_val[c] = r->Value(0);}
			hPulseHeight[c]->Fill(peak_val[c]);
		}
		//find PulseArea
		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			TGraph *g = new TGraph(2*area_range+1);
			for (int m = peak_idx[c]-area_range; m <= peak_idx[c]+area_range; m++){
				g->SetPoint(m-(peak_idx[c]-area_range), time[c][m],waveform[c][m]);
			}
			pulse_area[c] = g->Integral(0,-1)/1E-12; //gives in unit pC -> scale by known capacitance of source
			//printf("Pulse Area: %f (V-ns/pC)\n",pulse_area[c]);
			hPulseArea[c]->Fill(fabs(pulse_area[c]));
		}


		// Interpolate to find base
		double frac_time[4] = { 1000, 1000, 1000, 1000 }; // Beyond 200ns window

		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			double vf = baseline[c]+(peak_val[c]-baseline[c]) * cfd_frac;
			double vf9 = baseline[c]+(peak_val[c]-baseline[c]) * .8;
			double vf1 = baseline[c]+(peak_val[c]-baseline[c]) * .2;
			for (int j=peak_idx[c]; j>skip 	; j--) {
				//printf("Peak value: %f, wvfm value: %f, j: %d,c: %d\n",fabs(vf),fabs(waveform[c][j]),j,c);
				if (fabs(waveform[c][j]) < vf) {
					double t = -1000;
						for (int p=-interp_pts_down; p<=interp_pts_up; p++) {
							int idx = j + p;
							double t = time[c][idx];
							double v = fabs(waveform[c][idx]);
							// swap x and y because we can't eval on y
							interp_graph->SetPoint(interp_pts_down+p, t, v); //To solve for "Zeros" later
							}

						if (interp_type == 1) {
							t=interp_graph->Eval(vf, 0, "S");

						} else {
							TF1 *r = new TF1("f1","pol3");
							interp_graph->Fit("f1","SQIF","ROB=.9");
							t = r->GetX(vf,time[c][j-interp_pts_down],time[c][j+interp_pts_up]);//(vf-b)/m;
							rise_time[c] = r->GetX(vf9,time[c][j-interp_pts_down],time[c][j+interp_pts_up])-r->GetX(vf1,time[c][j-interp_pts_down],time[c][j+interp_pts_up]);//(vf9-b)/m - (vf1-b)/m;
							hrise_time[c]->Fill(rise_time[c]);
							//printf("y = %f*x+%f; t = %f\n",m,b,(vf-b)/m);
						}
						frac_time[c] = t;
						if (t != -1000){
							break;
						}
					if (j == skip) {
						printf("WARNING: %d: Failed to find fraction (ch%d)\n", i, c+1);
					}
				}
			}
		}
					// delta t

		std::vector<double> dt;
		for (int j=0; j<nactive_channels; j++) {
			for (int k=j+1; k<nactive_channels; k++) {
				unsigned int c1 = active_channels[j] - 1;
				unsigned int c2 = active_channels[k] - 1;
				double t = (frac_time[c1] - frac_time[c2]);
				//TGraph *printed_graph = new TGraph(1024,time[c1],waveform[c1]);
				//printed_graph->Draw("AC*");
				//print_graph(time[c1],waveform[c1],peak_val[c1],peak_idx[c1]);
				//return;
				//double da = (baseline[c1] - baseline[c2]);
				//dtda->Fill(da,t);
				//printf("%d chs. %d and %d t1=%f t2=%f dt=%f\n", i, c1+1, c2+1, frac_time[c1], frac_time[c2], t);
				if (floateq(t, 0)) {
					//printf("%d dt%d%d = 0\n", i, c1+1, c2+1);
				}
				dt.push_back(t);
			}
		}
		for (int j=0; j<dt.size(); j++) {
			hdt[j]->Fill(dt[j]);
		}
	}
	puts("100%\n");

	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(0011);

	TCanvas *c1 = new TCanvas("c1", "Channel 1 Waveforms");
	c1->SetLogz(1);
	ch1waveforms->GetXaxis()->SetTitle("time [ns]");
	ch1waveforms->GetYaxis()->SetTitle("voltage [V]");
	ch1waveforms->Draw("COLZ");

	TCanvas *c2 = new TCanvas("c2", "Channel 2 Waveforms");
	c2->SetLogz(1);
	ch2waveforms->GetXaxis()->SetTitle("time [ns]");
	ch2waveforms->GetYaxis()->SetTitle("voltage [V]");
	ch2waveforms->Draw("COLZ");

	TCanvas *c3 = new TCanvas("c3", "Channel 3 Waveforms");
	c3->SetLogz(1);
	ch3waveforms->GetXaxis()->SetTitle("time [ns]");
	ch3waveforms->GetYaxis()->SetTitle("voltage [V]");
	ch3waveforms->Draw("COLZ");

	TCanvas *c4 = new TCanvas("c4", "Channel 4 Waveforms");
	c4->SetLogz(1);
	ch4waveforms->GetXaxis()->SetTitle("time [ns]");
	ch4waveforms->GetYaxis()->SetTitle("voltage [V]");
	ch4waveforms->Draw("COLZ");

	TCanvas *cdt = new TCanvas("cdt", "dt Between Channels");

	double hdtmax = 0;
	for (int i=0; i<hdt.size(); i++) {
		double m = hdt[i]->GetMaximum();
		if (m > hdtmax) {
			hdtmax = m;
		}
	}
	hdt[0]->SetLineColor(1);
	hdt[0]->SetMaximum(hdtmax+hdtmax*0.1f);
	hdt[0]->GetXaxis()->SetTitle("time [ns]");
	hdt[0]->GetYaxis()->SetTitle("frequency");
	hdt[0]->Draw();
	for (int i=1; i<hdt.size(); i++) {
		//hdt[i]->Fit("gaus");
		hdt[i]->SetLineColor(i+1);
		//hdt[i]->Fit("gaus","V","ROB=.9");
		//printf("**********************Delta T Distribution %d**********************",i);
		hdt[i]->Draw("same");
	}
	*cdt->BuildLegend();

	TCanvas *cph = new TCanvas("cph", "Pulse Height");
	double hphmax = 0;
	for (int i=0; i<4; i++) {
		double m = hPulseHeight[i]->GetMaximum();
		if (m > hphmax) {
			hphmax = m;
		}
	}
	hPulseHeight[0]->SetLineColor(1);
	hPulseHeight[0]->SetMaximum(hphmax+hphmax*0.1f);
	hPulseHeight[0]->GetXaxis()->SetTitle("Amplitude [V]");
	hPulseHeight[0]->GetYaxis()->SetTitle("frequency");
	hPulseHeight[0]->Draw();
	for (int i=1; i<4; i++) {
		hPulseHeight[i]->SetLineColor(i+1);
		hPulseHeight[i]->Draw("same");
	}
	cph->BuildLegend();
	TCanvas *cprt = new TCanvas("cprt", "Pulse Rise Time");
	hrise_time[0]->SetLineColor(1);
	hrise_time[0]->SetMaximum(hphmax+hphmax*0.1f);
	hrise_time[0]->GetXaxis()->SetTitle("Rise Time [ns]");
	hrise_time[0]->GetYaxis()->SetTitle("frequency");
	hrise_time[0]->Draw();
	for (int i=1; i<4; i++) {
		hrise_time[i]->SetLineColor(i+1);
		hrise_time[i]->Draw("same");
		//hPulseArea[i]->Fit("gaus","V");
	}
	cprt->BuildLegend();
	TCanvas *cpa = new TCanvas("cpa", "Pulse Rise Time");
	hPulseArea[0]->SetLineColor(1);
	hPulseArea[0]->SetMaximum(hphmax+hphmax*0.1f);
	hPulseArea[0]->GetXaxis()->SetTitle("Cummulative Charge [v*ns/pC]");
	hPulseArea[0]->GetYaxis()->SetTitle("frequency");
	hPulseArea[0]->Draw();
	for (int i=1; i<4; i++) {
		hPulseArea[i]->SetLineColor(i+1);
		hPulseArea[i]->Draw("same");
		//hPulseArea[i]->Fit("gaus","V");
	}
	cpa->BuildLegend();
	printf("Kept %lld/%lld events\n", (nentries-ndiscarded), nentries);

	// XXX THIS IS BAD
	std::stringstream ss;
	ss << "histograms_" << filename;
	TFile *outfile = new TFile(ss.str().c_str(), "RECREATE");

	ch1waveforms->Write();
	ch2waveforms->Write();
	ch3waveforms->Write();
	ch4waveforms->Write();
	for (int i=0; i<hdt.size(); i++) {
		hdt[i]->Write();
	}
	for (int i=0; i<4; i++) {
		hPulseArea[i]->Write();
		hPulseHeight[i]->Write();
	}
	outfile->Close();

	std::cout << "Saved to " << ss.str() << std::endl;
}
