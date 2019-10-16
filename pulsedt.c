#include <TROOT.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include "Kalman.h"
#include "sgfilter.h"
float q = .002;
float r = .002;
float P = .002;
float K = 0;
Kalman myFilter = Kalman(q,r,P,K);//(q,r,p,k) = (process noise variance,measurement noise variance,estimated error covariance,kalman gain)
// Histogram parameters
double hist_ylo[4] = { -0.5, -0.5, -0.5, -0.5 }; // Volts
double hist_yhi[4] = { 0.5, 0.5, 0.5, 0.5 }; // Volts
double pulsearea_max = 10;

// Settings
double cfd_frac = 0.5;// Fraction of height to measure to
int skip = 32; // Skip start and end of signal
double repeat_thresh = 0.1; // When checking for clipping, signal must be > this threshold
unsigned int repeat_max = 4; // Discard if (# repeated values) > repeat_max
//double signal_thresh[4] = { 0.025, 0.02, 0.025, 0.025 }; // Discard if whole signal < signal_thresh
double signal_thresh[4] = { 0.01, 0.01, 0.01	, 0.01 }; // Discard if whole signal < signal_thresh
int interp_type = 0	; // 0 is linear, 1 is cubic spline
int smoothingbool = 0; // true to smooth with a kalman filter
//SINC Filtering Settings and array declaration
const float N = 100;
float T = 40;
int L = 6;
int FIRFilterBool = 1; //Steven Edit - testing new methods
double upsampledpoints[1024][100] = {{0}};
double upsampledtimes[1024][100] = {{0}};
// Floating-point equality check
const double epsilon = 0.001;
#define floateq(A,B) (fabs(A-B) < epsilon)

unsigned int active_channels[] = { 1, 2};
size_t nactive_channels = 2; // Length of the above array

int stop_doing_stuff = 0;

double sinc(double i){
	double result;
	if(i == 0){
		result = 0;
	}
	if(i!=0){
		result = sin(i)/i;
	}
	return result;
}
void upsampletime(double times[4][1024],int c, int idx,int N){
	upsampledtimes[idx][0] = times[c][idx];
	//upsampledtimes[idx][N] = times[c][idx+1];
	double dT = (times[c][idx]-times[c][idx-1])/double(N);
	for (int p = 1;p < int(N);p++){
		upsampledtimes[idx][p] = upsampledtimes[idx][p-1]+dT;
	}
}
void upsampledata(double waveform[4][1024],int c,int idx,int N, int L, int T){
	upsampledpoints[idx][0] =waveform[c][idx];
	//upsampledpoints[idx][N] =waveform[c][idx+1];
	for (int m = 1;m < int(N);m++){
		upsampledpoints[idx][m] = 0;
		for (int p = 0;p <= L-1;p++){
			double temp = (p*N+m);
			double sinc1 = sinc(temp*M_PI/N)*exp(-(temp/T)*(temp/T));
			temp = (p+1)*N-m;
			double sinc2 = sinc(temp*M_PI/N)*exp(-(temp/T)*(temp/T));
			upsampledpoints[idx][m] += waveform[c][idx-p]*sinc1 + waveform[c][idx+1+p]*sinc2;
			//printf("Test: %d,%d,%d,%f,%f,%f\n",m,p,idx,upsampledpoints[idx][m],waveform[c][idx-p]*sinc1,waveform[c][idx+1+p]*sinc2);
		}
	}
}

void pulsedt(const char *filename)
{
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

	TH2D *ht12difference = new TH2D("t12difference", "Channnel 1 - Channel 2", 1024, 0, 1024, 512, -0.1, 0.1);
	TH1D *hdtime12 = new TH1D("hdtime12", "Channel 1 - Channel 2 time difference", 512, -0.1, 0.1);

	TH2D *ch1waveforms_left = new TH2D("ch1waveforms_left", "Channel 1 Waveforms (left lobe)", 1024, 0, 1024, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch1waveforms_center = new TH2D("ch1waveforms_center", "Channel 1 Waveforms (center lobe)", 1024, 0, 1024, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch1waveforms_right = new TH2D("ch1waveforms_right", "Channel 1 Waveforms (right lobe)", 1024, 0, 1024, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch2waveforms_left = new TH2D("ch2waveforms_left", "Channel 2 Waveforms (left lobe)", 1024, 0, 1024, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch2waveforms_center = new TH2D("ch2waveforms_center", "Channel 2 Waveforms (center lobe)", 1024, 0, 1024, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch2waveforms_right = new TH2D("ch2waveforms_right", "Channel 2 Waveforms (right lobe)", 1024, 0, 1024, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *dtda = new TH2D("dtda", "Delta Amplitude vs dt", 512, -.005, .005, 512,-.6,.6);

	//TH1D *hch1peak = new TH1D("hch1peak", "hch1peak", 1024, -0.1, 200);
	//TH1D *hch2peak = new TH1D("hch2peak", "hch2peak", 1024, -0.1, 200);

	TH1I *hDPeakIndex = new TH1I("hDPeakIndex", "hDPeakIndex", 1025, -512, 512);

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
		sslabel << "Pulse Area (Channel " << i+1 << ")";
		hPulseArea[i] = new TH1D(ssname.str().c_str(), sslabel.str().c_str(), 50, 0, pulsearea_max);
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
			TH1D *h = new TH1D(ss1.str().c_str(), ss2.str().c_str(), 1024, -1, 1);
			hdt.push_back(h);

			//printf("%d %d (%d %d) %s\n", i, j, c1, c2, name);
		}
	}

	unsigned int ndiscarded = 0;

	long long nentries = tree->GetEntries();
	printf("Processing %lld entries\n", nentries);

	// For interpolating to find the constant-fraction time
	int interp_pts_up = 2; // # points above the principle point
	int interp_pts_down = 2; // # points below the principle point
	int ninterp_pts = interp_pts_up + 1 + interp_pts_down;
	TGraph *interp_graph = new TGraph(ninterp_pts);

	// filtered
	double waveform[4][1024];

	for (int i=0; i<nentries && !stop_doing_stuff; i++) {
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
				for (int j=0; j<1024; j++) { waveform[c][j] = raw_waveform[c][j]; }
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

		// Find peak

		double peak_val[4] = { -100, -100, -100, -100 };
		unsigned int peak_idx[4] = { 0, 0, 0, 0 };
		double baseline[4];
		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			baseline[c] = 0;
			for (int j = skip;j<=3*skip;j++){
				baseline[c]+=waveform[c][j];
			}
			baseline[c] = baseline[c]/(2*skip);
		}
		for (int j=0+skip; j<1024-skip; j++) {
			ht12difference->Fill(j, time[0][j] - time[1][j]);
			hdtime12->Fill(time[0][j] - time[1][j]);
			for (int k=0; k<nactive_channels; k++) {
				unsigned int c = active_channels[k] - 1;
				hwaveforms[c]->Fill(time[c][j], waveform[c][j]);
				double w = fabs(waveform[c][j]-baseline[c]);
				//double w2 = fabs(waveform[c][j+1]-baseline[c]);
				//double w3 = fabs(waveform[c][j-1]-baseline[c]);
				if (w > peak_val[c]) {

					peak_val[c] = w;
					peak_idx[c] = j;
				}
			}
		}

		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			upsampledata( waveform,c,peak_idx[c]-2,N, L, T);
			upsampledata( waveform,c,peak_idx[c],N, L, T);
			upsampledata( waveform,c,peak_idx[c]-1,N, L, T);
			upsampledata( waveform,c,peak_idx[c]+1,N, L, T);
			upsampledata( waveform,c,peak_idx[c]+2,N, L, T);
			for (int j=peak_idx[c]-2;j<=peak_idx[c]+2;j++){
				for(int m=0;m<int(N);m++){
					if(fabs(upsampledpoints[j][m]-baseline[c])>fabs(peak_val[c])){
						peak_val[c] = fabs(upsampledpoints[j][m]-baseline[c]);
					}
				}
			}
			hPulseHeight[c]->Fill(peak_val[c]);

			// integrate the whole signal
			// FIXME: There should be a smarter way to do this
			//TGraph g(1024-2*skip, time[c], waveform[c]);

			//std::stringstream ss;
			//ss << "waveform_area_func_event" << i << "_channel" << c+1;
			//TF1 f(ss.str().c_str(), [&](double *x, double *){ return g.Eval(x[0]); }, time[c][skip], time[c][1024-skip], 0);

			//hPulseArea[c]->Fill(f.Integral(time[c][skip], time[c][1024-skip], 1e-3));
		}

		hDPeakIndex->Fill(peak_idx[1] - peak_idx[0]);

		double frac_time[4] = { 1000, 1000, 1000, 1000 }; // Beyond 200ns window

		// Interpolate to find base

		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			double vf = baseline[c]+(peak_val[c]-baseline[c]) * cfd_frac;

			for (int j=peak_idx[c] +	1; j>skip 	; j--) {
				//printf("Peak value: %f, wvfm value: %f, j: %d\n",fabs(peak_val[c]),fabs(waveform[c][peak_idx[c]]),j);

				if (fabs(waveform[c][j]) < vf) {
					if (FIRFilterBool == 0) {
						for (int p=-interp_pts_down; p<=interp_pts_up; p++) {
							int idx = j + p;
							double t = time[c][idx];
							double v = fabs(waveform[c][idx]);
							// swap x and y because we can't eval on y
							interp_graph->SetPoint(interp_pts_down+p, v, t);
							}
						double t = -1000;
						if (interp_type == 1) {
							t = interp_graph->Eval(vf, 0, "S");
						} else {
							t = interp_graph->Eval(vf);
						}

						frac_time[c] = t;
						break;
					}
					//Start of SINC Upsampling
					if (FIRFilterBool==1){
							upsampledata(waveform,c,j,N, L, T);
							upsampletime(time,c, j,N);
							double t = -1000;
							for (int n=int(N-1); n >=0; n--) {
								//printf("Point: %f,vf: %f\n",fabs(upsampledpoints[j][n]),vf);
								if (fabs(upsampledpoints[j][n]) < vf) {
									double Y1 = fabs(upsampledpoints[j][n]);
									double Y2 = fabs(upsampledpoints[j][n+1]);
									double X1 = upsampledtimes[j][n];
									double X2 = upsampledtimes[j][n+1];
									double slope = (Y2 - Y1) / (X2 - X1);
									t = (vf - Y1) / slope + X1;
									// for (int p=-interp_pts_down; p<=interp_pts_up; p++) {
									// 	int idx = n + p;
									// 	double t = upsampledtimes[j][idx];
									// 	double v = fabs(upsampledpoints[j][idx]);
									// 	interp_graph->SetPoint(interp_pts_down+p, v, t);
									// }
									//
									// if (interp_type == 1) {
									// 	t = interp_graph->Eval(vf, 0, "S");
									// } else {
									// 	t = interp_graph->Eval(vf);
									// }
									frac_time[c] = t;
									//printf("Time=%f\n",t);
									break;
								}
							}
							if (t != -1000){
								//printf("Time=%f\n",t);
								//frac_time[c] = t;
								break;
							}

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
				double da = (baseline[c1] - baseline[c2]);
				dtda->Fill(da,t);
				printf("%d chs. %d and %d t1=%f t2=%f dt=%f\n", i, c1+1, c2+1, frac_time[c1], frac_time[c2], t);
				if (floateq(t, 0)) {
					//printf("%d dt%d%d = 0\n", i, c1+1, c2+1);
				}
				dt.push_back(t);
			}
		}
		for (int j=0; j<dt.size(); j++) {
			hdt[j]->Fill(dt[j]);
		}

		double deltat = frac_time[1] - frac_time[0];
		if (deltat > 0.4) {
			TGraph *hpulse1 = new TGraph(1024,time[0],waveform[0]);
			TGraph *hpulse2 = new TGraph(1024,time[1],waveform[1]);
			//TGraph *hpulse3 = new TGraph(1,time[1],waveform[1]);
			TLine *line1 = new TLine(0,cfd_frac*peak_val[0],400,cfd_frac*peak_val[0]);
			TLine *line2 = new TLine(0,cfd_frac*peak_val[1],400,cfd_frac*peak_val[1]);
			TLine *line3 = new TLine(frac_time[0],-1,frac_time[0],1);
			TLine *line4 = new TLine(frac_time[1],-1,frac_time[1],1);
			// for (int idx=0; idx<1024; idx++) {
			// 	hpulse1->Fill(time[0][idx], waveform[0][idx]);
			// 	hpulse2->Fill(time[1][idx], waveform[1][idx]);
			// }
			printf("%d tf1=%f tf2=%f tp1=%f tp2=%f dt=%f\n", i, frac_time[0], frac_time[1], peak_val[0], peak_val[1], deltat);
			stop_doing_stuff = 1;
			TCanvas *c = new TCanvas("c", "c");
			hpulse1->Draw("AC*");
			hpulse1->SetLineColor(4);
			///hpulse2->SetLineColor(1);
			hpulse2->SetLineWidth(3);
	    hpulse2->SetMarkerStyle(21);
	    hpulse2->SetLineColor(2);
	    hpulse2->Draw("CP");
			line1->Draw();
			line2->Draw();
			line2->SetLineColor(2);
			line1->SetLineColor(4);
			line3->Draw();
			line4->Draw();
			line3->SetLineColor(2);
			line4->SetLineColor(4);
			for (int i=0; i<1024; i++) {
				ch1waveforms_right->Fill(i, waveform[0][i]);
				ch2waveforms_right->Fill(i, waveform[1][i]);
			}
		} else if (deltat < -0.4) {
			for (int i=0; i<1024; i++) {
				ch1waveforms_left->Fill(i, waveform[0][i]);
				ch2waveforms_left->Fill(i, waveform[1][i]);
			}
		} else {
			for (int i=0; i<1024; i++) {
				ch1waveforms_center->Fill(i, waveform[0][i]);
				ch2waveforms_center->Fill(i, waveform[1][i]);
			}
		}
	}
	puts("100%\n");

	if (stop_doing_stuff) {

		return;
	}

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

	TCanvas *cdiff = new TCanvas("cdiff", "Channel 1 - Channel 2");
	ht12difference->GetXaxis()->SetTitle("time [ns]");
	ht12difference->GetYaxis()->SetTitle("#Deltat");
	ht12difference->Draw("COLZ");

	TCanvas *cdiff2 = new TCanvas("cdiff2", "Channel 1 - Channel 2");
	hdtime12->Draw();
	TCanvas *cdtda = new TCanvas("cdtda", "dt v da");
	dtda->GetXaxis()->SetTitle("#DeltaAmplitude");
	dtda->GetYaxis()->SetTitle("#Deltat");
	dtda->Draw("COLZ");

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

	TCanvas *c1l = new TCanvas("c1l", "Channel 1 Waveforms (left lobe)");
	c1l->SetLogz(1);
	ch1waveforms_left->GetXaxis()->SetTitle("index");
	ch1waveforms_left->GetYaxis()->SetTitle("voltage [V]");
	ch1waveforms_left->Draw("COLZ");
	TCanvas *c1c = new TCanvas("c1c", "Channel 1 Waveforms (center lobe)");
	c1c->SetLogz(1);
	ch1waveforms_center->GetXaxis()->SetTitle("index");
	ch1waveforms_center->GetYaxis()->SetTitle("voltage [V]");
	ch1waveforms_center->Draw("COLZ");
	TCanvas *c1r = new TCanvas("c1r", "Channel 1 Waveforms (left lobe)");
	c1r->SetLogz(1);
	ch1waveforms_right->GetXaxis()->SetTitle("index");
	ch1waveforms_right->GetYaxis()->SetTitle("voltage [V]");
	ch1waveforms_right->Draw("COLZ");

	TCanvas *cdpi = new TCanvas("cdpi", "cdpi");
	hDPeakIndex->Draw();

	/*
	TCanvas *cpa = new TCanvas("cpa", "Pulse Area");
	double hpamax = 0;
	for (int i=0; i<4; i++) {
		double m = hPulseArea[i]->GetMaximum();
		if (m > hpamax) {
			hpamax = m;
		}
	}
	hPulseArea[0]->SetLineColor(1);
	hPulseArea[0]->SetMaximum(hpamax+hpamax*10);
	hPulseArea[0]->GetXaxis()->SetTitle("Area [Vns]");
	hPulseArea[0]->GetYaxis()->SetTitle("frequency");
	hPulseArea[0]->Draw();
	for (int i=1; i<4; i++) {
		hPulseArea[i]->SetLineColor(i+1);
		hPulseArea[i]->Draw("same");
	}
	cpa->BuildLegend();
	*/

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
