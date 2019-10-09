#include <TROOT.h>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "sgfilter.h"

// Histogram parameters
double hist_ylo[4] = { -0.6, -0.6, -0.6, -0.6 }; // Volts
double hist_yhi[4] = { 0.1, 0.1, 0.1, 0.1 }; // Volts
double pulsearea_max = 10;

// Settings
double cfd_frac = 0.50;// Fraction of height to measure to
int skip = 32; // Skip start and end of signal
double repeat_thresh = 0.049; // When checking for clipping, signal must be > this threshold
unsigned int repeat_max = 4; // Discard if (# repeated values) > repeat_max
//double signal_thresh[4] = { 0.025, 0.02, 0.025, 0.025 }; // Discard if whole signal < signal_thresh
double signal_thresh[4] = { 0.01, 0.01, 0.01, 0.01 }; // Discard if whole signal < signal_thresh
int interp_type = 0; // 0 is linear, 1 is cubic spline
int smoothing = 1; // true to smooth

// Floating-point equality check
const double epsilon = 0.0001;
#define floateq(A,B) (fabs(A-B) < epsilon)

unsigned int active_channels[] = { 1, 2, 3 };
size_t nactive_channels = 3; // Length of the above array

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

	TH2D *ch1waveforms = new TH2D("ch1waveforms", "Channel 1 Waveforms", 1024, -0.1, 200, 512, hist_ylo[0], hist_yhi[0]);
	TH2D *ch2waveforms = new TH2D("ch2waveforms", "Channel 2 Waveforms", 1024, -0.1, 200, 512, hist_ylo[1], hist_yhi[1]);
	TH2D *ch3waveforms = new TH2D("ch3waveforms", "Channel 3 Waveforms", 1024, -0.1, 200, 512, hist_ylo[2], hist_yhi[2]);
	TH2D *ch4waveforms = new TH2D("ch4waveforms", "Channel 4 Waveforms", 1024, -0.1, 200, 512, hist_ylo[3], hist_yhi[3]);
	TH2D *hwaveforms[] = { ch1waveforms, ch2waveforms, ch3waveforms, ch4waveforms };

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
			ss1 << "dt" << c1-1 << c2-1;
			std::stringstream ss2;
			ss2 << "#Deltat_{" << c1-1 << c2-1 << "}";
			//TH1D *h = new TH1D(&name[0], &name[0], 1024, 0, 200);
			//TH1D *h = new TH1D(&name[0], &name[0], 2048, -60, 60);
			TH1D *h = new TH1D(ss1.str().c_str(), ss2.str().c_str(), 2048, -60, 60);
			hdt.push_back(h);

			//printf("%d %d (%d %d) %s\n", i, j, c1, c2, name);
		}
	}

	unsigned int ndiscarded = 0;

	long long nentries = tree->GetEntries();
	printf("Processing %lld entries\n", nentries);

	// For interpolating to find the constant-fraction time
	int interp_pts_up = 1; // # points above the principle point
	int interp_pts_down = 1; // # points below the principle point
	int ninterp_pts = interp_pts_up + 1 + interp_pts_down;
	TGraph *interp_graph = new TGraph(ninterp_pts);

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
			if (smoothing) {
				SGSmoothing::Smooth(1024, &raw_waveform[c][0], &waveform[c][0], 5, 3);
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

		for (int j=0+skip; j<1024-skip; j++) {
			for (int k=0; k<nactive_channels; k++) {
				unsigned int c = active_channels[k] - 1;
				hwaveforms[c]->Fill(time[c][j], waveform[c][j]);
				double w = fabs(waveform[c][j]);
				if (w > peak_val[c]) {
					peak_val[c] = w;
					peak_idx[c] = j;
				}
			}
		}

		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;

			hPulseHeight[c]->Fill(peak_val[c]);

			// integrate the whole signal
			// FIXME: There should be a smarter way to do this
			//TGraph g(1024-2*skip, time[c], waveform[c]);

			//std::stringstream ss;
			//ss << "waveform_area_func_event" << i << "_channel" << c+1;
			//TF1 f(ss.str().c_str(), [&](double *x, double *){ return g.Eval(x[0]); }, time[c][skip], time[c][1024-skip], 0);

			//hPulseArea[c]->Fill(f.Integral(time[c][skip], time[c][1024-skip], 1e-3));
		}

		double frac_time[4] = { 1000, 1000, 1000, 1000 }; // Beyond 200ns window

		// Interpolate to find base

		for (int k=0; k<nactive_channels; k++) {
			unsigned int c = active_channels[k] - 1;
			double vf = peak_val[c] * cfd_frac;
			for (int j=peak_idx[c]; j>skip-1; j--) {
				if (fabs(waveform[c][j]) <= vf) {
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
				if (j == skip) {
					printf("WARNING: %d: Failed to find fraction (ch%d)\n", i, c+1);
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
				if (floateq(t, 0)) {
					printf("%d dt%d%d = 0\n", i, c1+1, c2+1);
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
		hdt[i]->SetLineColor(i+1);
		hdt[i]->Draw("same");
	}
	cdt->BuildLegend();

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

	TCanvas *cpa = new TCanvas("cpa", "Pulse Area");
	double hpamax = 0;
	for (int i=0; i<4; i++) {
		double m = hPulseArea[i]->GetMaximum();
		if (m > hpamax) {
			hpamax = m;
		}
	}
	hPulseArea[0]->SetLineColor(1);
	hPulseArea[0]->SetMaximum(hpamax+hpamax*0.1f);
	hPulseArea[0]->GetXaxis()->SetTitle("Area [Vns]");
	hPulseArea[0]->GetYaxis()->SetTitle("frequency");
	hPulseArea[0]->Draw();
	for (int i=1; i<4; i++) {
		hPulseArea[i]->SetLineColor(i+1);
		hPulseArea[i]->Draw("same");
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
