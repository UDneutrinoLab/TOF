/*
 * decode.c
 * 2018-01-07
 *
 * Extracts voltage and timing data from DRS4 .dat data, saves as root file
 * It was adapted from a macro written by Stefan Ritt of PSI which can be found here:
 *   https://midas.psi.ch/elogs/DRS4+Forum/361
 *
 * It was adapted from an adaptation of the above by:
 *   Primary Author:  Abhishek Rajput abhi_rajput5@utexas.edu
 *   Other Author(s): Will Flanagan   will.flanagan@utexas.edu
 * The above version can be found here:
 *   https://github.com/UTKLgroup/DRS4macros
 *
 * This version by Aidan Medcalf amedcalf@udallas.edu
 *
 * usage:
 *   $ root -l 'decode.c("file.dat")'
 *     or
 *   root [0] .L decode.c+
 *   root [1] decode("filename.dat")
 * filename.root will be written to the same directory as the binary file upon completion of decoding
*/

#include <cstring>
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include <string.h>
#include <stdio.h>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Getline.h"
#include <stdlib.h>
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TLine.h"
#include <sstream>
#include <vector>
#include "TLegend.h"

using namespace std;

typedef struct {
	char           time_header[4];
	char           bn[2];
	unsigned short board_serial_number;

} THEADER;

typedef struct {
	char           event_header[4];
	unsigned int   event_serial_number;
	unsigned short year;
	unsigned short month;
	unsigned short day;
	unsigned short hour;
	unsigned short minute;
	unsigned short second;
	unsigned short millisecond;
	unsigned short reserved1;
	char           bs[2];
	unsigned short board_serial_number;
	char           tc[2];
	unsigned short trigger_cell;
} EHEADER;

typedef struct {
	char hdr[4];
	char binary_version;
	unsigned short voltage[1024];
	double waveform[4][1024], time[4][1024];
	float bin_width[4][1024];
	vector<int> inpch;
	int ext;
} DATAINFO;


typedef struct {
	int ev;
	string query0;
} QUERYUSER;

/**********************************************************************************************************************************/

DATAINFO find_header(FILE *f);
DATAINFO output_dat(FILE *f, DATAINFO dat);
QUERYUSER queries();

void decode(const char *filename)
{
	THEADER th;
	EHEADER eh;
	int counter[4] = { 0, 0, 0, 0 };

	// Open the binary waveform file
	FILE *f = fopen(Form("%s", filename), "r");
	if (f == NULL) {
		printf("Cannot find file \'%s\'\n", filename);
		return;
	}

	// Pass file to read_header function
	DATAINFO datfh = find_header(f);

	// Query user for information
	QUERYUSER quer = queries();

	// Loop over events
	printf("\nDecoding binary data...\n");

	TTree *tree = new TTree("tree", "Decoded DRS4 data");
	double waveforms[4][1024];
	double time[4][1024];
	tree->Branch("ch0waveforms", &waveforms[0][0], "ch0waveforms[1024]/D");
	tree->Branch("ch1waveforms", &waveforms[1][0], "ch1waveforms[1024]/D");
	tree->Branch("ch2waveforms", &waveforms[2][0], "ch2waveforms[1024]/D");
	tree->Branch("ch3waveforms", &waveforms[3][0], "ch3waveforms[1024]/D");
	tree->Branch("ch0time", &time[0][0], "ch0time[1024]/D");
	tree->Branch("ch1time", &time[1][0], "ch1time[1024]/D");
	tree->Branch("ch2time", &time[2][0], "ch2time[1024]/D");
	tree->Branch("ch3time", &time[3][0], "ch3time[1024]/D");

	for (int n=0; n < quer.ev; n++) {
		// Use output_dat function to decode binary data for each event and output waveform and time arrays for ROOT purposes
		// There is probably a faster way to do this
		DATAINFO dat1 = output_dat(f, datfh);

		if (dat1.ext < 1) {
			break;
		}

		// Loop over all the channels for each event
		for (int chn = 0; chn < (datfh.inpch).size(); chn++) {
			int c = datfh.inpch[chn];
			for (int i=0; i<1024; i++) {
				waveforms[c][i] = dat1.waveform[c][i];
				time[c][i] = dat1.time[c][i];
			}
			counter[chn]++;
		}
		tree->Fill();
	}

	for (int chn = 0; chn < (datfh.inpch).size(); chn++) {
		printf("\n%d events in channel %d decoded out of %d events \n", counter[chn],datfh.inpch[chn]+1,quer.ev);
	}

	// Open the root file
	// This is ugly and fragile
	char rootfile[1024];
	strcpy(rootfile, filename);
	if (strchr(rootfile, '.'))
		*strchr(rootfile, '.') = 0;
	strcat(rootfile, ".root");
	TFile *outfile = new TFile(rootfile, "RECREATE");

	tree->Write();

	printf("\n%s written\n", rootfile);
	outfile->Close();
}

/*****************************************************************************************************************************************************/

DATAINFO find_header(FILE *f)
{
	THEADER th;
	EHEADER eh;
	DATAINFO dat;
	int version = 0; // Zero if binary is old format, version number if new format
	char newvhdr[4]; // 'D' 'R' 'S' 'N' if new format, where N is the version

	// WARNING: there is currently no error checking on any file operations

	// detect version
	fread(&newvhdr[0], 4*sizeof(char), 1, f);
	// not elegant, but it works
	// if new format, set version, else rewind file and proceed as old format
	if (newvhdr[0]=='D' && newvhdr[1]=='R' && newvhdr[2]=='S') {
		version = newvhdr[3] - '0';
	} else {
		version = 0;
		fseek(f, 0, SEEK_SET);
	}
	dat.binary_version = version;

	// read time header
	/*sizeof(th) is 8
	  sizeof(hdr) is 4
	  sizeof(float) is 4
	  sizeof(eh) is 32
	  sizeof(short) is 2 */

	int t = fread(&th, sizeof(th), 1, f);

	//printf("Found data for board #%d\n", th.board_serial_number);

	// read time bin widths

	memset(dat.bin_width, sizeof(dat.bin_width), 0);

	for (int ch=0; ch<5 ; ch++) {
		int d = fread(dat.hdr, sizeof(dat.hdr), 1, f);

		if (dat.hdr[0] != 'C') {
			// event header found
			fseek(f, -4, SEEK_CUR);
			break;
		}

		int i = dat.hdr[3] - '0' - 1;
		printf("\n Found timing calibration for channel #%d\n\n", i+1);
		int b = fread(&dat.bin_width[i][0], sizeof(float), 1024, f);

		(dat.inpch).push_back(i);
	}

	return dat;
}

/*****************************************************************************************************************************************************/

DATAINFO output_dat(FILE *f, DATAINFO dat)
{
	THEADER th;
	EHEADER eh;
	DATAINFO dat1;

	// read event header

	dat1.ext = fread(&eh, sizeof(eh), 1, f);

	//printf("Found event #%d\n", eh.event_serial_number);

	// reach channel data

	for (int ch=0; ch<5; ch++) {

		int i = fread(dat1.hdr, sizeof(dat1.hdr), 1, f);

		if (i < 1)
			break;

		if (dat1.hdr[0] != 'C') {
			// event header found
			fseek(f, -4, SEEK_CUR);
			break;
		}

		int chn_index = dat1.hdr[3] - '0' - 1;

		unsigned int scaler = 1;
		if (dat.binary_version > 0) {
			fread(&scaler, 4, 1, f);
		}

		int k = fread(dat1.voltage, 2, 1024, f);

		for (i=0; i<1024; i++) {
			// convert data to volts
			dat1.waveform[chn_index][i] = -1*(dat1.voltage[i] / 65536. - 0.5);

			// calculate time for this cell
			dat1.time[chn_index][i] = 0;
			for (int j=0; j<i; j++){
				dat1.time[chn_index][i] += dat.bin_width[chn_index][(j+eh.trigger_cell) % 1024];
			}
		}
	}

	// align cell #0 of all channels

	double t1 = dat1.time[dat.inpch[0]][(1024-eh.trigger_cell) % 1024];

	for (int ch=1; ch<(dat.inpch).size(); ch++){
		double t2 = dat1.time[dat.inpch[ch]][(1024-eh.trigger_cell) % 1024];
		double dt = t1 - t2;

		for (int i=0; i<1024; i++){
			dat1.time[dat.inpch[ch]][i] += dt;
		}
	}

	return dat1;

}

/*****************************************************************************************************************************************************/

QUERYUSER queries()
{
	QUERYUSER quer;

	while (true) {
		cout << "How many events would you like to analyze? ";
		getline(cin,quer.query0);
		stringstream ss(quer.query0);
		if (ss >> quer.ev)
			break;
		cout << "Invalid entry." << endl;
	}

	return quer;
}
