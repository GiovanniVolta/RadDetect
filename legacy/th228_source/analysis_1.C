#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFileMerger.h"
#include "TAxis.h"
#include "TText.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <ctime>


int analysis(){

	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	std::string file_name = "Rn09022015";


	TFile *f = TFile::Open((base_file_path+file_name+"/"+file_name+base_file_ending).c_str());
	TTree *t = (TTree*) f->Get("t");

	TCanvas *c1 = new TCanvas("c1");

	t->Draw("channel:runtime", "runtime>4 &&runtime< 12");

	TCanvas *c2 = new TCanvas("c2");

	t->Draw("channel>>htmp(800,0,800)", "runtime>4 &&runtime< 12");
	c2->SetLogy();

	return 0;

}



