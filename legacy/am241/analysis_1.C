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
#include "TGraph.h"
#include "TGraphErrors.h"

#include <string>
#include <vector>
#include <utility>
#include <sstream>


int analysis_1(){

	// define the run names for each run
	std::vector<std::string> run_labels;
	std::vector<int> colors;
	std::vector<int> n_bins;
	std::vector<std::string> filenames;
	std::vector<TTree*> trees;
	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	int n_bins_overview = 500;
	double reject_ingrowth = 0.2;

	//First define the runs of interest and put them in the vectors

	run_labels.push_back("normal");
	filenames.push_back("Ap21112019_normal");

        run_labels.push_back("improved");
        filenames.push_back("Ap21112019_impr");


        for(unsigned int i = 0; i < filenames.size(); ++i){

                std::string this_file_name = base_file_path + filenames[i] + "/" + filenames[i] + base_file_ending;

                TFile *this_file = TFile::Open(this_file_name.c_str());
                TTree* this_tree = (TTree*) this_file->Get("t");
                trees.push_back(this_tree);
	}

	// Now the TTrees are loaded and ready
	
	/*****************
 	 * OVERVIEW PLOT *
	 *****************/

	TCanvas *c_overview = new TCanvas("overview");
	std::vector<TH1D*> histos;
	
	for(unsigned int i = 0; i < trees.size(); ++i){

		
		trees[i]->Draw(("channel>>hist_"+std::to_string(i)+"(500,0,500)").c_str());
		histos.push_back((TH1D*) c_overview->GetPrimitive(("hist_"+std::to_string(i)).c_str()));
		histos[i]->Scale(1./histos[i]->Integral(200,500));
		//histos[i]->Draw("SAME");
	}

        TCanvas *c_comp= new TCanvas("comp");
	histos[0]->Draw();
	//histos[1]->Draw("SAME");
	//histos[1]->SetLineColor(2);

	TF1 *gaus = new TF1("fit_func", "gaus", 340, 420);
	histos[0]->Fit(gaus, "", "", 340,420);

	histos[0]->GetXaxis()->SetRangeUser(200,430);
	c_comp->SetLogy();	



	return 0;

}



