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

#include <string>
#include <vector>
#include <utility>
#include <sstream>


double IngrowFunc212Po(double* x, double* par){
// This function is used to fit the ingrow of 212Po due to emanation
// Code stolen from Guillaume's script

//      int toffs = par[1]*60;

        int toffs = par[1];
        double lam_b = 0.000018098*60.;
        double lam_c = 0.00019095*60.;
        double lam_d = 2318218.*60.;

        double A = ((lam_c*lam_d)/((lam_b-lam_c)*(lam_d-lam_b)))*(exp(-lam_b*(x[0]+toffs))-exp(-lam_d*(x[0]+toffs)));
        double B = ((lam_b*lam_d)/((lam_b-lam_c)*(lam_c-lam_d)))*(exp(-lam_c*(x[0]+toffs))-exp(-lam_d*(x[0]+toffs)));
        double C = (1-exp(-lam_d*(x[0]*60+toffs)));
        double X = par[0]*(A+B+C);

	return X;
}



double FitFunc212Po(double* x, double* par){
// This function is used to fit the ingrow of 212Po due to emanation
// Code stolen from Guillaume's script

//	int toffs = par[1]*60;	
	
//	int toffs = -par[2];
	int toffs = 0;

	double lam_b = 0.000018098*60.;
	double lam_c = 0.00019095*60.;
	double lam_d = 2318218.*60.;

	double A = ((lam_c*lam_d)/((lam_b-lam_c)*(lam_d-lam_b)))*(exp(-lam_b*(x[0]+toffs))-exp(-lam_d*(x[0]+toffs)));
	double B = ((lam_b*lam_d)/((lam_b-lam_c)*(lam_c-lam_d)))*(exp(-lam_c*(x[0]+toffs))-exp(-lam_d*(x[0]+toffs)));
	double C = (1-exp(-lam_d*(x[0]+toffs)));
	double X = par[0]*(A+B+C);

	double init_act = par[1]*exp(-lam_b*(x[0]+toffs));
	
// describe the 224Ra decay
	double lambRa = 2.22848e-06*60;
//	double lambRa = par[1];
	double dec = exp(-lambRa*(x[0]+toffs));

	return X*dec+init_act;
//	return X;
}

double FitFunc224Ra(double* x, double* par){
// This function is used to fit the ingrow of 212Po due to emanation
// Code stolen from Guillaume's script

//	int toffs = par[1]*60;	
        int toffs = par[1];
	
// describe the 224Ra decay
	double lambRa = 2.22848e-06*60;
//	double lambRa = par[1];
	double dec = par[0]*exp(-lambRa*(x[0]+toffs));

	return dec;
//	return X;
}



double FitFuc216Po(double* x, double* par){

	int toffs = -0;
	double lambRa = par[1];
	double dec = exp(-lambRa*(x[0]*60+toffs));

	return par[0]*dec;
}

int comparison_1(){

	// define the run names for each run
	std::vector<std::string> run_labels;
	std::vector<std::string> legend_entry;
	std::vector<int> colors;
	std::vector<int> n_bins;
	std::vector<std::pair<int, int>> channel_cuts;
        std::vector<std::pair<int, int>> channel_cuts_bi;
	std::vector<std::vector<std::string>> filenames;
	std::vector<TTree*> merged_trees;
	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	std::vector<std::string> irina_a;
	run_labels.push_back("irina_CF63_1");
	legend_entry.push_back("CF63_test_1");
	irina_a.push_back("Rn16112020");
	filenames.push_back(irina_a);
	channel_cuts.push_back(std::make_pair(220,310));
	channel_cuts_bi.push_back(std::make_pair(270,310));

	n_bins.push_back(981);
	colors.push_back(2);



	for(unsigned int i = 0; i < filenames.size(); ++i){

		std::cout<<"Running measurement: "<<i<<" "+run_labels[i]<<std::endl;
		std::vector<std::string> this_file_names = filenames[i];
		std::string merged_file_name = "/d12/lin/xenon/radon_monitor/rnmonanalysis/radium224_implantation/merged_"+run_labels[i]+".root";
	        TFileMerger *m = new TFileMerger();
	        m->OutputFile(merged_file_name.c_str());
		for(unsigned int j = 0; j < this_file_names.size(); ++j){
			
			std::string this_file_name = base_file_path + this_file_names[j] + "/" + this_file_names[j] + base_file_ending;
			std::cout<<"Running on file: "<<this_file_name<<std::endl;
			m->AddFile(this_file_name.c_str());
		}
		m->Merge();
		TFile *this_file = TFile::Open(merged_file_name.c_str());
		TTree* this_merged_tree = (TTree*) this_file->Get("t");
		merged_trees.push_back(this_merged_tree);

	}

	// Now the merged_trees array has the data available and ready to use
	
	TCanvas *c_comp = new TCanvas("c_comp");
	TLegend* leg = new TLegend(0.7,0.7,0.85,0.85);
	std::vector<long int> unix_starts;
	std::vector<TH1D*> histos;
	for(unsigned int i = 0; i < merged_trees.size(); ++i){
		Long64_t this_unix;
		int this_channel;
		merged_trees[i]->SetBranchAddress("timestamp", &this_unix);
		merged_trees[i]->SetBranchAddress("channel", &this_channel);
		merged_trees[i]->GetEvent(0);
		Long64_t min_unix = this_unix;
                merged_trees[i]->GetEvent(merged_trees[i]->GetEntries()-1);
                Long64_t max_unix = this_unix;
		double runtime = (max_unix-min_unix)/60.;
		//double weight = 1E3/(runtime*60./n_bins[i]);
		double weight = 1/(runtime*60./n_bins[i]);

                TH1D* this_act = new TH1D(run_labels[i].c_str(), "", n_bins[i], 0, runtime);
		histos.push_back(this_act);
		int min_channel = channel_cuts[i].first;
		int max_channel = channel_cuts[i].second;
		int min_channel_bi = channel_cuts_bi[i].first;
                int max_channel_bi = channel_cuts_bi[i].second;

		for(long int j = 0; j < merged_trees[i]->GetEntries(); ++j){
			merged_trees[i]->GetEvent(j);
			// Po212 AND Bi212 selection
//			if (((this_channel > min_channel)&&(this_channel < max_channel)) || (
//				((this_channel > min_channel_bi)&&(this_channel < max_channel_bi))) )

			// Po 212 selection only
			if (((this_channel > min_channel)&&(this_channel < max_channel)))

				this_act->Fill((this_unix-min_unix)/60., weight);
		}
		this_act->Sumw2();
		this_act->Draw("SAME");
		this_act->SetLineColor(colors[i]);
		leg->AddEntry(this_act, legend_entry[i].c_str());
		//break;




		TF1 *fitFunc;
		fitFunc = new TF1(("212PoFit"+std::to_string(i)).c_str(), FitFunc212Po,  -5E3, 20E3, 2);
		fitFunc->SetParLimits(1, 0, 5);
		fitFunc->SetLineColor(colors[i]);

		this_act->Fit(fitFunc, "wllmi", 0, 20E3);
/*
		TH1D* fitFuncConfidence = new TH1D(("conf_band"+std::to_string(i)).c_str(), "", 500, -5E3, 20E3);
		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(fitFuncConfidence, 0.68);
		fitFuncConfidence->SetFillColor(colors[i]);
		fitFuncConfidence->Draw("SAME");
		fitFunc->Draw("SAME");
*/	
/*
	TF1 *expFunc = new TF1("decay", "[0]*exp(-[1]*x)", 0, 6E3);
	expFunc -> SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1));
	expFunc -> SetLineStyle(2);
	expFunc -> Draw("SAME");

	TF1 *ingrFunc = new TF1("ingr", IngrowFunc212Po, 0, 20E3, 2);
	ingrFunc -> SetParameter(0, fitFunc->GetParameter(0));
        ingrFunc -> SetLineStyle(2);
        ingrFunc -> Draw("SAME");

	double a0 = fitFunc->GetParameter(0)/width/60.*1000;
	double a0Err = fitFunc -> GetParError(0) / width/60.*1000;
	double T_half = log(2)/fitFunc -> GetParameter(1)/60./24.;
	double T_halfErr = T_half - log(2)/(fitFunc -> GetParameter(1)-fitFunc->GetParError(1))/60./24.;
	double chi2 = fitFunc -> GetChisquare();
	double ndf = fitFunc -> GetNDF();
	double entries = hist->Integral();
	
	stringstream a0Textstrstr;
	a0Textstrstr<<"$A_0(Po212) = $\\si{("<<a0<<"\\pm"<<a0Err<<")}{mBq}";
	string a0Textstr = a0Textstrstr.str();
	TText *a0Text = new TText(0.5*hist->GetXaxis()->GetXmax(),0.9*fitFunc->GetParameter(0), a0Textstr.c_str());
//	a0Text -> SetTextColor(2);
	a0Text->Draw();

        stringstream ThTextstrstr;
        ThTextstrstr<<"$T_{1/2} = $\\si{("<<T_half<<"\\pm"<<T_halfErr<<")}{days}";
        string ThTextstr = ThTextstrstr.str();
        TText *ThText = new TText(0.6*hist->GetXaxis()->GetXmax(),1.1*fitFunc->Eval(0.6*hist->GetXaxis()->GetXmax()), ThTextstr.c_str());
//        ThText -> SetTextColor(2);
        ThText->Draw();


        stringstream entriesstrstr;
        entriesstrstr << "Entries = "<<entries;
        string entriesstr = entriesstrstr.str();
        ThText -> DrawTextNDC(0.6, 0.7, entriesstr.c_str());

	stringstream chi2strstr;
	chi2strstr << "$\\varchi^2/DoF = $"<<chi2<<" / "<<ndf;
	string chi2str = chi2strstr.str();
	ThText -> DrawTextNDC(0.6, 0.6, chi2str.c_str());
*/






		
	}
	histos[0]->SetTitle("Recoil implanted stainless steel samples (Ra224)");
	//histos[0]->GetXaxis()->SetRangeUser(0,30000);
        histos[0]->GetYaxis()->SetTitle("Po212 Activity (Bq)");
        //histos[0]->GetYaxis()->SetRangeUser(0,3.5);
        histos[0]->GetXaxis()->SetTitle("Runtime (min)");
	
	histos[0]->SetStats(0);
	leg->Draw();


	return 0;

}
