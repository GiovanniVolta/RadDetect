#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFileMerger.h"

#include <string>
#include <sstream>


double FitFunc212Po(double* x, double* par){
// This function is used to fit the ingrow of 212Po due to emanation
// Code stolen from Guillaume's script

//	int toffs = par[1]*60;	
	int toffs = -0;
	double lam_b = 0.000018098;
	double lam_c = 0.00019095;
	double lam_d = 2318218.;

	double A = ((lam_c*lam_d)/((lam_b-lam_c)*(lam_d-lam_b)))*(exp(-lam_b*(x[0]*60+toffs))-exp(-lam_d*(x[0]*60+toffs)));
	double B = ((lam_b*lam_d)/((lam_b-lam_c)*(lam_c-lam_d)))*(exp(-lam_c*(x[0]*60+toffs))-exp(-lam_d*(x[0]*60+toffs)));
	double C = (1-exp(-lam_d*(x[0]*60+toffs)));
	double X = par[0]*(A+B+C);
	
// describe the 224Ra decay
//	double lambRa = 2.22848e-06;
	double lambRa = par[1];
	double dec = exp(-lambRa*(x[0]*60+toffs));

	return X*dec;
//	return X;
}

double FitFuc216Po(double* x, double* par){

	int toffs = -0;
	double lambRa = par[1];
	double dec = exp(-lambRa*(x[0]*60+toffs));

	return par[0]*dec;
}

int analysis(){

	TFileMerger *m = new TFileMerger();
	m->OutputFile("/home/lin/fjoerg/rnmonanalysis/radium224_implantation/merged.root");
        m->AddFile("/home/lin/fjoerg/rnmondatabase/backup/Rn21032019/Rn21032019.root");
	m->AddFile("/home/lin/fjoerg/rnmondatabase/backup/Rn22032019/Rn22032019.root");

	m->Merge();
	

	TFile *fIn = TFile::Open("merged.root");
	TTree *t = (TTree*) fIn->Get("t");

	TCanvas *c1 = new TCanvas("canv");


//	t->Draw("(timestamp-1499.69E6)/60.>>htemp(1000, 0, 70000)", "channel>390&&channel<550");
	double unix_min = t->GetMinimum("timestamp");
	std::cout<< unix_min <<std::endl;
	std::string min_timestamp = "1553.15E6";
	
//	t->Draw(("(timestamp-"+min_timestamp+")/60>>htemp(50, 0, 16000)").c_str(), "channel>270&&channel<330");	//212Po
        t->Draw(("(timestamp-"+min_timestamp+")/60>>htemp(100,0,15000)").c_str(), "channel>270&&channel<330"); //212Po
	//return 0;
//	t->Draw("(timestamp-1507.25E6)/60>>htemp(100, 0, 16000)", "channel>310&&channel<340");	//216Po

	TH1D *hist = (TH1D*)c1->GetPrimitive("htemp");
	hist -> GetXaxis() -> SetTitle("Runtime (min)");
	double width = hist->GetBinWidth(2);

	std::stringstream yTitle_stream;	
	yTitle_stream<<"Counts/"<<width<<" (1/min)";
	std::string yTitle = yTitle_stream.str();
	std::cout<<yTitle<<std::endl;

	hist -> GetYaxis() -> SetTitle(yTitle.c_str());
	hist->SetTitle("Detected activity in 212Po");
	
	TF1 *fitFunc = new TF1("212PoFit", FitFunc212Po,  0, 20E3, 2);
//	TF1 *fitFunc = new TF1("216PoFit", FitFuc216Po,  0, 20E3, 2);
//	fitFunc->SetParameter(1, -200);
	fitFunc->SetParameter(1, 2.22848e-06);

	hist->Fit(fitFunc, "ll", 0, 4E3);
	fitFunc->Draw("SAME");

	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(1111);


	return 1;	

}
