#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFileMerger.h"
#include "TAxis.h"
#include "TText.h"
#include "TBox.h"

#include <string>
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
	
	int toffs = -par[2];
	double lam_b = 0.000018098*60.;
	double lam_c = 0.00019095*60.;
	double lam_d = 2318218.*60.;

	double A = ((lam_c*lam_d)/((lam_b-lam_c)*(lam_d-lam_b)))*(exp(-lam_b*(x[0]+toffs))-exp(-lam_d*(x[0]+toffs)));
	double B = ((lam_b*lam_d)/((lam_b-lam_c)*(lam_c-lam_d)))*(exp(-lam_c*(x[0]+toffs))-exp(-lam_d*(x[0]+toffs)));
	double C = (1-exp(-lam_d*(x[0]*60+toffs)));
	double X = par[0]*(A+B+C);
	
// describe the 224Ra decay
//	double lambRa = 2.22848e-06;
	double lambRa = par[1];
	double dec = exp(-lambRa*(x[0]+toffs));

	return X*dec;
//	return X;
}

double FitFuc216Po(double* x, double* par){

	int toffs = -0;
	double lambRa = par[1];
	double dec = exp(-lambRa*(x[0]*60+toffs));

	return par[0]*dec;
}

int analysis_2(){


        TFileMerger *m = new TFileMerger();
        m->OutputFile("/home/lin/fjoerg/rnmonanalysis/radium224_implantation/merged.root");
        m->AddFile("/home/lin/fjoerg/rnmondatabase/backup/Rn21032019/Rn21032019.root");
        m->AddFile("/home/lin/fjoerg/rnmondatabase/backup/Rn22032019/Rn22032019.root");

        m->Merge();


        TFile *fIn = TFile::Open("merged.root");
        TTree *t = (TTree*) fIn->Get("t");

        TCanvas *c1 = new TCanvas("canv");


	std::string min_timestamp = "1553.15E6";
	t->Draw(("(timestamp-"+min_timestamp+")/60>>htemp(100,0,15000)").c_str(), "channel>270&&channel<330"); //212Po



	TH1D *hist = (TH1D*)c1->GetPrimitive("htemp");
	hist -> GetXaxis() -> SetTitle("Runtime (min)");
	double width = hist->GetBinWidth(2);

	std::stringstream yTitle_stream;	
	yTitle_stream<<"detected Po-212 counts/"<<width<<"min";
	std::string yTitle = yTitle_stream.str();
	std::cout<<yTitle<<std::endl;

	hist -> GetYaxis() -> SetTitle(yTitle.c_str());
	hist->SetTitle("Sample: IRINA\\_Al\\_1");
	
	TF1 *fitFunc = new TF1("212PoFit", FitFunc212Po,  0, 20E3, 2);
//	TF1 *fitFunc = new TF1("216PoFit", FitFuc216Po,  0, 20E3, 2);
//	fitFunc->SetParameter(1, -200);
	fitFunc->SetParameter(1, 1.3709e-4);

	hist->Fit(fitFunc, "llim", 0, 7E3);
//        hist->Fit(fitFunc, "ll", 0, 7E3);

	fitFunc->Draw("SAME");
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


	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);

	c1 -> SaveAs("irina_ss_uncoated.tex");

	return 1;	

}
