#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCut.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFileMerger.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TText.h"
#include "TPaveText.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TBox.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TVector.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDecompChol.h"
#include "TRandom3.h"
#include "TMath.h"


#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <ctime>
#include <algorithm>

double RN_DECAY_CONST = log(2)/3.824;
TRandom3 *gen = new TRandom3(0);
unsigned int N_ISOTOPES_FIT;
unsigned int N_PAR_CB;
double PI = 3.141592654;
double TINY = 1E-5;

// Vector of ranges to keep blind for the fit
// Needs to be global parameter...
std::vector<std::pair<double, double>> blinding_ranges;


int Tableau20(int i){ 
        std::vector <std::vector<int> > tableau20 = 
                {{31, 119, 180}, {174, 199, 232}, {255, 127, 14}, {255, 187, 120},
             {44, 160, 44}, {152, 223, 138}, {214, 39, 40}, {255, 152, 150},
             {148, 103, 189}, {197, 176, 213}, {140, 86, 75}, {196, 156, 148},
             {227, 119, 194}, {247, 182, 210}, {127, 127, 127}, {199, 199, 199},
             {188, 189, 34}, {219, 219, 141}, {23, 190, 207}, {158, 218, 229}};
    std::vector<int> this_color = tableau20[i];
    Int_t ci = TColor::GetFreeColorIndex();
    TColor *col = new TColor(ci, this_color[0]/255., this_color[1]/255., this_color[2]/255.);
        return ci; 

}


double CrystalBallFunc(double* energy, double* pars){
	// This is the definition of a single crystalball function used to 
	// describe alpha peak of one isotope

	// Paramters:	pars[0]:	Area 
	//		pars[1]:	Alpha
	//		pars[2]:	Mean
	//		pars[3]:	Sigma
	//		pars[4]:	N

	double gaus = exp(-(energy[0]-pars[2])*(energy[0]-pars[2])/(2*pars[3]*pars[3]));
	double poly = pow(pars[4]/abs(pars[1]), pars[4])*exp(-pars[1]*pars[1]/2.);
	poly *= pow(pars[4]/abs(pars[1]) - abs(pars[1]) - (energy[0]-pars[2])/pars[3], -pars[4]);
	double norm = pars[4]/abs(pars[1]) * 1./(pars[4]-1)*exp(-pars[1]*pars[1]/2.);
	norm += sqrt(PI/2.)*(1+erf(abs(pars[1])/sqrt(2)));
	norm = 1./(norm*pars[3]);
	if((energy[0]-pars[2])/pars[3] > -pars[1]){
		return norm*pars[0]*gaus;
	}else{
		return norm*pars[0]*poly;
	}
}


double CompleteSpectrum(double *energy, double *allParams){
// This function builds a joint function in order to fit the complete spectrum of peaks

	const unsigned int nParams = N_PAR_CB;
	double result = 0;

	// is this energy value in any exlusion interval?
	// Then tell TF1 to reject this point, but continue computation
	for(unsigned int i = 0; i < blinding_ranges.size(); ++i){
		if(energy[0] > blinding_ranges[i].first && energy[0] < blinding_ranges[i].second){
			TF1::RejectPoint();
		}
	}

	for(unsigned int i = 0; i < N_ISOTOPES_FIT; ++i){
		double *thisParams = (double*) malloc(nParams*sizeof(double));
		for(unsigned int j = 0; j < nParams; ++j){
			thisParams[j] = allParams[i*nParams + j];
		}
		result += CrystalBallFunc(energy, thisParams);
	}
	return result;
}


TF1* GetRandomFunction(TMatrixD u_matrix, TVectorT<double> mean_vec){

	// This function will generate a random TF1 func, where the parameters are
	// Set according to a Multivariate Gaussian, based on the Cholesky decomposed 
	// Covariance Matrix Cov = U^T U, and the vector containing the mean parameter position
	unsigned int total_n_par = N_ISOTOPES_FIT * N_PAR_CB;
	
	TF1 *result = new TF1("", CompleteSpectrum, 0, 500, total_n_par);
	int n_pars = result->GetNpar();
	if(n_pars != mean_vec.GetNrows())
		return NULL;
	// Generate that many samples of a normal Gaussian (mu = 0, sigma = 1)
	TVectorT<double> normal_samples;
	normal_samples.ResizeTo(u_matrix.GetNrows());
	for(int i = 0; i < normal_samples.GetNrows(); ++i){
		normal_samples[i] = gen->Gaus();
	}
	// Now transform this sample to the desired multivariate distribution
	// new_par = mu_vec + U * normal_sample
	// If necessary (D(mean_vec) > D(normal_samples), fill it with zeros
	TVectorT<double> rand_vec = (u_matrix * normal_samples);
	rand_vec.ResizeTo(n_pars);
	TVectorT<double> new_par = mean_vec + rand_vec;
	// Set these as parameters for the TF1
	result->SetParameters(new_par.GetMatrixArray());

	return result;

}

bool grater_as (double i, double j) {
  return (i >= j);
}


int analysis(){


	TF1* cryst = new TF1("test_func", CrystalBallFunc, -1000, 1000, 5);

	// define names of the files etc.

	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	std::vector<std::string> filenames;
        std::vector<std::string> run_labels;
	std::vector<std::string> isotope_labels;
	std::vector<double> iso_energies;
	std::vector<double> iso_energies_err;
	std::vector<double> iso_energies_indicator;
	std::vector<bool> fit_isotope;
	N_PAR_CB = 5; // Number of parameters for a CB function
	std::vector<int> isotope_color_tab;
	std::vector<int> isotope_color;
	std::vector<int> line_style;
	std::vector<int> true_index;
	// 3d Matrix storing [run, isotope, fit_param]
	std::vector <double> *isotope_fit_vals = new std::vector<double>;

	std::vector<double> tmp_line_pos;

	std::vector<TTree*> trees;
	std::vector<TH1D*> spectra;
	TMatrixDSym cov_mat;
	std::vector<TF1*> fit_funcs;
	std::vector<TF1*> fit_funcs_cp;
	std::vector<TMatrixDSym> cov_matrices;
	// int spec_range[2] = {230, 540};
	int spec_range[2] = {255, 530};
	int fit_ranges[2] = {220, 505};
	double aspect_point = 0.25;
	// double aspect_point = 0.33;


	blinding_ranges.push_back(std::make_pair(295, 304));
	blinding_ranges.push_back(std::make_pair(363, 367));
	blinding_ranges.push_back(std::make_pair(386, 408));
	blinding_ranges.push_back(std::make_pair(430, 442));
	blinding_ranges.push_back(std::make_pair(502, 539));
	
	//int hist_colors[2] = {Tableau20(8),Tableau20(6)};
	int hist_colors[2] = {Tableau20(0),Tableau20(6)};

	//int fit_colors[2] = {Tableau20(8+1),Tableau20(6+1)};
	int fit_colors[2] = {Tableau20(1),Tableau20(6+1)};


	int n_bins = spec_range[1] - spec_range[0];
	std::vector<double> runtimes;
	std::vector<double> min_runtime_cut;

// 226Ra subdominant alpha line ~5%
	isotope_labels.push_back("{}^{226}Ra");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(0);
	line_style.push_back(1);
	true_index.push_back(0);
	tmp_line_pos.push_back(278);
	iso_energies.push_back(4610.0);
	iso_energies_err.push_back(1.0);


// 226Ra dominant alpha line 94%
	isotope_labels.push_back("{}^{226}Ra");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(0);
	line_style.push_back(1);
	true_index.push_back(1);
	tmp_line_pos.push_back(290);
        iso_energies.push_back(4784.34);
        iso_energies_err.push_back(0.25);

	isotope_labels.push_back("{}^{210}Po");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(14);
	line_style.push_back(1);
	true_index.push_back(2);
	tmp_line_pos.push_back(315);
        iso_energies.push_back(5304.38);
        iso_energies_err.push_back(0.07);

	isotope_labels.push_back("{}^{222}Rn");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(0);
	line_style.push_back(1);
	true_index.push_back(3);
	tmp_line_pos.push_back(325);
        iso_energies.push_back(5489.52);
        iso_energies_err.push_back(0.3);

//223 Ra 

//true_index.push_back(4);
//true_index.push_back(5);

// 225Ac dominant line 50%
	isotope_labels.push_back("{}^{225}Ac");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(4);
	line_style.push_back(1);
	true_index.push_back(6);
	tmp_line_pos.push_back(345);
        iso_energies.push_back(5830.0);
        iso_energies_err.push_back(2.0);

	isotope_labels.push_back("{}^{218}Po");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(0);
	line_style.push_back(1);
	true_index.push_back(7);
	tmp_line_pos.push_back(360);
        iso_energies.push_back(6002.35);
        iso_energies_err.push_back(0.09);

	isotope_labels.push_back("{}^{221}Fr");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(4);
	tmp_line_pos.push_back(378);
	line_style.push_back(1);
	true_index.push_back(9);
        iso_energies.push_back(6341.0);
        iso_energies_err.push_back(1.3);

	isotope_labels.push_back("{}^{217}At");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(4);
	line_style.push_back(1);
	true_index.push_back(12);
	tmp_line_pos.push_back(420);
        iso_energies.push_back(7066.9);
        iso_energies_err.push_back(1.5);

	isotope_labels.push_back("{}^{214}Po");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(0);
	line_style.push_back(1);
	true_index.push_back(14);
	tmp_line_pos.push_back(455);
        iso_energies.push_back(7686.82);
        iso_energies_err.push_back(0.07);

	isotope_labels.push_back("{}^{213}Po");
	fit_isotope.push_back(true);
	isotope_color_tab.push_back(4);
	line_style.push_back(1);
	true_index.push_back(15);
	tmp_line_pos.push_back(490);
        iso_energies.push_back(8375.9);
        iso_energies_err.push_back(2.5);

	N_ISOTOPES_FIT = fit_isotope.size();

	// Define isotopes for drawing only

	isotope_labels.push_back("{}^{212}Po");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(19);
	line_style.push_back(2);
	true_index.push_back(16);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(8784.37);

	// Subdominant peaks from 223Ra chain
//223Ra
//5606.73 30 	25.7 5 
//5716.23 29 	52.6 13 	
	isotope_labels.push_back("{}^{223}Ra");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(2);
	true_index.push_back(4);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(5606.73);

	isotope_labels.push_back("{}^{223}Ra");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(2);
	true_index.push_back(5);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(5716.23);

//219 Rn
//6552.6 10 	12.9 6 
//6819.1 3 	79.4 10 
/*
	isotope_labels.push_back("{}^{219}Rn");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(2);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(6552.6);
*/
	isotope_labels.push_back("{}^{219}Rn");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(9);
	true_index.push_back(11);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(6819.1);

//215 Po
//7386.2 8 	100 
	isotope_labels.push_back("{}^{215}Po");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(9);
	true_index.push_back(13);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(7386.2);
//211 Bi
//6278.2 7 	16.23 14 
//6622.9 6 	83.77 14 
/*
	isotope_labels.push_back("{}^{211}Bi");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(3);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(6278.2);
*/
	isotope_labels.push_back("{}^{211}Bi");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(9);
	true_index.push_back(10);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(6622.9);

// Subdominant 225Ac line 18%
/*
	isotope_labels.push_back("{}^{225}Ac");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(3);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(5792.5);
*/

// Subdominant 221Fr line 15%

	isotope_labels.push_back("{}^{221}Fr");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(4);
	line_style.push_back(9);
	true_index.push_back(8);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(6126.3);


// Sub subdominant 221Fr line 1.3%
/*
	isotope_labels.push_back("{}^{221}Fr");
	fit_isotope.push_back(false);
	isotope_color_tab.push_back(2);
	line_style.push_back(3);
	tmp_line_pos.push_back(520);
	iso_energies_indicator.push_back(6243.0);
*/

	// Now fill the isotope_color vector with the corresponding ROOT color indices
	for(unsigned int i = 0; i < isotope_color_tab.size(); ++i)
		isotope_color.push_back(Tableau20(isotope_color_tab[i]));

	//isotope_color.push_back(Tableau20(6)); Back up color (last of the rainbow)

	filenames.push_back("Rm21122017");
	run_labels.push_back("Sample b");
	min_runtime_cut.push_back(0);


        for(unsigned int i = 0; i < filenames.size(); ++i){

                std::string this_file_name = base_file_path + filenames[i] + "/" + filenames[i] + base_file_ending;

                TFile *this_file = TFile::Open(this_file_name.c_str());
                TTree* this_tree = (TTree*) this_file->Get("t");
                trees.push_back(this_tree);
	}

	//bool refit = true;
	bool refit = false;

	TFile* fit_values_file;

	if(refit){
		fit_values_file = TFile::Open("fit_vals.root", "RECREATE");
	}else{
		fit_values_file = TFile::Open("fit_vals.root", "READ");
	}



	// Now the TTrees are loaded and ready
	
	/*****************
 	 * OVERVIEW PLOT *
	 *****************/

	TCanvas *c_overview = new TCanvas("overview");
	
	for(unsigned int i = 0; i < trees.size(); ++i){
		
		TCut runtime_cut = ("(runtime>"+std::to_string(min_runtime_cut[i])+")").c_str();
		
		trees[i]->Draw(("channel>>h_spec_"+std::to_string(i)+"("
				+std::to_string(n_bins)+","
				+std::to_string(spec_range[0])+","
				+std::to_string(spec_range[1])+")").c_str(), runtime_cut);
		TH1D* this_spec = (TH1D*) c_overview->GetPrimitive(("h_spec_"+std::to_string(i)).c_str());
		spectra.push_back(this_spec);
		double this_runtime = trees[i]->GetMaximum("runtime") - min_runtime_cut[i];
		runtimes.push_back(this_runtime);
		std::cout<<"Runtime ("<<i<<"): "<<this_runtime<<std::endl;
	}

        int canvas_dim[2] = {1200, 900};
        double top_margin = 0.15;
        int axis_title_size = 40;
        int axis_tick_size = 35;

	TCanvas *c_norm = new TCanvas("c2", "",canvas_dim[0], canvas_dim[1]);
	c_norm->SetRightMargin(TINY);
	// gStyle->SetCanvasBorderMode(0);
	// gStyle->SetPadBorderMode(0);
	// gStyle->SetFrameFillColor(0);
	// gStyle->SetFrameBorderMode(0);
	c_norm->Divide(1,2,TINY,TINY);	

	c_norm->cd(1);

	gPad->SetPad(0, aspect_point, 1, 1);
	gPad->SetBottomMargin(TINY);
        gPad->SetTopMargin(top_margin);
	gPad->SetRightMargin(TINY);
	gPad->SetLogy();	

	for(unsigned int i = 0; i < spectra.size(); ++i){
		//spectra[i]->Scale(1000./(runtimes[i]*60.));
		if (i == 0)
			spectra[i]->Draw("HIST");
		else
			spectra[i]->Draw("HISTSAME");


		unsigned int total_n_par = N_ISOTOPES_FIT * N_PAR_CB; 
		TF1* fit_f = new TF1(("fit_"+std::to_string(i)).c_str(), CompleteSpectrum, fit_ranges[0], fit_ranges[1], total_n_par);
		// Initialize fit paramters of the function
		for(unsigned int j = 0; j < N_ISOTOPES_FIT; ++j){

			// Areas are par[0]
			int area_idx = j*N_PAR_CB;
                        fit_f->SetParameter(area_idx, 2E3);
			fit_f->SetParLimits(area_idx,0, 1E7);
			// Alphas are par[1]
			int alpha_idx = j*N_PAR_CB + 1;
                        fit_f->SetParameter(alpha_idx, 2);
			fit_f->SetParLimits(alpha_idx, 0.5, 2.5);
			// Means are par[2]:
			int mean_idx = j*N_PAR_CB + 2;
			fit_f->SetParameter(mean_idx, tmp_line_pos[j]);
			fit_f->SetParLimits(mean_idx, tmp_line_pos[j]-4, tmp_line_pos[j]+4);
			// Sigmas are par[3]:
			int sigma_idx = j*N_PAR_CB + 3;
                        fit_f->SetParameter(sigma_idx, 4);
			fit_f->SetParLimits(sigma_idx, 0.5, 5);
			// Ns are par[4]:
			int ns_idx = j*N_PAR_CB + 4;
                        fit_f->SetParameter(ns_idx, 1);
			fit_f->SetParLimits(ns_idx, 0.2, 6);
			
		}

		fit_f->SetNpx(1000);

		if(refit){
			TFitResultPtr f_ptr;
			TVirtualFitter::Fitter(spectra[i])->SetMaxIterations(700000);
			//f_ptr = spectra[i]->Fit(fit_f, "LLS", "",fit_ranges[0], fit_ranges[1]);
			f_ptr = spectra[i]->Fit(fit_f, "S", "0",fit_ranges[0], fit_ranges[1]);
			TMatrixDSym cov = f_ptr->GetCovarianceMatrix();
			cov_mat.ResizeTo(cov);
			cov_mat = cov;
			// seems to be easier to save the pointer only

			for(int j = 0; j < fit_f->GetNpar(); j++){
				isotope_fit_vals->push_back(fit_f->GetParameter(j));
			}
	
			fit_values_file->WriteObject(isotope_fit_vals, "fit_vals");
			fit_values_file->WriteTObject(&cov_mat, "cov_mat");

		}else{
			// Retrieve the value from the File
			fit_values_file->GetObject("fit_vals", isotope_fit_vals);
			fit_f->SetParameters(&(isotope_fit_vals->at(0)));
			// Set TF1 with the stored values
			
			TMatrixDSym* cov_m = (TMatrixDSym*) fit_values_file->Get("cov_mat");
			cov_mat.ResizeTo(*cov_m);
			cov_mat = *cov_m;

			for(int j = 0; j < fit_f->GetNpar(); j++){
				double this_err = sqrt( cov_mat(j,j));
				fit_f->SetParError(j, this_err); 
                        }
		}
		

		cov_mat.Print();
		//cov_matrices.push_back(cov_mat);

		fit_funcs.push_back(fit_f);

		double this_scale = 1000./(runtimes[i]*60.);
		TF1* fit_f_scale = (TF1*) fit_f->Clone(("fit_scaled"+std::to_string(i)).c_str());
		fit_funcs_cp.push_back(fit_f_scale);
		for(unsigned int j = 0; j < N_ISOTOPES_FIT; ++j){
			double old_area = fit_f->GetParameter(j*N_PAR_CB);
			double old_area_err = fit_f->GetParError(j*N_PAR_CB);
			fit_f->SetParameter(j*N_PAR_CB, old_area*this_scale);
			fit_f->SetParError(j*N_PAR_CB, old_area_err*this_scale);
			std::cout<<"Isotope: "<<isotope_labels[j]<<" Constant = "<<
			fit_f->GetParameter(j*N_PAR_CB)<<" +/- "<<fit_f->GetParError(j*N_PAR_CB)<<std::endl;
		}
	
		spectra[i]->Scale(this_scale);

		fit_f->Draw("SAME");


		// Evaluate integral of fit function inside ROI for 215Po analysis
		std::cout<<"Integral(431 ... 441): "<< fit_f->Integral(431, 441) <<std::endl;
		std::cout<<"Integral(445 ... 457): "<< fit_f->Integral(445, 457) <<std::endl;
		std::cout<<"Integral(386 ... 407): "<< fit_f->Integral(386, 407) <<std::endl;
		std::cout<<"Integral(412 ... 422): "<< fit_f->Integral(412, 422) <<std::endl;


		TCanvas* c_calib = new TCanvas("calib", "", 1200, 900);

		c_calib->SetLeftMargin(0.13);
		c_calib->SetRightMargin(0.01);

		std::vector<double> fit_means;
		std::vector<double> fit_means_err;
		for(unsigned int j = 0; j < N_ISOTOPES_FIT; ++j){
			double this_mean = isotope_fit_vals->at(j*N_PAR_CB+2);
			fit_means.push_back(this_mean);
			fit_means_err.push_back(sqrt(cov_mat(j,j)));
		}
		
		TGraphErrors* calib = new TGraphErrors(N_ISOTOPES_FIT, &fit_means[0], &iso_energies[0], 0, &iso_energies_err[0]);
	
		calib->Draw("AP");
		TF1* ch_2_e = new TF1("ch_2_e", "[0] + [1]*x", 0, 800);
		calib->Fit(ch_2_e, "IM");
		TF1* e_2_ch = new TF1("e_2_ch", "(x - [0])/[1]", 0, 800);
		e_2_ch->SetParameters(ch_2_e->GetParameters());
		
/*	
	// Print the fitted energies and their uncertainty
	
	for(unsigned int j = 0; j < N_ISOTOPES_FIT; ++j){
		std::cout<<isotope_labels[j]<<" Energy: ";
		
		double mean_energy = ch_2_e->Eval(fit_means[j]);
		double energy_stat_err = mean_energy * fit_means_err[j] / fit_means[j];
		double energy_syst_err = ch_2_e->GetParError(0)*ch_2_e->GetParError(0);
		energy_syst_err += ch_2_e->GetParError(1) * fit_means[j] * ch_2_e->GetParError(1) * fit_means[j];
                energy_syst_err = sqrt(energy_syst_err);
		std::cout<<mean_energy<<" +/- ";
		std::cout<<energy_stat_err<<"(stat) +/- ";
		std::cout<<energy_syst_err<<"(syst)"<<std::endl;
	}
*/

	calib->SetMarkerStyle(20);
	calib->SetMarkerSize(2.5);
	ch_2_e->SetLineWidth(3);
	ch_2_e->SetLineColor(Tableau20(6));

	double axis_title_size = 40;	
	calib->GetXaxis()->SetTitle("Energy channel (A.U.)");
	calib->GetXaxis()->CenterTitle();
	calib->GetXaxis()->SetTitleFont(43);
	calib->GetXaxis()->SetLabelFont(43);
	calib->GetXaxis()->SetTitleSize(axis_title_size);
	calib->GetXaxis()->SetLabelSize(axis_title_size);

	calib->GetYaxis()->SetTitle("Alpha energy (keV)");
	calib->GetYaxis()->CenterTitle();
	calib->GetYaxis()->SetTitleFont(43);
	calib->GetYaxis()->SetLabelFont(43);
	calib->GetYaxis()->SetTitleSize(axis_title_size);
	calib->GetYaxis()->SetLabelSize(axis_title_size);
	calib->GetYaxis()->SetTitleOffset(1.4);

	calib->GetXaxis()->SetRangeUser(250,520);
	calib->GetYaxis()->SetRangeUser(4300,8700);

	calib->SetTitle("");

	c_calib->SaveAs("isolde_short_lived_energy_calibration.pdf");
	c_calib->SaveAs("isolde_short_lived_energy_calibration.png");






		c_norm->cd(1);

/*
		TMatrixDSym cov_mat = f_ptr->GetCovarianceMatrix();
		TMatrixDSym old_cov = GetPreviousCovMatrix();
		old_cov(0,0) = 0;old_cov(0,1) = 0;old_cov(0,2) = 0;old_cov(0,3) = 0;old_cov(0,4) = 0;
		old_cov(1,0) = 0;old_cov(1,1) = 0;old_cov(1,2) = 0;old_cov(1,3) = 0;old_cov(1,4) = 0;
		old_cov(2,0) = 0;old_cov(2,1) = 0;
		old_cov(3,0) = 0;old_cov(3,1) = 0;
		old_cov(4,0) = 0;old_cov(4,1) = 0;

		old_cov.Print();		
		//cov_mat.ResizeTo(2,2);
		cov_mat += old_cov;		
	
		cov_mat.Print();
		cov_matrices.push_back(cov_mat);

		fit_f_cp->SetRange(0, 500);
		fit_f_cp->SetLineStyle(2);
		fit_funcs_cp.push_back(fit_f_cp);
		fit_f_cp->Draw("SAME");

		fit_funcs.push_back(fit_f);
		fit_f->Draw("SAME");
*/



	// Draw isotope lines
	TLine *tmp_line = new TLine();
	TLatex *tmp_text = new TLatex();
	gPad->Update();

	double pad_x_width_ndc = 1 - gPad->GetLeftMargin() - gPad->GetRightMargin();

	for(unsigned int j = 0; j < isotope_labels.size(); ++j){
		int color = isotope_color[j];
		std::string text = isotope_labels[j];
		double line_pos;
		if(fit_isotope[j]){
			line_pos = isotope_fit_vals->at(j*N_PAR_CB+2);
		}else{
			line_pos = e_2_ch->Eval(iso_energies_indicator[j-N_ISOTOPES_FIT]);
			std::cout<<line_pos<<std::endl;
		}
		double line_x_ndc = pad_x_width_ndc*((line_pos-gPad->GetUxmin())/(gPad->GetUxmax()-gPad->GetUxmin()))+0.1;
		TLine* this_line = tmp_line->DrawLineNDC(line_x_ndc, 0.0, line_x_ndc, 1-top_margin);
		this_line->SetLineColor(color);
		this_line->SetLineWidth(3);
		this_line->SetLineStyle(line_style[j]);

	}
	// Loop again to draw the label texts (needs to be done in separate loop
	// to prevent some lines being drawn over label boxes...) 
	for(unsigned int j = 0; j < isotope_labels.size(); ++j){
		int color = isotope_color[j];
		std::string text = isotope_labels[j];
		double line_pos;
		if(fit_isotope[j]){
			line_pos = isotope_fit_vals->at(j*N_PAR_CB+2);
		}else{
			line_pos = e_2_ch->Eval(iso_energies_indicator[j-N_ISOTOPES_FIT]);
			std::cout<<line_pos<<std::endl;
		}
		double line_x_ndc = pad_x_width_ndc*((line_pos-gPad->GetUxmin())/(gPad->GetUxmax()-gPad->GetUxmin()))+0.1;

			
		// For the second 226Ra label, skip since it got its label in the previous ieration
		if(true_index[j] == 1){
			continue;
		}
		// The same for 223Ra label
		if(true_index[j] == 5){
			continue;
		}
		int n_prev_skips = 0;
		if(true_index[j]>1)
			++n_prev_skips;
		if(true_index[j]>5)
			++n_prev_skips;
		// The y position of the text needs to alternate to avoid overlapping labels
		double label_y_ndc = 0.1 +(0.1*((true_index[j]-n_prev_skips)%3));

		// Add the text to the line
		// double text_width = 0.05;
		// double text_height = 0.08;
		double text_width = 0.085;
		double text_height = 0.09;

		// For the second 226Ra label, move x position
		double dom_line_pos = line_pos;
		double subdom_line_pos = line_pos;
		if(true_index[j] == 0){
			dom_line_pos = isotope_fit_vals->at(1*N_PAR_CB+2);
			subdom_line_pos = isotope_fit_vals->at(0*N_PAR_CB+2);
		}
		if(true_index[j] == 4){
                        dom_line_pos = e_2_ch->Eval(iso_energies_indicator[j-N_ISOTOPES_FIT+1]);
                        subdom_line_pos = e_2_ch->Eval(iso_energies_indicator[j-N_ISOTOPES_FIT]);
                }
		double sl_ndc = pad_x_width_ndc*((subdom_line_pos-gPad->GetUxmin())/(gPad->GetUxmax()-gPad->GetUxmin()))+0.1;
		double dl_ndc = pad_x_width_ndc*((dom_line_pos-gPad->GetUxmin())/(gPad->GetUxmax()-gPad->GetUxmin()))+0.1;
			
		line_x_ndc = 0.5*(sl_ndc+dl_ndc);
		
		TPaveText *pt = new TPaveText(	line_x_ndc-0.5*text_width, label_y_ndc-0.5*text_height,
						line_x_ndc+0.5*text_width, label_y_ndc+0.5*text_height, 
						"NB NDC");
		TText* this_tex = pt->AddText(text.c_str());
		pt->Draw();

		this_tex->SetTextColor(color);
	}

	
	spectra[0]->SetLineColor(hist_colors[0]);
	spectra[0]->SetLineWidth(3);
	//spectra[1]->SetLineColor(hist_colors[1]);
	//spectra[1]->SetLineWidth(2);
	// fit_funcs[0]->SetLineColor(fit_colors[0]);
	fit_funcs[0]->SetLineColor(Tableau20(6));
	fit_funcs[0]->SetLineWidth(3);
	//fit_funcs[1]->SetLineColor(fit_colors[1]);


        spectra[0]->GetYaxis()->SetTitle("Detected rate (mHz)");
        spectra[0]->GetYaxis()->CenterTitle();
        spectra[0]->GetYaxis()->SetTitleFont(43);
        spectra[0]->GetYaxis()->SetTitleOffset(1.0);
        spectra[0]->GetYaxis()->SetLabelFont(43);
        spectra[0]->GetYaxis()->SetTitleSize(axis_title_size);
        spectra[0]->GetYaxis()->SetLabelSize(axis_tick_size);
        spectra[0]->GetYaxis()->SetRangeUser(5E-4, 5);
        spectra[0]->GetXaxis()->SetTitle("Energy channel (A.U.)");
        spectra[0]->GetXaxis()->CenterTitle();
        spectra[0]->GetXaxis()->SetTitleFont(43);
        spectra[0]->GetXaxis()->SetLabelFont(43);
        spectra[0]->GetXaxis()->SetTitleSize(axis_title_size);
        spectra[0]->GetXaxis()->SetLabelSize(axis_tick_size);
	spectra[0]->GetXaxis()->SetTitleOffset(3);

        spectra[0]->SetTitle("");
        spectra[0]->SetLineColor(1);
	
	gStyle->SetOptStat(0);

	gPad->Update();

	
	TLegend *leg = new TLegend(0.12, 0.65, 0.25, 0.85);
	leg->AddEntry(spectra[0], run_labels[0].c_str());
        //leg->AddEntry(spectra[1], run_labels[1].c_str());
	
	//leg->Draw();
	double y_max = 5;
	double x_new_lim[2] = {ch_2_e->Eval(spec_range[0]), ch_2_e->Eval(spec_range[1])};
	TGaxis* second_x = new TGaxis(	spec_range[0], y_max, 
						spec_range[1], y_max,
						x_new_lim[0], x_new_lim[1], 1010, "-");


	second_x->SetName("sec_x");
	second_x->SetTitle("Energy (keV)");
        second_x->CenterTitle();
        second_x->SetTitleFont(43);
        second_x->SetLabelFont(43);
        second_x->SetTitleSize(axis_title_size);
        second_x->SetLabelSize(axis_tick_size);
        second_x->SetTitleOffset(1.5);
	second_x->Draw();

	
	c_norm->cd(2);
	gPad->SetTopMargin(TINY);
	gPad->SetRightMargin(TINY);
	gPad->SetBottomMargin(0.3);

	//gPad->SetBottomMargin(2);
	// Change aspect ratio of the lower pad
	gPad->SetPad(0, 0.0, 1, aspect_point);
	gPad->SetBorderMode(0);

	//TCanvas *residual = new TCanvas("residual");

	TH1D* h_residual = (TH1D*) spectra[0]->Clone("res");
	for(int j = 0; j < h_residual->GetNbinsX(); ++j){
		double bin_center = h_residual->GetBinCenter(j);
		double fit_value = fit_funcs[0]->Eval(bin_center);
		double scale = 1000./(runtimes[i]*60.);
		double fit_counts = fit_value / scale;
		double this_counts = h_residual->GetBinContent(j) / scale;
		double this_counts_err = sqrt(this_counts);
	
		if(this_counts > 0)	
			h_residual->SetBinContent(j, (this_counts - fit_counts)/sqrt(this_counts));
		else
			h_residual->SetBinContent(j, 0);
		h_residual->SetBinError(j, 0.);
		

	}
	//h_residual->Add(fit_funcs[0], -1);
	h_residual->Draw("AXIS");

	// h_residual->GetYaxis()->SetRangeUser(-3.75, 9);
	h_residual->GetYaxis()->SetRangeUser(-4, 9);
	h_residual->GetYaxis()->SetTitle("Residual");
        h_residual->GetYaxis()->CenterTitle();
        h_residual->GetYaxis()->SetTitleFont(43);
        h_residual->GetYaxis()->SetLabelFont(43);
        h_residual->GetYaxis()->SetTitleSize(axis_title_size);
        h_residual->GetYaxis()->SetLabelSize(axis_tick_size);
	// h_residual->GetYaxis()->SetTitleOffset(0.3);
	h_residual->GetXaxis()->SetTitleOffset(3);
        int n_labels_res_ax = 19;
        h_residual->GetYaxis()->SetNdivisions(n_labels_res_ax,0,0);
        //h_residual->GetYaxis()->SetNdivisions(-n_labels_res_ax);

        // std::vector<int> labels_to_keep {4, 10, 16, 22};
        std::vector<int> labels_to_keep {2, 5, 8, 11};

        // Modify the labels by hand
        for(int jdx = 1; jdx < n_labels_res_ax+1; ++jdx){
            if(std::find(labels_to_keep.begin(), labels_to_keep.end(), jdx) != labels_to_keep.end()){
                // Label should be kept
            }else{
                // Label should not be kept
                h_residual->GetYaxis()->ChangeLabel(jdx, -1, 0);
            }
        }

        // Now modify the labels on the axis
        TAxis *res_ax = h_residual->GetYaxis();
        for(int jdx = res_ax->GetFirst(); jdx < res_ax->GetLast(); ++jdx){
            std::cout<<res_ax->GetBinLabel(jdx)<<std::endl;
        }
/*
	h_residual->GetYaxis()->SetTitleSize(0.11);
	h_residual->GetYaxis()->SetTitleOffset(0.3);
	h_residual->GetYaxis()->SetLabelSize(0.11);
	h_residual->GetXaxis()->SetTitleSize(0.12);
	h_residual->GetXaxis()->SetTitleOffset(1.);
	h_residual->GetXaxis()->SetLabelSize(0.13);
*/
	//residual->Update();
	c_norm->Update();

	// Draw boxes for the blinding regions
	for(unsigned int j = 0; j < blinding_ranges.size(); ++j){
		TBox* this_range = new TBox(	blinding_ranges[j].first,
						gPad->GetUymin(),
						blinding_ranges[j].second,
						gPad->GetUymax());
		this_range->SetFillColorAlpha(Tableau20(15), 1.);
		this_range->Draw("SAME");
	}

        TBox *conf_band_2 = new TBox(spec_range[0], -2, spec_range[1], 2);
        conf_band_2->SetFillColorAlpha(Tableau20(7), 0.5);
        conf_band_2->Draw("SAME");

	TBox *conf_band_1 = new TBox(spec_range[0], -1, spec_range[1], 1);
	conf_band_1->SetFillColorAlpha(Tableau20(6), 0.5);
	conf_band_1->Draw("SAME");


	h_residual->Draw("SAME P");
	h_residual->SetMarkerStyle(8);
	h_residual->SetMarkerSize(0.9);

        gPad->Modified();
        gPad->Update();


	gPad->Update();

	//residual->SaveAs("residual.pdf");
	//residual->SaveAs("residual.png");
	
	c_norm->Update();


	c_norm->SaveAs("isolde_4_spectra_pot.pdf");
	c_norm->SaveAs("isolde_4_spectra_pot.png");

	}

/*
383, 440 AT217
465, 550 Po213
*/
	/************************
 	 * TIME DEPENDENCE PLOT *
	 ************************/

	TCanvas *c_tdep = new TCanvas("t_dep", "", 1200, 900);
	
		
	TCut runtime_cut = ("(runtime>"+std::to_string(min_runtime_cut[0])+")").c_str();
	TCut channel_selection ="(channel>470 && channel<510)";

	double this_runtime = trees[0]->GetMaximum("runtime") - min_runtime_cut[0];
	int time_dep_bins = 75;
	trees[0]->Draw(("runtime/(60*24)>>h_tdep_"+std::to_string(0)+"("
			+std::to_string(time_dep_bins)+","
			+std::to_string(0)+","
			+std::to_string(trees[0]->GetMaximum("runtime")/(60.*24.))+")").c_str(), runtime_cut+channel_selection);
	TH1D* this_spec = (TH1D*) c_tdep->GetPrimitive(("h_tdep_"+std::to_string(0)).c_str());
	std::cout<<"Runtime ("<<0<<"): "<<this_runtime<<std::endl;


	// Calculate the bin width in seconds
	double bin_width = this_runtime * 60. / time_dep_bins;
	double this_scale = 1000./bin_width;
	this_spec->Scale(this_scale);
	this_spec->Sumw2();


        this_spec->GetYaxis()->SetTitle("Detected rate (mHz)");
        this_spec->GetYaxis()->CenterTitle();
        this_spec->GetYaxis()->SetTitleFont(43);
        // this_spec->GetYaxis()->SetTitleOffset(1.0);
        this_spec->GetYaxis()->SetLabelFont(43);
        this_spec->GetYaxis()->SetTitleSize(axis_title_size);
        this_spec->GetYaxis()->SetLabelSize(axis_tick_size);
        this_spec->GetYaxis()->SetRangeUser(-0.5, 17);
        this_spec->GetXaxis()->SetTitle("Measurement time (days)");
        this_spec->GetXaxis()->CenterTitle();
        this_spec->GetXaxis()->SetTitleFont(43);
        this_spec->GetXaxis()->SetLabelFont(43);
        this_spec->GetXaxis()->SetTitleSize(axis_title_size);
        this_spec->GetXaxis()->SetLabelSize(axis_tick_size);
	// this_spec->GetXaxis()->SetTitleOffset(3);

        this_spec->SetTitle("");
        this_spec->SetLineColor(1);
	this_spec->SetLineWidth(3);

	TF1 *exp_fit = new TF1("exp_fit", "expo", 0.5, 11);
	this_spec->Fit(exp_fit, "IMER");
        
	TF1 *exp_fit_extrapol = new TF1("exp_fit_e", "expo", -1, 15);
        exp_fit_extrapol->SetParameters(exp_fit->GetParameters());

	exp_fit->SetLineWidth(3);
	exp_fit_extrapol->SetLineWidth(3);
	exp_fit_extrapol->SetLineStyle(2);

	exp_fit_extrapol->Draw("SAME");
 	

	gStyle->SetOptStat(0);

	c_tdep->SaveAs("isolde_4_time_dep.pdf");
	c_tdep->SaveAs("isolde_4_time_dep.png");


	return 0;
	
/*
	// Here the calculation of the limit starts
	for(unsigned int i = 0; i < trees.size(); ++i){
	
		double this_scale = 1000./(runtimes[i]*60.);

		// Generate n_rand_tf1 random TF1 crystal ball functions based on the 
		// covariance matrix of the fitted function
		// Then use this functions to estimate the uncertainty in the 
		// integration region

		int n_rand_tf1 = 100;
		// Number of points in x that should be stored for each function
		// These points are only for later display of confidence belt
		int n_pts_rand_tf1 = 1000;
		std::vector<std::vector <double>> mc_scan_points;
		std::vector<double> mc_scan_x_values;
		// we also need to create a vector to store the evaluated
		// Background integrals
		std::vector<double> mc_integrals;
		

		TMatrixDSym this_cov_m = cov_matrices[i];
		// Do cholesky decomposition in order to transform the 1 sigma parameter
		TDecompChol* this_decomposition = new TDecompChol(this_cov_m);
		TMatrixD u_mat;
		if(this_decomposition->Decompose()){
			TMatrixD u = this_decomposition->GetU();
			u_mat.ResizeTo(u);
			u_mat = u;
		}else{
			//Decomposition failed
			std::cout<<"Colesky decomposition failed"<<std::endl;
			return 42;
		}

		// Use U matrix to transform the mean to +/- 1 sigma of the fitted value
		// the new fit vals are then mu_vec +/- U * unit_vec

		TVectorT<double> fit_means(fit_funcs[i]->GetNpar(), fit_funcs[i]->GetParameters());

		for(int j = 0; j < n_rand_tf1; ++j){
			TF1* fit_func_mc = GetRandomFunction(u_mat,  fit_means);

			//fit_func_mc->SetNpx(1000);
		
			for(int k = 0; k < n_pts_rand_tf1; ++k){
				
				double this_x_value =  (spec_range[1] - spec_range[0]) * (1.*k/n_pts_rand_tf1) +spec_range[0];
				// since this allways refers to the left edge of the bin, add half the width
				this_x_value += (spec_range[1] - spec_range[0])/(n_pts_rand_tf1*2.0);
				// If this is the first iteration, fill this x value also to the array
				// Also initialize the vector holding all the points
				if(j == 0){
					mc_scan_x_values.push_back(this_x_value);
					std::vector<double> empty_vec;
					mc_scan_points.push_back(empty_vec);
				}
				double this_point = fit_func_mc->Eval(this_x_value);
				mc_scan_points[k].push_back(this_point);
			}
	
			//mc_integrals.push_back(fit_func_mc->Integral(integration_range[0], integration_range[1]));
			// Delete this random function, since it was created with new in the Method
			fit_func_mc->Draw("SAME");
			return 0;
			delete fit_func_mc;
		}

		// Sort the area vector to save time later
		std::sort(mc_integrals.begin(), mc_integrals.end());

		// Initialize the histogram that will hold the confidence belt of the TF1
		TGraphAsymmErrors* confidence_band_1 = new TGraphAsymmErrors(n_pts_rand_tf1);
		TGraphAsymmErrors* confidence_band_2 = new TGraphAsymmErrors(n_pts_rand_tf1);

		for(int j = 0; j < n_pts_rand_tf1; ++j){
			// the graph value will be the mean (should be similar to the Eval of the original TF1)
			// graph error will be the +/- 1 sigma quantiles (34% -> 66% quantile range)
			double this_x = mc_scan_x_values[j];
			std::vector<double> this_vector = mc_scan_points[j];
			double mean = TMath::Mean(this_vector.begin(), this_vector.end());
			double quantiles[4] = {0.15865, 0.84135, 0.02275, 0.997725 };
			double returns[4];
			TMath::Quantiles(this_vector.size(), 4, &this_vector[0], &returns[0], &quantiles[0], false);
			double std = TMath::RMS(this_vector.begin(), this_vector.end());
			confidence_band_1->SetPoint(j, this_x, mean);
			confidence_band_1->SetPointError(j, 0., 0, mean-returns[0], returns[1]-mean);
			confidence_band_2->SetPoint(j, this_x, mean);
			confidence_band_2->SetPointError(j, 0., 0, mean-returns[2], returns[3]-mean);

		}

		confidence_band_1->Draw("SAME E3");
		confidence_band_1->SetFillColorAlpha(hist_colors[i],0.35);
		// Make the graph markers invisible not not spoil the dashed fit line
		confidence_band_1->SetMarkerColorAlpha(hist_colors[i], 1.);
		//confidence_band_2->Draw("SAME E3");
		//confidence_band_2->SetFillColorAlpha(hist_colors[i],0.15);
		//
		}




		std::cout<<"Measurement "<<i<<":"<<std::endl;
		std::cout<<"Exposure: "<<runtimes[i]<<" (min)"<<std::endl;
		int integral_bin_range[2] = {spectra[i]->FindBin(integration_range[0]), spectra[i]->FindBin(integration_range[1])};
		double counts = spectra[i]->Integral(integral_bin_range[0], integral_bin_range[1], "width");
		double cts_err = sqrt(counts);
		std::cout<<"Counts = "<<counts<<" +/- "<<cts_err<<std::endl;

		// calculate the background prediction with its error
		double background = fit_funcs[i]->Integral(integration_range[0], integration_range[1]);
		double background_quantiles[3];
		
		double quantiles[3] = {0.15865, 0.84135, 0.9};
                TMath::Quantiles(mc_integrals.size(), 3, &mc_integrals[0], &background_quantiles[0], &quantiles[0], true);
		// Now find the significance of the observed counts given H0
		// i.e. How often did an integral >= counts occur in the MC sample?
		std::vector<double>::iterator match;
		double search_for[1] = {counts};
		match = std::search(mc_integrals.begin(), mc_integrals.end(), search_for, search_for+1, grater_as);
		int match_i = match-mc_integrals.begin();
		std::cout<<"MC areas: "<<mc_integrals[match_i-1]<<" < "<<counts<<" < "<<mc_integrals[match_i]<<" At index: "<<match_i<<std::endl;
		std::cout<<"=> likelihood Background only: "<<1 - (double)match_i/n_rand_tf1<<std::endl;
		std::cout<<"Background = "<<background<<" + "<<background_quantiles[1]-background<<" - "<<background-background_quantiles[0]<<std::endl;
		std::cout<<"Limit on background counts (@ alpha = "<<quantiles[2]<<") = "<< background_quantiles[2]<<std::endl;
		// Draw a histogram of the background distribution
		double range_backgr[2] = {	TMath::MinElement(n_rand_tf1, &mc_integrals[0]),
						TMath::MaxElement(n_rand_tf1, &mc_integrals[0])};
		TH1D* backgr_dis = new TH1D(("h_back_"+std::to_string(i)).c_str(), "", 200, range_backgr[0], range_backgr[1]);
		for(unsigned int j = 0; j < mc_integrals.size(); ++j){
			backgr_dis->Fill(mc_integrals[j]);
		}
		background_distributions.push_back(backgr_dis);

		double signal = counts - background;
		// calculate uncertainty on signal by error propagation (WRONG!)
		double signal_err[2] = {sqrt(cts_err*cts_err + (-background_quantiles[0]+background)*(-background_quantiles[0]+background)),
					sqrt(cts_err*cts_err + (+background_quantiles[1]-background)*(+background_quantiles[1]-background))};

		std::cout<<"Signal (cts - backgr) = "<<signal<<" + "<<signal_err[1]<<" - "<<signal_err[0]<<std::endl;

		std::cout<<"------------------------------------------------------------"<<std::endl;
		// Now output all the values after correcting for the exposure "this_scale"
		std::cout<<"Observed activity [mBq] = "<<counts*this_scale<<" +/- "<<cts_err*this_scale<<std::endl;
		std::cout<<"Background [mBq] = "<<background*this_scale<<" + "<<(background_quantiles[1]-background)*this_scale<<" - "<<(background-background_quantiles[0])*this_scale<<std::endl;
		std::cout<<"Limit on background activity (@ alpha = "<<quantiles[2]<<") = "<< background_quantiles[2]*this_scale<<std::endl;
		std::cout<<"Signal (cts - backgr) [mBq] = "<<signal*this_scale<<" + "<<signal_err[1]*this_scale<<" - "<<signal_err[0]*this_scale<<std::endl;

		std::cout<<"============================================================"<<std::endl<<std::endl;
	}



	// Draw the background distributions
	for(unsigned int i = 0; i < background_distributions.size(); ++i){
		TCanvas *c_back = new TCanvas(("back_"+std::to_string(i)).c_str());
		background_distributions[i]->Scale(1./background_distributions[i]->Integral(), "WIDTH");
		background_distributions[i]->Draw("HIST");
		background_distributions[i]->SetLineWidth(2);
		background_distributions[i]->GetXaxis()->SetTitle("Background [counts]");
		background_distributions[i]->GetYaxis()->SetTitle("Probability");
		c_back->SetLogy();

		c_back->SaveAs(("dlc_backgr_"+std::to_string(i)+".pdf").c_str());
		c_back->SaveAs(("dlc_backgr_"+std::to_string(i)+".png").c_str());

	}
*/
}



