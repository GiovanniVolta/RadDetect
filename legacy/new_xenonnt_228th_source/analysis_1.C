#include "TTree.h"
#include "TStyle.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TFileMerger.h"
#include "TAxis.h"
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
#include "TGaxis.h"
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
#include <sstream>
#include <ctime>

double RN_DECAY_CONST = log(2)/3.824;
TRandom3 *gen = new TRandom3(0);
unsigned int N_ISOTOPES_FIT;
unsigned int N_PAR_CB;
double PI = 3.141592654;

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

double ModCrystalBallFunc(double* energy, double* pars){
	// This is the definition of a modified crystalball function used to 
	// describe alpha peak of 212Po with the right hand side modified
	// to have an exponential dependence on detected energy

	// Paramters:	pars[0]:	Area 
	//		pars[1]:	Alpha
	//		pars[2]:	Mean
	//		pars[3]:	Sigma
	//		pars[4]:	N
	//		pars[5]:	Lambda

	double gaus = exp(-(energy[0]-pars[2])*(energy[0]-pars[2])/(2*pars[3]*pars[3]));
	double poly = pow(pars[4]/abs(pars[1]), pars[4])*exp(-pars[1]*pars[1]/2.);
	poly *= pow(pars[4]/abs(pars[1]) - abs(pars[1]) - (energy[0]-pars[2])/pars[3], -pars[4]);
	double exponential = exp(-pars[5]*(energy[0] - pars[2])); 

	double norm = pars[4]/abs(pars[1]) * 1./(pars[4]-1)*exp(-pars[1]*pars[1]/2.);
	norm += sqrt(PI/2.)*(1+erf(abs(pars[1])/sqrt(2)));
	norm = 1./(norm*pars[3]);

	if((energy[0] - pars[2])/pars[3] < -pars[1])
		 return norm*pars[0]*poly;
	if((energy[0] - pars[2])/pars[3] > pars[1])
		 return norm*pars[0]*exponential;

	return norm*pars[0]*gaus;
}

double CompleteSpectrum(double *energy, double *allParams){
// This function builds a joint function in order to fit the complete spectrum of peaks

	const unsigned int nParams = N_PAR_CB;
	double result = 0;

	for(unsigned int i = 0; i < N_ISOTOPES_FIT; ++i){
		double *thisParams = (double*) malloc(nParams*sizeof(double));
		for(unsigned int j = 0; j < nParams; ++j){
			thisParams[j] = allParams[i*nParams + j];
		}
		result += CrystalBallFunc(energy, thisParams);
	}
	return result;
}

double CompleteSpectrum_mod(double *energy, double *allParams){
// This function builds a joint function in order to fit the complete spectrum of peaks
// The last function should however be the modified CB function (see above)

	const unsigned int nParams = N_PAR_CB;
	double result = 0;

	for(unsigned int i = 0; i < N_ISOTOPES_FIT; ++i){
                if(i < N_ISOTOPES_FIT -1 ){
			double *thisParams = (double*) malloc(nParams*sizeof(double));
			for(unsigned int j = 0; j < nParams; ++j){
				thisParams[j] = allParams[i*nParams + j];
			}
			result += CrystalBallFunc(energy, thisParams);
		}else{
			double *thisParams = (double*) malloc((nParams+1)*sizeof(double));
			for(unsigned int j = 0; j < nParams+1; ++j){
				thisParams[j] = allParams[i*nParams + j];
			}
			result += ModCrystalBallFunc(energy, thisParams);
		}
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



std::string corrected_channel(std::vector<double> correction_pars){   

    // Time dependent correction of energy response using linear
    // function

    std::stringstream return_str;                                               
    return_str<<"channel / ("                                                  
        <<correction_pars[0]<<"+"                                               
        <<"(runtime)*"<<correction_pars[1]*1e6<<"*1e-6)";                     
    return return_str.str();                                                    
}                           


int analysis_1(){

        TCanvas* blasd = new TCanvas("test");

	TF1* cryst = new TF1("test_func", CrystalBallFunc, -1000, 1000, 5);

	// define names of the files etc.

	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	std::vector<std::string> filenames;
        std::vector<std::string> run_labels;
	std::vector<std::string> isotope_labels;
	std::vector<bool> fit_isotope;
	N_PAR_CB = 5; // Number of parameters for a CB function
	std::vector<int> isotope_color;
	// 3d Matrix storing [run, isotope, fit_param]
	std::vector<std::vector<std::vector <double>>> isotope_fit_vals;
	std::vector<std::vector<std::vector <double>>> isotope_fit_errs;

	std::vector<double> tmp_line_pos;
	std::vector<double> iso_energies;
	std::vector<double> iso_energies_err;
	std::vector<double> iso_energies_indicator;

	std::vector<TTree*> trees;
	std::vector<TH1D*> spectra;
	std::vector<TF1*> fit_funcs;
	std::vector<TF1*> fit_funcs_cp;
	std::vector<TMatrixDSym> cov_matrices;
	int spec_range[2] = {230, 540};
	int fit_ranges[2] = {231, 526};


	// int hist_colors[2] = {Tableau20(6),Tableau20(6)};
	// int fit_colors[2] = {Tableau20(6+1),Tableau20(6+1)};

	int hist_colors[2] = {1, 2};
	int fit_colors[2] = {Tableau20(6),Tableau20(6+1)};
	int ra226_line_cols = Tableau20(0);
	int ra224_line_cols = Tableau20(2);


	int n_bins = spec_range[1] - spec_range[0];
	std::vector<double> runtimes;
	std::vector<double> min_runtime_cut;


	isotope_labels.push_back("{}^{212}Bi");
	fit_isotope.push_back(true);
	// isotope_color.push_back(Tableau20(0));
	isotope_color.push_back(ra224_line_cols);
	tmp_line_pos.push_back(280);
        iso_energies.push_back(5606.60);
        iso_energies_err.push_back(0.25);

	isotope_labels.push_back("{}^{212}Bi");
	fit_isotope.push_back(true);
	// isotope_color.push_back(Tableau20(18));
	isotope_color.push_back(ra224_line_cols);
	tmp_line_pos.push_back(290);
        iso_energies.push_back(6050.78);
        iso_energies_err.push_back(0.3);

	isotope_labels.push_back("{}^{216}Po");
	fit_isotope.push_back(true);
	// isotope_color.push_back(Tableau20(4));
	isotope_color.push_back(ra226_line_cols);
	tmp_line_pos.push_back(325);
        iso_energies.push_back(6778.3);
        iso_energies_err.push_back(0.09);
	double po216_sel_box[2] = {304, 348};
	double po216_sel_iso[2] = {285, 329};


	// Define isotopes for drawing only

	isotope_labels.push_back("{}^{212}Bi");
	fit_isotope.push_back(true);
	// isotope_color.push_back(Tableau20(6));
	isotope_color.push_back(ra224_line_cols);
	tmp_line_pos.push_back(517);
        iso_energies.push_back(10552.1);
        iso_energies_err.push_back(0.02);

	// iso_energies_indicator.push_back(10552.1);


	isotope_labels.push_back("{}^{212}Po");
	fit_isotope.push_back(true);
	isotope_color.push_back(Tableau20(4));
	// isotope_color.push_back(ra226_line_cols);
	tmp_line_pos.push_back(428);
        iso_energies.push_back(8784.37);
        iso_energies_err.push_back(0.07);
	double po212_sel_box[2] = {348, 496};
	double po212_sel_iso[2] = {329, 476};

	
	N_ISOTOPES_FIT = fit_isotope.size();

	//isotope_color.push_back(Tableau20(6)); Back up color (last of the rainbow)

	filenames.push_back("Rn22042022");
	// Define wether the TH2D of the data should be re-computed or not
	bool compress_data = false;

	long int max_events = 100000000;
//	long int first_event = 6*50000000;
	long int first_event = 1*50000000;

	run_labels.push_back("w/ top ring");
	min_runtime_cut.push_back(0);

	/*
	filenames.push_back("Ap03062019_without");
	run_labels.push_back("w/o top ring");
	min_runtime_cut.push_back(0);
	*/

        for(unsigned int i = 0; i < filenames.size(); ++i){

                std::string this_file_name = base_file_path + filenames[i] + "/" + filenames[i] + base_file_ending;

                TFile *this_file = TFile::Open(this_file_name.c_str());
                TTree* this_tree = (TTree*) this_file->Get("t");
		//this_tree->SetMaxEntryLoop(max_events);
                trees.push_back(this_tree);
	}




	// Now the TTrees are loaded and ready

	// Derive a time correction for the energy response
	// Maybe due to aging, the response drops as function of time,
	// leading to worsening in resolution
	// The decrease is fit with a linear function E = p1*time + p0:
	/*
	****************************************
	Minimizer is Linear
	Chi2                      =      88306.6
	NDf                       =         5461
	p0                        =      327.919   +/-   0.000488358 
	p1                        = -0.000155324   +/-   1.66286e-07 
	*/

	std::vector<double> energy_corr_pars;
	energy_corr_pars.push_back(327.919);
	energy_corr_pars.push_back(-0.000155324);
	

	// Due to high event rate, construct a 2D Histogram of the data and save the TH2D
	if(compress_data){
		Float_t ev_time;
		Int_t ev_channel;
		trees[0]->SetBranchAddress("runtime", &ev_time);
		trees[0]->SetBranchAddress("channel", &ev_channel);

		trees[0]->GetEvent(0);
		double min_time = ev_time;
		trees[0]->GetEvent(trees[0]->GetEntries()-1);
		double max_time = ev_time;
		double n_minutes = max_time - min_time; 

		double corrected_channel;

		int n_channels = spec_range[1]-spec_range[0];
		TH2D* complete_data = new TH2D("all_data", "",
                                               n_minutes, min_time, max_time,
                                               n_channels, spec_range[0], spec_range[1]);
		TH2D* complete_data_corr = new TH2D("all_data_corr", "",
                                               n_minutes, min_time, max_time,
                                               n_channels, spec_range[0], spec_range[1]);

		for(long int idx = 0; idx < trees[0]->GetEntries(); ++idx){
		// for(long int idx = 0; idx < 50000000; ++idx){
			if(idx % 10000000 == 0)
				std::cout<<idx<< " of "<<trees[0]->GetEntries()<<" ("<<(double)idx/trees[0]->GetEntries()*100.<<")%"<<std::endl;
			trees[0]->GetEvent(idx);
			complete_data->Fill(ev_time, ev_channel);
			
			corrected_channel = 	ev_channel * 
						energy_corr_pars[0] / (energy_corr_pars[0] + ev_time * energy_corr_pars[1]);
			complete_data_corr->Fill(ev_time, corrected_channel);	
			
		}

		TCanvas *c_time_stab = new TCanvas("time");
		complete_data->Draw("COLZ");

		// Save this histogram to the .root file for later useage.
		complete_data->SaveAs("compressed.root");

		TCanvas *c_time_stab_corr = new TCanvas("time_corr");
		complete_data_corr->Draw("COLZ");

		// Save the area corrected histogram to the .root file for later useage.
		complete_data_corr->SaveAs("compressed_corr.root");

	}

	
	TFile* data = TFile::Open("compressed.root");
	TH2D* all_data = (TH2D*) data->Get("all_data");

	TCanvas *c_time_stab = new TCanvas("time");
	all_data->Draw("COLZ");

	TProfile* med_profile_216 = all_data->ProfileX("pfx_216", 305, 312);
	med_profile_216->Draw("SAME");


	/*****************
 	/ * OVERVIEW PLOT *
	 *****************/

	TCanvas *c_overview = new TCanvas("overview");
	
	TH1D* this_spec = all_data->ProjectionY("h_spec_chan", min_runtime_cut[0], -1);
	this_spec->Draw();
	spectra.push_back(this_spec);

	// Calculate exposure of the data
	Float_t ev_time;
	trees[0]->SetBranchAddress("runtime", &ev_time);
	first_event = 0;
	max_events = trees[0]->GetEntries()-1;
	trees[0]->GetEvent(first_event);
	double min_time = ev_time;
	trees[0]->GetEvent(first_event + max_events);
	double max_time = ev_time;
	// double this_runtime = trees[0]->GetMaximum("runtime") - min_runtime_cut[0];
	double this_runtime = max_time - min_time; 
	runtimes.push_back(this_runtime);
	std::cout<<"Runtime : "<<this_runtime<<std::endl;


	//TCanvas *c_norm = new TCanvas("c2", "", 1200, 900);
	//TCanvas *c_norm = new TCanvas("c2", "",1200, 900);
	TCanvas *c_norm = new TCanvas("c2", "",1400, 900);

	double left_margin = 0.125;
	double right_margin = 0.01;
	double top_margin = 0.01;
	c_norm->SetLeftMargin(left_margin);
	c_norm->SetRightMargin(right_margin);
	c_norm->SetTopMargin(top_margin);


	c_norm->SetLogy();

	double this_scale = 0;	


	for(unsigned int i = 0; i < spectra.size(); ++i){
		//spectra[i]->Scale(1000./(runtimes[i]*60.));
		if (i == 0)
			this_spec->Draw();


		unsigned int total_n_par = N_ISOTOPES_FIT * N_PAR_CB; 
		// TF1* fit_f = new TF1(("fit_"+std::to_string(i)).c_str(), CompleteSpectrum, fit_ranges[0], fit_ranges[1], total_n_par);
		TF1* fit_f = new TF1(("fit_"+std::to_string(i)).c_str(), CompleteSpectrum_mod, fit_ranges[0], fit_ranges[1], total_n_par+1);

		// Initialize fit paramters of the function
		for(unsigned int j = 0; j < N_ISOTOPES_FIT; ++j){

			// Areas are par[0]
			int area_idx = j*N_PAR_CB;
                        fit_f->SetParameter(area_idx, 1E8);
			//fit_f->SetParLimits(area_idx,0, 2E7);
			// Alphas are par[1]
			int alpha_idx = j*N_PAR_CB + 1;
                        fit_f->SetParameter(alpha_idx, 2);
			fit_f->SetParLimits(alpha_idx, 1, 3);
			// Means are par[2]:
			int mean_idx = j*N_PAR_CB + 2;
			fit_f->SetParameter(mean_idx, tmp_line_pos[j]);
			fit_f->SetParLimits(mean_idx, 0.975*tmp_line_pos[j], 1.025*tmp_line_pos[j]);
			// Sigmas are par[3]:
			int sigma_idx = j*N_PAR_CB + 3;
                        fit_f->SetParameter(sigma_idx, 4);
			fit_f->SetParLimits(sigma_idx, 1, 10);
			// Ns are par[4]:
			int ns_idx = j*N_PAR_CB + 4;
                        fit_f->SetParameter(ns_idx, 1.9);
			fit_f->SetParLimits(ns_idx, 1.6, 2.2);
			
		}
		// Set the lambda parameter for the 212Po peak
	        fit_f->SetParameter(N_ISOTOPES_FIT*N_PAR_CB, 0.25);
		fit_f->SetParLimits(N_ISOTOPES_FIT*N_PAR_CB, 0.05, 0.5);


		fit_f->SetNpx(1000);
		TFitResultPtr f_ptr;
		TVirtualFitter::Fitter(all_data)->SetMaxIterations(50000);
		f_ptr = this_spec->Fit(fit_f, "S0", "",fit_ranges[0], fit_ranges[1]);
		
		std::vector<std::vector<double>> this_run_fit_vals;
		std::vector<std::vector<double>> this_run_fit_vals_err;
		for(unsigned int k = 0; k < N_ISOTOPES_FIT; ++k){
			std::vector<double> this_iso_vals;
			std::vector<double> this_iso_errs;
			TF1* cryst = new TF1("test_func", CrystalBallFunc, -10, 10, 5);
			cryst->SetParameters(1.,1.35700e+00,2.87874e+02,4.61077e+00,4.99516e+00),
			cryst->SetNpx(1000);

			for(unsigned int l = 0; l < N_PAR_CB; ++l){
				double this_param = fit_f->GetParameter(k*N_PAR_CB+l);
				this_iso_vals.push_back(this_param);
				double this_err = fit_f->GetParError(k*N_PAR_CB+l);
				this_iso_errs.push_back(this_err);

			}
			this_run_fit_vals.push_back(this_iso_vals);
			this_run_fit_vals_err.push_back(this_iso_errs);
		}
		isotope_fit_vals.push_back(this_run_fit_vals);
		isotope_fit_errs.push_back(this_run_fit_vals_err);

		TMatrixDSym cov_mat = f_ptr->GetCovarianceMatrix();

		cov_mat.Print();
		cov_matrices.push_back(cov_mat);

		fit_funcs.push_back(fit_f);

		this_scale = 1./(runtimes[i]*60.);

		// We like to express the y-axis in terms of counts per second
		// per keV. So we need to know the width of each bin in terms
		// of keV. This corresponds to the slope of the linear
		// calibration function determined below. Its slope was found
		// to be 20.403
	
		this_scale = this_scale / 20.403;

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

		this_spec->Scale(this_scale);

		// fit_f->Draw("SAME");


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
	}


	TCanvas* c_calib = new TCanvas("calib", "", 1200, 900);

	c_calib->SetLeftMargin(0.13);
	c_calib->SetRightMargin(0.01);
//   calib->SetLeftMargin(0.1252087);
//   calib->SetRightMargin(0.008347246);
//   calib->SetTopMargin(0.01027397);



	std::vector<double> fit_means;
	std::vector<double> fit_means_err;
	TMatrixDSym cov_mat = cov_matrices[0];
	for(unsigned int j = 0; j < N_ISOTOPES_FIT; ++j){
	double this_mean = isotope_fit_vals[0][j][2];
	double this_error = isotope_fit_errs[0][j][2];

		fit_means.push_back(this_mean);
		fit_means_err.push_back(this_error);
	}
	
	TGraphErrors* calib = new TGraphErrors(N_ISOTOPES_FIT, &fit_means[0], &iso_energies[0], &fit_means_err[0], &iso_energies_err[0]);
	
	calib->Draw("AP");
	TF1* ch_2_e = new TF1("ch_2_e", "[0] + [1]*x", 0, 800);
	calib->Fit(ch_2_e, "IM");
	TF1* e_2_ch = new TF1("e_2_ch", "(x - [0])/[1]", 0, 800);
	e_2_ch->SetParameters(ch_2_e->GetParameters());
	
	
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

	std::cout<<"Slope of calibration: "<<ch_2_e->GetParameter(1)<<std::endl;

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

	calib->GetXaxis()->SetRangeUser(160,550);
	calib->GetYaxis()->SetRangeUser(5200,10800);

	calib->SetTitle("");

	c_calib->SaveAs("implanted_spectra_calibration.pdf");
	c_calib->SaveAs("implanted_spectra_calibration.png");


	c_norm->cd();

	// double y_lim[2] = {5e-3, 5e3};
	double y_lim[2] = {0.00015, 7};

	spectra[0]->GetXaxis()->SetRangeUser(spec_range[0], spec_range[1]);

	// Indicate the isotope selection regions

	gStyle->SetHatchesLineWidth(2);
	
	TH1D *po212_sel_spec = (TH1D*) spectra[0]->Clone("po212_sel");
	po212_sel_spec->SetLineWidth(0);
	po212_sel_spec->SetFillColor(Tableau20(4));
	po212_sel_spec->SetFillStyle(3454);
	po212_sel_spec->GetXaxis()->SetRangeUser(po212_sel_box[0], po212_sel_box[1]);
	
	
	TH1D *po216_sel_spec = (TH1D*) spectra[0]->Clone("po216_sel");
	po216_sel_spec->SetLineWidth(0);
	po216_sel_spec->SetFillColor(Tableau20(0));
	po216_sel_spec->SetFillStyle(3445);
	po216_sel_spec->GetXaxis()->SetRangeUser(po216_sel_box[0], po216_sel_box[1]);


/*
	// Draw the typical 212Po selection
	TBox *sel_box = new TBox();
	sel_box->SetFillColor(Tableau20(4));
	sel_box->SetFillStyle(3654);
	sel_box->DrawBox(po212_sel_box[0], y_lim[0], po212_sel_box[1], y_lim[1]);

	sel_box->SetFillColor(Tableau20(0));
	sel_box->SetFillStyle(3654);
	sel_box->DrawBox(po216_sel_box[0], y_lim[0], po216_sel_box[1], y_lim[1]);
*/


	po212_sel_spec->Draw("HIST SAME");
	po216_sel_spec->Draw("HIST SAME");

	// Re-Draw Histogram and Fit-function to have it in front of the box
	spectra[0]->Draw("HISTSAME");
	fit_funcs[0]->Draw("SAME");

	// Draw isotope lines
	TLine *tmp_line = new TLine();
	TLatex *tmp_text = new TLatex();
	c_norm->Update();
	for(unsigned int j = 0; j < isotope_labels.size(); ++j){
		int color = isotope_color[j];
		std::string text = isotope_labels[j];
		double line_pos;
		if(fit_isotope[j]){
			line_pos = isotope_fit_vals[0][j][2];
		}else{
			line_pos = e_2_ch->Eval(iso_energies_indicator[j-N_ISOTOPES_FIT]);
		}
		double line_x_ndc = (1-left_margin-right_margin)*((line_pos-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+left_margin;

		TLine* this_line = tmp_line->DrawLineNDC(line_x_ndc, 0.1, line_x_ndc, 1-top_margin);

		//TLatex* this_tex = tmp_text->DrawLatexNDC(line_x_ndc, 0.8, text.c_str());
		
		this_line->SetLineColor(color);
		this_line->SetLineWidth(3);
		//this_tex->SetTextColor(color);
	}

	// Loop again to draw the label texts (needs to be done in separate loop
	// to prevent some lines being drawn over label boxes...) 
	for(unsigned int j = 0; j < isotope_labels.size(); ++j){
		int color = isotope_color[j];
		std::string text = isotope_labels[j];
		double line_pos;
		if(fit_isotope[j]){
			line_pos = isotope_fit_vals[0][j][2];
		}else{
			line_pos = e_2_ch->Eval(iso_energies_indicator[j-N_ISOTOPES_FIT]);
		}

		double line_x_ndc = (1-left_margin-right_margin)*((line_pos-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+left_margin;

		// The y position of the text needs to alternate to avoid overlapping labels
		double label_y_ndc = 0.2 +(0.075*((j)%4));

		// Skip the label of the subdominant Bi212 line
		if (j == 0)
			continue;
		if (j==1){
			label_y_ndc = 0.2;
			line_x_ndc = line_x_ndc - 0.025;
	  	} 
		if(text == "{}^{212}Po")
			label_y_ndc = 0.2 +(0.075*((2)%4));
                if(text == "{}^{216}Po")
			label_y_ndc = 0.2 +(0.075*((1)%4));
		if(text == "{}^{214}Po")
			label_y_ndc = 0.2 +(0.075*((0)%4));


		// Add the text to the line
		double text_width = 0.1;
		double text_height = 0.08;
		
		TPaveText *pt = new TPaveText(	line_x_ndc-0.5*text_width, label_y_ndc-0.5*text_height,
						line_x_ndc+0.5*text_width, label_y_ndc+0.5*text_height, 
						"NB NDC");
		TText* this_tex = pt->AddText(text.c_str());
		pt->Draw();

		this_tex->SetTextColor(color);
	}
	
	spectra[0]->SetLineColor(hist_colors[0]);
	spectra[0]->SetLineWidth(3);
	// spectra[1]->SetLineColor(hist_colors[1]);
	// spectra[1]->SetLineWidth(2);
	fit_funcs[0]->SetLineColor(fit_colors[0]);
	fit_funcs[0]->SetLineWidth(3);
	fit_funcs[0]->SetLineWidth(0);
	// fit_funcs[1]->SetLineColor(fit_colors[1]);

	spectra[0]->GetXaxis()->SetTitle("Energy channel (A.U.)");
	spectra[0]->GetXaxis()->CenterTitle();
	spectra[0]->GetXaxis()->SetTitleFont(43);
	spectra[0]->GetXaxis()->SetLabelFont(43);
	spectra[0]->GetXaxis()->SetTitleSize(axis_title_size);
	spectra[0]->GetXaxis()->SetLabelSize(axis_title_size);

	// Since we normalized to the (energy) width of the bins
	spectra[0]->GetYaxis()->SetTitle("Detected rate (Hz/keV)");
	spectra[0]->GetYaxis()->SetTitleOffset(1.25);
	spectra[0]->GetYaxis()->CenterTitle();
	spectra[0]->GetYaxis()->SetTitleFont(43);
	spectra[0]->GetYaxis()->SetLabelFont(43);
	spectra[0]->GetYaxis()->SetTitleSize(axis_title_size);
	spectra[0]->GetYaxis()->SetLabelSize(axis_title_size);

	spectra[0]->GetYaxis()->SetRangeUser(y_lim[0], y_lim[1]);
	gStyle->SetOptStat(0);
	spectra[0]->SetTitle("");


	// Create a second x-Axis that reads the energy in MeV
	// Then replace the original (meaningless) x-axis that reads the ADC channels
	
	double x_new_lim[2] = {ch_2_e->Eval(spec_range[0]), ch_2_e->Eval(spec_range[1])};
	TGaxis* second_x = new TGaxis(	spec_range[0], y_lim[0], 
						spec_range[1], y_lim[0],
						x_new_lim[0]/1e3, x_new_lim[1]/1e3, 1010, "+");

	second_x->SetName("sec_x");
	second_x->SetTitle("Alpha energy (MeV)");
        second_x->CenterTitle(); 
        second_x->SetTitleFont(43); 
        second_x->SetLabelFont(43);
	second_x->SetTitleSize(axis_title_size); 
	second_x->SetLabelSize(axis_title_size);


	// Make the original axis invsible
	spectra[0]->GetXaxis()->SetLabelSize(0);
	spectra[0]->GetXaxis()->SetTickLength(0);
	spectra[0]->GetXaxis()->SetTitle("");

	second_x->Draw();

	c_norm->Update();

	std::cout<<"Integral"<<spectra[0]->Integral(250,310)<<std::endl;

	// We need to evaluate which fraction of the area is not taken into account by that selection
	// For this we take the part of the CB Function describing 212Po and integrate the are between
	// po212_sel_box[0] and po212_sel_box[1] and compare this with the integral from 0 to infinity

	std::cout<<"Computing fraction of area within the 212Po selection region"<<std::endl;
	// Initialize and set a TF1 for the 212Po fit using the ModCrystalBallFunc
	TF1 *po212_part_function = new TF1("212PoPeakFunction", ModCrystalBallFunc, 0, 900, 6);
	for(int jdx = 0; jdx < 6; ++jdx){
		double this_value = fit_funcs[0]->GetParameter((N_ISOTOPES_FIT-1)*N_PAR_CB+jdx);
		po212_part_function->SetParameter(jdx, this_value);
	}
	// Now compute the integral in the selection and in the full range
	double int_212po_box = po212_part_function->Integral(po212_sel_box[0], po212_sel_box[1]);
	double int_212po_tot = po212_part_function->Integral(0, 900);
	std::cout<<"Integral in Selection box: "<<int_212po_box<<" (A.U.)"<<std::endl;
	std::cout<<"Integral Total: "<<int_212po_tot<<" (A.U.)"<<std::endl;
	std::cout<<"Fraction of area outside box "<<(int_212po_tot - int_212po_box)/int_212po_tot<<std::endl;

	TCanvas *test_canv = new TCanvas("test_test");
	po212_part_function->Draw();

	c_norm->SaveAs("spectrum_radon_monitor.pdf");
	c_norm->SaveAs("spectrum_radon_monitor.png");



	TCanvas *po216_evolution = new TCanvas("216Po_evolution");

	po216_evolution->SetLogy();	

	TH1D* po216_time_evol = all_data->ProjectionX("h_216po_evol", po216_sel_box[0], po216_sel_box[1]);

	TH1D* po212_time_evol = all_data->ProjectionX("h_212po_evol", po212_sel_box[0], po212_sel_box[1]);


	po212_time_evol->Draw();
	po216_time_evol->Draw("SAME");

	// between 20 and 500 minutes the 216Po is steady and unaffected by 
	// the high rate in the detector. Then average the 216Po rate 
	// over the full run duration and derive a correction factor
	// to this first part of the measurement

	//double range_early[2] = {20,500};
	// double range_early[2] = {20,240};
	//double range_early[2] = {100,240};
	// double range_early[2] = {15,235};
	double range_early[2] = {25,245};

	double range_late[2] = {4800,5200};
	//double range_late[2] = {2500,5200};
	//double range_late[2] = {3600,3900};

	double po216_counts_early = po216_time_evol->Integral(range_early[0], range_early[1]);
	double po216_act_early = po216_counts_early / (range_early[1] - range_early[0]) / 60.;
	std::cout<<	"N Po216 early: "<<po216_counts_early<<
			" => A(216Po) early: "<<po216_act_early<<std::endl;

	double po216_counts_total = po216_time_evol->Integral();
	double po216_act_total = po216_counts_total / runtimes[0] / 60.;
	std::cout<<	"N Po216 total: "<<po216_counts_total<<
			" => A(216Po) total: "<<po216_act_total<<std::endl;
	
	double act_correction = po216_act_early/po216_act_total;
	std::cout<<	"Activity correction factor: "<<act_correction<<std::endl;

	// Now estimate the 212Po activity from time late in the measurement
	// when an equilibrium is achieved (range_late)
	// Then apply the correction factor to this value
	
	double po212_counts_late = po212_time_evol->Integral(range_late[0], range_late[1]);
	double po212_act_late = po212_counts_late / (range_late[1] - range_late[0]) / 60.;

	std::cout<<	"N Po212(late): "<<po212_counts_late<<
			" => A(212Po) late: "<<po212_act_late<<std::endl<<
			"Corrected 212Po activity: "<<po212_act_late * act_correction<<std::endl;


	
	// Evaluate the 212Po activity very early in the measurement in order to extrapolate 
	// this value to the equilibrium value

	double po212_counts_early = po212_time_evol->Integral(range_early[0], range_early[1]);
	double po212_act_early = po212_counts_early / (range_early[1] - range_early[0]) / 60.;

	std::cout<<	"N Po212(early): "<<po212_counts_early<<
			" => A(212Po) early: "<<po212_act_early<<std::endl;
	
	return 0;
}



