#include "TTree.h"
#include "TStyle.h"
#include "TCut.h"
#include "TH1.h"
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
#include "TBox.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
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
TRandom3 *gen = new TRandom3(42);

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


TMatrixDSym GetPreviousCovMatrix(){
/*****************************************
 * Reference Ffit values from 218Po Peak *
 *****************************************/
	TMatrixDSym result = TMatrixDSym(5);
	double l[25] = {	5900, 0.1915, -0.3832, -0.01766, 0.02963,
				0.1915, 0.0001288, -6.356e-05, -0.0001228, 0.00029,
				-0.3832, -6.356e-05, 7.695e-05, 8.415e-05, -0.0001782,
				-0.01766, -0.0001228, 8.415e-05, 0.0003229, -0.0009467,
				0.02963, 0.00029, -0.0001782, -0.0009467, 0.003323};
	for(int i = 0; i < result.GetNrows(); ++i){
		for(int j = 0; j < result.GetNcols(); ++j){
			result[i][j] = l[i*result.GetNrows() + j];
		}
	}
	return result;
}

TF1* GetRandomFunction(std::string func, TMatrixD u_matrix, TVectorT<double> mean_vec){

	// This function will generate a random TF1 func, where the parameters are
	// Set according to a Multivariate Gaussian, based on the Cholesky decomposed 
	// Covariance Matrix Cov = U^T U, and the vector containing the mean parameter position
	
	TF1 *result = new TF1("", func.c_str(), 0, 500);
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

	// define names of the files etc.

	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	std::vector<std::string> filenames;
        std::vector<std::string> run_labels;
	std::vector<TTree*> trees;
	std::vector<TH1D*> spectra;
	std::vector<TH1D*> background_distributions;
	std::vector<TF1*> fit_funcs;
	std::vector<TF1*> fit_funcs_cp;
	std::vector<TMatrixDSym> cov_matrices;
	int spec_range[2] = {140, 225};
	// To not bias the crystal ball fit, blind the range of 224Ra
	// To not run into problems with lying on the edge of the bins,
	// Define this as a double to be clearly inside a bin
	// double integration_range[2] = {173.5,176.5};
	double integration_range[2] = {173.01,177.99};


	//std::string fit_function = "expo(0)+expo(2)";
	std::string fit_function = "crystalball";

	//int fit_ranges[2] = {179,193}; // For Crystal ball function
	int fit_ranges[2] = {179,193}; // For Crystal ball function
	int fit_ranges_ref[2] = {170,193}; //For reference 218Po spectrum

	//int fit_ranges[2] = {179,188}; // For exponential tail

	int hist_colors[2] = {1,2};

	int n_bins = spec_range[1] - spec_range[0];
	std::vector<double> runtimes;
	std::vector<double> min_runtime_cut;

/*	
	filenames.push_back("Rn17022020");
	run_labels.push_back("reference");
	min_runtime_cut.push_back(0);
*/
	filenames.push_back("Rn04112019");
	run_labels.push_back("uncoated");
	min_runtime_cut.push_back(0);
	filenames.push_back("Rn04032020");
	run_labels.push_back("DLC coated");
	min_runtime_cut.push_back(300);


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

	TCanvas *c_norm = new TCanvas("c2", "", 1200, 900);

	double TINY = 0.02;
	gPad->SetTopMargin(TINY);
	gPad->SetRightMargin(TINY);

	for(unsigned int i = 0; i < spectra.size(); ++i){
		// spectra[i]->Scale(1000./(runtimes[i]*60.));
		// spectra[i]->Sumw2();
		if (i == 0)
			spectra[i]->Draw("HIST");
		else
			spectra[i]->Draw("HISTSAME");

		TF1* fit_f = new TF1(("fit_"+std::to_string(i)).c_str(), fit_function.c_str(), fit_ranges[0],fit_ranges[1]);

		fit_f->SetNpx(1000);


/*****************************************
 * Reference Ffit values from 218Po Peak *
 *****************************************/

/*

 FCN=36.9607 FROM MINOS     STATUS=SUCCESSFUL     72 CALLS         919 TOTAL
                     EDM=3.07931e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  Constant     1.42450e+04   7.68090e+01   2.68514e-06  -1.16834e-02
   2  Mean         1.87177e+02   1.13482e-02   3.34465e-05   3.85746e-03
   3  Sigma        1.65644e+00   8.77185e-03  -3.19943e-05  -1.05594e-02
   4  Alpha        1.10160e+00   1.79699e-02  -1.63664e-04   6.47491e-03
   5  N            1.91903e+00   5.76260e-02   5.76260e-02  -6.74317e-03
                               ERR DEF= 0.5

5x5 matrix is as follows

     |      0    |      1    |      2    |      3    |      4    |
----------------------------------------------------------------------
   0 |       5900      0.1915     -0.3832    -0.01766     0.02963 
   1 |     0.1915   0.0001288  -6.356e-05  -0.0001228     0.00029 
   2 |    -0.3832  -6.356e-05   7.695e-05   8.415e-05  -0.0001782 
   3 |   -0.01766  -0.0001228   8.415e-05   0.0003229  -0.0009467 
   4 |    0.02963     0.00029  -0.0001782  -0.0009467    0.003323 

*/

		if (run_labels[i] == "reference"){	
			fit_f->SetParLimits(0, 0, 5E4);
		}else{
			fit_f->SetParLimits(0, 0, 1E3);
		}

		fit_f->SetParameter(1, 188);
		fit_f->SetParLimits(1, 185, 191);
		fit_f->SetParLimits(2, 0, 5);
		fit_f->SetParLimits(3, 0, 5);
		fit_f->SetParLimits(4, 0, 5);

/*		
		fit_f->SetParameter(1, 188);
		fit_f->SetParLimits(1, 185, 191);
		fit_f->SetParLimits(2, 0, 3);
		fit_f->SetParLimits(3, 0, 4);
		fit_f->SetParLimits(4, 0, 4);
*/
		fit_f->SetParameter(2, 1.65644e+00);
		fit_f->SetParameter(3, 1.10160e+00);
		fit_f->SetParameter(4, 1.91903e+00);

	
		fit_f->FixParameter(2, 1.65644e+00);
		fit_f->FixParameter(3, 1.10160e+00);
		fit_f->FixParameter(4, 1.91903e+00);
		
		TFitResultPtr f_ptr;
		if (run_labels[i] == "reference"){
			f_ptr = spectra[i]->Fit(fit_f, "LLIMES", "",fit_ranges_ref[0], fit_ranges_ref[1]);
		}else{
			f_ptr = spectra[i]->Fit(fit_f, "LLIMES", "",fit_ranges[0], fit_ranges[1]);
		}
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


		// Clone the function twice for drawing
		TF1* fit_f_cp = (TF1*) fit_f->Clone(("fit_cp_"+std::to_string(i)).c_str());
		TF1* fit_f_cp_draw = (TF1*) fit_f->Clone(("fit_cp_draw_"+std::to_string(i)).c_str());

		// Scale their height of the function by the runtime
		fit_f_cp->SetParameter(0, fit_f->GetParameter(0) * 1000./(runtimes[i]*60.));
		fit_f_cp_draw->SetParameter(0, fit_f->GetParameter(0) * 1000./(runtimes[i]*60.));

		fit_f_cp->SetRange(0, 500);
		fit_f_cp_draw->SetRange(fit_ranges[0], fit_ranges[1]);

		fit_f_cp->SetLineStyle(2);
		fit_f_cp->SetLineWidth(2);
		fit_f_cp->SetLineColor(hist_colors[i]);

		fit_f_cp_draw->SetLineColor(hist_colors[i]);
		fit_f_cp_draw->SetLineWidth(2);

		fit_f_cp->Draw("SAME");
		fit_f_cp_draw->Draw("SAME");

		fit_funcs.push_back(fit_f);

	}

	// After the fit, normalize the histograms by the runtime
	spectra[0]->Scale(1000./(runtimes[0]*60.));
	spectra[1]->Scale(1000./(runtimes[1]*60.));
	
	spectra[0]->SetLineColor(hist_colors[0]);
	spectra[0]->SetLineWidth(3);
	spectra[1]->SetLineColor(hist_colors[1]);
	spectra[1]->SetLineWidth(3);

	spectra[0]->GetXaxis()->SetTitle("Energy channel [A.U.]");
	spectra[0]->GetYaxis()->SetTitle("Couts");

	c_norm->SetLogy();

	gStyle->SetOptStat(0);
	spectra[0]->SetTitle("");

	c_norm->Update();	

	double line_1_x_ndc = (0.9-TINY)*((integration_range[0]-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+0.1;
	double line_2_x_ndc = (0.9-TINY)*((integration_range[1]-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+0.1;

	double iso_position[3] = {211.5, 188.5, 163.7};
	std::vector<std::string> iso_labels;
	iso_labels.push_back("{}^{216}Po");
	iso_labels.push_back("{}^{212}Bi");
	iso_labels.push_back("{}^{210}Pb");

	for(unsigned int kdx = 0; kdx < iso_labels.size(); ++kdx){
		double this_line_x_ndx = (0.9-TINY)*((iso_position[kdx]-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+0.1;
		TLine *isotope_line = new TLine();
		isotope_line->SetLineColor(Tableau20(0));
		isotope_line->SetLineWidth(3);
		isotope_line->DrawLineNDC(this_line_x_ndx, 0.1, this_line_x_ndx, 1-TINY);

		double text_width = 0.085;
		double text_height = 0.09;
		double label_y_ndc = 0.175;

		TPaveText *pt = new TPaveText(this_line_x_ndx-0.5*text_width, label_y_ndc-0.5*text_height,
					this_line_x_ndx+0.5*text_width, label_y_ndc+0.5*text_height, 
					"NB NDC");
		TText* this_tex = pt->AddText(iso_labels[kdx].c_str());
		pt->Draw();
		this_tex->SetTextColor(Tableau20(0));
	
	}
	double po_216_channel = 211.5;
	double bi_212_channel = 188.5;
	double pb_210_channel = 163.7;
	double line_216po_x_ndc = (0.9-TINY)*((po_216_channel-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+0.1; 
	double line_212bi_x_ndc = (0.9-TINY)*((bi_212_channel-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+0.1; 
	double line_210pb_x_ndc = (0.9-TINY)*((pb_210_channel-c_norm->GetUxmin())/(c_norm->GetUxmax()-c_norm->GetUxmin()))+0.1; 

	double y_axis_range[2] = {1E-1, 1E5};

	TLine *integr_line = new TLine();
	integr_line->SetLineColor(15);
	integr_line->SetLineStyle(2);
	integr_line->SetLineWidth(3);
	integr_line->DrawLineNDC(line_1_x_ndc, 0.1, line_1_x_ndc, 1-TINY);
	integr_line->DrawLineNDC(line_2_x_ndc, 0.1, line_2_x_ndc, 1-TINY);



	// Here the calculation of the limit starts
	for(unsigned int i = 0; i < trees.size(); ++i){
	
		double this_scale = 1000./(runtimes[i]*60.);

		// Generate n_rand_tf1 random TF1 crystal ball functions based on the 
		// covariance matrix of the fitted function
		// Then use this functions to estimate the uncertainty in the 
		// integration region

		int n_rand_tf1 = 10000;
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
			TF1* fit_func_mc = GetRandomFunction(fit_function.c_str(), u_mat,  fit_means);

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
	
			mc_integrals.push_back(fit_func_mc->Integral(integration_range[0], integration_range[1]));
			// Delete this random function, since it was created with new in the Method
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

			// For the display, the values need to be normalized by the runtime
			double scale = 1000./(runtimes[i]*60.);

			double this_x = mc_scan_x_values[j];
			std::vector<double> this_vector = mc_scan_points[j];
			double mean = TMath::Mean(this_vector.begin(), this_vector.end());
			double quantiles[4] = {0.15865, 0.84135, 0.02275, 0.997725 };
			double returns[4];
			TMath::Quantiles(this_vector.size(), 4, &this_vector[0], &returns[0], &quantiles[0], false);
			double std = TMath::RMS(this_vector.begin(), this_vector.end());
			confidence_band_1->SetPoint(j, this_x, mean*scale);
			confidence_band_1->SetPointError(j, 0., 0, (mean-returns[0])*scale, (returns[1]-mean)*scale);
			confidence_band_2->SetPoint(j, this_x, mean*scale);
			confidence_band_2->SetPointError(j, 0., 0, (mean-returns[2])*scale, (returns[3]-mean)*scale);

		}

		confidence_band_1->Draw("SAME E3");
		confidence_band_1->SetFillColorAlpha(hist_colors[i],0.35);
		// Make the graph markers invisible not not spoil the dashed fit line
		confidence_band_1->SetMarkerColorAlpha(hist_colors[i], 1.);
		//confidence_band_2->Draw("SAME E3");
		//confidence_band_2->SetFillColorAlpha(hist_colors[i],0.15);

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


	double axis_title_size = 40;	
	spectra[0]->GetXaxis()->SetTitle("Energy channel (A.U.)");
	spectra[0]->GetXaxis()->CenterTitle();
	spectra[0]->GetXaxis()->SetTitleFont(43);
	spectra[0]->GetXaxis()->SetLabelFont(43);
	spectra[0]->GetXaxis()->SetTitleSize(axis_title_size);
	spectra[0]->GetXaxis()->SetLabelSize(axis_title_size);

	spectra[0]->GetYaxis()->SetTitle("Detected rate (mHz)");
	spectra[0]->GetYaxis()->CenterTitle();
	spectra[0]->GetYaxis()->SetTitleFont(43);
	spectra[0]->GetYaxis()->SetLabelFont(43);
	spectra[0]->GetYaxis()->SetTitleSize(axis_title_size);
	spectra[0]->GetYaxis()->SetLabelSize(axis_title_size);


	TLegend *leg = new TLegend(0.15, 0.8, 0.35, 0.93);

	//leg->SetTextSize(40);
	leg->AddEntry(spectra[0], "uncoated");
        leg->AddEntry(spectra[1], "coated");

	
	leg->Draw();

	c_norm->SaveAs("dlc_comp_spec.pdf");
	c_norm->SaveAs("dlc_comp_spec.png");

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

	return 0;

}



