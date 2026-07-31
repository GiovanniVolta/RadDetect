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
#include <fstream>

double RN_DECAY_CONST = log(2)/3.824;

double GetTimeDifferenceFromString(string time1, string time2){

	cout << "WARNING!: The time string must be entered as follows: YYYYMMDD_HHMMSS!" << endl;

	struct tm time1_tm;

	time1_tm.tm_year   = stoi(time1.substr (0,4));
	time1_tm.tm_mon    = stoi(time1.substr (4,2));
	time1_tm.tm_mday   = stoi(time1.substr (6,2));
	time1_tm.tm_hour   = stoi(time1.substr (9,2));
	time1_tm.tm_min    = stoi(time1.substr (11,2));
	time1_tm.tm_sec    = stoi(time1.substr (13,2));

	cout << "Entered time1 time: Day: " << time1_tm.tm_year << " " << time1_tm.tm_mon
	     << " " << time1_tm.tm_mday << " Time: " << time1_tm.tm_hour << ":" << time1_tm.tm_min
	     << ":" << time1_tm.tm_sec << endl;

	time1_tm.tm_isdst=0;

	// For proper counting of the year and month
	time1_tm.tm_mon-=1;
	time1_tm.tm_year-=1900;
 
	// -------------- // 
	// Get the time2  //
	// -------------- //
	
 	struct tm time2_tm;

	time2_tm.tm_year   = stoi(time2.substr (0,4));
	time2_tm.tm_mon    = stoi(time2.substr (4,2));
	time2_tm.tm_mday   = stoi(time2.substr (6,2));
	time2_tm.tm_hour   = stoi(time2.substr (9,2));
	time2_tm.tm_min    = stoi(time2.substr (11,2));
	time2_tm.tm_sec    = stoi(time2.substr (13,2));
 
	// For proper counting of the year and month
	time2_tm.tm_mon-=1;
	time2_tm.tm_year-=1900;
 
	double seconds_from_start = difftime(mktime(&time1_tm), mktime(&time2_tm));

	return seconds_from_start;
}




std::pair<double, double> GetEmanationYield(std::string emanation_start, std::string emanation_end, double start_uncertainty = 0, double end_uncertainty = 0){
	// Function that returns the expected emanation yield (and uncertainty) after a given emanation 
	// time interval start and end time uncertainty are given in minutes!
	
	double seconds_emanation = GetTimeDifferenceFromString(emanation_start, emanation_end);
	double days_emanation = seconds_emanation/3600./24.;

	// calculate the yield as Y = ( 1-exp(-t*lambda) )
	double yield = 1 - exp(days_emanation*RN_DECAY_CONST);
	// For now dont calc the uncertainty
	double time_uncert = sqrt(start_uncertainty*start_uncertainty + end_uncertainty*end_uncertainty);
	time_uncert = time_uncert /60./24.; 
	double yield_uncert = RN_DECAY_CONST * exp(RN_DECAY_CONST*days_emanation) * time_uncert;

	std::pair<double, double> result(yield, yield_uncert);

	return result;
}

int analysis(){

	// define the run names for each run
	std::vector<std::string> run_labels;
	std::vector<int> colors;
	std::vector<int> n_bins;
	std::vector<std::pair<int, int>> channel_cuts_218;
        std::vector<std::pair<int, int>> channel_cuts_214;
	std::vector<std::string> filenames;
	std::vector<std::string> emanation_starts;
	std::vector<std::string> emanation_stops;
	std::vector<double> start_offsets;
	std::vector<TTree*> trees;
	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	int n_bins_overview = 500;
	std::vector<double> reject_ingrowth_vec;
	double reject_ingrowth_default = 0.5;
	//double reject_ingrowth = 0.;

	// Detection efficiency of bluemchen and its error
	// std::pair<double, double> detection_eff = {0.315, 0.015};
	// std::pair<double, double> detection_eff = {0.4327213, 0.06};

	// Efficiency from the two reliable PTB standard determinations
	// std::pair<double, double> detection_eff = {0.3547415, 0.0206204};

	// On 2022-11-04 it was found that two of the three efficiency measurements
	// were carried out directly following the high-activity from the ISOLDE sample
	// 222Rn might still have been present in the detector and let to a too high
	// estimate of the detection efficiency! Therefore, using only second PTB measurement.
	std::pair<double, double> detection_eff = {0.3466,0.0198};


	std::vector<double> activities_214;
        std::vector<double> activities_214_stat_err;
        std::vector<double> activities_214_syst_err;
        std::vector<double> activities_218;
        std::vector<double> activities_218_err;

	std::vector<double> t_values;
	std::vector<double> t_errors;
	double std_t_err = 2;

	std::vector<double> p_values;
	std::vector<double> p_errors;
	double std_p_err = 5.;


	//First define the runs of interest and put them in the vectors

	run_labels.push_back("room temperature");
	filenames.push_back("Rn17022020");
	// Taking preliminary values from second R.T. Measurement
	t_values.push_back(22.095165204854016);
	t_errors.push_back(1.0887331000054044);
	p_values.push_back(202);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200212_192900");
	emanation_stops.push_back("20200217_141800");
	start_offsets.push_back( GetTimeDifferenceFromString(
				 "20200217_141800", "20200217_153832"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));

        run_labels.push_back("0 degC");
        filenames.push_back("Rn24022020");
        t_values.push_back(1.9798246863746094);
        t_errors.push_back(0.09541574360871456);
	p_values.push_back(202);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200217_160300");
	emanation_stops.push_back("20200224_141900");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200224_141900", "20200224_152332"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));
	
        run_labels.push_back("-20 degC");
        filenames.push_back("Rn02032020"); 
        t_values.push_back(-19.286789180642227);
        t_errors.push_back(1.0458367082902686);
	p_values.push_back(207);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200224_152800");
	emanation_stops.push_back("20200302_141600");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200302_141600", "20200302_150337"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));
	
	
	run_labels.push_back("-30 degC");
        filenames.push_back("Rn09032020"); 
        t_values.push_back(-29.73250321309943);
        t_errors.push_back(1.1759928447865482);
	p_values.push_back(207);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200302_151000");
	emanation_stops.push_back("20200309_151400");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200309_151400", "20200309_155629"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));


	// Additional -30degC measurement for cross check of systematics
	run_labels.push_back("-30 degC");
        filenames.push_back("Rn16032020"); 
        t_values.push_back(-29.7714711397152);
        t_errors.push_back(1.1578615834009196);
	p_values.push_back(198);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200309_160200");
	emanation_stops.push_back("20200316_141300");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200316_141300", "20200316_145804"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));

	//Additional room temperature measurement for cross check of systematics
	run_labels.push_back("22 degC");
        filenames.push_back("Rn02072020");
        t_values.push_back(23.292437278781772);
        t_errors.push_back(0.7974974945413555);
	p_values.push_back(210);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200316_150200");
	emanation_stops.push_back("20200702_151400");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200702_151400", "20200702_161328"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));

	/**************************************
 	 * MEASUREMENT AT DIFFERENT PRESSURES *
	 **************************************/
/*
	//room temperature measurement for cross check of systematics
	run_labels.push_back("22 degC");
        filenames.push_back("Rn09072020"); 
        t_values.push_back(22);
        t_errors.push_back(std_t_err);
	p_values.push_back(505);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200702_160300");
	emanation_stops.push_back("20200709_151100");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200709_151100", "20200709_155851"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));

	run_labels.push_back("22 degC");
        filenames.push_back("Rn20072020"); 
        t_values.push_back(22);
        t_errors.push_back(std_t_err);
	p_values.push_back(795);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200709_160400");
	emanation_stops.push_back("20200720_154300");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200720_154300", "20200720_162333"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));

	run_labels.push_back("22 degC");
        filenames.push_back("Rn30072020"); 
        t_values.push_back(22);
        t_errors.push_back(std_t_err);
	p_values.push_back(1074);
	p_errors.push_back(std_p_err);
	reject_ingrowth_vec.push_back(reject_ingrowth_default);
	emanation_starts.push_back("20200709_160400");
	emanation_stops.push_back("20200720_154300");
        start_offsets.push_back( GetTimeDifferenceFromString(
                                 "20200720_154300", "20200720_162333"));
	channel_cuts_218.push_back(std::make_pair(170,200));
	channel_cuts_214.push_back(std::make_pair(200,250));
*/
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
	std::vector<Long64_t> unix_starts;
	std::vector<Long64_t> unix_ends;
	
	for(unsigned int i = 0; i < trees.size(); ++i){

		Long64_t timestamp;

		trees[i]->SetBranchAddress("timestamp", &timestamp);
		trees[i]->GetEvent(0);
		
		unix_starts.push_back(timestamp);
		trees[i]->GetEvent(trees[i]->GetEntries() -1 );
                unix_ends.push_back(timestamp);
		std::cout<<"TTree: "<<i<<" Min: "<<unix_starts[i]<<" Max: "<<unix_ends[i]<<std::endl;
	}

	// FOLLOWING WILL ONLY WORK, IF RUNS ARE SORTED CHRONOLOGICALLY!!
        //Long64_t unix_min = std::min_element(std::begin(unix_starts), std::end(unix_starts));
	//Long64_t unix_max = std::max_element(std::begin(unix_ends), std::end(unix_ends));
	Long64_t unix_min = unix_starts[0];
	Long64_t unix_max = unix_ends[unix_ends.size() - 1];
	double duration = (unix_max-unix_min)/3600./24.;


	TH1D *h_po_214 = new TH1D("po_214", "", n_bins_overview, 0., duration);
        TH1D *h_po_218 = new TH1D("po_218", "", n_bins_overview, 0., duration);

	std::vector<int> counts_214(trees.size(),0);
	std::vector<int> counts_218(trees.size(),0);


	for(unsigned int i = 0; i < trees.size(); ++i){

		double reject_ingrowth = reject_ingrowth_vec[i]; 
		int channel;
		Long64_t timestamp;
		trees[i]->SetBranchAddress("timestamp", &timestamp);
                trees[i]->SetBranchAddress("channel", &channel);

		for(unsigned int j = 0; j < trees[i]->GetEntries(); ++j){
			trees[i]->GetEvent(j);

			if (channel > channel_cuts_218[i].first && channel < channel_cuts_218[i].second){
				h_po_218->Fill( (timestamp - unix_min)/3600./24 );
                                if ( (timestamp - unix_starts[i])/3600./24. < reject_ingrowth ){
					continue;
				}

				++counts_218[i];
			}
			if (channel > channel_cuts_214[i].first && channel < channel_cuts_214[i].second){
                                h_po_214->Fill( (timestamp - unix_min)/3600./24 );
				if ( (timestamp - unix_starts[i])/3600./24. < reject_ingrowth ){
					continue;
				}

				++counts_214[i];
			}
		}
	}

	h_po_214->Draw();
	h_po_218->Draw("SAME");

	h_po_214->GetXaxis()->SetTitle("Time (days)");
        h_po_214->GetYaxis()->SetTitle("Counts");


	h_po_214->SetLineColor(4);
	h_po_218->SetLineColor(2);

	c_overview->Update();

	// Draw line to indicate the ingrowth cut
	for(unsigned int i = 0; i < trees.size(); ++i){
		double reject_ingrowth = reject_ingrowth_vec[i];
		TLine *indicator_ingrowth = new TLine(reject_ingrowth+(unix_starts[i]-unix_min)/24./3600., c_overview->GetUymin(), reject_ingrowth+(unix_starts[i]-unix_min)/24./3600., c_overview->GetUymax());
		indicator_ingrowth->Draw("SAME");
		indicator_ingrowth->SetLineWidth(2);
		indicator_ingrowth->SetLineStyle(2);
		indicator_ingrowth->SetLineColor(2);
	}
	c_overview->Update();

	// Draw lines to indicate the different running intervals
	
	for(unsigned int i = 0; i < unix_starts.size(); ++i){
	
		double this_time = (unix_starts[i]-unix_min)/3600./24.;
		double reject_ingrowth = reject_ingrowth_vec[i]; 
		TLine *indicator_line = new TLine(this_time,c_overview->GetUymin(),this_time,c_overview->GetUymax());
		indicator_line->Draw("SAME");
		indicator_line->SetLineStyle(2);
		indicator_line->SetLineWidth(2);
		indicator_line->SetLineColor(kGray);

		// let calculate the emanation yield
		std::pair<double,double> yield = GetEmanationYield(emanation_starts[i], emanation_stops[i], 30, 45);
		
		std::cout<<"interval: "<<i<<" Po218: "<<counts_218[i]<<" Po214: "<<counts_214[i]<<std::endl;
		std::cout<<"label: "<<run_labels[i]<<std::endl;
		std::cout<<"Emanation yield: "<<yield.first<<" +/- "<<yield.second<<std::endl;

		// calculate the emanation using 214Po 

		// This is the fraction to which the radon present in the monitor has decayed while it was measured (i.e. moment DAQ was started - rejection -> moment it was stopped)
		double detected_214_fraction = exp(-RN_DECAY_CONST * (start_offsets[i]/24./3600. + reject_ingrowth)) - exp(-RN_DECAY_CONST * (start_offsets[i] + unix_ends[i] - unix_starts[i])/24./3600.);
		double detected_214_fraction_err = 0;
std::cout<<"detected fraction: "<<detected_214_fraction<<std::endl;

		// This is the extrapolated value of A0 as it has been detected (converted in Bq)
		double detected_214_activity = counts_214[i] * RN_DECAY_CONST * 1/detected_214_fraction / 24./3600.;
		// INCORECT
		double detected_214_activity_err = sqrt(counts_214[i]) * RN_DECAY_CONST * 1/detected_214_fraction / 24./3600.;
std::cout<<"detected activity: "<<detected_214_activity<<" +/- "<<detected_214_activity_err<<std::endl;

		// This is the samples activity, corrected for the emanation yield
		double yield_corrected_activity_214 = detected_214_activity * 1./yield.first;
		double yield_corrected_activity_err_214 = yield_corrected_activity_214 * sqrt(  pow(detected_214_activity_err/detected_214_activity,2) + 
								pow(yield.second/yield.first,2));
std::cout<<"yield corrected activity: "<<yield_corrected_activity_214<<" +/- "<<yield_corrected_activity_err_214<<std::endl;

		// This is the samples activity corrected by the detection efficiency
		double real_activity_214 = yield_corrected_activity_214 * 1./detection_eff.first;
		double real_activity_214_err = real_activity_214 * 
							sqrt(	pow(yield_corrected_activity_err_214/yield_corrected_activity_214, 2) + 
								pow(detection_eff.second/detection_eff.first, 2));
std::cout<<"eff_corr activity: "<<real_activity_214<<" +/- "<<real_activity_214_err<<std::endl;



std::cout<<"cross-check: avg act: "<<counts_214[i] / (unix_ends[i] - unix_starts[i] - reject_ingrowth*24.*3600.)<<std::endl;


		activities_214.push_back(real_activity_214);
		activities_214_stat_err.push_back(yield_corrected_activity_err_214);
		activities_214_syst_err.push_back(real_activity_214_err);

		// From estimation: 2% (1% due to Rn yield +/- 10 minutes, 1% due to statistic)
		//activities_214_err.push_back(activity_214* sqrt( 1./sqrt(counts_214[i]) * 1./sqrt(counts_214[i]) + 0.01*0.01));

	}


	/**********************
	 * Plotting activites *
	 * ********************/

	// Pressure
	TCanvas *c_t_dep = new TCanvas("t_dependence");
	Int_t num_t_runs = (int) t_values.size();
	TLegend *t_legend = new TLegend(0.7, 0.7, 0.9,0.9);

	TGraphErrors *t_dep_activity_214 = new TGraphErrors(num_t_runs, &t_values[0], &activities_214[0], &t_errors[0], &activities_214_syst_err[0]);
	// Export the reconstructed points and errors to a file
	ofstream out_file ("data_points_out.dat");
	if(out_file.is_open()){
		out_file<<"Temp (degC)\tTemp err (degC)\tAct (Bq)\t Act Err [stat] (Bq)\t Act Err [tot] (Bq)\n";
		for(int idx = 0; idx < num_t_runs; ++idx){
			out_file<<t_values[idx]<<"\t"<<t_errors[idx]<<"\t"<<activities_214[idx]<<"\t"<<activities_214_stat_err[idx]<<"\t"<<activities_214_syst_err[idx]<<"\n";
		}
		out_file.close();
	}
	t_dep_activity_214->Draw("AP");
	t_dep_activity_214->SetLineColor(4);
	t_legend->AddEntry(t_dep_activity_214, "Rn-222 emanation", "l");

	TF1 *const_f = new TF1("const_f", "pol0", -40, 30);
	t_dep_activity_214->Fit(const_f, "", "");
	double fit_act = const_f->GetParameter(0);
	double fit_act_uncert = const_f->GetParError(0);
	c_t_dep->Update();
	TBox *uncert_box = new TBox(gPad->GetUxmin(), fit_act - fit_act_uncert, gPad->GetUxmax(), fit_act + fit_act_uncert);
	uncert_box->SetFillColorAlpha(kGray, 0.7);
	uncert_box->Draw("SAME");

	// I don't know how to better move the data points in front of the box...
	t_dep_activity_214->Draw("SAME P");

/*
        TGraphErrors *t_dep_efficiency_218 = new TGraphErrors(num_t_runs, &t_values[0], &effiecies_218[0], &t_errors[0], &effiecies_218_err[0]);
        t_dep_efficiency_218->Draw("SAMEP");
        t_dep_efficiency_218->SetLineColor(2);
        t_legend->AddEntry(t_dep_efficiency_218, "Po-218", "l");
*/

	t_dep_activity_214->GetXaxis()->SetTitle("Temperature (degC)");
        t_dep_activity_214->GetYaxis()->SetTitle("Activity (Bq)");
        t_dep_activity_214->SetTitle("");
	t_dep_activity_214->GetYaxis()->SetRangeUser(1.5,2.2);

	t_legend->Draw();
	
	c_t_dep->SaveAs("emanation_t_dependence.pdf");

	return 0;

}



