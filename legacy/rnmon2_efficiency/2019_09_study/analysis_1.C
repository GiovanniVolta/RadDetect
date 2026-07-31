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


std::pair<double, double> GetStandardCounts(std::pair<double,double> init_act, double start, double end){
	// Function that returns the expected amount of counts (and uncertainty) within a 
	// given time interval 
	
	// half life of 222Rn
	double t1_2_rn222 = 3.824;
	double rn_decay_const = log(2)/3.824;
	// Number of atoms N(t) = N_0 * exp(-t*lambda)
	// N_0 = init_act * 1/lambda
	double init_num_rn = init_act.first / 1000. * 3600*24. /rn_decay_const;
	double init_num_rn_err = init_act.second / 1000. * 3600*24. /rn_decay_const;

	// number of atoms present at interval start
	double decay_start = init_num_rn * exp(-start*rn_decay_const);
	double decay_start_err = init_num_rn_err * exp(-start*rn_decay_const);

	// number of atoms present at interval end
	double decay_end = init_num_rn * exp(-end*rn_decay_const);
        double decay_end_err = init_num_rn_err * exp(-end*rn_decay_const);

	double num_decays = decay_start - decay_end;
	double num_decays_err = decay_start_err - decay_end_err;
	std::pair<double, double> result(num_decays, num_decays_err);

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
	std::vector<TTree*> trees;
	std::string base_file_path = "/d12/lin/xenon/radon_monitor/database/backup/";
	std::string base_file_ending = ".root";

	int n_bins_overview = 500;
	double reject_ingrowth = 0.2;
	// Activites of the 222Rn standard as from the database
	std::pair<double, double> standard_act_err(192.51, 3.46);
	std::vector<double> effiecies_214;
        std::vector<double> effiecies_214_err;
        std::vector<double> effiecies_218;
        std::vector<double> effiecies_218_err;

	std::vector<double> p_values;
	std::vector<double> p_errors;
	double std_p_err = 2;
	std::vector<double> u_values;
        std::vector<double> u_errors;
	double std_u_err = 1;


	//First define the runs of interest and put them in the vectors

	run_labels.push_back("850mbar");
	filenames.push_back("Ap19092019");
	p_values.push_back(850);
	p_errors.push_back(std_p_err);
	channel_cuts_218.push_back(std::make_pair(120,200));
	channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("900mbar");
        filenames.push_back("Ap20092019_900mbar");
        p_values.push_back(900);
        p_errors.push_back(std_p_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("950mbar");
        filenames.push_back("Ap20092019");
        p_values.push_back(950);
        p_errors.push_back(std_p_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("1000mbar");
        filenames.push_back("Ap23092019");
        p_values.push_back(1000);
        p_errors.push_back(std_p_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("1050mbar");
        filenames.push_back("Ap24092019");
        p_values.push_back(1050);
        p_errors.push_back(std_p_err);
	// This run is also the first for the voltage dependence
        u_values.push_back(500);
        u_errors.push_back(std_u_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("750V");
        filenames.push_back("Ap25092019");
        u_values.push_back(750);
        u_errors.push_back(std_u_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("1000V");
        filenames.push_back("Ap27092019");
        u_values.push_back(1000);
        u_errors.push_back(std_u_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("1250V");
        filenames.push_back("Ap30092019");
        u_values.push_back(1250);
        u_errors.push_back(std_u_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));

        run_labels.push_back("1000V");
        filenames.push_back("Ap04102019");
        u_values.push_back(1000);
        u_errors.push_back(std_u_err);
        channel_cuts_218.push_back(std::make_pair(120,200));
        channel_cuts_214.push_back(std::make_pair(250,320));



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

		int channel;
		Long64_t timestamp;
		trees[i]->SetBranchAddress("timestamp", &timestamp);
                trees[i]->SetBranchAddress("channel", &channel);

		for(unsigned int j = 0; j < trees[i]->GetEntries(); ++j){
			trees[i]->GetEvent(j);

			if (channel > channel_cuts_218[i].first && channel < channel_cuts_218[i].second){
				h_po_218->Fill( (timestamp - unix_min)/3600./24 );
                                if (i == 0){
                                        if ( (timestamp - unix_min)/3600./24. < reject_ingrowth ){
                                                continue;
					}
				}

				++counts_218[i];
			}
			if (channel > channel_cuts_214[i].first && channel < channel_cuts_214[i].second){
                                h_po_214->Fill( (timestamp - unix_min)/3600./24 );
				if (i == 0){
					if ( (timestamp - unix_min)/3600./24. < reject_ingrowth ){
						continue;
					}
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
	TLine *indicator_ingrowth = new TLine(reject_ingrowth, c_overview->GetUymin(), reject_ingrowth, c_overview->GetUymax());
	indicator_ingrowth->Draw("SAME");
	indicator_ingrowth->SetLineWidth(2);
	indicator_ingrowth->SetLineStyle(2);
	indicator_ingrowth->SetLineColor(2);

	// Draw lines to indicate the different running intervals
	
	for(unsigned int i = 0; i < unix_starts.size(); ++i){
	
		double this_time = (unix_starts[i]-unix_min)/3600./24.;
		TLine *indicator_line = new TLine(this_time,c_overview->GetUymin(),this_time,c_overview->GetUymax());
		indicator_line->Draw("SAME");
		indicator_line->SetLineStyle(2);
		indicator_line->SetLineWidth(2);
		indicator_line->SetLineColor(kGray);

		std::cout<<"interval: "<<i<<" Po218: "<<counts_218[i]<<" Po214: "<<counts_214[i]<<std::endl;
		std::cout<<"label: "<<run_labels[i]<<std::endl;
		std::pair<double, double> standard_counts;
	       
		if(i == 0){
			standard_counts = GetStandardCounts(standard_act_err, (unix_starts[i]-unix_min)/3600./24. + reject_ingrowth, (unix_ends[i]-unix_min)/3600./24.);
		}else{
                        standard_counts = GetStandardCounts(standard_act_err, (unix_starts[i]-unix_min)/3600./24., (unix_ends[i]-unix_min)/3600./24.);
		}

		std::cout<<"expected counts: "<<standard_counts.first<<" +/- "<<standard_counts.second<<std::endl;

		// calculate the 214Po efficiency and error
		double efficiency_214 = counts_214[i]/standard_counts.first;
		double efficiency_err_214 = sqrt(pow(efficiency_214,2) * pow(sqrt(counts_214[i])/counts_214[i],2) + pow(standard_counts.second/standard_counts.first,2));
	        effiecies_214.push_back(efficiency_214);
	        effiecies_214_err.push_back(efficiency_err_214);
		std::cout<<"Efficiency(Po214): "<<efficiency_214<<" +/- "<<efficiency_err_214<<std::endl;

                // calculate the 218Po efficiency and error
                double efficiency_218 = counts_218[i]/standard_counts.first;
                double efficiency_err_218 = sqrt(pow(efficiency_218,2) * pow(sqrt(counts_218[i])/counts_218[i],2) + pow(standard_counts.second/standard_counts.first,2));;
                effiecies_218.push_back(efficiency_218);
                effiecies_218_err.push_back(efficiency_err_218);
                std::cout<<"Efficiency(Po214): "<<efficiency_218<<" +/- "<<efficiency_err_218<<std::endl;

	}


	/*************************
	 * Plotting efficiencies *
	 * ***********************/

	// Pressure
	TCanvas *c_p_dep = new TCanvas("p_dependence");
	Int_t num_p_runs = (int) p_values.size();
	TLegend *p_legend = new TLegend(0.7, 0.7, 0.9,0.9);

	TGraphErrors *p_dep_efficiency_214 = new TGraphErrors(num_p_runs, &p_values[0], &effiecies_214[0], &p_errors[0], &effiecies_214_err[0]);
	p_dep_efficiency_214->Draw("AP");
	p_dep_efficiency_214->SetLineColor(4);
	p_legend->AddEntry(p_dep_efficiency_214, "Po-214", "l");

        TGraphErrors *p_dep_efficiency_218 = new TGraphErrors(num_p_runs, &p_values[0], &effiecies_218[0], &p_errors[0], &effiecies_218_err[0]);
        p_dep_efficiency_218->Draw("SAMEP");
        p_dep_efficiency_218->SetLineColor(2);
        p_legend->AddEntry(p_dep_efficiency_218, "Po-218", "l");


	p_dep_efficiency_214->GetXaxis()->SetTitle("Pressure (mbar)");
        p_dep_efficiency_214->GetYaxis()->SetTitle("Efficiency");
        p_dep_efficiency_214->SetTitle("Monitor Pressure");
	p_dep_efficiency_214->GetYaxis()->SetRangeUser(.24,.38);

	p_legend->Draw();

	// High Voltage
        TCanvas *c_u_dep = new TCanvas("u_dependence");
        Int_t num_u_runs = (int) u_values.size();
        TLegend *u_legend = new TLegend(0.7, 0.7, 0.9,0.9);


        TGraphErrors *u_dep_efficiency_214 = new TGraphErrors(num_u_runs, &u_values[0], &effiecies_214[num_p_runs-1], &u_errors[0], &effiecies_214_err[num_p_runs-1]);
        u_dep_efficiency_214->Draw("AP");
        u_dep_efficiency_214->SetLineColor(4);
        u_legend->AddEntry(u_dep_efficiency_214, "Po-214", "l");


        TGraphErrors *u_dep_efficiency_218 = new TGraphErrors(num_u_runs, &u_values[0], &effiecies_218[num_p_runs-1], &u_errors[0], &effiecies_218_err[num_p_runs-1]);
        u_dep_efficiency_218->Draw("SAMEP");
        u_dep_efficiency_218->SetLineColor(2);
        u_legend->AddEntry(u_dep_efficiency_218, "Po-218", "l");


        u_dep_efficiency_214->GetXaxis()->SetTitle("Voltage (V)");
        u_dep_efficiency_214->GetYaxis()->SetTitle("Efficiency");
	u_dep_efficiency_214->SetTitle("High Voltage");
        u_dep_efficiency_214->GetYaxis()->SetRangeUser(.24,.38);

	u_legend->Draw();


	return 0;

}



