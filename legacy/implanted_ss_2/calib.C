void calib()
{
//=========Macro generated from canvas: calib/
//=========  (Tue Nov  2 10:16:48 2021) by ROOT version6.08/06
   TCanvas *calib = new TCanvas("calib", "",318,382,1200,900);
   gStyle->SetOptStat(0);
   calib->Range(132.1818,4567.394,352.0041,10864.7);
   calib->SetFillColor(0);
   calib->SetBorderMode(0);
   calib->SetBorderSize(2);
   calib->SetLeftMargin(0.1252087);
   calib->SetRightMargin(0.008347246);
   calib->SetTopMargin(0.01027397);
   calib->SetFrameBorderMode(0);
   calib->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1005[5] = {
   173.875,
   187.6106,
   210.763,
   334.0044,
   276.2706};
   Double_t Graph0_fy1005[5] = {
   5606.6,
   6050.78,
   6778.3,
   10552.1,
   8784.37};
   Double_t Graph0_fex1005[5] = {
   0.06110683,
   0.002215218,
   0.003361294,
   0.1324653,
   0.002196981};
   Double_t Graph0_fey1005[5] = {
   0.25,
   0.3,
   0.09,
   0.02,
   0.07};
   TGraphErrors *gre = new TGraphErrors(5,Graph0_fx1005,Graph0_fy1005,Graph0_fex1005,Graph0_fey1005);
   gre->SetName("Graph0");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_Graph_Graph100110031005 = new TH1F("Graph_Graph_Graph_Graph100110031005","",100,157.7816,350.1692);
   Graph_Graph_Graph_Graph100110031005->SetMinimum(5200);
   Graph_Graph_Graph_Graph100110031005->SetMaximum(10800);
   Graph_Graph_Graph_Graph100110031005->SetDirectory(0);
   Graph_Graph_Graph_Graph100110031005->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph_Graph100110031005->SetLineColor(ci);
   Graph_Graph_Graph_Graph100110031005->GetXaxis()->SetTitle("Energy channel (A.U.)");
   Graph_Graph_Graph_Graph100110031005->GetXaxis()->SetRange(2,100);
   Graph_Graph_Graph_Graph100110031005->GetXaxis()->CenterTitle(true);
   Graph_Graph_Graph_Graph100110031005->GetXaxis()->SetLabelFont(43);
   Graph_Graph_Graph_Graph100110031005->GetXaxis()->SetLabelSize(40);
   Graph_Graph_Graph_Graph100110031005->GetXaxis()->SetTitleSize(40);
   Graph_Graph_Graph_Graph100110031005->GetXaxis()->SetTitleFont(43);
   Graph_Graph_Graph_Graph100110031005->GetYaxis()->SetTitle("Alpha energy (keV)");
   Graph_Graph_Graph_Graph100110031005->GetYaxis()->CenterTitle(true);
   Graph_Graph_Graph_Graph100110031005->GetYaxis()->SetLabelFont(43);
   Graph_Graph_Graph_Graph100110031005->GetYaxis()->SetLabelSize(40);
   Graph_Graph_Graph_Graph100110031005->GetYaxis()->SetTitleSize(40);
   Graph_Graph_Graph_Graph100110031005->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph_Graph_Graph100110031005->GetYaxis()->SetTitleFont(43);
   Graph_Graph_Graph_Graph100110031005->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph_Graph100110031005->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph_Graph100110031005->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph_Graph100110031005->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph_Graph100110031005);
   
   
   TF1 *ch_2_e1006 = new TF1("ch_2_e","[0] + [1]*x",157.7816,350.1692);
   ch_2_e1006->SetFillColor(19);
   ch_2_e1006->SetFillStyle(0);
   ch_2_e1006->SetLineColor(2);
   ch_2_e1006->SetLineWidth(2);
   ch_2_e1006->SetChisquare(2996.288);
   ch_2_e1006->SetNDF(3);
   ch_2_e1006->GetXaxis()->SetLabelFont(42);
   ch_2_e1006->GetXaxis()->SetLabelSize(0.035);
   ch_2_e1006->GetXaxis()->SetTitleSize(0.035);
   ch_2_e1006->GetXaxis()->SetTitleFont(42);
   ch_2_e1006->GetYaxis()->SetLabelFont(42);
   ch_2_e1006->GetYaxis()->SetLabelSize(0.035);
   ch_2_e1006->GetYaxis()->SetTitleSize(0.035);
   ch_2_e1006->GetYaxis()->SetTitleFont(42);
   ch_2_e1006->SetParameter(0,306.2213);
   ch_2_e1006->SetParError(0,0.5480183);
   ch_2_e1006->SetParLimits(0,0,0);
   ch_2_e1006->SetParameter(1,30.68971);
   ch_2_e1006->SetParError(1,0.002172998);
   ch_2_e1006->SetParLimits(1,0,0);
   gre->GetListOfFunctions()->Add(ch_2_e1006);
   gre->Draw("ap");
   calib->Modified();
   calib->cd();
   calib->SetSelected(calib);
}
