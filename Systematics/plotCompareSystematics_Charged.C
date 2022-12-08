/*
 * plotCompareSystematics.C: Plot the relative systematic uncetainties.
 * Caitie Beattie <caitie.beattie@yale.edu>
 *
 */

void plotCompareSystematics_Charged(const char* tag = "hi", int R = 4) {
 

  bool high = false;
  bool low = false;
  bool highlow = false;
  bool highinout = false;
  bool lowinout = false;
  bool lowin = false;
  bool lowout = false;
  bool highin = false;
  bool highout = false;
  bool epup = false;

  const char* hl = "hilo";
  const char* li = "loin";
  const char* lo = "loout";
  const char* hi = "hiin";
  const char* ho = "hiout";
  const char* hio = "hiinout";
  const char* lio = "loinout";
  const char* oup = "outup";
  const char* iup = "inup";
  if (strncmp(tag, hl, 4) == 0)   highlow = true;
  if (strncmp(tag, li, 4) == 0)   lowin = true;
  if (strncmp(tag, lo, 5) == 0)   lowout = true;
  if (strncmp(tag, hi, 5) == 0)   highin = true;
  if (strncmp(tag, ho, 5) == 0)   highout = true;
  if (strncmp(tag, hio, 7) == 0)  highinout = true;
  if (strncmp(tag, lio, 7) == 0)  lowinout = true;
  if (strncmp(tag, oup, 5) == 0 || strncmp(tag, iup, 4) == 0)   epup = true;



  int sysnum = 5;
  if (highlow == true || epup == true) sysnum = 4;

  stringstream ssinput;
  //ssinput << "preSystematics_R020_3050_" << tag << "_Jul8.root";
  if (R == 2) ssinput << "preSystematics_R020_3050_20_" << tag << "_Dec1.root";
  if (R == 4) ssinput << "preSystematics_R040_3050_30_" << tag << "_Dec1.root";

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // Get the file for the systematic uncertanties
  TFile* _file0     = TFile::Open(ssinput.str().c_str());

  //build the legends
  TLegend* leg = new TLegend(0.55, 0.3, 0.8, 0.65);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  TLegend* lat = new TLegend(0.15, 0.7, 0.4, 0.85);
  //lat->SetTextFont(42);
  lat->SetBorderSize(0);
  lat->SetTextSize(0.035);
  lat->SetFillColorAlpha(kWhite, 0.0);
  lat->AddEntry((TObject*)0, "ALICE Pb--Pb 5.02 TeV, 30-50%", "");
  if (R == 2) lat->AddEntry((TObject*)0, "Charged jets, #it{R} = 0.2, |#eta_{jet}| < 0.7", "");
  if (R == 4) lat->AddEntry((TObject*)0, "Charged jets, #it{R} = 0.4, |#eta_{jet}| < 0.5", "");
  if (highlow == true) lat->AddEntry((TObject*)0, "high/low #it{q}_{2}^{V0C}", "");
  if (highinout == true)  lat->AddEntry((TObject*)0, "out-of-plane/in-plane, highest 30\%", "");
  if (lowinout == true)   lat->AddEntry((TObject*)0, "out-of-plane/in-plane, lowest 30\%", "");
  if (lowin == true && lowinout == false)  lat->AddEntry((TObject*)0, "in-plane, lowest 30\%", "");
  if (lowout == true)  lat->AddEntry((TObject*)0, "out-of-plane, lowest 30\%", "");
  if (highin == true && highinout == false)  lat->AddEntry((TObject*)0, "in-plane, highest 30\%", "");
  if (highout == true)  lat->AddEntry((TObject*)0, "out-of-plane, highest 30\%", "");
  //lat->AddEntry((TObject*)0, "highest 30\%, in-plane", "");

  //retrieve errors from file
  TGraphAsymmErrors* g_nom_err = (TGraphAsymmErrors*)_file0->Get("h1_nom_err");
  TH1D* hnom = (TH1D*)_file0->Get("h1_nom_copy");
  TH1D* hnom_err = (TH1D*)hnom->Clone("h_err");
  hnom_err->Reset();
  
   stringstream ssoutput;
   if (R == 2) ssoutput << "Systematics3050_R020_" << tag << "_Dec1.root";
   if (R == 4) ssoutput << "Systematics3050_R040_" << tag << "_Dec1.root";

   TFile *fout = new TFile (ssoutput.str().c_str(),"RECREATE");
   fout->cd();

  cout << "Total Systematics per Bin: \n";
  for (int i = 1; i <= hnom_err->GetNbinsX(); i++)   {
    double x = 0;
    double y = 0;
    //cout << "hnom(i): " << hnom->GetBinContent(i) <<"\n";
    //cout << "g_nom_err(i): " << g_nom_err->GetPoint(i-1, x, y) <<"\n";
    g_nom_err->GetPoint(i-1, x, y);
    //cout << "x: " << x << ", y: " << y <<"\n";
    //double y_low = g_nom_err->GetErrorYlow(i-1)/y;
    //get error from previos file and convert to fractional error
    float y_high = g_nom_err->GetErrorYhigh(i-1)/y;
    float err = std::sqrt(y_high*y_high);
    hnom_err->SetBinContent(i, err);
    if (i >= 5+(R/2) && i < 12) cout << "	Bin: " << i<< ", " << hnom_err->GetXaxis()->GetBinLowEdge(i) << " - " << hnom_err->GetXaxis()->GetBinLowEdge(i+1) << ", Error: " << err <<"\n";
  }
  hnom_err->Write();
  fout->Close();

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  c1->cd();

  hnom_err->SetLineColor(kBlack);
  leg->AddEntry(hnom_err, "Total", "l");
  hnom_err->GetYaxis()->SetRangeUser(0., 1.);
  hnom_err->GetXaxis()->SetRangeUser(25.01,120);
  hnom_err->GetYaxis()->SetTitle("Relative Systematic Uncertainty");
  hnom_err->GetYaxis()->SetTitleOffset(1.0);
  hnom_err->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hnom_err->SetLineWidth(4);
  hnom_err->SetMaximum(0.3);
  hnom_err->Draw("hist");


  //int col[20] = {kRed+2, kRed-4, kOrange+7, kOrange, kYellow-4, kSpring+10, kSpring, kGreen-3, kGreen+3, kTeal-7, kTeal, kAzure+10, kAzure-4, kBlue+2, kViolet+8, kViolet-1,kMagenta+1, kMagenta-4, kPink+7, kPink-4};
  int col[20] = {kBlue, kBlue, kCyan, kCyan, kMagenta, kMagenta, kViolet+7, kViolet+7, kOrange+7, kOrange+7, kGreen+3, kTeal-7, kTeal, kAzure+10, kAzure-4, kBlue+2, kViolet+8, kViolet-1,kMagenta+1, kMagenta-4};

  std::string label[10] = {"Truncation", "Iterations", "Reweighting Prior", "Tracking Efficiency", "EP Resolution"};
  std::string sys[10] = {"trunc", "iter", "rw", "eff", "gluons"};
  std::stringstream ss;

  //loop through systematics
  for (int i = 0; i  < sysnum; i++){
      TH1D* h1;
      TH1D* hratio;
      if ((sys[i] != "trunc") && (sys[i] != "iter"))
	{
	  ss << "h1_" << sys[i] << "_copy";
	  //cout <<ss.str()<<endl;
	  h1 = (TH1D*)_file0->Get(ss.str().c_str());
	  ss << "_ratiio";
	  hratio =(TH1D*)h1->Clone(ss.str().c_str());
	  ss.str("");
	  hratio->Add(hnom, -1);
	  for (int j = 1; j <= hratio->GetNbinsX(); j++){
	      double cont = hratio->GetBinContent(j);
	      hratio->SetBinContent(j, std::sqrt(cont*cont));
	  }
	}

      else{

	std::string ty1;
	std::string ty2;
	if (sys[i] == "trunc") { ty1 = "plus5"; ty2 = "min5";}
	else {ty1 = "iterhigh"; ty2 = "iterlow";}
	ss << "h1_" << ty1 << "_copy";
	//cout << ss.str() <<endl;
	h1 = (TH1D*)_file0->Get(ss.str().c_str());
	ss << "_ratiio";
	hratio =(TH1D*)h1->Clone(ss.str().c_str());
	ss.str("");
	hratio->Add(hnom, -1);
	ss <<"h1_" << ty2 << "_copy";
	//cout << ss.str() <<endl;

	TH1D* h11 = (TH1D*)_file0->Get(ss.str().c_str());
	ss << "_ratiio";
	TH1D* hratio1 =(TH1D*)h11->Clone(ss.str().c_str());
	ss.str("");
	hratio1->Add(hnom, -1);
	for (int i = 1; i <= hratio->GetNbinsX(); i++)
	  {
	    double cont =  hratio->GetBinContent(i);
	    double cont1 =  hratio1->GetBinContent(i);
	    hratio->SetBinContent(i, std::sqrt(cont*cont+cont1*cont1));
	  }
      }
      hratio->Divide(hnom);
      cout <<"\n";
      cout << label[i] <<"\n";
      for (int z = 5+(R/2); z <= 11; z++) {
          cout << "   Bin: " << z << ", " << hratio->GetXaxis()->GetBinLowEdge(z) << " - " 
               << hratio->GetXaxis()->GetBinLowEdge(z+1) << ", Error: " << hratio->GetBinContent(z) <<"\n";} 
      hratio->SetLineColor(col[i*2]);
      hratio->SetLineWidth(4);
      hratio->Draw("same Hist");
      leg->AddEntry(hratio, label[i].c_str(), "l");
    }
  leg->Draw("same");
  lat->Draw("same");
  hnom_err->Draw("hist same");





  /*
  TH1D* iterHighErr = (TH1D*)_file0->Get("h1_iterhigh_copy_delta_up");
  TH1D* iterLowErr  = (TH1D*)_file0->Get("h1_iterlow_copy_delta_down");
  TH1D* plus5Err    = (TH1D*)_file0->Get("h1_plus5_copy_delta_up");
  TH1D* min5Err     = (TH1D*)_file0->Get("h1_min5_copy_delta_down");

  //plot the errors as percent on the total
  iterHighErr->Divide(nominal);
  iterLowErr->Divide(nominal);
  plus5Err->Divide(nominal);
  min5Err->Divide(nominal);

  plus5Err->Draw("HIST");
  iterHighErr->Draw("HIST same");
  iterLowErr->Draw("HIST same");
  min5Err->Draw("HIST same");
  */






}
