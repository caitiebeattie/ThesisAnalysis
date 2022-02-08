/*
 * plotCompareSystematics.C: Plot the relative systematic uncetainties.
 * Hannah Bossi <hannah.bossi@yale.edu>
 *
 */
void plotCompareSystematics_Charged(const char* tag = "hi") {
 

  stringstream ssinput;
  ssinput << "preSystematics_R020_3050_" << tag << "_Feb4.root";


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // Get the file for the systematic uncertanties
  TFile* _file0     = TFile::Open(ssinput.str().c_str());

  //build the legends
  TLegend* leg = new TLegend(0.15, 0.6, 0.3, 0.8);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  TLatex* lat = new TLatex();
  lat->SetTextFont(42);
  lat->SetNDC(true);
  lat->SetTextSize(0.035);
  lat->SetTextAlign(33);


  //retrieve errors from file
  TGraphAsymmErrors* g_nom_err = (TGraphAsymmErrors*)_file0->Get("h1_nom_err");
  TH1D* hnom = (TH1D*)_file0->Get("h1_nom_copy");
  TH1D* hnom_err = (TH1D*)hnom->Clone("h_err");
  hnom_err->Reset();
  
   stringstream ssoutput;
   ssoutput << "Systematics3050_" << tag << ".root";

   TFile *fout = new TFile (ssoutput.str().c_str(),"RECREATE");
   fout->cd();

  //
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
    cout << "g_nom_err->GetErrorYhigh(i-1): " << g_nom_err->GetErrorYhigh(i-1)<<"\n";
    float err = std::sqrt(y_high*y_high);
    hnom_err->SetBinContent(i, err);
    cout << "Bin: " << i<< ", " << hnom_err->GetXaxis()->GetBinLowEdge(i) << " - " << hnom_err->GetXaxis()->GetBinLowEdge(i+1) << ", Error: " << err <<"\n";
  }
  hnom_err->Write();
  fout->Close();

  hnom_err->SetLineColor(kBlack);
  leg->AddEntry(hnom_err, "Total", "l");
  hnom_err->GetYaxis()->SetRangeUser(0., 1.);
  hnom_err->GetXaxis()->SetRangeUser(40.01,100);
  hnom_err->GetYaxis()->SetTitle("Relative Systematic Uncertainty");
  hnom_err->GetYaxis()->SetTitleOffset(1.0);
  hnom_err->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hnom_err->SetLineWidth(4);
  hnom_err->Draw("hist");


  //int col[20] = {kRed+2, kRed-4, kOrange+7, kOrange, kYellow-4, kSpring+10, kSpring, kGreen-3, kGreen+3, kTeal-7, kTeal, kAzure+10, kAzure-4, kBlue+2, kViolet+8, kViolet-1,kMagenta+1, kMagenta-4, kPink+7, kPink-4};
  int col[20] = {kBlue, kBlue, kCyan, kCyan, kMagenta, kMagenta, kSpring, kGreen-3, kGreen+3, kTeal-7, kTeal, kAzure+10, kAzure-4, kBlue+2, kViolet+8, kViolet-1,kMagenta+1, kMagenta-4, kPink+7, kPink-4};

  std::string label[10] = {"Truncation", "Iterations", "Reweighting Prior", "Tracking Efficiency", "Quark/Gluon"};
  std::string sys[10] = {"trunc", "iter", "rw", "eff", "gluons"};
  std::stringstream ss;
  for (int i = 0; i  < 4; i++){
      TH1D* h1;
      TH1D* hratio;
      if ((sys[i] != "trunc") && (sys[i] != "iter"))
	{
	  ss << "h1_" << sys[i] << "_copy";
	  cout <<ss.str()<<endl;
	  h1 = (TH1D*)_file0->Get(ss.str().c_str());
	  ss << "_ratiio";
	  hratio =(TH1D*)h1->Clone(ss.str().c_str());
	  ss.str("");
	  hratio->Add(hnom, -1);
	  for (int i = 1; i <= hratio->GetNbinsX(); i++){
	      double cont = hratio->GetBinContent(i);
	      hratio->SetBinContent(i, std::sqrt(cont*cont));
	  }
	}

      else{

	std::string ty1;
	std::string ty2;
	if (sys[i] == "trunc") { ty1 = "plus5"; ty2 = "min5";}
	else {ty1 = "iterhigh"; ty2 = "iterlow";}
	ss << "h1_" << ty1 << "_copy";
	cout << ss.str() <<endl;
	h1 = (TH1D*)_file0->Get(ss.str().c_str());
	ss << "_ratiio";
	hratio =(TH1D*)h1->Clone(ss.str().c_str());
	ss.str("");
	hratio->Add(hnom, -1);
	ss <<"h1_" << ty2 << "_copy";
	cout << ss.str() <<endl;

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
      hratio->SetLineColor(col[i*2]);
      hratio->SetLineWidth(4);
      hratio->Draw("same Hist");
      leg->AddEntry(hratio, label[i].c_str(), "l");
    }
  leg->Draw("same");
  hnom_err->Draw("hist same");

  lat->DrawLatex(0.88, 0.83, "ALICE Pb--Pb 5.02 TeV, 30-50%");
  lat->DrawLatex(0.88, 0.78, "Charged jets, #it{R} = 0.2,|#eta_{jet}| < 0.7");
  lat->DrawLatex(0.88, 0.73, "Low q_{2}, out-of-plane/in-plane");






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
