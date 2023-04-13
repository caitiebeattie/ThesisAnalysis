using namespace std;

void setColor(TH1F* h, int kcolor, int kMarkerShape = 20);
void setErrorColor(TGraphErrors* h, int kLineColor, int kFillColor, double trans, int kMarkerShape = 20, int kLineWidth = 1);
void setAsymmErrorColor(TGraphAsymmErrors* h, int kLineColor, int kFillColor, double trans, int kMarkerShape = 20, int kLineWidth = 1);
int centl = 30;
int centr = 50;
int R = 2;

void plotESE() {


  //set custom colors
  int kLilac, kLightPink, kOcean, kFuschia, kCaitiePink, kCaitieLightPink, kCaitieBlue, kCaitieDarkBlue, kSunPurple, kSunBlue, kSunOrange, kSunPink;
      kOcean = TColor::GetColor("#240e8f");   
      kFuschia = TColor::GetColor("#d62443");
      kLilac = TColor::GetColor("#3918d9");
      kLightPink = TColor::GetColor("#de3e86");
        kCaitiePink = TColor::GetColor("#CA267A");
        kCaitieLightPink = TColor::GetColor("#F42ECF");
      kCaitieBlue = TColor::GetColor("#00c6c9");
        kCaitieDarkBlue = TColor::GetColor("#00576A");
      kSunPurple = TColor::GetColor("#2100a3");
      kSunPink = TColor::GetColor("#b8239c");
      kSunBlue = TColor::GetColor("#006a80"); 
      kSunOrange = TColor::GetColor("#e86400");

  //create spacer histogram
  TH1D *spacer1 = new TH1D("spacer1", "spacer1", 10, 30, 120);
  TH1D *spacer2 = new TH1D("spacer2", "spacer2", 10, 30, 120);
  TH1D *spacer3 = new TH1D("spacer3", "spacer3", 10, 30, 120);




 //--------------------------------------------------------------------------------------------
 //------------------------ Get Results from Root Files ---------------------------------------
 //--------------------------------------------------------------------------------------------

  //q2-large/q2-small plots
  TFile *file1 = new TFile("ESEplots_R04_Dec1.root");
  TH1D *yieldratio1 = (TH1D*)file1->Get("yieldratio"); 
  TGraphErrors *syst1 = (TGraphErrors*)file1->Get("systematics"); 
    yieldratio1->SetLineColor(kFuschia);
    yieldratio1->SetMarkerColor(kFuschia);
    yieldratio1->SetMarkerSize(0.8);
    syst1->SetLineColor(kFuschia);
    syst1->SetMarkerColor(kFuschia);
    syst1->SetMarkerStyle(21);
    syst1->SetFillColorAlpha(kFuschia, 0.2);

  TFile *file2 = new TFile("ESEplots_R02_Dec1.root");
  TH1D *yieldratio2 = (TH1D*)file2->Get("yieldratio"); 
  TGraphErrors *syst2 = (TGraphErrors*)file2->Get("systematics"); 
    yieldratio2->SetLineColor(kOcean);
    yieldratio2->SetMarkerColor(kOcean);
    yieldratio2->SetMarkerSize(0.8);
    syst2->SetLineColor(kOcean);
    syst2->SetMarkerColor(kOcean);
    syst2->SetFillColorAlpha(kLilac, 0.2);


  //out-of-plane/in-plane plots
  TFile *file3 = new TFile("ESEEPplots_R02_Dec1.root");
  TH1F *ratioLo1 = (TH1F*)file3->Get("ratioOILow"); 
  TH1F *ratioHi1 = (TH1F*)file3->Get("ratioOIHigh"); 
  TGraphErrors *systLo1 = (TGraphErrors*)file3->Get("systematicsLo"); 
  TGraphErrors *systHi1 = (TGraphErrors*)file3->Get("systematicsHi"); 
    setColor(ratioLo1, kCaitiePink);
    setColor(ratioHi1, kCaitieDarkBlue, 21);    
    setErrorColor(systLo1, kCaitiePink, kCaitieLightPink, 0.2);
    setErrorColor(systHi1, kCaitieDarkBlue, kCaitieBlue, 0.3, 21);


  TFile *file4 = new TFile("ESEEPplots_R04_Dec1.root");
  TH1F *ratioLo2 = (TH1F*)file4->Get("ratioOILow"); 
  TH1F *ratioHi2 = (TH1F*)file4->Get("ratioOIHigh"); 
  TGraphErrors *systLo2 = (TGraphErrors*)file4->Get("systematicsLo"); 
  TGraphErrors *systHi2 = (TGraphErrors*)file4->Get("systematicsHi"); 
    setColor(ratioLo2, kCaitiePink);
    setColor(ratioHi2, kCaitieDarkBlue, 21);    
    setErrorColor(systLo2, kCaitiePink, kCaitieLightPink, 0.2);
    setErrorColor(systHi2, kCaitieDarkBlue, kCaitieBlue, 0.3, 21);


  //spectra plots
  TFile *file5 = new TFile("spectraplots_R02_Dec1.root");
  TH1F* spectraR02_loin = (TH1F*)file5->Get("yieldLowInPlane");
  TH1F* spectraR02_loout = (TH1F*)file5->Get("yieldLowOutPlane");
  TH1F* spectraR02_hiin = (TH1F*)file5->Get("yieldHighInPlane");
  TH1F* spectraR02_hiout = (TH1F*)file5->Get("yieldHighOutPlane");
  TGraphAsymmErrors *systR02_loin = (TGraphAsymmErrors*)file5->Get("sysLowInPlane");
  TGraphAsymmErrors *systR02_loout = (TGraphAsymmErrors*)file5->Get("sysLowOutPlane");
  TGraphAsymmErrors *systR02_hiin = (TGraphAsymmErrors*)file5->Get("sysHighInPlane");
  TGraphAsymmErrors *systR02_hiout = (TGraphAsymmErrors*)file5->Get("sysHighOutPlane");
    setColor(spectraR02_loin, kSunPurple);
    setColor(spectraR02_loout, kSunPink);
    setColor(spectraR02_hiin, kSunBlue, 21);
    setColor(spectraR02_hiout, kSunOrange, 21);
    setAsymmErrorColor(systR02_loin, kSunPurple, kSunPurple, 0.0, 20, 2);
    setAsymmErrorColor(systR02_loout, kSunPink, kSunPink, 0.0, 20, 2);
    setAsymmErrorColor(systR02_hiout, kSunOrange, kSunOrange, 0.4);
    setAsymmErrorColor(systR02_hiin, kCaitieDarkBlue, kCaitieBlue, 0.3);

  TFile *file6 = new TFile("spectraplots_R04_Dec1.root");
  TH1F* spectraR04_loin = (TH1F*)file6->Get("yieldLowInPlane");
  TH1F* spectraR04_loout = (TH1F*)file6->Get("yieldLowOutPlane");
  TH1F* spectraR04_hiin = (TH1F*)file6->Get("yieldHighInPlane");
  TH1F* spectraR04_hiout = (TH1F*)file6->Get("yieldHighOutPlane");
  TGraphAsymmErrors *systR04_loin = (TGraphAsymmErrors*)file6->Get("sysLowInPlane");
  TGraphAsymmErrors *systR04_loout = (TGraphAsymmErrors*)file6->Get("sysLowOutPlane");
  TGraphAsymmErrors *systR04_hiin = (TGraphAsymmErrors*)file6->Get("sysHighInPlane");
  TGraphAsymmErrors *systR04_hiout = (TGraphAsymmErrors*)file6->Get("sysHighOutPlane");
    setColor(spectraR04_loin, kSunPurple);
    setColor(spectraR04_loout, kSunPink);
    setColor(spectraR04_hiin, kSunBlue, 21);
    setColor(spectraR04_hiout, kSunOrange, 21);
    setAsymmErrorColor(systR04_loin, kSunPurple, kSunPurple, 0.0, 20, 2);
    setAsymmErrorColor(systR04_loout, kSunPink, kSunPink, 0.0, 20, 2);
    setAsymmErrorColor(systR04_hiout, kSunOrange, kSunOrange, 0.4);
    setAsymmErrorColor(systR04_hiin, kCaitieDarkBlue, kCaitieBlue, 0.3);



 //--------------------------------------------------------------------------------------------
 //------------------------ Generate Legends --------------------------------------------------
 //--------------------------------------------------------------------------------------------

  auto gen = new TLegend(0.2, 0.725, 0.3, 0.9);
        gen->SetTextSize(0.04);
        gen->SetBorderSize(0);
        gen->SetFillColorAlpha(kWhite, 0.0);
        gen->AddEntry((TObject*)0, "ALICE, 30#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
        gen->AddEntry((TObject*)0, "Charged-particle jets, anti-#it{k}_{T}", "");
        gen->AddEntry((TObject*)0, "#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.9-#it{R}", "");
   auto gen2 = new TLegend(0.175, 0.65, 0.3, 0.9);
        gen2->SetTextSize(0.045);
        gen2->SetBorderSize(0);
        gen2->SetFillColorAlpha(kWhite, 0.0);
        gen2->AddEntry((TObject*)0, "ALICE", "");
        gen2->AddEntry((TObject*)0, Form("%.0d#font[122]{-}%.0d%% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", centl, centr), "");
        gen2->AddEntry((TObject*)0, Form("Charged-particle jets, anti-#it{k}_{T}, #it{R} = 0.%.0d", 2), "");
        gen2->AddEntry((TObject*)0, Form("#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.%.0d", 7), "");
   auto gen4 = new TLegend(0.175, 0.65, 0.3, 0.9);
        gen4->SetTextSize(0.045);
        gen4->SetBorderSize(0);
        gen4->SetFillColorAlpha(kWhite, 0.0);
        gen4->AddEntry((TObject*)0, "ALICE", "");
        gen4->AddEntry((TObject*)0, Form("%.0d#font[122]{-}%.0d%% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", centl, centr), "");
        gen4->AddEntry((TObject*)0, Form("Charged-particle jets, anti-#it{k}_{T}, #it{R} = 0.%.0d", 4), "");
        gen4->AddEntry((TObject*)0, Form("#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.%.0d", 5), "");
  auto key = new TLegend(0.55, 0.25, 0.8, 0.4);
       key->SetTextSize(0.045);
       key->SetBorderSize(0);
       key->SetFillColorAlpha(kWhite, 0.0);
       key->AddEntry(syst2, "#it{R} = 0.2", "plf");
       key->AddEntry(syst1, "#it{R} = 0.4", "plf");
   auto datadetails = new TLegend(0.725, 0.725, 0.9, 0.875);
        datadetails->SetTextSize(0.04);
        datadetails->SetBorderSize(0);
        datadetails->AddEntry(systHi1, "#it{q}_{2}-large", "plf");
        datadetails->AddEntry(systLo1, "#it{q}_{2}-small", "plf");
   auto yieldLeg = new TLegend(0.25, 0.15, 0.45, 0.35);
        yieldLeg->SetTextSize(0.03);
        yieldLeg->SetBorderSize(0);
        yieldLeg->AddEntry(spectraR02_hiin, "#it{q}_{2}-large, in-plane", "plf");
        yieldLeg->AddEntry(spectraR02_hiout, "#it{q}_{2}-large, out-of-plane", "plf");
        yieldLeg->AddEntry(spectraR02_loin, "#it{q}_{2}-small, in-plane", "plf");
        yieldLeg->AddEntry(spectraR02_loout, "#it{q}_{2}-small, out-of-plane", "plf");
   auto gen5 = new TLegend(0.3, 0.7, 0.8, 0.9);
        gen5->SetTextSize(0.0315);
        gen5->SetBorderSize(0);
        gen5->SetFillColorAlpha(kWhite, 0.0);
        gen5->AddEntry((TObject*)0, "ALICE", "");
        gen5->AddEntry((TObject*)0, "30#font[122]{-}50\% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
        gen5->AddEntry((TObject*)0, "Charged-particle jets, anti-#it{k}_{T}, #it{R} = 0.2", "");
        gen5->AddEntry((TObject*)0, "#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.7", "");
   auto gen6 = new TLegend(0.3, 0.7, 0.8, 0.9);
        gen6->SetTextSize(0.0315);
        gen6->SetBorderSize(0);
        gen6->SetFillColorAlpha(kWhite, 0.0);
        gen6->AddEntry((TObject*)0, "ALICE", "");
        gen6->AddEntry((TObject*)0, "30#font[122]{-}50\% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
        gen6->AddEntry((TObject*)0, "Charged-particle jets, anti-#it{k}_{T}, #it{R} = 0.4", "");
        gen6->AddEntry((TObject*)0, "#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.5", "");




 //--------------------------------------------------------------------------------------------
 //------------------------ Draw Results ------------------------------------------------------
 //--------------------------------------------------------------------------------------------



  TLine *line = new TLine(30.0, 1, 120.0, 1);
  line->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("Plot1", "Plot1", 800, 600);
  c1->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
          pad1->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad1->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad1->SetRightMargin(0.05); // Upper and lower plot are joined
          pad1->SetTopMargin(0.05); // Upper and lower plot are joined
          pad1->Draw();
          pad1->cd();               // pad1 becomes the current pad   
          spacer1->SetTitle("");
          spacer1->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
          spacer1->SetYTitle("#frac{d^{2}#it{N}}{d#it{p}_{T, ch jet}d#it{#eta}_{jet}}#cbar_{#it{q}_{2}-large} / #frac{d^{2}#it{N}}{d#it{p}_{T, ch jet}d#it{#eta}_{jet}}#cbar_{#it{q}_{2}-small}");
          spacer1->GetXaxis()->SetTitleSize(0.05);
          spacer1->GetXaxis()->SetTitleOffset(1.2);
          spacer1->GetYaxis()->SetTitleSize(0.05);
          spacer1->GetYaxis()->SetTitleOffset(1.2);
          spacer1->GetXaxis()->SetLabelSize(0.04);
          spacer1->GetYaxis()->SetLabelSize(0.04);
          spacer1->SetMinimum(0.51);
          spacer1->SetMaximum(1.49);
          spacer1->Draw("same");
       syst1->Draw("ezp 5 same");
       yieldratio1->Draw("same");
       syst2->Draw("ezp 5 same");
       yieldratio2->Draw("same");
       line->Draw("same");
       gen->Draw("same");
       key->Draw("same");

     /*
       cout << "\nR = 0.2 Values\n";
       for (int i = 1; i <= yieldratio2->GetNbinsX(); i++) {
           cout << "Bin " << yieldratio2->GetBinLowEdge(i) << "-" << yieldratio2->GetBinLowEdge(i+1) << ": " << yieldratio2->GetBinContent(i) <<"\n";}
       cout << "R = 0.2 Stat Errors\n";
       for (int i = 1; i <= yieldratio2->GetNbinsX(); i++) {
           cout << "Bin " << yieldratio2->GetBinLowEdge(i) << "-" << yieldratio2->GetBinLowEdge(i+1) << ": " << yieldratio2->GetBinError(i) <<"\n";
           cout << "	" << 100*yieldratio2->GetBinError(i)/yieldratio2->GetBinContent(i) <<"\n";}
       cout << "R = 0.4\n";
       for (int i = 1; i <= yieldratio1->GetNbinsX(); i++) {
           cout << "Bin " << yieldratio1->GetBinLowEdge(i) << "-" << yieldratio1->GetBinLowEdge(i+1) << ": " << yieldratio1->GetBinContent(i) <<"\n";}
       cout << "R = 0.4 Stat Errors\n";
       for (int i = 1; i <= yieldratio1->GetNbinsX(); i++) {
           cout << "Bin " << yieldratio1->GetBinLowEdge(i) << "-" << yieldratio1->GetBinLowEdge(i+1) << ": " << yieldratio1->GetBinError(i) <<"\n";
           cout << "	" << 100*yieldratio1->GetBinError(i)/yieldratio1->GetBinContent(i) <<"\n";}
     */
     

  TCanvas *c2 = new TCanvas("Plot2", "Plot2", 800, 600);
  c2->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 1.0, 1.0);
          pad2->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad2->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad2->SetRightMargin(0.05); // Upper and lower plot are joined
          pad2->SetTopMargin(0.05); // Upper and lower plot are joined
          pad2->Draw();
          pad2->cd();               // pad1 becomes the current pad   
         //pad1->SetLogx();
         //pad1->SetLogy();
             spacer2->SetTitle("");
             spacer2->SetMinimum(0.01);
             spacer2->SetMaximum(1.99);
             spacer2->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
             spacer2->GetXaxis()->SetTitleSize(0.05);
             spacer2->GetXaxis()->SetLabelSize(0.045);
             spacer2->GetXaxis()->SetTitleOffset(1.3);
             spacer2->GetYaxis()->SetTitleSize(0.05);
             spacer2->GetYaxis()->SetLabelSize(0.045);
             spacer2->GetYaxis()->SetTitleOffset(1.2);
             spacer2->GetXaxis()->SetRangeUser(35.0, 120.0);
             spacer2->SetYTitle("#frac{d^{2}#it{N}}{d#it{p}_{T, ch jet}d#it{#eta}_{jet}}#cbar_{out-of-plane} / #frac{d^{2}#it{N}}{d#it{p}_{T, ch jet}d#it{#eta}_{jet}}#cbar_{in-plane}");
         spacer2->Draw("same");
         //ratioOILow->GetXaxis()->SetRangeUser(ratioOILow->GetBinLowEdge(startbin), 100.0);
         //ratioHi1->GetXaxis()->SetRangeUser(ratioHi1->GetBinLowEdge(startbin), 120.0);
         ratioLo1->Draw("same");
         ratioHi1->Draw("same");
         systLo1->GetXaxis()->SetRangeUser(35.0, 120.0);
         systLo1->Draw("ezp 5 same");
         systHi1->Draw("ezp 5 same");
         ratioLo1->Draw("same");
         ratioHi1->Draw("same");
     //caitie->Draw("ezp 2 same");
   gen2->Draw("same");
   datadetails->Draw("same");
   line->Draw("same");

  TCanvas *c3 = new TCanvas("Plot3", "Plot3", 800, 600);
  c3->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad3 = new TPad("pad3", "", 0.0, 0.05, 1.0, 1.0);
          pad3->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad3->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad3->SetRightMargin(0.05); // Upper and lower plot are joined
          pad3->SetTopMargin(0.05); // Upper and lower plot are joined
          pad3->Draw();
          pad3->cd();               // pad1 becomes the current pad   
         //pad1->SetLogx();
         //pad1->SetLogy();
         spacer2->Draw("same");
         //ratioOILow->GetXaxis()->SetRangeUser(ratioOILow->GetBinLowEdge(startbin), 100.0);
         //ratioHi1->GetXaxis()->SetRangeUser(ratioHi1->GetBinLowEdge(startbin), 120.0);
         ratioLo2->Draw("same");
         ratioHi2->Draw("same");
         systLo2->GetXaxis()->SetRangeUser(35.0, 120.0);
         systLo2->Draw("ezp 5 same");
         systHi2->Draw("ezp 5 same");
         ratioLo2->Draw("same");
         ratioHi2->Draw("same");
     //caitie->Draw("ezp 2 same");
   gen4->Draw("same");
   datadetails->Draw("same");
   line->Draw("same");


  TCanvas *c4 = new TCanvas("Plot4", "Plot4", 600, 700);
  c4->cd();
  TPad *pad4 = new TPad("pad4", "", 0.0, 0.05, 1.0, 1.0);
        pad4->SetBottomMargin(0.1); // Upper and lower plot are joined
        pad4->SetLeftMargin(0.19); // Upper and lower plot are joined
        pad4->SetRightMargin(0.05); // Upper and lower plot are joined
        pad4->SetTopMargin(0.05); // Upper and lower plot are joined
        pad4->Draw();
        pad4->cd();
        pad4->SetLogy();
       spacer3->GetXaxis()->SetTitleSize(0.04);
       spacer3->GetXaxis()->SetRangeUser(35.0, 120.0);
       spacer3->SetTitle("");
       spacer3->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
       spacer3->GetXaxis()->SetTitleOffset(1.2);
       spacer3->SetYTitle("#frac{1}{#it{N}_{event}} #frac{d^{2}#it{N}}{d#it{p}_{T, ch jet}d#it{#eta}_{jet}}");
       spacer3->GetYaxis()->SetTitleOffset(2.1);
       spacer3->SetMaximum(2*pow(10,-4));
       spacer3->SetMinimum(8*pow(10,-8));
   spacer3->Draw("same");
   systR02_hiin->Draw("ezp 2 same");  
   systR02_hiout->Draw("ezp 2 same");  
   systR02_loin->Draw("ezp 5 same");  
   systR02_loout->Draw("ezp 5 same");  
   spectraR02_hiin->Draw("same ex0");
   spectraR02_hiout->Draw("same ex0");
   spectraR02_loin->Draw("same ex0");
   spectraR02_loout->Draw("same ex0");
   yieldLeg->Draw("same");
   gen5->Draw("same");
  
  TCanvas *c5 = new TCanvas("Plot5", "Plot5", 600, 700);
  c5->cd();
  TPad *pad5 = new TPad("pad5", "", 0.0, 0.05, 1.0, 1.0);
        pad5->SetBottomMargin(0.1); // Upper and lower plot are joined
        pad5->SetLeftMargin(0.19); // Upper and lower plot are joined
        pad5->SetRightMargin(0.05); // Upper and lower plot are joined
        pad5->SetTopMargin(0.05); // Upper and lower plot are joined
        pad5->Draw();
        pad5->cd();
        pad5->SetLogy();
       spacer3->GetXaxis()->SetTitleSize(0.04);
       spacer3->GetXaxis()->SetRangeUser(35.0, 120.0);
       spacer3->SetTitle("");
       spacer3->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
       spacer3->GetXaxis()->SetTitleOffset(1.2);
       spacer3->SetYTitle("#frac{1}{#it{N}_{event}} #frac{d^{2}#it{N}}{d#it{p}_{T, ch jet}d#it{#eta}_{jet}}");
       spacer3->GetYaxis()->SetTitleOffset(2.1);
       spacer3->SetMaximum(2*pow(10,-4));
       spacer3->SetMinimum(8*pow(10,-8));
   spacer3->Draw("same");
   systR04_hiin->Draw("ezp 2 same");  
   systR04_hiout->Draw("ezp 2 same");  
   systR04_loin->Draw("ezp 5 same");  
   systR04_loout->Draw("ezp 5 same");  
   spectraR04_hiin->Draw("same ex0");
   spectraR04_hiout->Draw("same ex0");
   spectraR04_loin->Draw("same ex0");
   spectraR04_loout->Draw("same ex0");
   yieldLeg->Draw("same");
   gen6->Draw("same");
}









//------------------------------------------------------------------------------------------------------------




void setColor(TH1F* h, int kcolor, int kMarkerShape = 20) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(2);
  h->SetMarkerStyle(kMarkerShape);
  h->SetMarkerSize(0.9);
}

void setErrorColor(TGraphErrors* h, int kLineColor, int kFillColor, double trans, int kMarkerShape = 21, int kLineWidth = 1) {
  h->SetFillColorAlpha(kFillColor, trans);
  h->SetMarkerColor(kLineColor);
  h->SetLineColor(kLineColor);
  h->SetLineWidth(kLineWidth);
  h->SetMarkerStyle(kMarkerShape);
  h->SetMarkerSize(1.0);
}

void setAsymmErrorColor(TGraphAsymmErrors* h, int kLineColor, int kFillColor, double trans, int kMarkerShape = 21, int kLineWidth = 1) {
  h->SetFillColorAlpha(kFillColor, trans);
  h->SetMarkerColor(kLineColor);
  h->SetLineColor(kLineColor);
  h->SetLineWidth(kLineWidth);
  h->SetMarkerStyle(kMarkerShape);
  h->SetMarkerSize(1.0);
}
