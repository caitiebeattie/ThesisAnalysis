using namespace std;


void plotESE() {

  //set custom colors
  int kLilac, kPistachio, kMint, kForest, kLightPink, kOcean, kFuschia, kChocolate;
      kPistachio = TColor::GetColor("#caf0c1");
      kMint = TColor::GetColor("#87e4db");
      //kForest = TColor::GetColor("#3b5e33");
      kForest = TColor::GetColor("#167D0D");
      kOcean = TColor::GetColor("#240e8f");   
      kChocolate = TColor::GetColor("#2b1d07");
      //kFuschia = TColor::GetColor("#c23072");
      kFuschia = TColor::GetColor("#d62443");
      kLilac = TColor::GetColor("#3918d9");
      kLightPink = TColor::GetColor("#de3e86");

  //create spacer histogram
  TH1D *spacer = new TH1D("spacer", "spacer", 10, 30, 120);

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

 
  //Legend
  auto gen = new TLegend(0.12, 0.65, 0.3, 0.9);
        gen->SetTextSize(0.04);
        gen->SetBorderSize(0);
        gen->SetFillColorAlpha(kWhite, 0.0);
        gen->AddEntry((TObject*)0, "ALICE, 30#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
        gen->AddEntry((TObject*)0, "Charged-particle jets, anti-#it{k}_{T}", "");
        gen->AddEntry((TObject*)0, "#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.9-R", "");
  auto key = new TLegend(0.6, 0.2, 0.8, 0.4);
       key->SetTextSize(0.04);
       key->SetBorderSize(0);
       key->SetFillColorAlpha(kWhite, 0.0);
       key->AddEntry(syst2, "R = 0.2, #it{q}_{2}^{V0C}", "plf");
       key->AddEntry(syst1, "R = 0.4, #it{q}_{2}^{V0C}", "plf");

  TLine *line = new TLine(30.0, 1, 120.0, 1);
  line->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("QMPlot1", "QMPlot1", 800, 600);
  c1->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
          pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
          pad1->SetLeftMargin(0.1); // Upper and lower plot are joined
          pad1->SetRightMargin(0.05); // Upper and lower plot are joined
          pad1->SetTopMargin(0.05); // Upper and lower plot are joined
          pad1->Draw();
          pad1->cd();               // pad1 becomes the current pad   
          spacer->SetTitle("");
          spacer->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          spacer->SetYTitle("Ratio (large/small #it{q}_{2})");
          spacer->SetMinimum(0.01);
          spacer->SetMaximum(2.39);
          spacer->Draw("same");
       syst1->Draw("ezp 5 same");
       yieldratio1->Draw("same");
       syst2->Draw("ezp 5 same");
       yieldratio2->Draw("same");
       line->Draw("same");
       gen->Draw("same");
       key->Draw("same");

}
