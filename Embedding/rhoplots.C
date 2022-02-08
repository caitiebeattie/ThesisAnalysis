
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;

void getbin(const char *analysisfile);

int lbin = 0;
int rbin = 500;
int binres = 100;
    
TH1D *rho0_10 = new TH1D("rho0_10", "central", binres, lbin, rbin);
TH1D *rho10_30 = new TH1D("rho10_30", "other", binres, lbin, rbin);
TH1D *rho30_50 = new TH1D("rho30_50", "semi-central", binres, lbin, rbin);
TH1D *rho50_90 = new TH1D("rho70_90", "peripheral", binres, lbin, rbin);

long int totevents = 0;
    

  vector<const char*> filenames = {"TrainOutputs/AnalysisResults7007.root", "TrainOutputs/AnalysisResults7008.root", "TrainOutputs/AnalysisResults7009.root", "TrainOutputs/AnalysisResults7010.root", "TrainOutputs/AnalysisResults7011.root", "TrainOutputs/AnalysisResults7012.root", "TrainOutputs/AnalysisResults7013.root", "TrainOutputs/AnalysisResults7014.root", "TrainOutputs/AnalysisResults7015.root", "TrainOutputs/AnalysisResults7016.root", "TrainOutputs/AnalysisResults7017.root", "TrainOutputs/AnalysisResults7018.root", "TrainOutputs/AnalysisResults7022.root", "TrainOutputs/AnalysisResults7019.root", "TrainOutputs/AnalysisResults7020.root", "TrainOutputs/AnalysisResults7021.root", "TrainOutputs/AnalysisResults7023.root", "TrainOutputs/AnalysisResults7024.root", "TrainOutputs/AnalysisResults7025.root", "TrainOutputs/AnalysisResults7026.root"};

  vector<const char*> cents = {"cent1", "cent2", "cent3", "cent4", "cent5", "cent6", "cent7", "cent8", "cent9", "cent10", "cent11", "cent12", "cent13", "cent14", "cent15", "cent16", "cent17", "cent18", "cent19", "cent20"};
  vector<const char*> othrs = {"othr1", "othr2", "othr3", "othr4", "othr5", "othr6", "othr7", "othr8", "othr9", "othr10", "othr11", "othr12", "othr13", "othr14", "othr15", "othr16", "othr17", "othr18", "othr19", "othr20"};
  vector<const char*> semis = {"semi1", "semi2", "semi3", "semi4", "semi5", "semi6", "semi7", "semi8", "semi9", "semi10", "semi11", "semi12", "semi13", "semi14", "semi15", "semi16", "semi17", "semi18", "semi19", "semi20"};
  vector<const char*> peris = {"peri1", "peri2", "peri3", "peri4", "peri5", "peri6", "peri7", "peri8", "peri9", "peri10", "peri11", "peri12", "peri13", "peri14", "peri15", "peri16", "peri17", "peri18", "peri19", "peri20"};


void rhoplots () {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);



vector<TH1D*> cent(0);
vector<TH1D*> othr(0);
vector<TH1D*> semi(0);
vector<TH1D*> peri(0);

  //Fill Vectors with Hists
  for (int i = 0; i < filenames.size(); i++)    { 
      getbin(filenames[i]);
        cent.push_back((TH1D*)rho0_10->Clone(cents[i]));
        othr.push_back((TH1D*)rho10_30->Clone(othrs[i]));
        semi.push_back((TH1D*)rho30_50->Clone(semis[i]));
        peri.push_back((TH1D*)rho50_90->Clone(peris[i]));
      }


  //Get Cumulative Hists
  TH1D *centsum = (TH1D*)cent[0]->Clone("centsum");
  for (int i = 1; i < 20; i++)       centsum->Add(cent[i]);
  TH1D *othrsum = (TH1D*)othr[0]->Clone("othrsum");
  for (int i = 1; i < 20; i++)       othrsum->Add(othr[i]);
  TH1D *semisum = (TH1D*)cent[0]->Clone("semisum");
  for (int i = 1; i < 20; i++)       semisum->Add(semi[i]);
  TH1D *perisum = (TH1D*)cent[0]->Clone("perisum");
  for (int i = 1; i < 20; i++)       perisum->Add(peri[i]);

  cent.push_back(centsum);
  othr.push_back(othrsum);
  semi.push_back(semisum);
  peri.push_back(perisum);
 
  
  cout << "Total events: " << totevents <<"\n";
  double avgevents = (double)totevents/20.0;

  //Scale by Average Events
  for (int i = 0; i <= filenames.size(); i++) {
     cent[i]->Scale(avgevents);
     othr[i]->Scale(avgevents);
     semi[i]->Scale(avgevents);
     peri[i]->Scale(avgevents);
  }

  vector<int> colors = {kRed+2, kRed, kOrange+10, kOrange-4, kYellow, kSpring+10, kSpring, kGreen-3, kGreen+2, kTeal+6, kTeal, kCyan-3, kAzure+7, kBlue+1, kBlue-2, kViolet+2, kViolet, kMagenta-4, kPink+6, kPink+10, kBlack};

  //Set Plot Aesthetics
  for (int i = 0; i < filenames.size(); i++) {
      cent[i]->SetMarkerStyle(20);
        cent[i]->SetMarkerSize(0.5);
        cent[i]->SetLineWidth(2);
        cent[i]->SetLineColor(colors[i]);
        cent[i]->SetMarkerColor(colors[i]);
      othr[i]->SetMarkerStyle(20);
        othr[i]->SetMarkerSize(0.5);
        othr[i]->SetLineWidth(2);
        othr[i]->SetLineColor(colors[i]);
        othr[i]->SetMarkerColor(colors[i]);
      semi[i]->SetMarkerStyle(20);
        semi[i]->SetMarkerSize(0.5);
        semi[i]->SetLineWidth(2);
        semi[i]->SetLineColor(colors[i]);
        semi[i]->SetMarkerColor(colors[i]);
      peri[i]->SetMarkerStyle(20);
        peri[i]->SetMarkerSize(0.5);
        peri[i]->SetLineWidth(2);
        peri[i]->SetLineColor(colors[i]);
        peri[i]->SetMarkerColor(colors[i]);
  }


  centsum->Scale(avgevents);
  othrsum->Scale(avgevents);
  semisum->Scale(avgevents);
  perisum->Scale(avgevents);

  centsum->SetLineColor(kBlack);
  centsum->SetMarkerColor(kBlack);
  centsum->SetMarkerStyle(20);
  centsum->SetMarkerSize(0.5);
  centsum->SetLineWidth(2);
  othrsum->SetLineColor(kBlack);
  othrsum->SetMarkerColor(kBlack);
  othrsum->SetMarkerStyle(20);
  othrsum->SetMarkerSize(0.5);
  othrsum->SetLineWidth(2);
  semisum->SetLineColor(kBlack);
  semisum->SetMarkerColor(kBlack);
  semisum->SetMarkerStyle(20);
  semisum->SetMarkerSize(0.5);
  semisum->SetLineWidth(2);
  perisum->SetLineColor(kBlack);
  perisum->SetMarkerColor(kBlack);
  perisum->SetMarkerStyle(20);
  perisum->SetMarkerSize(0.5);
  perisum->SetLineWidth(2);




/*
//Differential Scaling
y0->Scale(1.0/y0->Integral(), "width");
y5->Scale(1.0/y5->Integral(), "width");
y7->Scale(1.0/y7->Integral(), "width");
*/
  
  double yd = 0.2;
  double yu = 0.75;  
  double x1 = 0.55;
  double x2 = 0.7;
  double x3 = 0.85;
  

  auto datarun1 = new TLegend(x1, yu, x3, 0.85);
    datarun1->SetTextSize(0.035);
    datarun1->SetBorderSize(0);
    datarun1->AddEntry((TObject*)0, "0%-10%, R = 0.4", "");
  auto datarun2 = new TLegend(x1, yu, x3, 0.85);
    datarun2->SetTextSize(0.035);
    datarun2->SetBorderSize(0);
    datarun2->AddEntry((TObject*)0, "10%-30%, R = 0.4", "");
  auto datarun3 = new TLegend(x1, yu, x3, 0.85);
    datarun3->SetTextSize(0.035);
    datarun3->SetBorderSize(0);
    datarun3->AddEntry((TObject*)0, "30%-50%, R = 0.4", "");
  auto datarun4 = new TLegend(x1, yu, x3, 0.85);
    datarun4->SetTextSize(0.035);
    datarun4->SetBorderSize(0);
    datarun4->AddEntry((TObject*)0, "50%-90%, R = 0.4", "");


  //Key: Central Mean
  auto cleg1 = new TLegend(x1, yd, x2, yu);
    cleg1->SetTextSize(0.025);
    cleg1->SetBorderSize(0);
    for (int i = 0; i < 10; i++)   cleg1->AddEntry(cent[i], Form("%d, #rho: %.2f", i, cent[i]->GetMean()), "pl");
    auto cleg2 = new TLegend(x2, yd, x3, yu);
    cleg2->SetTextSize(0.025);
    cleg2->SetBorderSize(0);
    for (int i = 10; i < 20; i++)  cleg2->AddEntry(cent[i], Form("%d, #rho: %.2f", i, cent[i]->GetMean()), "pl");

  //Key: Other Mean
  auto oleg1 = new TLegend(x1, yd, x2, yu);
    oleg1->SetTextSize(0.025);
    oleg1->SetBorderSize(0);
    for (int i = 0; i < 10; i++)   oleg1->AddEntry(othr[i], Form("%d, #rho: %.2f", i, othr[i]->GetMean()), "pl");
    auto oleg2 = new TLegend(x2, yd, x3, yu);
    oleg2->SetTextSize(0.025);
    oleg2->SetBorderSize(0);
    for (int i = 10; i < 20; i++)  oleg2->AddEntry(othr[i], Form("%d, #rho: %.2f", i, othr[i]->GetMean()), "pl");

  //Key: Semi-Central Mean
  auto sleg1 = new TLegend(x1, yd, x2, yu);
    sleg1->SetTextSize(0.025);
    sleg1->SetBorderSize(0);
    for (int i = 0; i < 10; i++)   sleg1->AddEntry(semi[i], Form("%d, #rho: %.2f", i, semi[i]->GetMean()), "pl");
    auto sleg2 = new TLegend(x2, yd, x3, yu);
    sleg2->SetTextSize(0.025);
    sleg2->SetBorderSize(0);
    for (int i = 10; i < 20; i++)  sleg2->AddEntry(semi[i], Form("%d, #rho: %.2f", i, semi[i]->GetMean()), "pl");

  //Key: Pripheral Mean
  auto pleg1 = new TLegend(x1, yd, x2, yu);
    pleg1->SetTextSize(0.025);
    pleg1->SetBorderSize(0);
    for (int i = 0; i < 10; i++)   pleg1->AddEntry(peri[i], Form("%d, #rho: %.2f", i, peri[i]->GetMean()), "pl"); 
    auto pleg2 = new TLegend(x2, yd, x3, yu);
    pleg2->SetTextSize(0.025);
    pleg2->SetBorderSize(0);
    for (int i = 10; i < 20; i++)  pleg2->AddEntry(peri[i], Form("%d, #rho: %.2f", i, peri[i]->GetMean()), "pl");

  //Key: Overall Comp
  auto comp = new TLegend(x1, 0.5, x3, yu);
    comp->SetTextSize(0.035);
    comp->SetBorderSize(0);
    comp->AddEntry(centsum, Form("0%%-10%%, #rho: %.2f", centsum->GetMean()), "pl"); 
    comp->AddEntry(othrsum, Form("10%%-30%%, #rho: %.2f", othrsum->GetMean()), "pl"); 
    comp->AddEntry(semisum, Form("30%%-50%%, #rho: %.2f", semisum->GetMean()), "pl"); 
    comp->AddEntry(perisum, Form("50%%-90%%, #rho: %.2f", perisum->GetMean()), "pl"); 
  


  
  TCanvas *c1 = new TCanvas("c1", "R = 0.4", 1200, 1200);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.5, 0.5, 1.0);
   pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   pad1->SetLogy();
   cent[0]->SetXTitle("#rho");
   cent[0]->SetTitle("");
   cent[0]->SetMinimum(pow(10,-5));
   cent[0]->SetMaximum(pow(10,10));
   cent[0]->Draw("same");
   for (int i = 1; i < 20; i++)  cent[i]->Draw("same");
   cleg1->Draw("same");
   cleg2->Draw("same");
   datarun1->Draw("same");
     c1->cd();
   TPad *pad2 = new TPad("pad2", "pad2", 0.5, 0.5, 1.0, 1.0);
   pad2->SetBottomMargin(0.1);
   pad2->Draw();
   pad2->cd();
   pad2->SetLogy();
   othr[0]->SetXTitle("#rho");
   othr[0]->SetTitle("");
   othr[0]->SetMinimum(pow(10,-5));
   othr[0]->SetMaximum(pow(10,10));
   othr[0]->Draw("same");
   for (int i = 1; i < 20; i++)  othr[i]->Draw("same");
   oleg1->Draw("same");
   oleg2->Draw("same");
   datarun2->Draw("same");
     c1->cd();
   TPad *pad3 = new TPad("pad3", "pad3", 0.0, 0.0, 0.5, 0.5);
   pad3->SetBottomMargin(0.1);
   pad3->Draw();
   pad3->cd();
   pad3->SetLogy();
   semi[0]->SetXTitle("#rho");
   semi[0]->SetTitle("");
   semi[0]->SetMinimum(pow(10,-5));
   semi[0]->SetMaximum(pow(10,10));
   semi[0]->Draw("same");
   for (int i = 1; i < 20; i++)  semi[i]->Draw("same");
   sleg1->Draw("same");
   sleg2->Draw("same");
   datarun3->Draw("same");
     c1->cd();
   TPad *pad4 = new TPad("pad4", "pad4", 0.5, 0.0, 1.0, 0.5);
   pad4->SetBottomMargin(0.1);
   pad4->Draw();
   pad4->cd();
   pad4->SetLogy();
   peri[0]->SetXTitle("#rho");
   peri[0]->SetTitle("");
   peri[0]->SetMinimum(pow(10, -5));
   peri[0]->SetMaximum(pow(10,10));
   peri[0]->Draw("same");
   for (int i = 1; i < 20; i++)  peri[i]->Draw("same");
   pleg1->Draw("same");
   pleg2->Draw("same");
   datarun4->Draw("same");
  
  TH1D *totsum = (TH1D*)centsum->Clone("totsum");
     totsum->Add(othrsum);
     totsum->Add(semisum);
     totsum->Add(perisum);
     totsum->SetLineColor(kBlack);
     totsum->SetLineWidth(2);
  TCanvas *c5 = new TCanvas("c5", "", 500, 500);
   c5->cd();
   gStyle->SetOptStat(0);
   TPad *pad5 = new TPad("pad5", "pad5", 0.0, 0.05, 1.0, 1.0);
   pad5->SetBottomMargin(0.1);
   pad5->SetLeftMargin(0.15);
   pad5->Draw();
   pad5->cd();
   pad5->SetLogy();
     centsum->SetLineColor(kRed);
     centsum->SetMarkerColor(kRed);
     othrsum->SetLineColor(kOrange);
     othrsum->SetMarkerColor(kOrange);
     semisum->SetLineColor(kGreen);
     semisum->SetMarkerColor(kGreen);
     perisum->SetLineColor(kBlue);
     perisum->SetMarkerColor(kBlue);
   centsum->SetXTitle("#rho");
   centsum->SetYTitle("#frac{N_{jet}}{d#rho}");
   centsum->SetTitle("");
   centsum->SetMinimum(pow(10,-5));
   centsum->SetMaximum(pow(10,20));
   centsum->Draw("hist same");
   othrsum->Draw("hist same");
   semisum->Draw("hist same");
   perisum->Draw("hist same");
   comp->Draw("same");

}
  

void getbin(const char *analysisfile) {

    rho0_10->Reset();
    rho10_30->Reset();
    rho30_50->Reset(); 
    rho50_90->Reset(); 


    string rootfilename = analysisfile;
    const char *lerootfile = rootfilename.c_str(); 
    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);
 
    AliEmcalList *EmbeddingHelper = (AliEmcalList*)lefile->Get("AliAnalysisTaskEmcalEmbeddingHelper_histos");    
    TProfile *pThardscale = (TProfile*)EmbeddingHelper->FindObject("fHistXsection");
    TH1D *pThardbin = (TH1D*)EmbeddingHelper->FindObject("fHistTrials");
    TH1D *eventcount = (TH1D*)EmbeddingHelper->FindObject("fHistEventCount");

    TTree *T = nullptr;
    lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_hybJet", T);
    float hybridpT, detpT, partpT, centrality, rho;
    T->SetBranchAddress("Jet_Pt", &hybridpT);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &detpT);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &partpT);
    T->SetBranchAddress("Event_Centrality", &centrality);
    T->SetBranchAddress("Event_BackgroundDensity", &rho);

    int entries = T->GetEntries();

    //Determine pT hard scaling
    double xsect = pThardscale->GetBinContent(pThardbin->GetMean() + 1);
    double xscale = xsect*pThardscale->GetEntries();
    double trials = pThardbin->GetBinContent(pThardbin->GetMean() + 1);
 
    double n_event = eventcount->GetBinContent(1);
    double pThardFrac = xscale/(n_event*trials);


    //Determine Extractor Bin Scaling
    const int extractBins = 8;
    double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05, 0.05};
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0};
    
    int c = 0;
    int o = 0;
    int s = 0;
    int p = 0;

    for (int i = 0; i < entries; i++) {
       T->GetEntry(i);
       int scaleindex = -1;
         if (hybridpT >= extractEdge[0] && hybridpT < extractEdge[1])  scaleindex = 0;
         if (hybridpT >= extractEdge[1] && hybridpT < extractEdge[2])  scaleindex = 1;
         if (hybridpT >= extractEdge[2] && hybridpT < extractEdge[3])  scaleindex = 2;
         if (hybridpT >= extractEdge[3] && hybridpT < extractEdge[4])  scaleindex = 3;
         if (hybridpT >= extractEdge[4] && hybridpT < extractEdge[5])  scaleindex = 4;
         if (hybridpT >= extractEdge[5] && hybridpT < extractEdge[6])  scaleindex = 5;
         if (hybridpT >= extractEdge[6] && hybridpT < extractEdge[7])  scaleindex = 6;
         if (hybridpT >= extractEdge[7] && hybridpT < extractEdge[8])  scaleindex = 7;
       if (scaleindex == 0)   continue;
       if (hybridpT < 10.0 || detpT < 10.0 || partpT < 10.0 )  continue;
          
          if (centrality > 0.0 && centrality < 10.0)   rho0_10->Fill(rho, pThardFrac/extractFrac[scaleindex]);
          if (centrality > 10.0 && centrality < 30.0)  rho10_30->Fill(rho, pThardFrac/extractFrac[scaleindex]);
          if (centrality > 30.0 && centrality < 50.0)  rho30_50->Fill(rho, pThardFrac/extractFrac[scaleindex]);
          if (centrality > 50.0 && centrality < 90.0)  rho50_90->Fill(rho, pThardFrac/extractFrac[scaleindex]);
  }

    cout << "pT hard bin: " << pThardbin->GetMean() <<"\n";
   // cout << "	c: " << c << ", o: " << o << ", s: " << s << ", p: " << p <<"\n";

    totevents += n_event;

    lefile->Close();

}
