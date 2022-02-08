
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

long int totevents = 0;

using namespace std;
    
    int xres = 140;
    int yres = 140;

  std::vector<double> kBinsTru = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0/*, 190.0, 250.0*/};
  std::vector<double> kBinsRaw = {30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0};   //semi-central
  
    //Hyb <-> Det
    TH2D *HybDet = new TH2D("HybDet", "Background Response", kBinsRaw.size()-1, kBinsRaw.data(), kBinsTru.size()-1, kBinsTru.data()); 
    //Det <-> Truth
    TH2D *DetTru = new TH2D("DetTru", "Detector Response", kBinsRaw.size()-1, kBinsRaw.data(), kBinsTru.size()-1, kBinsTru.data());
    //Hyb <-> Truth
    TH2D *HybTru = new TH2D("HybTru", "Overall Response", kBinsRaw.size()-1, kBinsRaw.data(), kBinsTru.size()-1, kBinsTru.data());
    
    //Hyb <-> Det
    TH2D *HD = new TH2D("HD", "Background Response", kBinsRaw.size()-1, kBinsRaw.data(), kBinsTru.size()-1, kBinsTru.data()); 
    //Det <-> Truth
    TH2D *DT = new TH2D("DT", "Detector Response", kBinsRaw.size()-1, kBinsRaw.data(), kBinsTru.size()-1, kBinsTru.data());
    //Hyb <-> Truth
    TH2D *HT = new TH2D("HT", "Overall Response", kBinsRaw.size()-1, kBinsRaw.data(), kBinsTru.size()-1, kBinsTru.data());

void getbin(const char *analysisfile);


void response () {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

   getbin("TrainOutputs/AnalysisResults7440.root");  //pThard bin 20
   getbin("TrainOutputs/AnalysisResults7441.root");  //pThard bin 1
   getbin("TrainOutputs/AnalysisResults7442.root");  //pThard bin 2
   getbin("TrainOutputs/AnalysisResults7443.root");  //pThard bin 3
   getbin("TrainOutputs/AnalysisResults7444.root");  //pThard bin 4
   getbin("TrainOutputs/AnalysisResults7445.root");  //pThard bin 5
   getbin("TrainOutputs/AnalysisResults7446.root");  //pThard bin 6
   getbin("TrainOutputs/AnalysisResults7447.root");  //pThard bin 7
   getbin("TrainOutputs/AnalysisResults7448.root");  //pThard bin 8
   getbin("TrainOutputs/AnalysisResults7449.root");  //pThard bin 9
   getbin("TrainOutputs/AnalysisResults7450.root");  //pThard bin 10
   getbin("TrainOutputs/AnalysisResults7451.root");  //pThard bin 11
   getbin("TrainOutputs/AnalysisResults7452.root");  //pThard bin 12
   getbin("TrainOutputs/AnalysisResults7453.root");  //pThard bin 13
   getbin("TrainOutputs/AnalysisResults7454.root");  //pThard bin 14
   getbin("TrainOutputs/AnalysisResults7455.root");  //pThard bin 15
   getbin("TrainOutputs/AnalysisResults7456.root");  //pThard bin 16
   getbin("TrainOutputs/AnalysisResults7457.root");  //pThard bin 17
   getbin("TrainOutputs/AnalysisResults7458.root");  //pThard bin 18
   getbin("TrainOutputs/AnalysisResults7459.root");  //pThard bin 19

  
   double avgevents = (double)totevents/20.0;
   HT->Scale(avgevents);
   HD->Scale(avgevents);
   DT->Scale(avgevents);

   
   TH1D* truth = HT->ProjectionX();
   TH1D* raw = HT->ProjectionY();

   truth->SetLineColor(kBlue);
   truth->SetLineWidth(2);
   raw->SetLineColor(kRed);
   raw->SetLineWidth(2);

  TCanvas *c1 = new TCanvas("c1", "", 1500, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.05, 0.33, 1.0);
   pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
   pad1->SetLeftMargin(0.15);
   pad1->SetRightMargin(0.15);
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   pad1->SetLogz();
   HT->SetTitle("");
   HT->SetXTitle("p_{T}^{hybrid}");
   HT->SetYTitle("p_{T}^{truth}");
   HT->Draw("same colz");
     c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad2 = new TPad("pad2", "pad2", 0.33, 0.05, 0.67, 1.0);
   pad2->SetBottomMargin(0.1);
   pad2->SetLeftMargin(0.15);
   pad2->SetRightMargin(0.15);
   pad2->Draw();
   pad2->cd();
   pad2->SetLogz();
   HD->SetTitle("");
   HD->SetXTitle("p_{T}^{hybrid}");
   HD->SetYTitle("p_{T}^{detector}");
   HD->Draw("same colz");
     c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad3 = new TPad("pad3", "pad3", 0.67, 0.05, 1.0, 1.0);
   pad3->SetBottomMargin(0.1);
   pad3->SetLeftMargin(0.15);
   pad3->SetRightMargin(0.15);
   pad3->Draw();
   pad3->cd();
   pad3->SetLogz();
   DT->SetTitle("");
   DT->SetXTitle("p_{T}^{detector}");
   DT->SetYTitle("p_{T}^{truth}");
   DT->Draw("same colz");


  TCanvas *c4 = new TCanvas("c4", "", 500, 500);
   c4->cd();
   TPad *pad4 = new TPad("pad4", "pad4", 0.0, 0.0, 1.0, 1.0);
   pad4->Draw();
   pad4->cd();
   pad4->SetLogy();
   truth->Draw("same");
   raw->Draw("same");

}
  
void getbin (const char* analysisfile) {
 
    HybTru->Reset();
    DetTru->Reset();
    HybDet->Reset();   

    string rootfilename = analysisfile;
    const char *lerootfile = rootfilename.c_str(); 
    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);

    float hybridpT, detpT, partpT;
    float centrality;

    TTree* T = nullptr;
    lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_hybJet",T);
    T->SetBranchAddress("Jet_Pt", &hybridpT);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &detpT);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &partpT);
    T->SetBranchAddress("Event_Centrality", &centrality);

    int entries	= T->GetEntries();
    

    AliEmcalList *EmbeddingHelper = (AliEmcalList*)lefile->Get("AliAnalysisTaskEmcalEmbeddingHelper_histos");    
    TProfile *pThardscale = (TProfile*)EmbeddingHelper->FindObject("fHistXsection");
    TH1D *pThardbin = (TH1D*)EmbeddingHelper->FindObject("fHistTrials");
    TH1D *eventcount = (TH1D*)EmbeddingHelper->FindObject("fHistEventCount");
    
    //Determine pT hard scaling
    double xsect = pThardscale->GetBinContent(pThardbin->GetMean() + 1);
    double xscale = xsect*pThardscale->GetEntries();
    double trials = pThardbin->GetBinContent(pThardbin->GetMean() + 1);
    double pThardFrac = xscale/trials;
 
    double n_event = eventcount->GetBinContent(1);
    totevents += n_event;

    //Determine Extractor Bin Scaling
    const int extractBins = 8;
    double extractFrac[extractBins] = {0.0, 0.03, 0.03, 0.05, 0.05, 0.09, 0.09, 0.10};
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0};



    int lowhyb = 0;
    int lowdet = 0;
    int lowpart = 0;
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
        if (centrality < 30.0 || centrality > 50.0)   continue;
        if (hybridpT < 10.0)  {
           lowhyb++;
           continue;}
      /*  if (detpT < 10.0) {
           lowdet++;
          continue;} */
      if (partpT < 10.0) {
           lowpart++;
           continue;}   
       HybDet->Fill(hybridpT, detpT, pThardFrac/(/*n_event*/extractFrac[scaleindex]));
        DetTru->Fill(detpT, partpT, pThardFrac/(/*n_event*/extractFrac[scaleindex])); 
        HybTru->Fill(hybridpT, partpT, pThardFrac/(/*n_event*/extractFrac[scaleindex]));
    }
   
   HT->Add(HybTru);
   HD->Add(HybDet);
   DT->Add(DetTru);
    
  cout << "pT hard bin: " << pThardbin->GetMean() <<"\n";

  lefile->Close();
}

