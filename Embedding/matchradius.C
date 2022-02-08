
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;


    int xres = 250;
    int yres = 250;
    int xl = 10.0;
    int xr = 100.0;


  vector<const char*> filenames1 = {"TrainOutputs/old/AnalysisResults7007.root", "TrainOutputs/old/AnalysisResults7008.root", "TrainOutputs/old/AnalysisResults7009.root", "TrainOutputs/old/AnalysisResults7010.root", "TrainOutputs/old/AnalysisResults7011.root", "TrainOutputs/old/AnalysisResults7012.root", "TrainOutputs/old/AnalysisResults7013.root", "TrainOutputs/old/AnalysisResults7014.root", "TrainOutputs/old/AnalysisResults7015.root", "TrainOutputs/old/AnalysisResults7016.root", "TrainOutputs/old/AnalysisResults7017.root", "TrainOutputs/old/AnalysisResults7018.root", "TrainOutputs/old/AnalysisResults7022.root", "TrainOutputs/old/AnalysisResults7019.root", "TrainOutputs/old/AnalysisResults7020.root", "TrainOutputs/old/AnalysisResults7021.root", "TrainOutputs/old/AnalysisResults7023.root", "TrainOutputs/old/AnalysisResults7024.root", "TrainOutputs/old/AnalysisResults7025.root", "TrainOutputs/old/AnalysisResults7026.root"};

  vector<const char*> old = {"TrainOutputs/old/AnalysisResults6853.root"};

  
 
    TH1D *A_wtemp = new TH1D("A_wtemp", "", xres, -0.1, 0.4); 
    TH1D *A_utemp = new TH1D("A_utemp", "", xres, -0.1, 0.4); 
    TH1D *A_whist = new TH1D("A_whist", "", xres, -0.1, 0.4);
    TH1D *A_uhist = new TH1D("A_uhist", "", xres, -0.1, 0.4);



long int totevents = 0;
void getarea(const char *analysisfile, double cl, double cr);
void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3);
void setcolor(TH1D* h, int kcolor);

void matchradius (double centl, double centr) {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  for (int i = 0; i < filenames1.size(); i++)   getarea(filenames1[i], centl, centr);
      TH1D *A_w1 = (TH1D*)A_whist->Clone("A_w1");
      TH1D *A_u1 = (TH1D*)A_uhist->Clone("A_u1");
      A_whist->Reset();
      A_uhist->Reset();  

   A_w1->Scale(1.0/A_w1->Integral());
   A_u1->Scale(1.0/A_u1->Integral());
 

   setcolor(A_u1, kRed);
   setcolor(A_w1, kBlue);
    
   stringstream ss;
   stringstream ssw;
   stringstream ssu;
   ss << "R = 0.4, pthbcut 0, tru10, " << Form("%.0f-%.0f", centl, centr) << "%\n";
   ssw << "weighted, " << Form("%.2f", 100*A_w1->Integral(A_w1->GetXaxis()->FindBin(0.25), A_w1->GetXaxis()->FindBin(0.35))/A_w1->Integral()) << "%" << " #geq 0.25\n"; 
   ssu << "unweighted, " << Form("%.2f", 100*A_u1->Integral(A_u1->GetXaxis()->FindBin(0.25), A_u1->GetXaxis()->FindBin(0.35))/A_u1->Integral()) << "%" << " #geq 0.25\n"; 

   auto *leg = new TLegend(0.4, 0.6, 0.85, 0.85);
   leg->SetTextSize(0.03);
   leg->SetBorderSize(0);
   leg->AddEntry((TObject*)0, ss.str().c_str(), "");
   leg->AddEntry(A_u1, ssu.str().c_str(), "pl");
   leg->AddEntry(A_w1, ssw.str().c_str(), "pl");

   gStyle->SetOptStat(0);

   cout << "Weighted Integral: " << A_w1->Integral() << "\n";
   cout << "#geq 0.25: " << A_w1->Integral(A_w1->GetXaxis()->FindBin(0.25), A_w1->GetXaxis()->FindBin(0.35)) << "\n";


 TCanvas *c1 = new TCanvas("c1", "Match Radius", 800, 500);
   c1->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.05, 1.0, 1.0);
   pad1->Draw();
   pad1->cd();
   //pad1->SetLogy();
   //A_old->Draw("hist same");
   A_u1->Draw("hist same");
   A_w1->Draw("hist same");
   leg->Draw("same");

 }
 


 
void getarea(const char *analysisfile, double cl, double cr) {

    A_wtemp->Reset();
    A_utemp->Reset();

    string rootfilename = analysisfile;
    const char *lerootfile = rootfilename.c_str(); 
    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);

    AliEmcalList *EmbeddingHelper = (AliEmcalList*)lefile->Get("AliAnalysisTaskEmcalEmbeddingHelper_histos");    
    TProfile *pThardscale = (TProfile*)EmbeddingHelper->FindObject("fHistXsection");
    TH1D *pThardbin = (TH1D*)EmbeddingHelper->FindObject("fHistTrials");
    TH1D *eventcount = (TH1D*)EmbeddingHelper->FindObject("fHistEventCount");
    //Determine pT hard scaling
    double xsect = pThardscale->GetBinContent(pThardbin->GetMean() + 1);
    double xscale = xsect*pThardscale->GetEntries();
    double trials = pThardbin->GetBinContent(pThardbin->GetMean() + 1);
 
    double n_event = eventcount->GetBinContent(1);
    totevents += n_event;
    double pThardFrac = xscale/(n_event*trials);

    cout << "File: " << lerootfile <<"\n"; 

    float hybpT, detpT, trupT;
    float centrality, matchradius, area;

    TTree* T = nullptr;
    lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_hybJet",T);
    T->SetBranchAddress("Jet_Pt", &hybpT);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &detpT);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &trupT);
    T->SetBranchAddress("Event_Centrality", &centrality);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Distance", &matchradius);
    T->SetBranchAddress("Jet_Area", &area);

    int entries	= T->GetEntries();
    int passnum = 0;
    //cout << "Entries: " << entries <<"\n";
    
    const int extractBins = 8;
    double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05, 0.05};
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0};

    for (int i = 0; i < entries; i++) {
        T->GetEntry(i);
           int scaleindex = -1;
           if (hybpT >= extractEdge[0] && hybpT < extractEdge[1])  scaleindex = 0;
           if (hybpT >= extractEdge[1] && hybpT < extractEdge[2])  scaleindex = 1;
           if (hybpT >= extractEdge[2] && hybpT < extractEdge[3])  scaleindex = 2;
           if (hybpT >= extractEdge[3] && hybpT < extractEdge[4])  scaleindex = 3;
           if (hybpT >= extractEdge[4] && hybpT < extractEdge[5])  scaleindex = 4;
           if (hybpT >= extractEdge[5] && hybpT < extractEdge[6])  scaleindex = 5;
           if (hybpT >= extractEdge[6] && hybpT < extractEdge[7])  scaleindex = 6;
           if (hybpT >= extractEdge[7] && hybpT < extractEdge[8])  scaleindex = 7;
        if (scaleindex == 0)   continue;
        if (hybpT < 10.0 || trupT < 10.0 || detpT < 10.0)   continue;
        if (centrality < cl || centrality > cr)  continue;
          passnum++;
 
        A_wtemp->Fill(matchradius, pThardFrac/extractFrac[scaleindex]);
        A_utemp->Fill(matchradius);
   }

  A_whist->Add(A_wtemp);
  A_uhist->Add(A_utemp);

  lefile->Close();

}


void setcolor(TH1D* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}



