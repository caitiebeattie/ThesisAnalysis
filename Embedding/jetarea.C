
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;


    int xres = 100;
    int yres = 100;
    int xl = 10.0;
    int xr = 100.0;

  vector<const char*> proper = {"TrainOutputs/AnalysisResults7007.root", "TrainOutputs/AnalysisResults7008.root", "TrainOutputs/AnalysisResults7009.root", "TrainOutputs/AnalysisResults7010.root", "TrainOutputs/AnalysisResults7011.root", "TrainOutputs/AnalysisResults7012.root", "TrainOutputs/AnalysisResults7013.root", "TrainOutputs/AnalysisResults7014.root", "TrainOutputs/AnalysisResults7015.root", "TrainOutputs/AnalysisResults7016.root", "TrainOutputs/AnalysisResults7017.root", "TrainOutputs/AnalysisResults7018.root", "TrainOutputs/AnalysisResults7022.root", "TrainOutputs/AnalysisResults7019.root", "TrainOutputs/AnalysisResults7020.root", "TrainOutputs/AnalysisResults7021.root", "TrainOutputs/AnalysisResults7023.root", "TrainOutputs/AnalysisResults7024.root", "TrainOutputs/AnalysisResults7025.root", "TrainOutputs/AnalysisResults7026.root"};

  vector<const char*> old = {"TrainOutputs/old/AnalysisResults6853.root"};

  
 
    TH1D *A_temp = new TH1D("A_temp", "", xres, 0.3, 0.8); 
    TH1D *A_hist = new TH1D("A_hist", "", xres, 0.3, 0.8);



long int totevents = 0;
void getarea(const char *analysisfile);
void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3);
void setcolor(TH1D* h, int kcolor);

void jetarea () {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  for (int i = 0; i < proper.size(); i++)   getarea(proper[i]);
      TH1D *A_proper = (TH1D*)A_hist->Clone("A_proper");
      A_hist->Reset();
  for (int i = 0; i < old.size(); i++)      getarea(old[i]);
      TH1D *A_old = (TH1D*)A_hist->Clone("A_old");
  
   A_proper->Scale(1.0/A_proper->Integral());
   A_old->Scale(1.0/A_old->Integral());
 
   double belowcut = A_proper->Integral(A_proper->FindBin(0.3), A_proper->FindBin(0.4));
   double abovecut = A_proper->Integral(A_proper->FindBin(0.4), A_proper->FindBin(0.8));

   cout << "Below cut: " << belowcut << ", Above cut: " << abovecut <<"\n";
  
   setcolor(A_old, kBlue);
   setcolor(A_proper, kRed);

    


 TCanvas *c1 = new TCanvas("c1", "Jet Area", 500, 500);
   c1->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.05, 1.0, 1.0);
   pad1->Draw();
   pad1->cd();
   //pad1->SetLogy();
   A_old->Draw("hist same");
   A_proper->Draw("hist same");


 }
 


 
void getarea(const char *analysisfile) {

    A_temp->Reset();

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
        if (centrality < 30.0 || centrality > 50.0)  continue;
          passnum++;
 
        A_temp->Fill(area, pThardFrac/extractFrac[scaleindex]);
        
   }

  A_hist->Add(A_temp);
 

  lefile->Close();

}


void setcolor(TH1D* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}



