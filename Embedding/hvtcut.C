
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;


    int xres = 50;
    int yres = 250;
    int xl = 10.0;
    int xr = 100.0;


  //vector<const char*> filenames1 = {"TrainOutputs/old/AnalysisResults7007.root", "TrainOutputs/old/AnalysisResults7008.root", "TrainOutputs/old/AnalysisResults7009.root", "TrainOutputs/old/AnalysisResults7010.root", "TrainOutputs/old/AnalysisResults7011.root", "TrainOutputs/old/AnalysisResults7012.root", "TrainOutputs/old/AnalysisResults7013.root", "TrainOutputs/old/AnalysisResults7014.root", "TrainOutputs/old/AnalysisResults7015.root", "TrainOutputs/old/AnalysisResults7016.root", "TrainOutputs/old/AnalysisResults7017.root", "TrainOutputs/old/AnalysisResults7018.root", "TrainOutputs/old/AnalysisResults7022.root", "TrainOutputs/old/AnalysisResults7019.root", "TrainOutputs/old/AnalysisResults7020.root", "TrainOutputs/old/AnalysisResults7021.root", "TrainOutputs/old/AnalysisResults7023.root", "TrainOutputs/old/AnalysisResults7024.root", "TrainOutputs/old/AnalysisResults7025.root", "TrainOutputs/old/AnalysisResults7026.root"};
//vector<const char*> filenames1 = {"TrainOutputs/AnalysisResults7158.root", "TrainOutputs/AnalysisResults7121.root", "TrainOutputs/AnalysisResults7122.root", "TrainOutputs/AnalysisResults7123.root", "TrainOutputs/AnalysisResults7124.root", "TrainOutputs/AnalysisResults7126.root", "TrainOutputs/AnalysisResults7127.root", "TrainOutputs/AnalysisResults7128.root", "TrainOutputs/AnalysisResults7129.root", "TrainOutputs/AnalysisResults7131.root", "TrainOutputs/AnalysisResults7132.root", "TrainOutputs/AnalysisResults7133.root", "TrainOutputs/AnalysisResults7134.root", "TrainOutputs/AnalysisResults7135.root", "TrainOutputs/AnalysisResults7136.root", "TrainOutputs/AnalysisResults7159.root", "TrainOutputs/AnalysisResults7160.root"};
vector<const char*> filenames1 = {"TrainOutputs/AnalysisResults7208.root", "TrainOutputs/AnalysisResults7209.root", "TrainOutputs/AnalysisResults7210.root", "TrainOutputs/AnalysisResults7211.root", "TrainOutputs/AnalysisResults7212.root", "TrainOutputs/AnalysisResults7213.root", "TrainOutputs/AnalysisResults7214.root", "TrainOutputs/AnalysisResults7215.root", "TrainOutputs/AnalysisResults7216.root", "TrainOutputs/AnalysisResults7217.root", "TrainOutputs/AnalysisResults7218.root", "TrainOutputs/AnalysisResults7219.root", "TrainOutputs/AnalysisResults7220.root", "TrainOutputs/AnalysisResults7221.root", "TrainOutputs/AnalysisResults7222.root", "TrainOutputs/AnalysisResults7223.root", "TrainOutputs/AnalysisResults7224.root"};

  vector<const char*> old = {"TrainOutputs/old/AnalysisResults6853.root"};

  
 
    TH1D *A_wtemp = new TH1D("A_wtemp", "", xres, 0.0, 5.0); 
    TH1D *A_utemp = new TH1D("A_utemp", "", xres, 0.0, 5.0); 
    TH1D *A_etemp = new TH1D("A_etemp", "", xres, 0.0, 5.0); 
    TH1D *A_ptemp = new TH1D("A_ptemp", "", xres, 0.0, 5.0); 
    TH1D *A_whist = new TH1D("A_whist", "", xres, 0.0, 5.0);
    TH1D *A_uhist = new TH1D("A_uhist", "", xres, 0.0, 5.0);
    TH1D *A_ehist = new TH1D("A_ehist", "", xres, 0.0, 5.0);
    TH1D *A_phist = new TH1D("A_phist", "", xres, 0.0, 5.0);



long int totevents = 0;
void getarea(const char *analysisfile, double cl, double cr);
void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3);
void setcolor(TH1D* h, int kcolor);

void hvtcut (double centl, double centr) {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  for (int i = 0; i < filenames1.size(); i++)   getarea(filenames1[i], centl, centr);
      TH1D *A_w1 = (TH1D*)A_whist->Clone("A_w1");
      TH1D *A_u1 = (TH1D*)A_uhist->Clone("A_u1");
      TH1D *A_e1 = (TH1D*)A_ehist->Clone("A_e1");
      TH1D *A_p1 = (TH1D*)A_phist->Clone("A_p1");
      A_whist->Reset();
      A_uhist->Reset();  
      A_ehist->Reset();  
      A_phist->Reset();  

   TH1D *comp = (TH1D*)A_w1->Clone("comp");
   comp->Divide(A_p1);

   A_w1->Scale(1.0/A_w1->Integral());
   A_u1->Scale(1.0/A_u1->Integral());
   A_e1->Scale(1.0/A_e1->Integral());
   A_p1->Scale(1.0/A_p1->Integral());
 

   setcolor(A_u1, kRed);
   setcolor(A_w1, kGreen+1);
   setcolor(A_e1, kBlue);
   setcolor(A_p1, kViolet+1);
    
   stringstream ss;
   stringstream ssw1, ssw2, ssw3;
   stringstream ssu1, ssu2, ssu3;
   stringstream sse1, sse2, sse3;
   stringstream ssp1, ssp2, ssp3;
   ss << "R = 0.4, pthbcut 123, tru10, " << Form("%.0f-%.0f", centl, centr) << "%\n";
   ssw1 << Form("%.2f", 100*A_w1->Integral(A_w1->GetXaxis()->FindBin(2.5), A_w1->GetXaxis()->FindBin(5.0))/A_w1->Integral()) << "%" << " #geq 2.5\n"; 
   ssw2 << Form("%.2f", 100*A_w1->Integral(A_w1->GetXaxis()->FindBin(3.0), A_w1->GetXaxis()->FindBin(5.0))/A_w1->Integral()) << "%" << " #geq 3.0\n"; 
   ssw3 << Form("%.2f", 100*A_w1->Integral(A_w1->GetXaxis()->FindBin(3.5), A_w1->GetXaxis()->FindBin(5.0))/A_w1->Integral()) << "%" << " #geq 3.5\n"; 
   ssu1 << Form("%.2f", 100*A_u1->Integral(A_u1->GetXaxis()->FindBin(2.5), A_u1->GetXaxis()->FindBin(5.0))/A_u1->Integral()) << "%" << " #geq 2.5\n"; 
   ssu2 << Form("%.2f", 100*A_u1->Integral(A_u1->GetXaxis()->FindBin(3.0), A_u1->GetXaxis()->FindBin(5.0))/A_u1->Integral()) << "%" << " #geq 3.0\n"; 
   ssu3 << Form("%.2f", 100*A_u1->Integral(A_u1->GetXaxis()->FindBin(3.5), A_u1->GetXaxis()->FindBin(5.0))/A_u1->Integral()) << "%" << " #geq 3.5\n"; 
   sse1 << Form("%.2f", 100*A_e1->Integral(A_e1->GetXaxis()->FindBin(2.5), A_e1->GetXaxis()->FindBin(5.0))/A_e1->Integral()) << "%" << " #geq 2.5\n"; 
   sse2 << Form("%.2f", 100*A_e1->Integral(A_e1->GetXaxis()->FindBin(3.0), A_e1->GetXaxis()->FindBin(5.0))/A_e1->Integral()) << "%" << " #geq 3.0\n"; 
   sse3 << Form("%.2f", 100*A_e1->Integral(A_e1->GetXaxis()->FindBin(3.5), A_e1->GetXaxis()->FindBin(5.0))/A_e1->Integral()) << "%" << " #geq 3.5\n"; 
   ssp1 << Form("%.2f", 100*A_p1->Integral(A_p1->GetXaxis()->FindBin(2.5), A_p1->GetXaxis()->FindBin(5.0))/A_p1->Integral()) << "%" << " #geq 2.5\n"; 
   ssp2 << Form("%.2f", 100*A_p1->Integral(A_p1->GetXaxis()->FindBin(3.0), A_p1->GetXaxis()->FindBin(5.0))/A_p1->Integral()) << "%" << " #geq 3.0\n"; 
   ssp3 << Form("%.2f", 100*A_p1->Integral(A_p1->GetXaxis()->FindBin(3.5), A_p1->GetXaxis()->FindBin(5.0))/A_p1->Integral()) << "%" << " #geq 3.5\n"; 

   auto *leg = new TLegend(0.4, 0.7, 0.7, 0.85);
   leg->SetTextSize(0.03);
   leg->SetBorderSize(0);
   leg->AddEntry((TObject*)0, ss.str().c_str(), "");

   auto *uleg = new TLegend (0.4, 0.55, 0.6, 0.7);
   uleg->SetTextSize(0.03);
   uleg->SetBorderSize(0);
   uleg->AddEntry(A_u1, "unweighted", "pl");
   uleg->AddEntry((TObject*)0, ssu1.str().c_str(), "");
   uleg->AddEntry((TObject*)0, ssu2.str().c_str(), "");
   uleg->AddEntry((TObject*)0, ssu3.str().c_str(), "");

   auto *wleg = new TLegend (0.6, 0.55, 0.8, 0.7);
   wleg->SetTextSize(0.03);
   wleg->SetBorderSize(0);
   wleg->AddEntry(A_w1, "weighted", "pl");
   wleg->AddEntry((TObject*)0, ssw1.str().c_str(), "");
   wleg->AddEntry((TObject*)0, ssw2.str().c_str(), "");
   wleg->AddEntry((TObject*)0, ssw3.str().c_str(), "");
   
   auto *eleg = new TLegend(0.4, 0.4, 0.6, 0.55);
   eleg->SetTextSize(0.03);
   eleg->SetBorderSize(0);
   eleg->AddEntry(A_e1, "EB weighted", "pl");
   eleg->AddEntry((TObject*)0, sse1.str().c_str(), "");
   eleg->AddEntry((TObject*)0, sse2.str().c_str(), "");
   eleg->AddEntry((TObject*)0, sse3.str().c_str(), "");
 
   auto *pleg = new TLegend(0.6, 0.4, 0.8, 0.55);
   pleg->SetTextSize(0.03);
   pleg->SetBorderSize(0);
   pleg->AddEntry(A_p1, "pThard weighted", "pl");
   pleg->AddEntry((TObject*)0, ssp1.str().c_str(), "");
   pleg->AddEntry((TObject*)0, ssp2.str().c_str(), "");
   pleg->AddEntry((TObject*)0, ssp3.str().c_str(), "");

   gStyle->SetOptStat(0);

   cout << "Weighted Integral: " << A_w1->Integral() << "\n";
   cout << "#geq 0.25: " << A_w1->Integral(A_w1->GetXaxis()->FindBin(0.25), A_w1->GetXaxis()->FindBin(0.35)) << "\n";


 TCanvas *c1 = new TCanvas("c1", "Match Radius", 800, 500);
   c1->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.05, 1.0, 1.0);
   pad1->Draw();
   pad1->cd();
   //pad1->SetLogy();
   A_u1->Draw("hist same");
   A_w1->Draw("hist same");
   A_e1->Draw("hist same");
   A_p1->Draw("hist same");
   leg->Draw("same");
   uleg->Draw("same");
   wleg->Draw("same");
   eleg->Draw("same");
   pleg->Draw("same");

 TCanvas *c2 = new TCanvas("c2", "", 800, 500);
   c2->cd();
   TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.05, 1.0, 1.0);
   pad2->Draw();
   pad2->cd();
   comp->Draw("hist same");


 }
 


 
void getarea(const char *analysisfile, double cl, double cr) {

    A_wtemp->Reset();
    A_utemp->Reset();
    A_etemp->Reset();
    A_ptemp->Reset();

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
    double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.03, 0.03, 0.07, 0.07, 0.05};
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

        double hvt;
        hvt = hybpT/trupT; 

        A_wtemp->Fill(hvt, pThardFrac/extractFrac[scaleindex]);
        A_utemp->Fill(hvt);
        A_etemp->Fill(hvt, 1.0/extractFrac[scaleindex]);
        A_ptemp->Fill(hvt, pThardFrac);
   }

  A_ehist->Add(A_etemp);
  A_whist->Add(A_wtemp);
  A_uhist->Add(A_utemp);
  A_phist->Add(A_ptemp);

  lefile->Close();

}


void setcolor(TH1D* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}



