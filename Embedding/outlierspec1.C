
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;


void getbin(const char *analysisfile, double centl, double centr);
void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3);

int lbin = 10;
int rbin = 200;
int binres = 95;
    
TH1D *binspec = new TH1D("bin", "", binres, lbin, rbin);
TH1D *genspec = new TH1D("gen", "", binres, lbin, rbin);
TH1D *detspec = new TH1D("det", "", binres, lbin, rbin);

TH1D *binsum1 = new TH1D("bin1sum", "", binres, lbin, rbin);
TH1D *gensum1 = new TH1D("gen1sum", "", binres, lbin, rbin);
TH1D *detsum1 = new TH1D("det1sum", "", binres, lbin, rbin);
TH1D *binsum2 = new TH1D("bin2sum", "", binres, lbin, rbin);
TH1D *gensum2 = new TH1D("gen2sum", "", binres, lbin, rbin);
TH1D *detsum2= new TH1D("det2sum", "", binres, lbin, rbin);

long int totevents = 0;

  vector<const char*> filenames1 = {"TrainOutputs/old/AnalysisResults7007.root", "TrainOutputs/old/AnalysisResults7008.root", "TrainOutputs/old/AnalysisResults7009.root", "TrainOutputs/old/AnalysisResults7010.root", "TrainOutputs/old/AnalysisResults7011.root", "TrainOutputs/old/AnalysisResults7012.root", "TrainOutputs/old/AnalysisResults7013.root", "TrainOutputs/old/AnalysisResults7014.root", "TrainOutputs/old/AnalysisResults7015.root", "TrainOutputs/old/AnalysisResults7016.root", "TrainOutputs/old/AnalysisResults7017.root", "TrainOutputs/old/AnalysisResults7018.root", "TrainOutputs/old/AnalysisResults7022.root", "TrainOutputs/old/AnalysisResults7019.root", "TrainOutputs/old/AnalysisResults7020.root", "TrainOutputs/old/AnalysisResults7021.root", "TrainOutputs/old/AnalysisResults7023.root", "TrainOutputs/old/AnalysisResults7024.root", "TrainOutputs/old/AnalysisResults7025.root", "TrainOutputs/old/AnalysisResults7026.root"};
  vector<const char*> filenames2 = {/*"TrainOutputs/old/AnalysisResults6846.root", "TrainOutputs/old/AnalysisResults6847.root", "TrainOutputs/old/AnalysisResults6848.root", */"TrainOutputs/old/AnalysisResults6849.root", "TrainOutputs/old/AnalysisResults6850.root", "TrainOutputs/old/AnalysisResults6851.root", "TrainOutputs/old/AnalysisResults6852.root", "TrainOutputs/old/AnalysisResults6853.root", "TrainOutputs/old/AnalysisResults6854.root", "TrainOutputs/old/AnalysisResults6855.root", "TrainOutputs/old/AnalysisResults6856.root", "TrainOutputs/old/AnalysisResults6857.root", "TrainOutputs/old/AnalysisResults6858.root", "TrainOutputs/old/AnalysisResults6859.root", "TrainOutputs/old/AnalysisResults6860.root", "TrainOutputs/old/AnalysisResults6861.root", "TrainOutputs/old/AnalysisResults6862.root", "TrainOutputs/old/AnalysisResults6863.root", "TrainOutputs/old/AnalysisResults6864.root", "TrainOutputs/old/AnalysisResults6865.root"};
vector<const char*> filenames4 = {"TrainOutputs/AnalysisResults7158.root", "TrainOutputs/AnalysisResults7121.root", "TrainOutputs/AnalysisResults7122.root", "TrainOutputs/AnalysisResults7123.root", "TrainOutputs/AnalysisResults7124.root", "TrainOutputs/AnalysisResults7126.root", "TrainOutputs/AnalysisResults7127.root", "TrainOutputs/AnalysisResults7128.root", "TrainOutputs/AnalysisResults7129.root", "TrainOutputs/AnalysisResults7131.root", "TrainOutputs/AnalysisResults7132.root", "TrainOutputs/AnalysisResults7133.root", "TrainOutputs/AnalysisResults7134.root", "TrainOutputs/AnalysisResults7135.root", "TrainOutputs/AnalysisResults7136.root", "TrainOutputs/AnalysisResults7159.root", "TrainOutputs/AnalysisResults7160.root"};

vector<const char*> filenames3 = {"TrainOutputs/AnalysisResults7208.root", "TrainOutputs/AnalysisResults7209.root", "TrainOutputs/AnalysisResults7210.root", "TrainOutputs/AnalysisResults7211.root", "TrainOutputs/AnalysisResults7212.root", "TrainOutputs/AnalysisResults7213.root", "TrainOutputs/AnalysisResults7214.root", "TrainOutputs/AnalysisResults7215.root", "TrainOutputs/AnalysisResults7216.root", "TrainOutputs/AnalysisResults7217.root", "TrainOutputs/AnalysisResults7218.root", "TrainOutputs/AnalysisResults7219.root", "TrainOutputs/AnalysisResults7220.root", "TrainOutputs/AnalysisResults7221.root", "TrainOutputs/AnalysisResults7222.root", "TrainOutputs/AnalysisResults7223.root", "TrainOutputs/AnalysisResults7224.root"};

int j = 0;
int start = 0;

  vector<TH1D*> bin(0);
  vector<TH1D*> gen(0);
  vector<TH1D*> det(0);

  vector<const char*> bins = {"bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20"};
  vector<const char*> gens = {"gen1", "gen2", "gen3", "gen4", "gen5", "gen6", "gen7", "gen8", "gen9", "gen10", "gen11", "gen12", "gen13", "gen14", "gen15", "gen16", "gen17", "gen18", "gen19", "gen20"};
  vector<const char*> dets = {"det1", "det2", "det3", "det4", "det5", "det6", "det7", "det8", "det9", "det10", "det11", "det12", "det13", "det14", "det15", "det16", "det17", "det18", "det19", "det20"};

void outlierspec1 (double cl, double cr) {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  cout << "line 50\n";
  for (int i = 0; i < filenames3.size(); i++)    { 
      getbin(filenames3[i], cl, cr);
        bin.push_back((TH1D*)binspec->Clone(bins[i]));
        gen.push_back((TH1D*)genspec->Clone(gens[i]));
        det.push_back((TH1D*)detspec->Clone(dets[i]));
          binsum1->Add(bin[i]);
          gensum1->Add(gen[i]);
          detsum1->Add(det[i]);
      }
  cout << "line60\n";  
  bin.push_back(binsum1);
  gen.push_back(gensum1);
  det.push_back(detsum1);


  cout << "Total events: " << totevents <<"\n";
  double avgevents = (double)totevents/20.0;
  cout << "Average Events: " << avgevents <<"\n";
  
  //Scale by Average Events
  /*for (int i = 0; i <= filenames3.size(); i++) {
      bin[i]->Scale(avgevents);  
      gen[i]->Scale(avgevents);
      det[i]->Scale(avgevents);}
  */
  setstyle(bin, gen, det);


  stringstream ss;
  ss << Form("R = 0.4, %.0f-%.0f", cl, cr) << "%";
  auto datarun = new TLegend(0.15, 0.75, 0.35, 0.85);
  datarun->SetTextSize(0.03);
  datarun->SetBorderSize(0);
  //datarun->AddEntry((TObject*)0, "LHC18q + LHC18r", "");
  //datarun->AddEntry((TObject*)0, "Pb-Pb #sqrt{s_{NN} = 5.02 TeV}", "");
  datarun->AddEntry((TObject*)0, ss.str().c_str(), "");

vector<string> newleg(0);
for (int i = 0; i < filenames3.size(); i++) {
    int index = i+start;
    string nl = Form("bin %d", index);      //bin key
    newleg.push_back(nl);
  }

cout << "start: " << start <<"\n";

vector<const char*> newlegchar(0);
for (int i = 0; i < newleg.size(); i++) newlegchar.push_back(newleg[i].c_str());

  auto leg1 = new TLegend(0.45, 0.5, 0.65, 0.85);
  leg1->SetTextSize(0.035);
  leg1->SetBorderSize(0);
  auto leg2 = new TLegend(0.65, 0.5, 0.85, 0.85);
  leg2->SetTextSize(0.035);
  leg2->SetBorderSize(0);
  for (int i = 0; i < filenames3.size(); i++) {
      if (i < 10)      leg1->AddEntry(bin[i], newlegchar[i], "pl");
      if (i >= 10)     leg2->AddEntry(bin[i], newlegchar[i], "pl");
  }


  TLine *line = new TLine(10.0, pow(10,4), 10.0, pow(10,16));
  line->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1", "outlier1", 1500, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad11 = new TPad("pad11", "pad11", 0.0, 0.05, 0.33, 1.0);
     pad11->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad11->Draw();             // Draw the upper pad: pad1
     pad11->cd();               // pad1 becomes the current pad
     pad11->SetLogy();
     binsum1->SetXTitle("p_{T}^{hybrid}");
     double ymax = binsum1->GetBinContent(2);
     double ymin = binsum1->GetBinContent(40);
       cout << "ymax: " << ymax <<"\n";
       cout << "ymin: " << ymin <<"\n";
     binsum1->SetMinimum(pow(10,-4)*ymin);
     binsum1->SetMaximum(pow(10,5)*ymax);
     binsum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   bin[i]->Draw("same");  
       datarun->Draw("same");
       leg1->Draw("same");
       if (filenames3.size() > 10) leg2->Draw("same");
       //line->Draw("same");
     c1->cd();
   TPad *pad12 = new TPad("pad12", "pad12", 0.33, 0.05, 0.67, 1.0);
     pad12->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad12->Draw();             // Draw the upper pad: pad1
     pad12->cd();               // pad1 becomes the current pad
     pad12->SetLogy();
     detsum1->SetXTitle("p_{T}^{detector}");
     detsum1->SetMinimum(pow(10,-4)*ymin);
     detsum1->SetMaximum(pow(10,5)*ymax);
     detsum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   det[i]->Draw("same");  
       datarun->Draw("same");
       leg1->Draw("same");
       if (filenames3.size() > 10) leg2->Draw("same");
       //line->Draw("same");
     c1->cd();
   TPad *pad13 = new TPad("pad13", "pad13", 0.67, 0.05, 1.0, 1.0);
     pad13->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad13->Draw();             // Draw the upper pad: pad1
     pad13->cd();               // pad1 becomes the current pad
     pad13->SetLogy();
     gensum1->SetXTitle("p_{T}^{truth}");
     gensum1->SetMinimum(pow(10,-4)*ymin);
     gensum1->SetMaximum(pow(10,5)*ymax);
     gensum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   gen[i]->Draw("same");  
       datarun->Draw("same");
       leg1->Draw("same");
       if (filenames3.size() > 10) leg2->Draw("same");
       //line->Draw("same");
 

  
}








void getbin(const char *analysisfile, double centl, double centr) {

    j++;    

    binspec->Reset();
    genspec->Reset();
    detspec->Reset(); 


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
    float hybridpT, detpT, partpT, centrality, area, matchradius;
    T->SetBranchAddress("Jet_Pt", &hybridpT);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &detpT);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &partpT);
    T->SetBranchAddress("Event_Centrality", &centrality);
    T->SetBranchAddress("Jet_Area", &area);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Distance", &matchradius); 

    int entries = T->GetEntries();

    //Determine pT hard scaling
    double xsect = pThardscale->GetBinContent(pThardbin->GetMean() + 1);
    double xscale = xsect*pThardscale->GetEntries();
    double trials = pThardbin->GetBinContent(pThardbin->GetMean() + 1);
 
    double n_event = eventcount->GetBinContent(1);
    double pThardFrac = xscale/(/*n_event*/trials);

    if (j == 1)  start = pThardbin->GetMean();  


    //Determine Extractor Bin Scaling
    const int extractBins = 8;
    //double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05, 0.05};
    double extractFrac[extractBins] = {0.0, 0.03, 0.03, 0.05, 0.05, 0.09, 0.09, 0.1};
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0};
    

    int fillcount = 0;
    for (int i = 0; i < entries; i++) {
       T->GetEntry(i);
       if (centrality < centl || centrality > centr)  continue;
       //if (area < 0.4)  continue;
       //if (matchradius > 0.299) continue;
       if (hybridpT < 10.0 || detpT < 10.0 || partpT < 10.0 )  continue;
       //if (hybridpT > 3.0*partpT)   continue;
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
       binspec->Fill(hybridpT);
       genspec->Fill(partpT);
       detspec->Fill(detpT);
   }

    vector<int> manset(0);
    for (int b = 0; b < binres; b++) {
        if (binspec->GetBinContent(b+1) != 1) continue;
        if (binspec->GetBinContent(b) == 0 && binspec->GetBinContent(b+2) ==0)  manset.push_back(b+1);
    }
   
   binspec->Reset();
   genspec->Reset();
   detspec->Reset();

    for (int i = 0; i < entries; i++) {
       T->GetEntry(i);
       if (centrality < centl || centrality > centr)  continue;
       //if (area < 0.4)  continue;
       //if (matchradius > 0.299) continue;
       if (hybridpT < 10.0 || detpT < 10.0 || partpT < 10.0 )  continue;
       //if (hybridpT > 3.0*partpT)   continue;
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
       binspec->Fill(hybridpT, pThardFrac/extractFrac[scaleindex]);
       genspec->Fill(partpT, pThardFrac/extractFrac[scaleindex]);
       detspec->Fill(detpT, pThardFrac/extractFrac[scaleindex]);
   } 
  
   for (int i = 0; i < manset.size(); i++)   {
     binspec->SetBinContent(manset[i], 0);
   }

   manset.resize(0);

    cout << "pT hard bin: " << pThardbin->GetMean() <<"\n";
    //cout << "	Scale Factor: " << pThardFrac << "*avg_events" << "\n";    

    totevents += n_event;
    int pTindex = pThardbin->GetMean() - 1;
   
    lefile->Close();

}

void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3) {
  vector<int> colors = {kRed+2, kRed, kOrange+8, kOrange-4, kYellow, kSpring+10, kSpring, kGreen-3, kGreen+2, kTeal+6, kTeal, kCyan-3, kAzure+7, kBlue+1, kBlue-2, kViolet+2, kViolet, kMagenta-4, kPink+6, kPink+10};
  int entryremoval = colors.size() - filenames3.size();
  cout << "Entry removal" << entryremoval <<"\n"; 
  for (int i = 0; i < entryremoval; i++)   colors.pop_back();
  colors.push_back(kBlack);

  for (int i = 0; i < vec1.size(); i++) {
      vec1[i]->SetLineColor(colors[i]);
      vec1[i]->SetMarkerColor(colors[i]);
        vec1[i]->SetLineWidth(2);
        vec1[i]->SetMarkerSize(0.5);
        vec1[i]->SetMarkerStyle(20);
      vec2[i]->SetLineColor(colors[i]);
      vec2[i]->SetMarkerColor(colors[i]);
        vec2[i]->SetLineWidth(2);
        vec2[i]->SetMarkerSize(0.5);
        vec2[i]->SetMarkerStyle(20);
      vec3[i]->SetLineColor(colors[i]);
      vec3[i]->SetMarkerColor(colors[i]);
        vec3[i]->SetLineWidth(2);
        vec3[i]->SetMarkerSize(0.5);
        vec3[i]->SetMarkerStyle(20);
  }



}
