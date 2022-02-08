#include "TH1D.h"
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

int lbin = 0;
int rbin = 200;
int binres = 200;
    
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
vector<const char*> filenames3 = {"TrainOutputs/AnalysisResults7158.root", "TrainOutputs/AnalysisResults7121.root", "TrainOutputs/AnalysisResults7122.root", "TrainOutputs/AnalysisResults7123.root", "TrainOutputs/AnalysisResults7124.root", "TrainOutputs/AnalysisResults7126.root", "TrainOutputs/AnalysisResults7127.root", "TrainOutputs/AnalysisResults7128.root", "TrainOutputs/AnalysisResults7129.root", "TrainOutputs/AnalysisResults7131.root", "TrainOutputs/AnalysisResults7132.root", "TrainOutputs/AnalysisResults7133.root", "TrainOutputs/AnalysisResults7134.root", "TrainOutputs/AnalysisResults7135.root", "TrainOutputs/AnalysisResults7136.root", "TrainOutputs/AnalysisResults7159.root", "TrainOutputs/AnalysisResults7160.root"};


  vector<TH1D*> bin(0);
  vector<TH1D*> gen(0);
  vector<TH1D*> det(0);

  vector<const char*> bins = {"bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20"};
  vector<const char*> gens = {"gen1", "gen2", "gen3", "gen4", "gen5", "gen6", "gen7", "gen8", "gen9", "gen10", "gen11", "gen12", "gen13", "gen14", "gen15", "gen16", "gen17", "gen18", "gen19", "gen20"};
  vector<const char*> dets = {"det1", "det2", "det3", "det4", "det5", "det6", "det7", "det8", "det9", "det10", "det11", "det12", "det13", "det14", "det15", "det16", "det17", "det18", "det19", "det20"};



void laura(double cl, double cr, const char* fileo) {

  /*
  TFile *newfile = new TFile(fileo);
  TH2D* h2_pT = (TH2D*)newfile->Get("h2_pT"); 
  TH1D* h1_tot = (TH1D*)h2_pT->ProjectionX("h1", 5, 100);
  h1_tot->GetXaxis()->SetRangeUser(0.0, 200.0);
  */

  TFile* newfile = new TFile(fileo);
  std::vector<TH1D*> histvec_rg;
  std::stringstream qq;
  for (int i = 0; i < 20; i++)
    {
      qq << "h2_pthard" << i+1;
      TH2D* h2 = (TH2D*)newfile->Get(qq.str().c_str());
      qq.str("");
      qq << "h1_pthard" << i+1;
      TH1D* h1 = (TH1D*)h2->ProjectionY(qq.str().c_str());
         //h1->Scale(1./h1->Integral(), "width");
         histvec_rg.push_back(h1);
         qq.str("");
     }
      TH1D* h1_tot = (TH1D*)histvec_rg.at(0)->Clone("h1_tot");
      for (int i = 1; i < 20; i++) h1_tot->Add(histvec_rg.at(i));

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  cout << "h1_tot->GetXaxis()->GetNbins(): " << h1_tot->GetXaxis()->GetNbins() <<"\n";
  cout << "	Low(0): " << h1_tot->GetXaxis()->GetBinLowEdge(0) <<"\n";
  cout << "	Low(1): " << h1_tot->GetXaxis()->GetBinLowEdge(1) <<"\n";
  cout << "	Low(199): " << h1_tot->GetXaxis()->GetBinLowEdge(199) <<"\n";

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
  for (int i = 0; i <= filenames3.size(); i++) {
      bin[i]->Scale(1.0/avgevents);  
      gen[i]->Scale(1.0/avgevents);
      det[i]->Scale(1.0/avgevents);}

  setstyle(bin, gen, det);
  h1_tot->SetLineColor(kBlack);
  h1_tot->SetLineWidth(2);
  h1_tot->SetMarkerColor(kBlack);
  detsum1->SetLineColor(kBlue+1);
  detsum1->SetMarkerColor(kBlue+1);
  detsum1->SetLineWidth(2);
  binsum1->SetLineColor(kGreen+1);
  binsum1->SetLineWidth(2);
  binsum1->SetMarkerColor(kGreen+1);

  stringstream ss;
  ss << Form("R = 0.4, %.0f-%.0f", cl, cr) << "%";
  auto datarun = new TLegend(0.55, 0.7, 0.75, 0.85);
  datarun->SetTextSize(0.03);
  datarun->SetBorderSize(0);
  //datarun->AddEntry((TObject*)0, "LHC18q + LHC18r", "");
  //datarun->AddEntry((TObject*)0, "Pb-Pb #sqrt{s_{NN} = 5.02 TeV}", "");
  datarun->AddEntry((TObject*)0, ss.str().c_str(), "");
  datarun->AddEntry(h1_tot, "Laura", "pl");
  datarun->AddEntry(binsum1, "Caitie pre-cut", "pl");
  datarun->AddEntry(detsum1, "Caitie post-cut", "pl");


  TLine *line = new TLine(40.0, 1.0, 120.0, 1.0);
  line->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1", "new", 800, 400);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad11 = new TPad("pad11", "pad11", 0.0, 0.05, 0.5, 1.0);
     pad11->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad11->Draw();             // Draw the upper pad: pad1
     pad11->cd();               // pad1 becomes the current pad
     pad11->SetLogy();
     binsum1->SetXTitle("p_{T}^{hybrid}");
     double ymax = binsum1->GetBinContent(100);
     double ymin = binsum1->GetBinContent(100);
     h1_tot->SetTitle("");
     h1_tot->GetXaxis()->SetRangeUser(40.0, 120.0);
     h1_tot->SetMinimum(pow(10,-1)*ymin);
     h1_tot->SetMaximum(pow(10,4)*ymax);
     h1_tot->Draw("same");
     binsum1->Draw("same");
     detsum1->Draw("same");
       datarun->Draw("same");
   c1->cd();
   TPad *pad12 = new TPad("pad12", "pad12", 0.5, 0.05, 1.0, 1.0);
     pad12->SetBottomMargin(0.1);
     pad12->Draw();
     pad12->cd();
     //pad12->SetLogy();
     TH1D* divi = (TH1D*)h1_tot->Clone("divi");
     divi->Divide(binsum1);
       TH1D* divi2 = (TH1D*)h1_tot->Clone("divi2");
       divi2->Divide(detsum1);
       divi2->SetLineColor(kBlue);
     divi->SetTitle("");
     divi->GetXaxis()->SetRangeUser(40.0, 120.0);
     divi->SetMinimum(0.01);
     divi->SetMaximum(10.0);
     divi->Draw("same");
       divi2->Draw("same");
  auto leg = new TLegend(0.55, 0.7, 0.75, 0.85);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->AddEntry(divi, "Laura/Caitie (pre)", "pl");
  leg->AddEntry(divi2, "Laura/Caitie (post)", "pl");
     leg->Draw("same");
     line->Draw("same");

}



void getbin(const char *analysisfile, double centl, double centr) {

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


    //Determine Extractor Bin Scaling
    const int extractBins = 8;
    //double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05, 0.05};
    double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.03, 0.03, 0.07, 0.07, 0.05};
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0};
    

    int fillcount = 0;
    for (int i = 0; i < entries; i++) {
       T->GetEntry(i);
       if (centrality < centl || centrality > centr)  continue;
       //if (area < 0.4)  continue;
       //if (matchradius > 0.299) continue;
       if (hybridpT < 10.0 || detpT < 10.0 || partpT < 10.0 )  continue;
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
       if (hybridpT > 3.0*partpT)   continue;
       genspec->Fill(partpT, pThardFrac/extractFrac[scaleindex]);
       detspec->Fill(detpT, pThardFrac/extractFrac[scaleindex]);
   }

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

