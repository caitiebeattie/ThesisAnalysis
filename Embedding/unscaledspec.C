
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;


void getq2(const char *anfile, double cel, double cer);
void getbin(const char *analysisfile, double centl, double centr, vector<vector<long double>> q2vec, int ese);
void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3);

int lbin = 10;
int rbin = 200;
int binres = 95;
    
vector<vector<long double>> q2percentiles(0);

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

//R = 0.2
/*vector<const char*> filenames3 = {"../ESE/LHC18/R02/AnalysisResults8145.root", "../ESE/LHC18/R02/AnalysisResults8146.root",
  "../ESE/LHC18/R02/AnalysisResults8147.root", "../ESE/LHC18/R02/AnalysisResults8148.root", "../ESE/LHC18/R02/AnalysisResults8149.root", 
  "../ESE/LHC18/R02/AnalysisResults8150.root", "../ESE/LHC18/R02/AnalysisResults8151.root", "../ESE/LHC18/R02/AnalysisResults8152.root", 
  "../ESE/LHC18/R02/AnalysisResults8153.root", "../ESE/LHC18/R02/AnalysisResults8154.root", "../ESE/LHC18/R02/AnalysisResults8155.root",
  "../ESE/LHC18/R02/AnalysisResults8156.root", "../ESE/LHC18/R02/AnalysisResults8157.root", "../ESE/LHC18/R02/AnalysisResults8158.root",
  "../ESE/LHC18/R02/AnalysisResults8159.root", "../ESE/LHC18/R02/AnalysisResults8160.root", "../ESE/LHC18/R02/AnalysisResults8161.root",
  "../ESE/LHC18/R02/AnalysisResults8162.root", "../ESE/LHC18/R02/AnalysisResults8163.root", "../ESE/LHC18/R02/AnalysisResults8164.root"};
*/
//R = 0.4
    vector<const char*> filenames3 = {"../ESE/LHC18/R04/AnalysisResults8124.root",
       "../ESE/LHC18/R04/AnalysisResults8299.root", "../ESE/LHC18/R04/AnalysisResults8126.root", "../ESE/LHC18/R04/AnalysisResults8298.root", 
       "../ESE/LHC18/R04/AnalysisResults8301.root", "../ESE/LHC18/R04/AnalysisResults8129.root", "../ESE/LHC18/R04/AnalysisResults8130.root", 
       "../ESE/LHC18/R04/AnalysisResults8131.root", "../ESE/LHC18/R04/AnalysisResults8132.root", "../ESE/LHC18/R04/AnalysisResults8133.root",
       "../ESE/LHC18/R04/AnalysisResults8134.root", "../ESE/LHC18/R04/AnalysisResults8135.root", "../ESE/LHC18/R04/AnalysisResults8136.root",
       "../ESE/LHC18/R04/AnalysisResults8137.root", "../ESE/LHC18/R04/AnalysisResults8138.root", "../ESE/LHC18/R04/AnalysisResults8140.root",
       "../ESE/LHC18/R04/AnalysisResults8141.root", "../ESE/LHC18/R04/AnalysisResults8142.root", "../ESE/LHC18/R04/AnalysisResults8143.root"};
//R = 0.4
//vector<const char*> filenames3 = {"../ESE/LHC18/R04/AnalysisResults8001.root", "../ESE/LHC18/R04/AnalysisResults8002.root", 
//    "../ESE/LHC18/R04/AnalysisResults8003.root", "../ESE/LHC18/R04/AnalysisResults8004.root", "../ESE/LHC18/R04/AnalysisResults8005.root",
//    "../ESE/LHC18/R04/AnalysisResults8006.root", "../ESE/LHC18/R04/AnalysisResults8007.root", "../ESE/LHC18/R04/AnalysisResults8008.root"};


//R = 0.2, 30-50% with pileup cuts
//vector<const char*> filenames3 = {"../ESE/LHC18/AnalysisResults7909.root", "../ESE/LHC18/AnalysisResults7930.root", 
//       "../ESE/LHC18/AnalysisResults7911.root", "../ESE/LHC18/AnalysisResults7912.root", "../ESE/LHC18/AnalysisResults7913.root",
//       "../ESE/LHC18/AnalysisResults7914.root", "../ESE/LHC18/AnalysisResults7915.root", "../ESE/LHC18/AnalysisResults7916.root",
//       "../ESE/LHC18/AnalysisResults7917.root", "../ESE/LHC18/AnalysisResults7918.root", "../ESE/LHC18/AnalysisResults7919.root", 
//       "../ESE/LHC18/AnalysisResults7920.root", "../ESE/LHC18/AnalysisResults7921.root", "../ESE/LHC18/AnalysisResults7922.root",
//       "../ESE/LHC18/AnalysisResults7923.root", "../ESE/LHC18/AnalysisResults7924.root", "../ESE/LHC18/AnalysisResults7925.root",
//       "../ESE/LHC18/AnalysisResults7926.root", "../ESE/LHC18/AnalysisResults7927.root", "../ESE/LHC18/AnalysisResults7928.root"};

int j = 0;
int start = 0;

  vector<TH1D*> bin(0);
  vector<TH1D*> gen(0);
  vector<TH1D*> det(0);

  vector<const char*> bins = {"bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20"};
  vector<const char*> gens = {"gen1", "gen2", "gen3", "gen4", "gen5", "gen6", "gen7", "gen8", "gen9", "gen10", "gen11", "gen12", "gen13", "gen14", "gen15", "gen16", "gen17", "gen18", "gen19", "gen20"};
  vector<const char*> dets = {"det1", "det2", "det3", "det4", "det5", "det6", "det7", "det8", "det9", "det10", "det11", "det12", "det13", "det14", "det15", "det16", "det17", "det18", "det19", "det20"};





//========================================
//========================================
//=========== main function ==============
//========================================
//========================================

void unscaledspec (double cl, double cr, int ese = 2) {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  getq2("../ESE/AnalysisResults8119.root", cl, cr);

  cout << "line 50\n";
  for (int i = 0; i < filenames3.size(); i++)    { 
      getbin(filenames3[i], cl, cr, q2percentiles, ese);
        bin.push_back((TH1D*)binspec->Clone(bins[i]));
        gen.push_back((TH1D*)genspec->Clone(gens[i]));
        det.push_back((TH1D*)detspec->Clone(dets[i]));
          binsum1->Add(bin[i]);
          gensum1->Add(gen[i]);
          detsum1->Add(det[i]);
      }
  bin.push_back(binsum1);
  gen.push_back(gensum1);
  det.push_back(detsum1);


  cout << "Total events: " << totevents <<"\n";
  double avgevents = (double)totevents/20.0;
  cout << "Average Events: " << avgevents <<"\n";
  
  //Scale by Average Events
  for (int i = 0; i <= filenames3.size(); i++) {
      //bin[i]->Scale(avgevents);  
      //gen[i]->Scale(avgevents);
      //det[i]->Scale(avgevents);}
     }
  setstyle(bin, gen, det);

  auto info = new TLegend(0.2, 0.7, 0.35, 0.85);
    info->SetTextSize(0.035);
    info->SetBorderSize(0);
    info->AddEntry((TObject*)0, "ALICE, Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 30#font[122]{-}50%", "");
    info->AddEntry((TObject*)0, "Unscaled Embedding Spectra, LHC18 pass3", "");

  auto datarun = new TLegend(0.15, 0.75, 0.35, 0.85);
  datarun->SetTextSize(0.035);
  datarun->SetBorderSize(0);
  datarun->AddEntry((TObject*)0, "Charged Jets, R = 0.4, anti-#it{k}_{T}, |#it{#eta}_{jet}| < 0.5", "");
  if (ese == 0)  datarun->AddEntry((TObject*)0, "10-40\% smallest #it{q}_{2}, #it{p}_{T, true} > 10 GeV/#it{c}", ""); 
  if (ese == 1)  datarun->AddEntry((TObject*)0, "30\% largest #it{q}_{2}, #it{p}_{T, true} > 10 GeV/#it{c}", ""); 

vector<string> newleg(0);
for (int i = 0; i < filenames3.size(); i++) {
    int index = i+start;
    string nl = Form("bin %d", index);      //bin key
    newleg.push_back(nl);
  }

cout << "start: " << start <<"\n";

vector<const char*> newlegchar(0);
for (int i = 0; i < newleg.size(); i++) newlegchar.push_back(newleg[i].c_str());

  auto leg1 = new TLegend(0.45, 0.55, 0.65, 0.85);
  leg1->SetTextSize(0.035);
  leg1->SetBorderSize(0);
  auto leg2 = new TLegend(0.65, 0.55, 0.85, 0.85);
  leg2->SetTextSize(0.035);
  leg2->SetBorderSize(0);
  for (int i = 0; i < filenames3.size(); i++) {
      if (i < 10)      leg1->AddEntry(bin[i], newlegchar[i], "pl");
      if (i >= 10)     leg2->AddEntry(bin[i], newlegchar[i], "pl");
  }


  TLine *line = new TLine(10.0, 1.0, 200.0, 1.0);
  line->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1", "new", 1500, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad11 = new TPad("pad11", "pad11", 0.0, 0.05, 0.33, 1.0);
     pad11->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad11->SetLeftMargin(0.15);
     pad11->Draw();             // Draw the upper pad: pad1
     pad11->cd();               // pad1 becomes the current pad
     pad11->SetLogy();
     gensum1->SetXTitle("#it{p}_{T}^{truth} (GeV/#it{c})");
     gensum1->SetYTitle("N_{jet}");
     double ymax = gensum1->GetBinContent(2);
     double ymin = gensum1->GetBinContent(40);
       cout << "ymax: " << ymax <<"\n";
       cout << "ymin: " << ymin <<"\n";
     gensum1->SetMinimum(pow(10,-1));
     gensum1->SetMaximum(pow(10,5)*ymax);
     gensum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   gen[i]->Draw("same");  
       info->Draw("same");
       line->Draw("same");
     c1->cd();
   TPad *pad12 = new TPad("pad12", "pad12", 0.33, 0.05, 0.67, 1.0);
     pad12->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad12->Draw();             // Draw the upper pad: pad1
     pad12->cd();               // pad1 becomes the current pad
     pad12->SetLogy();
     detsum1->SetXTitle("#it{p}_{T}^{detector} (GeV/#it{c})");
     detsum1->SetMinimum(pow(10,-1));
     detsum1->SetMaximum(pow(10,5)*ymax);
     detsum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   det[i]->Draw("same");  
       datarun->Draw("same");
       line->Draw("same");
     c1->cd();
   TPad *pad13 = new TPad("pad13", "pad13", 0.67, 0.05, 1.0, 1.0);
     pad13->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad13->Draw();             // Draw the upper pad: pad1
     pad13->cd();               // pad1 becomes the current pad
     pad13->SetLogy();
     binsum1->SetXTitle("#it{p}_{T}^{hybrid} (GeV/#it{c})");
     binsum1->SetMinimum(pow(10,-1));
     binsum1->SetMaximum(pow(10,5)*ymax);
     binsum1->Draw("same"); 
     for (int i = 0; i <= filenames3.size(); i++)   bin[filenames3.size()-i]->Draw("same");  
       binsum1->Draw("same");
       leg1->Draw("same");
       if (filenames3.size() > 10) leg2->Draw("same");
       line->Draw("same");

  
}


//------------------------------------------------------------------------------------------------

void getq2(const char *anfile, double cel, double cer)  {

      TFile *fin = new TFile(anfile);
      //Calculate q2 percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)fin->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2_EP");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrFullV0");
      //initialize variables used to calculate percentiles
      vector<long double> percentiles(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercent(0);                      //vector that stores percentiles vectors for all centralities
      int centdiff = (int)cer-(int)cel;
      int leftc, rightc;
      //get percentiles in bins of 1% centrality
      for (int j = 0; j < centdiff; j++) {
        long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
        int percentilen = 1;                                             //number percentile we're on
        leftc = (int)cel + j; 
        rightc = (int)cel + j + 1;    
        q2hist2D->GetXaxis()->SetRangeUser(leftc, rightc);
        TH1F *q2hist = (TH1F*)q2hist2D->ProjectionY();
        //q2hist->Rebin(100);
        long double percentilemarker = (long double)q2hist->GetEntries()/10.0;      //number of entries that compose 10% of sample
        for (int i = 0; i <= q2hist->GetNbinsX(); i++)   {
          percentileticker += q2hist->GetBinContent(i);         //increment the percentileticker
          if (percentileticker >= percentilen*percentilemarker)  {
              percentiles.push_back(q2hist->GetBinCenter(i));
              percentilen++;
          }
        if (i == q2hist->GetNbinsX())  {
            allPercent.push_back(percentiles);
            percentiles.resize(0);
            break;}             
        }
      }
  q2percentiles = allPercent;

}





void getbin(const char *analysisfile, double centl, double centr, vector<vector<long double>> q2vec, int ese) {

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
    float hybridpT, detpT, partpT, centrality, area, matchradius, q2;
    T->SetBranchAddress("Jet_Pt", &hybridpT);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &detpT);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &partpT);
    T->SetBranchAddress("Event_Centrality", &centrality);
    T->SetBranchAddress("Jet_Area", &area);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Distance", &matchradius); 
    T->SetBranchAddress("Event_Q2VectorV0M", &q2); 

    int entries = T->GetEntries();

    //Determine pT hard scaling
    double xsect = pThardscale->GetBinContent(pThardbin->GetMean() + 1);
    double xscale = xsect*pThardscale->GetEntries();
    double trials = pThardbin->GetBinContent(pThardbin->GetMean() + 1);
 
    double n_event = eventcount->GetBinContent(1);
    double pThardFrac = xscale/(/*n_event*/trials);

    if (j == 1)  start = pThardbin->GetMean();  

    

    int fillcount = 0;
    for (int i = 0; i < entries; i++) {
       T->GetEntry(i);
       if (centrality < centl || centrality > centr)  continue;
       if (partpT < 10.0 )  continue;
       if (ese == 0)  {
          if (q2 < q2vec[(int)centrality-30][0])   continue;
          if (q2 > q2vec[(int)centrality-30][3])   continue;}
       if (ese == 1)  {if (q2 < q2vec[(int)centrality-30][6])   continue;}
       binspec->Fill(hybridpT);
       genspec->Fill(partpT);
       detspec->Fill(detpT);
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
