
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
 
long int prematch = 0;
long int postmatch = 0;
   
TH1D *det1spec = new TH1D("det1spec", "", binres, lbin, rbin);
TH1D *det2spec = new TH1D("det2spec", "", binres, lbin, rbin);
TH1D *gen1spec = new TH1D("gen1spec", "", binres, lbin, rbin);
TH1D *gen2spec = new TH1D("gen2spec", "", binres, lbin, rbin);
TH1D *leftoversh = new TH1D("leftoversh", "", binres, lbin, rbin);
TH1D *leftoverst = new TH1D("leftoverst", "", binres, lbin, rbin);
TH1D * dpt = new TH1D("dpt", "", 100, -100, 100);

TH1D *det1sum = new TH1D("det1sum", "", binres, lbin, rbin);
TH1D *det2sum = new TH1D("det2sum", "", binres, lbin, rbin);

long int totevents = 0;

  //vector<const char*> filenames1 = {"TrainOutputs/old/AnalysisResults7007.root", "TrainOutputs/old/AnalysisResults7008.root", "TrainOutputs/old/AnalysisResults7009.root", "TrainOutputs/old/AnalysisResults7010.root", "TrainOutputs/old/AnalysisResults7011.root", "TrainOutputs/old/AnalysisResults7012.root", "TrainOutputs/old/AnalysisResults7013.root", "TrainOutputs/old/AnalysisResults7014.root", "TrainOutputs/old/AnalysisResults7015.root", "TrainOutputs/old/AnalysisResults7016.root", "TrainOutputs/old/AnalysisResults7017.root", "TrainOutputs/old/AnalysisResults7018.root", "TrainOutputs/old/AnalysisResults7022.root", "TrainOutputs/old/AnalysisResults7019.root", "TrainOutputs/old/AnalysisResults7020.root", "TrainOutputs/old/AnalysisResults7021.root", "TrainOutputs/old/AnalysisResults7023.root", "TrainOutputs/old/AnalysisResults7024.root", "TrainOutputs/old/AnalysisResults7025.root", "TrainOutputs/old/AnalysisResults7026.root"};
  vector<const char*> filenames2 = {/*"TrainOutputs/old/AnalysisResults6846.root", "TrainOutputs/old/AnalysisResults6847.root", "TrainOutputs/old/AnalysisResults6848.root", */"TrainOutputs/old/AnalysisResults6849.root", "TrainOutputs/old/AnalysisResults6850.root", "TrainOutputs/old/AnalysisResults6851.root", "TrainOutputs/old/AnalysisResults6852.root", "TrainOutputs/old/AnalysisResults6853.root", "TrainOutputs/old/AnalysisResults6854.root", "TrainOutputs/old/AnalysisResults6855.root", "TrainOutputs/old/AnalysisResults6856.root", "TrainOutputs/old/AnalysisResults6857.root", "TrainOutputs/old/AnalysisResults6858.root", "TrainOutputs/old/AnalysisResults6859.root", "TrainOutputs/old/AnalysisResults6860.root", "TrainOutputs/old/AnalysisResults6861.root", "TrainOutputs/old/AnalysisResults6862.root", "TrainOutputs/old/AnalysisResults6863.root", "TrainOutputs/old/AnalysisResults6864.root", "TrainOutputs/old/AnalysisResults6865.root"};
vector<const char*> filenames1 = {"TrainOutputs/AnalysisResults7158.root", "TrainOutputs/AnalysisResults7121.root", "TrainOutputs/AnalysisResults7122.root", "TrainOutputs/AnalysisResults7123.root", "TrainOutputs/AnalysisResults7124.root", "TrainOutputs/AnalysisResults7126.root", "TrainOutputs/AnalysisResults7127.root", "TrainOutputs/AnalysisResults7128.root", "TrainOutputs/AnalysisResults7129.root", "TrainOutputs/AnalysisResults7131.root", "TrainOutputs/AnalysisResults7132.root", "TrainOutputs/AnalysisResults7133.root", "TrainOutputs/AnalysisResults7134.root", "TrainOutputs/AnalysisResults7135.root", "TrainOutputs/AnalysisResults7136.root", "TrainOutputs/AnalysisResults7159.root", "TrainOutputs/AnalysisResults7160.root"};

  vector<TH1D*> d1(0);
  vector<TH1D*> d2(0);

  vector<const char*> det1 = {"gen1", "gen2", "gen3", "gen4", "gen5", "gen6", "gen7", "gen8", "gen9", "gen10", "gen11", "gen12", "gen13", "gen14", "gen15", "gen16", "gen17", "gen18", "gen19", "gen20"};
  vector<const char*> det2 = {"det1", "det2", "det3", "det4", "det5", "det6", "det7", "det8", "det9", "det10", "det11", "det12", "det13", "det14", "det15", "det16", "det17", "det18", "det19", "det20"};

void ptspecchanges () {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  cout << "line 50\n";
  for (int i = 0; i < filenames1.size(); i++)    { 
      getbin(filenames1[i], 0.0, 10.0);
        d1.push_back((TH1D*)det1spec->Clone(det1[i]));
        d2.push_back((TH1D*)det2spec->Clone(det2[i]));
          det1sum->Add(d1[i]);
          det2sum->Add(d2[i]);
      }
  cout << "line60\n";  
  d1.push_back(det1sum);
  d2.push_back(det2sum);


  cout << "Total events: " << totevents <<"\n";
  double avgevents = (double)totevents/20.0;
  cout << "Average Events: " << avgevents <<"\n";
  
  //Scale by Average Events
  for (int i = 0; i <= filenames1.size(); i++) {
      d1[i]->Scale(avgevents);  
      d2[i]->Scale(avgevents);}

  TH1D* diff = (TH1D*)leftoversh->Clone("diff");
  diff->Add(leftoverst, -1);

  //det1sum->Scale(1/det1sum->Integral());
  //det2sum->Scale(1/det2sum->Integral());
  //gen1spec->Scale(1/gen1spec->Integral());
  //gen2spec->Scale(1/gen2spec->Integral());

  det1sum->SetLineColor(kRed);
  det1sum->SetLineWidth(2);
  det1sum->SetMarkerStyle(20);
  det1sum->SetMarkerSize(0.7);
  det1sum->SetMarkerColor(kRed+2);
  det2sum->SetLineColor(kBlack);
  det2sum->SetLineWidth(2);
  gen1spec->SetLineColor(kRed);
  gen1spec->SetLineWidth(2);
  gen2spec->SetLineColor(kBlack);
  gen2spec->SetLineWidth(2);
  leftoversh->SetLineColor(kRed);
  leftoverst->SetLineColor(kBlue);

  auto datarun = new TLegend(0.45, 0.7, 0.85, 0.85);
  datarun->SetTextSize(0.03);
  datarun->SetBorderSize(0);
  //datarun->AddEntry((TObject*)0, "LHC18q + LHC18r", "");
  //datarun->AddEntry((TObject*)0, "Pb-Pb #sqrt{s_{NN} = 5.02 TeV}", "");
  datarun->AddEntry((TObject*)0, "R = 0.4, 0-10%", "");

  auto lentries = new TLegend(0.5, 0.6, 0.85, 0.85);
  lentries->SetTextSize(0.03);
  lentries->SetBorderSize(0);
  datarun->AddEntry(det1sum, Form("Integral: %.2f", det1sum->Integral()), "pl");
  datarun->AddEntry(det2sum, Form("Integral: %.2f", det2sum->Integral()), "pl");


  cout << "Prematch: " << prematch <<"\n";
  cout << "Postmatch: " << postmatch <<"\n";



  TLine *line = new TLine(10.0, pow(10,4), 10.0, pow(10,16));
  line->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1", "new", 1200, 600);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad11 = new TPad("pad11", "pad11", 0.0, 0.05, 0.5, 1.0);
     pad11->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad11->Draw();             // Draw the upper pad: pad1
     pad11->cd();               // pad1 becomes the current pad
     pad11->SetLogy();
     det1sum->SetXTitle("p_{T}^{hybrid}");
     //det1sum->SetMinimum(pow(10,12));
     //det1sum->SetMaximum(5*pow(10,13));
       det1sum->Draw("hist same");
       det2sum->Draw("hist same");
       datarun->Draw("same");
       //line->Draw("same");
     c1->cd();
   TPad *pad12 = new TPad("pad12", "pad12", 0.5, 0.05, 1.0, 1.0);
      pad12->Draw();
      pad12->cd();
      pad12->SetLogy();
        gen1spec->Draw("hist same");
        gen2spec->Draw("hist same");
 
  TCanvas *c2 = new TCanvas("c2", "leftovers", 500, 500);
    c2->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.05, 1.0, 1.0);
    pad2->Draw();
    pad2->cd();
    pad2->SetLogy();
    leftoverst->Draw("hist same");
    leftoversh->Draw("hist same");

  TCanvas *c3 = new TCanvas("c3", "diff", 500, 500);
    c3->cd();
    diff->Draw();

  TCanvas *c4 = new TCanvas("c4", "dpt", 500, 500);
    c4->cd();
    dpt->Draw();


  TCanvas *c5 = new TCanvas("c5", "", 500, 500);
    c5->cd();
    TH1D* divi = (TH1D*)det1sum->Clone("divi");
    divi->Divide(det2sum);
    divi->Draw("hist same");
    auto *ldivi = new TLegend(0.6, 0.8, 0.6, 0.8);
    ldivi->SetTextSize(0.03);
    ldivi->SetBorderSize(0);
    ldivi->AddEntry((TObject*)0, "pre-hvt/post-hvt", "");

}
  



void getbin(const char *analysisfile, double centl, double centr) {

    det1spec->Reset();
    det2spec->Reset();
    //dpt->Reset();
    leftoversh->Reset();
    leftoverst->Reset();  

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
    double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05, 0.05};
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0};
    

    int fillcount = 0;
    for (int i = 0; i < entries; i++) {
       T->GetEntry(i);
       if (centrality < centl || centrality > centr)  continue;
       //if (area < 0.4)  continue;
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
       det1spec->Fill(hybridpT, pThardFrac/extractFrac[scaleindex]);
       gen1spec->Fill(partpT, pThardFrac/extractFrac[scaleindex]);
       prematch++;
       if (hybridpT/partpT > 3.0)    {
          leftoversh->Fill(hybridpT, pThardFrac/extractFrac[scaleindex]);
          leftoverst->Fill(partpT, pThardFrac/extractFrac[scaleindex]);
          dpt->Fill(hybridpT-partpT, pThardFrac/extractFrac[scaleindex]);
          continue;}
       det2spec->Fill(hybridpT, pThardFrac/extractFrac[scaleindex]);
       gen2spec->Fill(partpT, pThardFrac/extractFrac[scaleindex]);
       postmatch++;

   }

    cout << "pT hard bin: " << pThardbin->GetMean() <<"\n";
    //cout << "	Scale Factor: " << pThardFrac << "*avg_events" << "\n";    

    totevents += n_event;
    int pTindex = pThardbin->GetMean() - 1;
   
    lefile->Close();

}


void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3) {

  vector<int> colors = {kRed+2, kRed, kOrange+10, kOrange-4, kYellow, kSpring+10, kSpring, kGreen-3, kGreen+2, kTeal+6, kTeal, kCyan-3, kAzure+7, kBlue+1, kBlue-2, kViolet+2, kViolet, kMagenta-4, kPink+6, kPink+10, kBlack};

 // colors[filenames1.size()] = kBlack;

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
