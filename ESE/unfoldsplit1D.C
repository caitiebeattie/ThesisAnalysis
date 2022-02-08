/*#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"
*/
using namespace std;

double pi = 3.14159265359;
double lbound = 30.0;
double rbound = 120.0;

void setcolor(TH1F* h, int kcolor);

vector<const char*> names = {"i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12","i13", "i14", "i15"};


bool cent = false;
bool semi = false; 
bool r020 = false;
bool r040 = false;

void unfoldsplit1D (int iteration, string rootfilename, bool zoom = false) {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

   /*TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(15);
   axis->Draw("same");
   */

    const char *lerootfile = rootfilename.c_str(); 
    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);

   // RooUnfoldResponse* responsefolder = (RooUnfoldResponse*)lefile->Get("smeared_true");
   // TH1D* responseraw = (TH1D*)responsefolder->Hmeasured(); 
   // TH2F* responsesplit = (TH2F*)responsefolder->Hresponse();

vector<TH1F*> unfoldHists1D(0);
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_1"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_2"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_3"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_4"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_5"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_6"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_7"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_8"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_9"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_10"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_11"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_12"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_13"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_14"));
unfoldHists1D.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_15"));
 

    for (int i = 0; i < unfoldHists1D.size(); i++)  {
        unfoldHists1D[i]->Scale(1.0, "width");
     }

     TH1F* truth = (TH1F*)lefile->Get("true1");
     truth->Scale(1.0, "width");
     TH1F* raw = (TH1F*)lefile->Get("rawdata1");
     raw->Scale(1.0, "width");


vector<TH1F*> ratiohists(0);
for (int i = 0; i < unfoldHists1D.size(); i++)   {
      ratiohists.push_back((TH1F*)unfoldHists1D[i]->Clone(names[i]));
      ratiohists[i]->Divide(unfoldHists1D[iteration]);
     // ratiohists[i]->GetXaxis()->SetRangeUser(lbound, rbound); 
      }

TH1F* trivial  = (TH1F*)unfoldHists1D[iteration]->Clone("trivial");
trivial->Divide(truth);
//trivial->Scale(1.0, "width");

vector<int> colorlist = {kRed-7, kOrange+1, kOrange-2, kYellow-7, kSpring+5, kGreen-6, kGreen-2, kTeal+2, kTeal-5, kCyan-3, kCyan+2, kAzure-9, kAzure-4, kBlue-4, kBlue-6, kViolet+6, kViolet-8, kMagenta-8};
for (int i = 0; i < unfoldHists1D.size(); i++)    {
     setcolor(unfoldHists1D[i], colorlist[i]);
     setcolor(ratiohists[i], colorlist[i]);}
     unfoldHists1D[iteration]->SetLineWidth(3);
  
  TLine *line = new TLine(lbound, 1, rbound, 1);
  line->SetLineStyle(2);

   //Get File tag for legend
   vector<const char*> letters(0);
   stringstream ss;
   const char* delim = "_";
   for (int i = 0; i < rootfilename.length(); i++)    {
       const char *a = &rootfilename[i];
       letters.push_back(a);
       if (strncmp(letters[i], delim, 1) == 0)    ss << " ";
       else ss << rootfilename[i];
       }
    string tags;
    vector<string> words(0);
    while (ss >> tags)  {
       words.push_back(tags);
    }  
    stringstream s2;
    const char* ccent = "010";
    const char* scent = "3050";
    const char* r02 = "R02";
    const char* r04 = "R04";
    for (int i = 2; i < words.size() - 1; i++) {
        s2 << words[i] << " ";
        const char *newword = words[i].c_str();
        if (strncmp(newword, ccent, 3) == 0)   cent = true;
        if (strncmp(newword, scent, 3) == 0)   semi = true;
        if (strncmp(newword, r02, 3) == 0)     r020 = true;
        if (strncmp(newword, r04, 3) == 0)     r040 = true;
        }


    auto fileleg = new TLegend(0.4, 0.2, 0.85, 0.45);
         fileleg->SetTextSize(0.055);
         fileleg->SetBorderSize(0);
         fileleg->AddEntry((TObject*)0, s2.str().c_str(), "");


  auto itleg = new TLegend(0.6, 0.3, 0.85, 0.85);
    itleg->SetTextSize(0.045);
    itleg->SetBorderSize(0); 
    for (int i = 0; i < unfoldHists1D.size(); i++)    itleg->AddEntry(unfoldHists1D[i], Form("%d", i+1), "pl");

  auto stab = new TLegend(0.15, 0.65, 0.4, 0.85);
    stab->SetTextSize(0.05);
    stab->SetBorderSize(0); 
    stab->AddEntry((TObject*)0, "Stability", "");
    if (cent == true)  stab->AddEntry((TObject*)0, "R = 0.2, 0-10%", "");
    if (semi == true)  stab->AddEntry((TObject*)0, "R = 0.2, 30-50%", "");
  auto split = new TLegend(0.15, 0.65, 0.4, 0.85);
    split->SetTextSize(0.05);
    split->SetBorderSize(0); 
    if (cent == true)  split->AddEntry((TObject*)0, "R = 0.2, 0-10%", "");
    if (semi == true)  split->AddEntry((TObject*)0, "R = 0.2, 30-50%", "");

  float ymin, ymax;
  if (zoom == false)  {
      ymin = 0.5;
      ymax = 1.5;}
  if (zoom == true)   {
      ymin = 0.8;
      ymax = 1.2;}


  TCanvas *c1 = new TCanvas("Stability Test", lerootfile, 800, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad1 = new TPad("pad1", "", 0.0, 0.4, 0.5, 1.0);
   TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 0.5, 0.4);
     pad1->SetBottomMargin(0.0); // Upper and lower plot are joined
     pad1->SetRightMargin(0.0);
     pad2->SetTopMargin(0.0);
     pad2->SetRightMargin(0.0);
     pad2->SetBottomMargin(0.2);
     pad1->Draw();
     pad2->Draw();
   pad1->cd();               // pad1 becomes the current pad   
   //pad1->SetLogx();
   pad1->SetLogy();
   double topmax = unfoldHists1D[0]->GetBinContent(3)*pow(10,3);
   double topmin = unfoldHists1D[0]->GetBinContent(9)*pow(10,-2);
   unfoldHists1D[0]->SetMaximum(topmax);
   unfoldHists1D[0]->SetMinimum(topmin);
   //unfoldHists1D[0]->SetMinimum(pow(10,-9));
   unfoldHists1D[0]->SetTitle("");
   for (int i = 0; i < 15; i++)         unfoldHists1D[i]->Draw("same");
   itleg->Draw("same");
   stab->Draw("same");
   pad2->cd();
   ratiohists[0]->GetXaxis()->SetLabelSize(0.05);
   ratiohists[0]->GetYaxis()->SetLabelSize(0.05);
   ratiohists[0]->SetTitle("");
   ratiohists[0]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
   ratiohists[0]->GetXaxis()->SetTitleSize(0.09);
   ratiohists[0]->SetMinimum(ymin);
   ratiohists[0]->SetMaximum(ymax);
   for (int i = 0; i < unfoldHists1D.size(); i++)         ratiohists[i]->Draw("same");
   line->Draw("same");  

  //TCanvas *c3 = new TCanvas("Split Test", "Split Test", 500, 500);
   //c3->cd();
   gStyle->SetOptStat(0);
   c1->cd();
   TPad *pad3 = new TPad("pad3", "", 0.5, 0.4, 1.0, 1.0);
   TPad *pad4 = new TPad("pad4", "", 0.5, 0.05, 1.0, 0.4);
   TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.85);
     setcolor(unfoldHists1D[iteration], kBlue);
     setcolor(truth, kGreen);
     setcolor(raw, kRed);
     setcolor(trivial, kBlack);
   leg->AddEntry(unfoldHists1D[iteration], "unfolded", "pl");
   leg->AddEntry(truth, "truth", "pl");
   leg->AddEntry(raw, "raw", "pl");
   leg->SetBorderSize(0);
     pad3->SetBottomMargin(0.0); // Upper and lower plot are joined
     pad4->SetTopMargin(0.0);
     pad4->SetBottomMargin(0.2);
     pad3->SetLeftMargin(0.0);
     pad4->SetLeftMargin(0.0);
     pad3->Draw();
     pad4->Draw();
   pad3->cd();               // pad1 becomes the current pad   
     //pad1->SetLogx();
     pad3->SetLogy();
     unfoldHists1D[iteration]->GetXaxis()->SetRangeUser(lbound, rbound);
     //unfoldHists1D[iteration]->SetMaximum(pow(10,-1));
     //unfoldHists1D[iteration]->SetMinimum(pow(10,-9));
     unfoldHists1D[iteration]->SetTitle("");
     unfoldHists1D[iteration]->SetMaximum(topmax);
     unfoldHists1D[iteration]->SetMinimum(topmin);
     unfoldHists1D[iteration]->Draw("same");      
     truth->GetXaxis()->SetRangeUser(lbound, rbound);
     truth->Draw("same");
     raw->GetXaxis()->SetRangeUser(lbound, rbound);
     raw->Draw("same");
     leg->Draw("same"); 
     split->Draw("same"); 
   pad4->cd();
     trivial->SetTitle("");
     trivial->GetXaxis()->SetRangeUser(lbound, rbound);
     trivial->SetMinimum(ymin);
     trivial->SetMaximum(ymax);
     trivial->Draw("same");
     trivial->SetXTitle("#it{p}_{T} (GeV/#it{c})");
     trivial->GetXaxis()->SetTitleSize(0.09);
     trivial->GetXaxis()->SetLabelSize(0.05);
     trivial->GetYaxis()->SetLabelSize(0.05);
     line->Draw("same");
     fileleg->Draw("same");


}


void setcolor(TH1F* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
}

