/*
 * EvaluateSys_RaaCaitie.C : Evaluates the systematic uncertainties for the Charged Jet RAA
 * Caitie Beattie <caitie.beattie@yale.edu>, adapted Hannah Bossi, adapted from Laura Havener
 * 08/03/2021
 */
TGraphAsymmErrors* MakeError(TGraphAsymmErrors* f_data, TGraphAsymmErrors* f_mc, TH1D* h1, int debug);
std::pair<TH1D*, TH1D*> CalculateDeltaCut(TH1D* h, TH1D* h_up, TH1D* h_down, bool sym);
TGraphAsymmErrors* MakeAsymmError(std::vector<TH1D*> histvec_up, std::vector<TH1D*> histvec_down, TH1D* h);
TH1D* ScaleSpec(TH1D* h1,  TH1D* heffnum, TH1D* hrecEff);
void applyR2(TH1D* NinPre, TH1D* NoutPre, double EPR);
void GetKinematicEfficiency(string fname, TH1D *KEout, TH1D *KEin, string fTag);
//void EvaluateSys_RaaCaitie();
//TH1D* getRuedigerRecEffML(Double_t radius);

double centl, centr;
double pi = 3.14159265359;
bool cent = false;
bool semi = false;
bool epup = false;
double R2 = 0.68;
vector<double> nomerge = {35, 40, 50, 65, 80, 100, 120}; 
std::vector<double> kBinsUnfolded = {10.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 190.0, 250.0};


//========================================
//========================================
//=========== main function ==============
//========================================
//========================================
void EvaluateSys_lotest(int lowbin = 0, int highbin = -1, const char* tag = "lo", int debug = 0, int R = 2)
{
   bool out = false;
   bool in = false;
   if (lowbin == 0) out = true;
   if (lowbin == 4) in = true;
   

  //=================================================================
  //============= List Files for Different Systematics ==============
  //=================================================================
      // nominal file
      std::string file_nom       = "UnfoldingData_2D_non_R02_qV0C_epV0A_3050_hi30_20_Nov14.root";
      std::string file_plus5     = "UnfoldingData_2D_non_R02_qV0C_epV0A_3050_hi30_hipT_20_Nov14.root";
      std::string file_min5      = "UnfoldingData_2D_non_R02_qV0C_epV0A_3050_hi30_lopT_20_Nov14.root";
      std::string file_eff       = "UnfoldingData_2D_non_R02_qV0C_epV0A_3050_hi30_trak_20_Nov14.root";
      std::string file_reweight  = "UnfoldingData_2D_non_R02_qV0C_epV0A_3050_hi30_rewe_20_Nov14.root";
  
      // file for the reconstruction efficiency
      //TH1D* hreceff = (TH1D*)getRuedigerRecEffML(0.4);
      TFile *receffile  = new TFile("ChargedJetRecEfficiencies.root");
      stringstream recname;
      recname << "RecEff_R0" << R << "0_5GeV";
      TH1D* hreceff0                = (TH1D*)receffile->Get(recname.str().c_str());
      hreceff0->SetName("rectemp");
      TH1D* hreceff = (TH1D*)hreceff0->Rebin(nomerge.size()-1, recname.str().c_str(), nomerge.data());
      std::string file_rec          = "pubDataCharged_ppR040.root";
      //std::string file_ruediger     = "../ResultsRuediger.root";
      // vary the number of iterations by +/-1
      std::string histname          = "Bayesian_Unfolded_6";
      std::string histnameplus1     = "Bayesian_Unfolded_7";
      std::string histnamemin1      = "Bayesian_Unfolded_5";
      std::string histnameKinNum    = "KinPos2";
      std::string histnameKinDenom  = "KinPre2";
      std::vector<TH1D*> histvec;
      std::vector<TH1D*> histvec_up;
      std::vector<TH1D*> histvec_down;

  

  //=================================================================
  //============= Determine Centrality Details ======================
  //=================================================================
     vector<const char*> letters(0);
     stringstream ss;
     const char* delim = "_";
     for (int i = 0; i < file_nom.length(); i++)    {
         const char *a = &file_nom[i];
         letters.push_back(a);
         if (strncmp(letters[i], delim, 1) == 0)    ss << " ";
         else ss << file_nom[i];
         }
        string tags;
           vector<string> words(0);
           while (ss >> tags)   words.push_back(tags);
           stringstream s2;
        const char* ccent = "010";
        const char* scent = "3050";
        for (int i = 2; i < words.size() - 1; i++) {
          s2 << words[i] << " ";
          const char *newword = words[i].c_str();
          if (strncmp(newword, ccent, 3) == 0)   cent = true;
          if (strncmp(newword, scent, 3) == 0)   semi = true;
        } 
       if (cent == true)  { 
         centl = 0.0;
         centr = 10.0;
        }
       if (semi == true)  {
         centl = 30.0;
         centr = 50.0;
        } 


  //===============================================================
  //============= Calculate Systematic Contributions ==============
  //===============================================================
      // ============ handle the nominal case =========================================
      TFile* fnom = new TFile(file_nom.c_str());
      //Get Kinematic efficiencies
      TH1D *hKinO_nom = new TH1D("hKinO_nom", "hKinO_nom", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_nom = new TH1D("hKinI_nom", "hKinI_nom", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_nom, hKinO_nom, hKinI_nom, "Nom");
      //retrieve 2D hist
      TH2D* h2_nom = (TH2D*)fnom->Get(histname.c_str()); 
      //project nominal in and out of plane
      TH1D* h1_in = (TH1D*)h2_nom->ProjectionX("h1_in", 4, -1); 
      TH1D* h1_out = (TH1D*)h2_nom->ProjectionX("h1_out", 0, 1); 
      TH1D* h1_icorr = (TH1D*)ScaleSpec(h1_in, hKinI_nom, hreceff);
      TH1D* h1_ocorr = (TH1D*)ScaleSpec(h1_out, hKinO_nom, hreceff);
      applyR2(h1_icorr, h1_ocorr, R2);
      //asign nominal and KE hists to in or out
      TH1D *h1_nom;
      if (lowbin == 0)    h1_nom = (TH1D*)h1_ocorr->Clone("h1_nom"); 
      if (lowbin == 4)    h1_nom = (TH1D*)h1_icorr->Clone("h1_nom"); 
      TH1D* nom = (TH1D*)h1_nom->Clone("h1_nom_copy");
      histvec.push_back(nom);
 
      // ============= handle the systematic for number of iterations +/- 1 ===========
      // +1 Iterations
      TH2D* h2_iterhigh = (TH2D*)fnom->Get(histnameplus1.c_str());
      //project nominal in and out of plane
      TH1D* h1_iterhigh_in = (TH1D*)h2_iterhigh->ProjectionX("h1_iterhigh_in", 4, -1); 
      TH1D* h1_iterhigh_out = (TH1D*)h2_iterhigh->ProjectionX("h1_iterhigh_out", 0, 1); 
      TH1D* h1_iterhigh_icorr = (TH1D*)ScaleSpec(h1_iterhigh_in, hKinI_nom, hreceff);
      TH1D* h1_iterhigh_ocorr = (TH1D*)ScaleSpec(h1_iterhigh_out, hKinO_nom, hreceff);
      applyR2(h1_iterhigh_icorr, h1_iterhigh_ocorr, R2);
      //asign nominal and KE hists to in or out
      TH1D *h1_iterhigh;
      if (lowbin == 0)    h1_iterhigh = (TH1D*)h1_iterhigh_ocorr->Clone("h1_iterhigh"); 
      if (lowbin == 4)    h1_iterhigh = (TH1D*)h1_iterhigh_icorr->Clone("h1_iterhigh"); 
      TH1D* iterhigh = (TH1D*)h1_iterhigh->Clone("h1_iterhigh_copy");
      histvec.push_back(iterhigh);
      // -1 Iterations
      TH2D* h2_iterlow = (TH2D*)fnom->Get(histnamemin1.c_str());
      //project nominal in and out of plane
      TH1D* h1_iterlow_in = (TH1D*)h2_iterlow->ProjectionX("h1_iterlow_in", 4, -1); 
      TH1D* h1_iterlow_out = (TH1D*)h2_iterlow->ProjectionX("h1_iterlow_out", 0, 1); 
      TH1D* h1_iterlow_icorr = (TH1D*)ScaleSpec(h1_iterlow_in, hKinI_nom, hreceff);
      TH1D* h1_iterlow_ocorr = (TH1D*)ScaleSpec(h1_iterlow_out, hKinO_nom, hreceff);
      applyR2(h1_iterlow_icorr, h1_iterlow_ocorr, R2);
      //asign nominal and KE hists to in or out
      TH1D *h1_iterlow;
      if (lowbin == 0)    h1_iterlow = (TH1D*)h1_iterlow_ocorr->Clone("h1_iterlow"); 
      if (lowbin == 4)    h1_iterlow = (TH1D*)h1_iterlow_icorr->Clone("h1_iterlow"); 
      TH1D* iterlow = (TH1D*)h1_iterlow->Clone("h1_iterlow_copy");
      histvec.push_back(iterlow);

      std::pair<TH1D*, TH1D*> iterUp = CalculateDeltaCut(nom, iterhigh, iterhigh, true);
      std::pair<TH1D*, TH1D*> iterDown = CalculateDeltaCut(nom, iterlow, iterlow, true);
      histvec_up.push_back(iterUp.first);
      histvec_down.push_back(iterUp.second);
      // these lines are commented out as we want to make this symmetric so we instead take the largest variation...
      histvec_up.push_back(iterDown.first);
      histvec_down.push_back(iterDown.second);
  
      // ============ handle the systematic for the tracking efficiency ================
      TFile* feff = new TFile(file_eff.c_str());
      //Get Kinematic efficiencies
      TH1D *hKinO_eff = new TH1D("hKinO_eff", "hKinO_eff", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_eff = new TH1D("hKinI_eff", "hKinI_eff", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_eff, hKinO_eff, hKinI_eff, "Eff");
      //retrieve 2D hist
      TH2D* h2_eff = (TH2D*)feff->Get(histname.c_str());
      //project nominal in and out of plane
      TH1D* h1_eff_in = (TH1D*)h2_eff->ProjectionX("h1_eff_in", 4, -1); 
      TH1D* h1_eff_out = (TH1D*)h2_eff->ProjectionX("h1_eff_out", 0, 1); 
      TH1D* h1_eff_icorr = (TH1D*)ScaleSpec(h1_eff_in, hKinI_eff, hreceff);
      TH1D* h1_eff_ocorr = (TH1D*)ScaleSpec(h1_eff_out, hKinO_eff, hreceff);
      applyR2(h1_eff_icorr, h1_eff_ocorr, R2);
      //asign nominal and KE hists to in or out
      TH1D *h1_eff;
      if (lowbin == 0)    h1_eff = (TH1D*)h1_eff_ocorr->Clone("h1_eff"); 
      if (lowbin == 4)    h1_eff = (TH1D*)h1_eff_icorr->Clone("h1_eff"); 
      TH1D* eff = (TH1D*)h1_eff->Clone("h1_eff_copy");
      histvec.push_back(eff);
      std::pair<TH1D*, TH1D*> efficiency = CalculateDeltaCut(nom, eff, eff, true);
      histvec_up.push_back(efficiency.first);
      histvec_down.push_back(efficiency.second);
  
      //========== handle the systematic for the changing the pt rec range +/5 ========
      //+5 GeV/c
      TFile* fplus5 = new TFile(file_plus5.c_str());
      //Get Kinematic efficiencies
      TH1D *hKinO_plus5 = new TH1D("hKinO_plus5", "hKinO_plus5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_plus5 = new TH1D("hKinI_plus5", "hKinI_plus5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_plus5, hKinO_plus5, hKinI_plus5, "Plus5");
      //retrieve 2D hist
      TH2D* h2_plus5 = (TH2D*)fplus5->Get(histname.c_str());
      //project nominal in and out of plane
      TH1D* h1_plus5_in = (TH1D*)h2_plus5->ProjectionX("h1_plus5_in", 4, -1); 
      TH1D* h1_plus5_out = (TH1D*)h2_plus5->ProjectionX("h1_plus5_out", 0, 1); 
      TH1D* h1_plus5_icorr = (TH1D*)ScaleSpec(h1_plus5_in, hKinI_plus5, hreceff);
      TH1D* h1_plus5_ocorr = (TH1D*)ScaleSpec(h1_plus5_out, hKinO_plus5, hreceff);
      applyR2(h1_plus5_icorr, h1_plus5_ocorr, R2);
      //asign nominal and KE hists to in or out
      TH1D *h1_plus5;
      if (lowbin == 0)    h1_plus5 = (TH1D*)h1_plus5_ocorr->Clone("h1_plus5"); 
      if (lowbin == 4)    h1_plus5 = (TH1D*)h1_plus5_icorr->Clone("h1_plus5"); 
      TH1D* plus5 = (TH1D*)h1_plus5->Clone("h1_plus5_copy");
      histvec.push_back(plus5);
      //-5 GeV/c
      TFile* fmin5 = new TFile(file_min5.c_str());
      //Get Kinematic efficiencies
      TH1D *hKinO_min5 = new TH1D("hKinO_min5", "hKinO_min5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_min5 = new TH1D("hKinI_min5", "hKinI_min5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_min5, hKinO_min5, hKinI_min5, "Min5");
      //retrieve 2D hist
      TH2D* h2_min5 = (TH2D*)fmin5->Get(histname.c_str());
      //project nominal in and out of plane
      TH1D* h1_min5_in = (TH1D*)h2_min5->ProjectionX("h1_min5_in", 4, -1); 
      TH1D* h1_min5_out = (TH1D*)h2_min5->ProjectionX("h1_min5_out", 0, 1); 
      TH1D* h1_min5_icorr = (TH1D*)ScaleSpec(h1_min5_in, hKinI_min5, hreceff);
      TH1D* h1_min5_ocorr = (TH1D*)ScaleSpec(h1_min5_out, hKinO_min5, hreceff);
      applyR2(h1_min5_icorr, h1_min5_ocorr, R2);
      //asign nominal and KE hists to in or out
      TH1D *h1_min5;
      if (lowbin == 0)    h1_min5 = (TH1D*)h1_min5_ocorr->Clone("h1_min5"); 
      if (lowbin == 4)    h1_min5 = (TH1D*)h1_min5_icorr->Clone("h1_min5"); 
      TH1D* min5 = (TH1D*)h1_min5->Clone("h1_min5_copy");
      histvec.push_back(min5);
      std::pair<TH1D*, TH1D*> recRangeUp = CalculateDeltaCut(nom, plus5, plus5, true);
      std::pair<TH1D*, TH1D*> recRangeDown = CalculateDeltaCut(nom, min5, min5, true);
      histvec_up.push_back(recRangeUp.first);
      histvec_down.push_back(recRangeUp.second);
      histvec_up.push_back(recRangeDown.first);
      histvec_down.push_back(recRangeDown.second);
  
      // ============== handle the systematic for reweighting the prior ===========                                                               
      TFile* frw = new TFile(file_reweight.c_str());
      //Get Kinematic efficiencies
      TH1D *hKinO_rw = new TH1D("hKinO_rw", "hKinO_rw", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_rw = new TH1D("hKinI_rw", "hKinI_rw", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_reweight, hKinO_rw, hKinI_rw, "Reweight");
      //retrieve 2D hist
      TH2D* h2_rw = (TH2D*)frw->Get(histname.c_str());
      //project nominal in and out of plane
      TH1D* h1_rw_in = (TH1D*)h2_rw->ProjectionX("h1_rw_in", 4, -1); 
      TH1D* h1_rw_out = (TH1D*)h2_rw->ProjectionX("h1_rw_out", 0, 1); 
      TH1D* h1_rw_icorr = (TH1D*)ScaleSpec(h1_rw_in, hKinI_rw, hreceff);
      TH1D* h1_rw_ocorr = (TH1D*)ScaleSpec(h1_rw_out, hKinO_rw, hreceff);
      applyR2(h1_rw_icorr, h1_rw_ocorr, R2);
      //asign nominal and KE hists to in or out
      TH1D *h1_rw;
      if (lowbin == 0)    h1_rw = (TH1D*)h1_rw_ocorr->Clone("h1_rw"); 
      if (lowbin == 4)    h1_rw = (TH1D*)h1_rw_icorr->Clone("h1_rw"); 
      TH1D* rw = (TH1D*)h1_rw->Clone("h1_rw_copy");
      histvec.push_back(rw);
      std::pair<TH1D*, TH1D*> weight = CalculateDeltaCut(nom, rw, rw, true);
      histvec_up.push_back(weight.first);
      histvec_down.push_back(weight.second);
  
      //========== handle the systematic for ep resolution ========                                       
      //retrieve 2D hist
      if (epup == false)  {
      TH2D* h2_gluon = (TH2D*)fnom->Get(histname.c_str()); 
      //project nominal in and out of plane
      TH1D* h1_gluon_in = (TH1D*)h2_gluon->ProjectionX("h1_gluon_in", 4, -1); 
      TH1D* h1_gluon_out = (TH1D*)h2_gluon->ProjectionX("h1_gluon_out", 0, 1); 
      TH1D* h1_gluon_icorr = (TH1D*)ScaleSpec(h1_gluon_in, hKinI_nom, hreceff);
      TH1D* h1_gluon_ocorr = (TH1D*)ScaleSpec(h1_gluon_out, hKinO_nom, hreceff);
      applyR2(h1_gluon_icorr, h1_gluon_ocorr, (R2-(0.02*R2)));
      //asign nominal and KE hists to in or out
      TH1D *h1_gluons;
      if (lowbin == 0)    h1_gluons = (TH1D*)h1_gluon_ocorr->Clone("h1_gluons"); 
      if (lowbin == 4)    h1_gluons = (TH1D*)h1_gluon_icorr->Clone("h1_gluons"); 
      TH1D* gluons = (TH1D*)h1_gluons->Clone("h1_gluons_copy");
      histvec.push_back(gluons);
      std::pair<TH1D*, TH1D*> fragmentationG = CalculateDeltaCut(nom, gluons, gluons, true);
      histvec_up.push_back(fragmentationG.first);
      histvec_down.push_back(fragmentationG.second);
      }


      TFile* frec = new TFile(file_rec.c_str());

      TGraphAsymmErrors* h1_nom_err = MakeAsymmError(histvec_up, histvec_down, nom);
      h1_nom_err->SetName("h1_nom_err");





  //================================================================
  //============= Write Hists to Root File =========================
  //================================================================
      std::stringstream name;
      if (epup == false) name << "preSystematics_R0" << R << "0_3050_" << (R/2) + 1 << "0_" << tag << "_test_Nov21.root";
      if (epup == true)  name << "preSystematics_R0" << R << "0_3050_" << (R/2) + 1 << "0_" << tag << "sm_Nov21.root";

      std::cout << histvec_down.size() << std::endl;
      TFile* fout = new TFile(name.str().c_str(), "recreate");
      fout->cd("/");
      for (int i = 0; i < histvec.size(); i++) histvec.at(i)->Write();
      for (int i = 0; i < histvec_up.size(); i++) histvec_up.at(i)->Write();
      for (int i = 0; i < histvec_down.size(); i++) histvec_down.at(i)->Write();
      h1_nom_err->Write();
      fout->Close();
      name.str("");

  //======================================================================
  //============= Plot Hists for Quality Control =========================
  //======================================================================
  TCanvas *c1 = new TCanvas("hists", "hists", 500, 1000);
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
     pad1->SetBottomMargin(0.05); // Upper and lower plot are joined
     pad1->Draw();
     pad1->SetLogy(); 
  pad1->cd();               // pad1 becomes the current pad   
   nom->GetXaxis()->SetRangeUser(40.0, 120.0);
   nom->Draw("same");
     nom->SetLineColor(kRed+1);
     nom->SetMarkerColor(kRed+1);
   iterhigh->Draw("same");
     iterhigh->SetLineColor(kOrange+7);
     iterhigh->SetMarkerColor(kOrange+7);
   iterlow->Draw("same");
     iterlow->SetLineColor(kYellow);
     iterlow->SetMarkerColor(kYellow);
   eff->Draw("same");
     eff->SetLineColor(kSpring);
     eff->SetMarkerColor(kSpring);
   rw->Draw("same");
     rw->SetLineColor(kGreen+3);
     rw->SetMarkerColor(kGreen+3);
   plus5->Draw("same");
     plus5->SetLineColor(kTeal);
     plus5->SetMarkerColor(kTeal);
     plus5->SetLineWidth(2);
   min5->Draw("same");
     min5->SetLineColor(kBlue);
     min5->SetMarkerColor(kBlue);
   nom->Draw("same");
  /* if (epup == false) {gluons->Draw("same");
     gluons->SetLineColor(kViolet);
     gluons->SetMarkerColor(kViolet);
   }  

  TLegend *leg = new TLegend(0.65, 0.65, 0.8, 0.8);
           leg->SetTextSize(0.035);
           leg->SetBorderSize(0);
           leg->AddEntry(nom, "nom", "pl");
           leg->AddEntry(iterhigh, "iterhigh", "pl");
           leg->AddEntry(iterlow, "iterlow", "pl");
           leg->AddEntry(eff, "eff", "pl");
           leg->AddEntry(rw, "rw", "pl");
           leg->AddEntry(plus5, "plus5", "pl");
           leg->AddEntry(min5, "min5", "pl");
           leg->AddEntry(gluons, "gluons", "pl");
           leg->Draw("same");

*/
}









//========================================================
//========================================================
//=========== ScaleSpec Function Definition ==============
//========================================================
//========================================================
TH1D* ScaleSpec(TH1D* h1_Orig,  TH1D* heffnum, TH1D* hrecEff){
  TH1D* h1 = (TH1D*) h1_Orig->Clone(Form("%s_copy", h1_Orig->GetName()));

  // correct the Pb--Pb Spectra for the kinematic efficiency and rec eff at once
  for(int i = 1; i < h1->GetNbinsX()+1; i++){
    float pT = h1->GetXaxis()->GetBinCenter(i);
    int binRec = hrecEff->GetXaxis()->FindBin(pT);
    h1->SetBinContent(i, h1->GetBinContent(i)*(1./(heffnum->GetBinContent(i)*hrecEff->GetBinContent(binRec))));
  }
  return h1;
}





//=============================================================
//=============================================================
//=========== MakeAsymmError Function Definition ==============
//=============================================================
//=============================================================
TGraphAsymmErrors* MakeAsymmError(std::vector<TH1D*> histvec_up, std::vector<TH1D*> histvec_down, TH1D* h)
{

  std::stringstream ss;
  ss << h->GetName() << "_tot_up";
  TH1D* h_up = (TH1D*)histvec_up.at(0)->Clone(ss.str().c_str());
  ss.str("");
  cout << h_up->GetName() << "\n";

  for (int i = 0; i < histvec_up.size(); i++)  {
      cout <<  histvec_up.at(i)->GetName() <<"\n";
      h_up->Add(histvec_up.at(i));
    }

  ss << h->GetName() << "_tot_down";
  TH1D* h_down = (TH1D*)histvec_down.at(0)->Clone(ss.str().c_str());
  ss.str("");
  for (int i = 0; i < histvec_down.size(); i++)
    {
      h_down->Add(histvec_down.at(i));
    }

  ss << h->GetName() << "_sys";
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  g->SetName(ss.str().c_str());
  ss.str("");

  for (int i = 1; i <= h_up->GetNbinsX(); i++)
    {
      double y = h->GetBinContent(i);
      double x = h->GetBinCenter(i);
      double dx = h->GetBinWidth(i)*0.5;

      double dy_low = sqrt(h_down->GetBinContent(i));
      double dy_high = sqrt(h_up->GetBinContent(i));
      g->SetPoint(i-1, x, y);
      g->SetPointError(i-1, dx, dx, dy_low, dy_high);
    }
  return g;
}





//=================================================================
//=================================================================
//=========== CalaculateDeltaCut Function Definition ==============
//=================================================================
//=================================================================
std::pair<TH1D*, TH1D*> CalculateDeltaCut(TH1D* h, TH1D* h_up, TH1D* h_down, bool sym)
{
  std::stringstream name;
  name << h_up->GetName() << "_delta_up";
  TH1D* h_delta_up = (TH1D*)h_up->Clone(name.str().c_str());
  h_delta_up->Reset();
  name.str("");
  if (!h_down) name << h_up->GetName() << "_delta_down";
  else name << h_down->GetName() << "_delta_down";
  TH1D* h_delta_down = (TH1D*)h_up->Clone(name.str().c_str());
  h_delta_down->Reset();
  cout << "Calculate Delta Cut for: " << name.str() <<"\n";
  
  for (int i = 1; i <= h->GetNbinsX(); i++)
    {
      float cont = h->GetBinContent(i);
      float delta_up = cont - h_up->GetBinContent(i);
      float delta_down = 0.;
      if (h_down) delta_down = cont - h_down->GetBinContent(i);
      float delta_low =  0.;
      float delta_high =  0.;
      if (!sym)   {
	  if ((delta_up < 0) &&  (delta_down < 0)){
	    delta_high = delta_up*delta_up + delta_down*delta_down;
	  }
	  else if ((delta_up < 0) && (delta_down > 0)) {
	    delta_high = delta_up*delta_up;
	    delta_low = delta_down*delta_down;
	  }
	  else if ((delta_up > 0) && (delta_down < 0)) {
	    delta_high = delta_down*delta_down;
	    delta_low = delta_up*delta_up;
	  }
	  else if ((delta_up > 0) && (delta_down > 0)) {
	    delta_low = delta_up*delta_up + delta_down*delta_down;
	  }
	  else;                                                                                                         
	}
      else {
	  if (h_up && h_down)    {
	      delta_low = delta_up*delta_up;
	      delta_high = delta_down*delta_down;
	  }
	  else  {
	      if (delta_up > 0) delta_low = delta_up*delta_up;
	      else             delta_high = delta_up*delta_up;
	  }  
	}
        
      h_delta_down->SetBinContent(i, delta_low);
      h_delta_up->SetBinContent(i, delta_high);
      //std::cout << "	spec bin: " << i << " : " << cont <<  std::endl;
      //std::cout << "	delta_low: " << delta_low << std::endl;
      //std::cout << "	delta_high: " << delta_high << std::endl;
    }

  return std::pair<TH1D*, TH1D*> (h_delta_up, h_delta_down);
}





//========================================================
//========================================================
//=========== MakeError Function Definition ==============
//========================================================
//========================================================
TGraphAsymmErrors* MakeError(TGraphAsymmErrors* g_PbPb, TGraphAsymmErrors* g_pp, TH1D* h1, int debug = 0)
{
  std::stringstream ss;
  ss << h1->GetName() << "_err";
  TGraphAsymmErrors* g = new TGraphAsymmErrors();
  g->SetName(ss.str().c_str());
  ss.str("");       
  cout << g->GetName() << ": combined errors for " << g_PbPb->GetName() << " + " << g_pp->GetName() <<"\n"; 
  for (int i = 1; i <= h1->GetNbinsX(); i++)      {
      //get baseline points from histogram
      float x = h1->GetBinCenter(i);
      float y = h1->GetBinContent(i);
      float dx = h1->GetBinWidth(i)*0.5;

      double y_PbPb = -1;
      double x_PbPb = -1;
      g_PbPb->GetPoint(i, x_PbPb, y_PbPb);
      double y_pp = -1;
      double x_pp = -1;
      if (g_pp) g_pp->GetPoint(i-1, x_pp, y_pp);
      //cout << "x: " << x << " y: "<< y <<" x_PbPb: " << x_PbPb <<" y_PbPb: " <<y_PbPb << " y_pp: " << y_pp << " x_pp: " << x_pp <<"\n";

      //get fractional error for PbPb
      double y_low_PbPb = g_PbPb->GetErrorYlow(i)/y_PbPb;
      double y_high_PbPb = g_PbPb->GetErrorYhigh(i)/y_PbPb;
      if (y_PbPb == 0.) {y_low_PbPb = 0.; y_high_PbPb = 0.;}

      //get fractional error for pp
      double y_low_pp = 0.;
      double y_high_pp = 0.;
      if (g_pp)    {
	  y_low_pp = g_pp->GetErrorYlow(i-1)/y_pp;
	  y_high_pp = g_pp->GetErrorYhigh(i-1)/y_pp;
	}
      if (y_pp == 0.) {y_low_pp = 0.; y_high_pp = 0.;}

      //combine and symmetrize your errors
      double dy_low = 0;
      double dy_high = 0;
      if (!g_pp) {
	dy_low = sqrt(y_low_PbPb*y_low_PbPb)*y;
	dy_high = sqrt(y_high_PbPb*y_high_PbPb)*y;
      }
      else {
	dy_low = sqrt(y_low_PbPb*y_low_PbPb + y_high_pp*y_high_pp)*y;
	dy_high = sqrt(y_high_PbPb*y_high_PbPb + y_low_pp*y_low_pp)*y;
      }

      //print out and set relevant info
      if (debug == 1) {
          cout << "	x: " << x << ", y: " << y << "\n";
          cout << "		dy_low:  " << dy_low << ", dy_high: " << dy_high <<"\n"; 
         }
      g->SetPoint(i-1, x, y);
      g->SetPointError(i-1, dx, dx, dy_low, dy_high);
    }
  return g;
}

//========================================================
//========================================================
//=========== EP Resolution Correction ===================
//========================================================
//========================================================

void applyR2(TH1D* NinPre, TH1D* NoutPre, double EPR) {
        TH1D* sumIO = (TH1D*)NinPre->Clone("sumIO");
        sumIO->Add(NoutPre);
        //calculate the v2 bin-by-bin
        for (int i = 1; i <= NinPre->GetNbinsX(); i++)   {
          double Nin = NinPre->GetBinContent(i);
          double Nout = NoutPre->GetBinContent(i);
          double v2 = (pi/(3*sqrt(3)))*(1./EPR)*(Nin-Nout)/(Nin+Nout);
          double Nrat = ((2*pi)-(6*sqrt(3)*v2))/((2*pi)+(6*sqrt(3)*v2));

        NinPre->SetBinContent(i, sumIO->GetBinContent(i)/(Nrat+1.0));
        NoutPre->SetBinContent(i, sumIO->GetBinContent(i)/(1.0+(1.0/Nrat)));

        }
 }



//========================================================
//========================================================
//=========== Get Kinematic Efficiencies =================
//========================================================
//========================================================


  void GetKinematicEfficiency(string fname, TH1D *KEout, TH1D *KEin, string fTag) {
      TFile* f = new TFile(fname.c_str());
      //Get Kinematic efficiencies
      TH2D* hKinNum2 = (TH2D*)f->Get("KinPos2");
      TH2D* hKinDenom2 = (TH2D*)f->Get("KinPre2");
      stringstream numName, denomName;
      numName << "Num" << fTag;
      denomName << "Den" << fTag;
      hKinNum2->SetName(numName.str().c_str());
      hKinDenom2->SetName(denomName.str().c_str());
         //out-of-plane
         numName << "1";
         denomName << "1";
         hKinNum2->GetYaxis()->SetRangeUser(0.0, sqrt(1)/2.0); //out-of-plane
         TH1D* hKinNumO = (TH1D*)hKinNum2->ProjectionX(numName.str().c_str());
         hKinDenom2->GetYaxis()->SetRangeUser(0.0, sqrt(1)/2.0);
         TH1D* hKinDenomO = (TH1D*)hKinDenom2->ProjectionX(denomName.str().c_str());
         hKinNumO->Divide(hKinDenomO);
         for (int i = 1; i <= KEout->GetNbinsX(); i++)  KEout->SetBinContent(i, hKinNumO->GetBinContent(i));
         //in-plane
         numName << "2";
         denomName << "2";
         hKinNum2->GetYaxis()->SetRangeUser(sqrt(3)/2.0, 1.0); //in-plane
         TH1D* hKinNumI = (TH1D*)hKinNum2->ProjectionX(numName.str().c_str());
         hKinDenom2->GetYaxis()->SetRangeUser(sqrt(3)/2.0, 1.0);
         TH1D* hKinDenomI = (TH1D*)hKinDenom2->ProjectionX(denomName.str().c_str());
         hKinNumI->Divide(hKinDenomI);
         for (int i = 1; i <= hKinNumI->GetNbinsX(); i++)  KEin->SetBinContent(i, hKinNumI->GetBinContent(i));
   }

