/*
 * EvaluateSys_RaaCaitie.C : Evaluates the systematic uncertainties for the Charged Jet RAA
 * Caitie Beattie <caitie.beattie@yale.edu>, adapted Hannah Bossi, adapted from Laura Havener
 * 08/03/2021
 */
TGraphAsymmErrors* MakeError(TGraphAsymmErrors* f_data, TGraphAsymmErrors* f_mc, TH1D* h1, int debug);
std::pair<TH1D*, TH1D*> CalculateDeltaCut(TH1D* h, TH1D* h_up, TH1D* h_down, bool sym);
TGraphAsymmErrors* MakeAsymmError(std::vector<TH1D*> histvec_up, std::vector<TH1D*> histvec_down, TH1D* h);
TH1D* ScaleSpec(TH1D* h1,  TH1D* heffnum, TH1D* hrecEff);
//void EvaluateSys_RaaCaitie();
//TH1D* getRuedigerRecEffML(Double_t radius);

double centl, centr;
bool cent = false;
bool semi = false;
//========================================
//========================================
//=========== main function ==============
//========================================
//========================================
void EvaluateSys_2DCaitie(int debug = 0)
{
  //=================================================================
  //============= List Files for Different Systematics ==============
  //=================================================================
      // nominal file
      std::string file_nom       = "UnfoldingData_2D_non_R02_3050_lo30_Feb4.root";
      std::string file_plus5     = "UnfoldingData_2D_non_R02_3050_lo30_hipT_Feb4.root";
      std::string file_min5      = "UnfoldingData_2D_non_R02_3050_lo30_lopT_Feb4.root";
      std::string file_eff       = "UnfoldingData_2D_non_R02_3050_lo30_trak_Feb4.root";
      std::string file_reweight  = "UnfoldingData_2D_non_R02_3050_lo30_rewe_Feb1.root";
   
      // file for the reconstruction efficiency
      //TH1D* hreceff = (TH1D*)getRuedigerRecEffML(0.4);
      TFile *receffile  = new TFile("ChargedJetRecEfficiencies.root");
      TH1D* hreceff                 = (TH1D*)receffile->Get("RecEff_R020_5GeV");
      std::string file_rec          = "pubDataCharged_ppR020.root";
      //std::string file_ruediger     = "../ResultsRuediger.root";
      // vary the number of iterations by +/-1
      std::string histname          = "Bayesian_Unfolded_6";
      std::string histnameplus1     = "Bayesian_Unfolded_7";
      std::string histnamemin1      = "Bayesian_Unfolded_5";
      std::string histnameKinNum    = "KinPos";
      std::string histnameKinDenom  = "KinPre";
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
      //fnom->ls();
      //retrieve 2D hist
      TH2D* h2_nom = (TH2D*)fnom->Get(histname.c_str()); 
      //project in-plane
      TH1D* h1_nomi = (TH1D*)h2_nom->ProjectionX("h1_nomi", 4, -1);
      //project out-of-plane
      TH1D* h1_nomo = (TH1D*)h2_nom->ProjectionX("h1_nomo", 0, 1);
      //Divide hists to get nominal
      TH1D* h1_nom = (TH1D*)h1_nomo->Clone("h1_nom");
      h1_nom->Divide(h1_nomi);
      //Get Kinematic efficiencies
      TH1D* hKinNum_nom = (TH1D*)fnom->Get(histnameKinNum.c_str());
      hKinNum_nom->SetName("hKinNum_nom");
      TH1D* hKinDenom_nom = (TH1D*)fnom->Get(histnameKinDenom.c_str());
      hKinDenom_nom->SetName("hKinDenom_nom");
      hKinNum_nom->Divide(hKinDenom_nom);
      TH1D* nom = (TH1D*)ScaleSpec(h1_nom, hKinNum_nom, hreceff);
      histvec.push_back(nom);
 
      // ============= handle the systematic for number of iterations +/- 1 ===========
      // +1 Iterations
      TH2D* h2_iterhigh = (TH2D*)fnom->Get(histnameplus1.c_str());
      //project in-plane
      TH1D* h1_iterhighi = (TH1D*)h2_iterhigh->ProjectionX("h1_iterhighi", 4, -1);
      //project out-of-plane
      TH1D* h1_iterhigho = (TH1D*)h2_iterhigh->ProjectionX("h1_iterhigho", 0, 1);
      //Divide hists to get nominal
      TH1D* h1_iterhigh = (TH1D*)h1_iterhigho->Clone("h1_iterhigh");
      h1_iterhigh->Divide(h1_iterhighi);
      TH1D* iterhigh = (TH1D*)ScaleSpec(h1_iterhigh, hKinNum_nom, hreceff);
      histvec.push_back(iterhigh);
      // -1 Iterations
      TH2D* h2_iterlow = (TH2D*)fnom->Get(histnamemin1.c_str());
      TH1D* h1_iterlowi = (TH1D*)h2_iterlow->ProjectionX("h1_iterlowi", 4, -1);
      TH1D* h1_iterlowo = (TH1D*)h2_iterlow->ProjectionX("h1_iterlowo", 0, 1);
      TH1D* h1_iterlow = (TH1D*)h1_iterlowo->Clone("h1_iterlow");
      h1_iterlow->Divide(h1_iterlowi);
      TH1D* iterlow = (TH1D*)ScaleSpec(h1_iterlow, hKinNum_nom, hreceff);
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
      TH2D* h2_eff = (TH2D*)feff->Get(histname.c_str());
      //project in-plane
      TH1D* h1_effi = (TH1D*)h2_eff->ProjectionX("h1_effi", 4, -1);
      //project out-of-plane
      TH1D* h1_effo = (TH1D*)h2_eff->ProjectionX("h1_effo", 0, 1);
      TH1D* h1_eff = (TH1D*)h1_effo->Clone("h1_eff");
      h1_eff->Divide(h1_effi);
      TH1D* hKinNum_eff = (TH1D*)feff->Get(histnameKinNum.c_str());
      hKinNum_eff->SetName("hKinNum_eff");
      TH1D* hKinDenom_eff = (TH1D*)feff->Get(histnameKinDenom.c_str());
      hKinDenom_eff->SetName("hKinDenom_eff");
      hKinNum_eff->Divide(hKinDenom_eff);
      //TFile* fr = new TFile(file_ruediger.c_str());
      //TDirectoryFile* refFiles = (TDirectoryFile*)fr->Get("JetRecEfficiencies");
      //TH1D* hreceffRed = (TH1D*)refFiles->Get("hResponse_DET_LHC15o_R040_LeadingTrackBias0GeV_TrackEff096");
      TH1D* eff = (TH1D*)ScaleSpec(h1_eff, hKinNum_eff, hreceff);
      histvec.push_back(eff);
      std::pair<TH1D*, TH1D*> efficiency = CalculateDeltaCut(nom, eff, eff, true);
      histvec_up.push_back(efficiency.first);
      histvec_down.push_back(efficiency.second);
  
      //========== handle the systematic for the changing the pt rec range +/5 ========
      //+5 GeV/c
      TFile* fplus5 = new TFile(file_plus5.c_str());
      TH2D* h2_plus5 = (TH2D*)fplus5->Get(histname.c_str());
      //project in-plane
      TH1D* h1_plus5i = (TH1D*)h2_plus5->ProjectionX("h1_plus5i", 4, -1);
      //project out-of-plane
      TH1D* h1_plus5o = (TH1D*)h2_plus5->ProjectionX("h1_plus5o", 0, 1);
      //Divide hists to get nominal
      TH1D* h1_plus5 = (TH1D*)h1_plus5o->Clone("h1_plus5");
      h1_plus5->Divide(h1_plus5i);
      TH1D* hKinNum_plus5 = (TH1D*)fplus5->Get(histnameKinNum.c_str());
      hKinNum_plus5->SetName("hKinNum_plus5");
      TH1D* hKinDenom_plus5 = (TH1D*)fplus5->Get(histnameKinDenom.c_str());
      hKinDenom_plus5->SetName("hKinDenom_plus5");
      hKinNum_plus5->Divide(hKinDenom_plus5);
      TH1D* plus5 = (TH1D*)ScaleSpec(h1_plus5, hKinNum_plus5, hreceff);
      histvec.push_back(plus5);
      //-5 GeV/c
      TFile* fmin5 = new TFile(file_min5.c_str());
      TH2D* h2_min5 = (TH2D*)fmin5->Get(histname.c_str());
      //project in-plane
      TH1D* h1_min5i = (TH1D*)h2_min5->ProjectionX("h1_min5i", 4, -1);
      //project out-of-plane
      TH1D* h1_min5o = (TH1D*)h2_min5->ProjectionX("h1_min5o", 0, 1);
      TH1D* h1_min5 = (TH1D*)h1_min5o->Clone("h1_min5");
      h1_min5->Divide(h1_min5i);
      TH1D* hKinNum_min5 = (TH1D*)fmin5->Get(histnameKinNum.c_str());
      hKinNum_min5->SetName("hKinNum_min5");
      TH1D* hKinDenom_min5 = (TH1D*)fmin5->Get(histnameKinDenom.c_str());
      hKinDenom_min5->SetName("hKinDenom_min5");
      hKinNum_min5->Divide(hKinDenom_min5);
      TH1D* min5 = (TH1D*)ScaleSpec(h1_min5, hKinNum_min5, hreceff);
      histvec.push_back(min5);
      std::pair<TH1D*, TH1D*> recRangeUp = CalculateDeltaCut(nom, plus5, plus5, true);
      std::pair<TH1D*, TH1D*> recRangeDown = CalculateDeltaCut(nom, min5, min5, true);
      histvec_up.push_back(recRangeUp.first);
      histvec_down.push_back(recRangeUp.second);
      histvec_up.push_back(recRangeDown.first);
      histvec_down.push_back(recRangeDown.second);
  
      // ============== handle the systematic for reweighting the prior ===========                                                               
      TFile* frw = new TFile(file_reweight.c_str());
      TH2D* h2_rw = (TH2D*)frw->Get(histname.c_str());
      TH1D* h1_rwi = (TH1D*)h2_rw->ProjectionX("h1_rwi", 4, -1);
      TH1D* h1_rwo = (TH1D*)h2_rw->ProjectionX("h1_rwo", 0, 1);
      TH1D* h1_rw = (TH1D*)h1_rwo->Clone("h1_rw");
      h1_rw->Divide(h1_rwi);
      TH1D* hKinNum_rw = (TH1D*)frw->Get(histnameKinNum.c_str());
      hKinNum_rw->SetName("hKinNum_rw");
      TH1D* hKinDenom_rw = (TH1D*)frw->Get(histnameKinDenom.c_str());
      hKinDenom_rw->SetName("hKinDenom_rw");
      hKinNum_rw->Divide(hKinDenom_rw);
      TH1D* rw = (TH1D*)ScaleSpec(h1_rw, hKinNum_rw, hreceff);
      histvec.push_back(rw);
      std::pair<TH1D*, TH1D*> weight = CalculateDeltaCut(nom, rw, rw, true);
      histvec_up.push_back(weight.first);
      histvec_down.push_back(weight.second);
  
      /*
      //========== handle the systematic for q/g ========                                                                                                    TFile* fgluons = new TFile(file_gluons.c_str());
      TH1D* h1_gluons = (TH1D*)fgluons->Get(histname.c_str());
      h1_gluons->SetName("h1_gluons");
      TH1D* hKinNum_gluons = (TH1D*)fgluons->Get(histnameKinNum.c_str());
      hKinNum_gluons->SetName("hKinNum_gluons");
      TH1D* hKinDenom_gluons = (TH1D*)fgluons->Get(histnameKinDenom.c_str());
      hKinDenom_gluons->SetName("hKinDenom_gluons");
      hKinNum_gluons->Divide(hKinDenom_gluons);
      TH1D* gluons = (TH1D*)ScaleSpec(h1_gluons, hKinNum_gluons, hreceff);
      histvec.push_back(gluons);
      TFile* fquarks = new TFile(file_quarks.c_str());
      TH1D* h1_quarks = (TH1D*)fquarks->Get(histname.c_str());
      h1_quarks->SetName("h1_quarks");
      TH1D* hKinNum_quarks = (TH1D*)fquarks->Get(histnameKinNum.c_str());
      hKinNum_quarks->SetName("hKinNum_quarks");
      TH1D* hKinDenom_quarks = (TH1D*)fquarks->Get(histnameKinDenom.c_str());
      hKinDenom_quarks->SetName("hKinDenom_quarks");
      hKinNum_quarks->Divide(hKinDenom_quarks);
      TH1D* quarks = (TH1D*)ScaleSpec(h1_quarks, hKinNum_quarks, hreceff);
      //histvec.push_back(quarks);                                                                                                                           std::pair<TH1D*, TH1D*> fragmentationG = CalculateDeltaCut(nom, gluons, gluons, true);
      histvec_up.push_back(fragmentationG.first);
      histvec_down.push_back(fragmentationG.second);
      std::pair<TH1D*, TH1D*> fragmentationQ = CalculateDeltaCut(nom, quarks, quarks, true);
      //histvec_up.push_back(fragmentationQ.first);                                                                                                 
      //histvec_down.push_back(fragmentationQ.second);  
      */

      TFile* frec = new TFile(file_rec.c_str());

      TGraphAsymmErrors* h1_nom_err = MakeAsymmError(histvec_up, histvec_down, nom);
      h1_nom_err->SetName("h1_nom_err");



  //===============================================================
  //============= Retrieve pp Systematics =========================
  //===============================================================
      std::vector<double> binning = {25, 30, 40, 50, 60, 70, 85, 100};  

      //Values from Publication (https://www.hepdata.net/record/ins1733689) 
      double pp_centers[7] = {27.5, 35, 45, 55, 65, 77.5, 92.5};
      double pp_values[7] = {0.000412768, 0.000139049, 3.81726*pow(10, -5), 1.33766*pow(10, -5), 5.4954*pow(10, -6), 2.18345*pow(10, -6), 8.16514*pow(10, -7)};
      double pp_nul[7] = {2.5, 5.0, 5.0, 5.0, 5.0, 7.5, 7.5};
      double pp_err[7] = {2.51789*pow(10,-5), 9.59436*pow(10,-6), 2.97746*pow(10,-6), 1.17714*pow(10,-6), 5.33054*pow(10,-7), 2.37996*pow(10,-7), 9.96147*pow(10,-8)};

      //Generate Histogram
      double pTedgepp[8]  = {25, 30, 40, 50, 60, 70, 85, 100};
      TH1D* hist_pp = new TH1D("refpp_R040_ltb5", "", 7, pTedgepp);
      for (int i = 0; i < 7; i++)    hist_pp->Fill(pTedgepp[i] + 0.01, pp_values[i]); 

      //Generate TGraphAsymmErrors
      auto h1_pp_err = new TGraphAsymmErrors(7, pp_centers, pp_values, pp_nul, pp_nul, pp_err, pp_err);
      
      //calculate RAA
      TH1D* numerator = (TH1D*)nom->Rebin(binning.size()-1, "raa", binning.data());
      numerator->Divide(hist_pp);

      //Write raa_err to file
      TGraphAsymmErrors* raa_err = MakeError(h1_nom_err, h1_pp_err, numerator, debug);



  //================================================================
  //============= Write Hists to Root File =========================
  //================================================================
      std::stringstream name;
      name << "preSystematics_R020_3050_lo_Feb1.root";

      std::cout << histvec_down.size() << std::endl;
      TFile* fout = new TFile(name.str().c_str(), "recreate");
      fout->cd("/");
      for (int i = 0; i < histvec.size(); i++) histvec.at(i)->Write();
      for (int i = 0; i < histvec_up.size(); i++) histvec_up.at(i)->Write();
      for (int i = 0; i < histvec_down.size(); i++) histvec_down.at(i)->Write();
      h1_nom_err->Write();
      numerator->Write();
      raa_err->Write();
      fout->Close();
      name.str("");

  //======================================================================
  //============= Plot Hists for Quality Control =========================
  //======================================================================
  TCanvas *c1 = new TCanvas("hists", "hists", 1000, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
     pad1->SetBottomMargin(0.15); // Upper and lower plot are joined
     pad1->Draw();
   pad1->cd();               // pad1 becomes the current pad   
   h1_nom->Draw("same");
     h1_nom->SetLineColor(kRed);
   h1_iterhigh->Draw("same");
     h1_iterhigh->SetLineColor(kOrange);
   h1_iterlow->Draw("same");
     h1_iterlow->SetLineColor(kYellow);
   h1_eff->Draw("same");
     h1_eff->SetLineColor(kGreen);
   h1_rw->Draw("same");
     h1_rw->SetLineColor(kPink);
   h1_plus5->Draw("same");
     h1_plus5->SetLineColor(kBlue);
     h1_plus5->SetLineWidth(2);
   h1_min5->Draw("same");
     h1_min5->SetLineColor(kViolet);
  


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

  // now perform the numeric scaling
    Double_t Taa; //mb^{-1}
    if (cent == true)   Taa = 23.3;
    if (semi == true)   Taa = 3.9; 
    // Number of selected events in the desired centrality range (30-50%)
    //Double_t Nevents    = 4619963.0; 
    string specfilename = "../ESE/AnalysisResults7812.root";   //central
    //string specfilename = "../Spectra/AnalysisResults7283.root";   //semi-central


    const char *lespecfile = specfilename.c_str(); 
    //Load File/Tree into system
    TFile *specfile = new TFile(lespecfile);

      TDirectoryFile* tdf = (TDirectoryFile*)specfile->Get("ChargedJetsHadronCF");
      AliEmcalList *tlist = (AliEmcalList*)tdf->Get("AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet_histos"); 
      TH1F *centhist = (TH1F*)tlist->FindObject("fHistCentrality");
      centhist->GetXaxis()->SetRangeUser(centl, centr);
      double Nevents = centhist->Integral();

  //perform the final scaling
  h1->Scale(1./(2*(0.9-0.4)));
  h1->Scale(1./Taa); 
  h1->Scale(1./Nevents, "width");
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
      std::cout << "	spec bin: " << i << " : " << cont <<  std::endl;
      std::cout << "	delta_low: " << delta_low << std::endl;
      std::cout << "	delta_high: " << delta_high << std::endl;
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


//==================================================================================                                                                                                    
/*
TH1D* getRuedigerRecEffML(Double_t radius){
  TFile* prevResults = TFile::Open("~/ResultsRuediger.root");
  TDirectoryFile* respDirect = (TDirectoryFile*)prevResults->Get("ResponseMatrices");
  TDirectoryFile* direct;
  if (radius == 0.2){
    direct = (TDirectoryFile*)prevResults->Get("LHC15o_R020_NeuralNetwork_TrainedOnLHC18b8_DetLevel_PbPb_FlatMult_Simplified_000_010_LeadingTrackBias0GeV");
  }
  else if (radius == 0.4){
    direct = (TDirectoryFile*)prevResults->Get("LHC15o_R040_NeuralNetwork_TrainedOnLHC18b8_DetLevel_PbPb_FlatMult_Simplified_000_010_LeadingTrackBias0GeV");
  }
  else if (radius == 0.6){
    direct = (TDirectoryFile*)prevResults->Get("LHC15o_R060_NeuralNetwork_TrainedOnLHC18b8_DetLevel_PbPb_FlatMult_Simplified_000_010_LeadingTrackBias0GeV");
  }
  else{
    std::cout << "Error: Did not recognize jet radius." << std::endl;
  }
  TDirectoryFile* unfolding  = (TDirectoryFile*)direct->Get("UnfoldingPlots");
  TH1D* receffFullML = (TH1D*)unfolding->Get("hJetRecEfficiency");
  return receffFullML;
}*/
//=================================================================================  
