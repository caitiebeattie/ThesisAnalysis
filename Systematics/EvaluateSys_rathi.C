/*
 * EvaluateSys_RaaCaitie.C : Evaluates the systematic uncertainties for the Charged Jet RAA
 * Caitie Beattie <caitie.beattie@yale.edu>, adapted Hannah Bossi, adapted from Laura Havener
 * 08/03/2021
 */

//declare functions
TGraphAsymmErrors* MakeError(TGraphAsymmErrors* f_data, TGraphAsymmErrors* f_mc, TH1D* h1, int debug);
std::pair<TH1D*, TH1D*> CalculateDeltaCut(TH1D* h, TH1D* h_up, TH1D* h_down, bool sym);
TGraphAsymmErrors* MakeAsymmError(std::vector<TH1D*> histvec_up, std::vector<TH1D*> histvec_down, TH1D* h);
TH1D* ScaleSpec(TH1D* h1,  TH1D* heffnum, TH1D* hrecEff);
void applyR2(TH1D* Nratio, TH1D* NinPre, TH1D* NoutPre, double EPR);
void GetKinematicEfficiency(string fname, TH1D *KEout, TH1D *KEin, string fTag);
//void EvaluateSys_RaaCaitie();
//TH1D* getRuedigerRecEffML(Double_t radius);


//useful variables/constants
double centl, centr;
bool cent = false;
bool semi = false;
bool r02 = false;
bool r04 = false;
double pi = 3.14159265;
double R2 = 0.55;
int R;

bool sigma = false;
//========================================
//========================================
//=========== main function ==============
//========================================
//========================================
void EvaluateSys_rathi(bool epup = false, int debug = 0)
{
  //=================================================================
  //============= List Files for Different Systematics ==============
  //=================================================================
      // nominal file
      std::string file_nom       = "UnfoldingData_2D_non_R04_qV0C_epV0A_3050_hi30_30_Nov14.root";
      std::string file_plus5     = "UnfoldingData_2D_non_R04_qV0C_epV0A_3050_hi30_hipT_30_Nov14.root";
      std::string file_min5      = "UnfoldingData_2D_non_R04_qV0C_epV0A_3050_hi30_lopT_30_Nov14.root";
      std::string file_eff       = "UnfoldingData_2D_non_R04_qV0C_epV0A_3050_hi30_trak_30_Nov14.root";
      std::string file_reweight  = "UnfoldingData_2D_non_R04_qV0C_epV0A_3050_hi30_rewe_30_Nov14.root";
   
     vector<double> merge = {35, 50, 60, 80, 100, 120};  
     vector<double> nomerge = {35, 40, 50, 60, 80, 100, 120};  
     std::vector<double> kBinsUnfolded = {10.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 190.0, 250.0};

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
        const char* r020 = "R02";
        const char* r040 = "R04";
        for (int i = 2; i < words.size() - 1; i++) {
          s2 << words[i] << " ";
          const char *newword = words[i].c_str();
          if (strncmp(newword, ccent, 3) == 0)   cent = true;
          if (strncmp(newword, scent, 3) == 0)   semi = true;
          if (strncmp(newword, r020, 3) == 0)     r02 = true;
          if (strncmp(newword, r040, 3) == 0)     r04 = true;
        } 
       if (cent == true)  { 
         centl = 0.0;
         centr = 10.0;
        }
       if (semi == true)  {
         centl = 30.0;
         centr = 50.0;
        } 
       if (r02 == true)   R = 2;
       if (r04 == true)   R = 4;

      // file for the reconstruction efficiency
      //TH1D* hreceff = (TH1D*)getRuedigerRecEffML(0.4);
      TFile *receffile  = new TFile("ChargedJetRecEfficiencies.root");
      stringstream recname;
      recname << "RecEff_R0" << R << "0_5GeV";
      TH1D* hreceff0                = (TH1D*)receffile->Get(recname.str().c_str());
      hreceff0->SetName("rectemp");
      TH1D* hreceff;
      if (sigma == true)  hreceff   = (TH1D*)hreceff0->Rebin(merge.size()-1, recname.str().c_str(), merge.data());
      if (sigma == false) hreceff   = (TH1D*)hreceff0->Rebin(nomerge.size()-1, recname.str().c_str(), nomerge.data());
      cout << "recname: " << recname.str().c_str() << "\n";
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
      TH1D* nomi;
      if (sigma == true)  nomi = (TH1D*)h1_nomi->Rebin(merge.size()-1, "h1_nomi_copy", merge.data());
      if (sigma == false) nomi = (TH1D*)h1_nomi->Clone("h1_nomi_copy");
      //project out-of-plane
      TH1D* h1_nomo = (TH1D*)h2_nom->ProjectionX("h1_nomo", 0, 1);
      TH1D* nomo;
      if (sigma == true)  nomo = (TH1D*)h1_nomo->Rebin(merge.size()-1, "h1_nomo_copy", merge.data());
      if (sigma == false) nomo = (TH1D*)h1_nomo->Clone("h1_nomo_copy");
      //GetKinematic Efficiencies
      TH1D *hKinO_nom = new TH1D("hKinO_nom", "hKinO_nom", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_nom = new TH1D("hKinI_nom", "hKinI_nom", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_nom, hKinO_nom, hKinI_nom, "Nom");
      TH1D* nomoKE = (TH1D*)ScaleSpec(nomo, hKinO_nom, hreceff);
      nomoKE->SetName("nomoKE");
      TH1D* nomiKE = (TH1D*)ScaleSpec(nomi, hKinI_nom, hreceff);
      nomiKE->SetName("nomiKE");
      //Divide hists to get nominal
      TH1D* h1_nom = (TH1D*)nomoKE->Clone("h1_nom");
      h1_nom->Divide(nomiKE);
      applyR2(h1_nom, nomiKE, nomoKE, R2);
      TH1D* nom = (TH1D*)h1_nom->Clone("h1_nom_copy");
      histvec.push_back(nom);
 
      // ============= handle the systematic for number of iterations +/- 1 ===========
      // +1 Iterations
      TH2D* h2_iterhigh = (TH2D*)fnom->Get(histnameplus1.c_str());
      //project in-plane
      TH1D* h1_iterhighi = (TH1D*)h2_iterhigh->ProjectionX("h1_iterhighi", 4, -1);
      TH1D* iterhighi;
      if (sigma == true)  iterhighi = (TH1D*)h1_iterhighi->Rebin(merge.size()-1, "h1_iterhighi_copy", merge.data());
      if (sigma == false) iterhighi = (TH1D*)h1_iterhighi->Clone("h1_iterhighi_copy");
      //project out-of-plane
      TH1D* h1_iterhigho = (TH1D*)h2_iterhigh->ProjectionX("h1_iterhigho", 0, 1);
      TH1D* iterhigho;
      if (sigma == true)  iterhigho = (TH1D*)h1_iterhigho->Rebin(merge.size()-1, "h1_iterhigho_copy", merge.data());
      if (sigma == false) iterhigho = (TH1D*)h1_iterhigho->Clone("h1_iterhigho_copy");
      //Get Kinematic Efficiencies
      TH1D *hKinO_iterhigh = new TH1D("hKinO_iterhigh", "hKinO_iterhigh", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_iterhigh = new TH1D("hKinI_iterhigh", "hKinI_iterhigh", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_nom, hKinO_iterhigh, hKinI_iterhigh, "Iterhigh");
      TH1D* iterhighoKE = (TH1D*)ScaleSpec(iterhigho, hKinO_iterhigh, hreceff);
      iterhighoKE->SetName("iterhighoKE");
      TH1D* iterhighiKE = (TH1D*)ScaleSpec(iterhighi, hKinI_iterhigh, hreceff);
      iterhighiKE->SetName("iterhighiKE");
      //Divide hists to get nominal
      TH1D* h1_iterhigh = (TH1D*)iterhighoKE->Clone("h1_iterhigh");
      h1_iterhigh->Divide(iterhighiKE);
      //apply EPR correction
      applyR2(h1_iterhigh, iterhighiKE, iterhighoKE, R2);
      TH1D* iterhigh = (TH1D*)h1_iterhigh->Clone("h1_iterhigh_copy");
      histvec.push_back(iterhigh);

      // -1 Iterations
      TH2D* h2_iterlow = (TH2D*)fnom->Get(histnamemin1.c_str());
      TH1D* h1_iterlowi = (TH1D*)h2_iterlow->ProjectionX("h1_iterlowi", 4, -1);
      TH1D* iterlowi;
      if (sigma == true)  iterlowi = (TH1D*)h1_iterlowi->Rebin(merge.size()-1, "h1_iterlowi_copy", merge.data());
      if (sigma == false) iterlowi = (TH1D*)h1_iterlowi->Clone("h1_iterlowi_copy");
      TH1D* h1_iterlowo = (TH1D*)h2_iterlow->ProjectionX("h1_iterlowo", 0, 1);
      TH1D* iterlowo;
      if (sigma == true)  iterlowo = (TH1D*)h1_iterlowo->Rebin(merge.size()-1, "h1_iterlowo_copy", merge.data());
      if (sigma == false) iterlowo = (TH1D*)h1_iterlowo->Clone("h1_iterlowo_copy");
      //Get Kinematic efficiencies
      TH1D *hKinO_iterlow = new TH1D("hKinO_iterlow", "hKinO_iterlow", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_iterlow = new TH1D("hKinI_iterlow", "hKinI_iterlow", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_nom, hKinO_iterlow, hKinI_iterlow, "Iterlow");
      TH1D* iterlowoKE = (TH1D*)ScaleSpec(iterlowo, hKinO_iterlow, hreceff);
      iterlowoKE->SetName("iterlowoKE");
      TH1D* iterlowiKE = (TH1D*)ScaleSpec(iterlowi, hKinI_iterlow, hreceff);
      iterlowiKE->SetName("iterlowiKE");
      //TH1D* iterlow = (TH1D*)ScaleSpec(h1_iterlow, hKinNum_nom, hreceff);
      TH1D* h1_iterlow = (TH1D*)iterlowoKE->Clone("h1_iterlow");
      h1_iterlow->Divide(iterlowiKE);
      //apply EPR correction
      applyR2(h1_iterlow, iterlowiKE, iterlowoKE, R2);
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
      TH2D* h2_eff = (TH2D*)feff->Get(histname.c_str());
      //project in-plane
      TH1D* h1_effi = (TH1D*)h2_eff->ProjectionX("h1_effi", 4, -1);
      TH1D* effi;
      if (sigma == true)  effi = (TH1D*)h1_effi->Rebin(merge.size()-1, "h1_effi_copy", merge.data());
      if (sigma == false) effi = (TH1D*)h1_effi->Clone("h1_effi_copy");
      //project out-of-plane
      TH1D* h1_effo = (TH1D*)h2_eff->ProjectionX("h1_effo", 0, 1);
      TH1D* effo;
      if (sigma == true)  effo = (TH1D*)h1_effo->Rebin(merge.size()-1, "h1_effo_copy", merge.data());
      if (sigma == false) effo = (TH1D*)h1_effo->Clone("h1_effo_copy");
      //Get Kinematic efficiencies
      TH1D *hKinO_eff = new TH1D("hKinO_eff", "hKinO_eff", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_eff = new TH1D("hKinI_eff", "hKinI_eff", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_eff, hKinO_eff, hKinI_eff, "Eff");
      TH1D* effoKE = (TH1D*)ScaleSpec(effo, hKinO_eff, hreceff);
      effoKE->SetName("effoKE");
      TH1D* effiKE = (TH1D*)ScaleSpec(effi, hKinI_eff, hreceff);
      effiKE->SetName("effiKE");
      //Make Ratio
      TH1D* h1_eff = (TH1D*)effoKE->Clone("h1_eff");
      h1_eff->Divide(effiKE);
      //apply EPR correction
      applyR2(h1_eff, effiKE, effoKE, R2);
      TH1D* eff = (TH1D*)h1_eff->Clone("h1_eff_copy");
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
      TH1D* plus5i;
      if (sigma == true)  plus5i = (TH1D*)h1_plus5i->Rebin(merge.size()-1, "h1_plus5i_copy", merge.data());
      if (sigma == false) plus5i = (TH1D*)h1_plus5i->Clone("h1_plus5i_copy");
      //project out-of-plane
      TH1D* h1_plus5o = (TH1D*)h2_plus5->ProjectionX("h1_plus5o", 0, 1);
      TH1D* plus5o;
      if (sigma == true)  plus5o = (TH1D*)h1_plus5o->Rebin(merge.size()-1, "h1_plus5o_copy", merge.data());
      if (sigma == false) plus5o = (TH1D*)h1_plus5o->Clone("h1_plus5o_copy");
      //Get Kinematic efficiencies
      TH1D *hKinO_plus5 = new TH1D("hKinO_plus5", "hKinO_plus5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_plus5 = new TH1D("hKinI_plus5", "hKinI_plus5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_plus5, hKinO_plus5, hKinI_plus5, "plus5");
      TH1D* plus5oKE = (TH1D*)ScaleSpec(plus5o, hKinO_plus5, hreceff);
      plus5oKE->SetName("plus5oKE");
      TH1D* plus5iKE = (TH1D*)ScaleSpec(plus5i, hKinI_plus5, hreceff);
      plus5iKE->SetName("plus5iKE");

      //Divide hists to get nominal
      TH1D* h1_plus5 = (TH1D*)plus5oKE->Clone("h1_plus5");
      h1_plus5->Divide(plus5iKE);
      //apply EPR correction
      applyR2(h1_plus5, plus5iKE, plus5oKE, R2);
      TH1D* plus5 = (TH1D*)h1_plus5->Clone("h1_plus5_copy");
      histvec.push_back(plus5);

      //-5 GeV/c
      TFile* fmin5 = new TFile(file_min5.c_str());
      TH2D* h2_min5 = (TH2D*)fmin5->Get(histname.c_str());
      //project in-plane
      TH1D* h1_min5i = (TH1D*)h2_min5->ProjectionX("h1_min5i", 4, -1);
      TH1D* min5i;
      if (sigma == true)  min5i = (TH1D*)h1_min5i->Rebin(merge.size()-1, "h1_min5i_copy", merge.data());
      if (sigma == false) min5i = (TH1D*)h1_min5i->Clone("h1_min5i_copy");
      //project out-of-plane
      TH1D* h1_min5o = (TH1D*)h2_min5->ProjectionX("h1_min5o", 0, 1);
      TH1D* min5o;
      if (sigma == true)  min5o = (TH1D*)h1_min5o->Rebin(merge.size()-1, "h1_min5o_copy", merge.data());
      if (sigma == false) min5o = (TH1D*)h1_min5o->Clone("h1_min5o_copy");
      //Get Kinematic efficiencies
      TH1D *hKinO_min5 = new TH1D("hKinO_min5", "hKinO_min5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_min5 = new TH1D("hKinI_min5", "hKinI_min5", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_min5, hKinO_min5, hKinI_min5, "min5");
      TH1D* min5oKE = (TH1D*)ScaleSpec(min5o, hKinO_min5, hreceff);
      min5oKE->SetName("min5oKE");
      TH1D* min5iKE = (TH1D*)ScaleSpec(min5i, hKinI_min5, hreceff);
      min5iKE->SetName("min5iKE");
      TH1D* h1_min5 = (TH1D*)min5oKE->Clone("h1_min5");
      h1_min5->Divide(min5iKE);
      //apply EPR correction
      applyR2(h1_min5, min5iKE, min5oKE, R2);
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
      TH2D* h2_rw = (TH2D*)frw->Get(histname.c_str());
      TH1D* h1_rwi = (TH1D*)h2_rw->ProjectionX("h1_rwi", 4, -1);
      TH1D* rwi;
      if (sigma == true)  rwi = (TH1D*)h1_rwi->Rebin(merge.size()-1, "h1_rwi_copy", merge.data());
      if (sigma == false) rwi = (TH1D*)h1_rwi->Clone("h1_rwi_copy");
      TH1D* h1_rwo = (TH1D*)h2_rw->ProjectionX("h1_rwo", 0, 1);
      TH1D* rwo;
      if (sigma == true)  rwo = (TH1D*)h1_rwo->Rebin(merge.size()-1, "h1_rwo_copy", merge.data());
      if (sigma == false) rwo = (TH1D*)h1_rwo->Clone("h1_rwo_copy");
      //Get Kinematic efficiencies
      TH1D *hKinO_rw = new TH1D("hKinO_rw", "hKinO_rw", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_rw = new TH1D("hKinI_rw", "hKinI_rw", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_reweight, hKinO_rw, hKinI_rw, "rw");
      TH1D* rwoKE = (TH1D*)ScaleSpec(rwo, hKinO_rw, hreceff);
      rwoKE->SetName("rwoKE");
      TH1D* rwiKE = (TH1D*)ScaleSpec(rwi, hKinI_rw, hreceff);
      rwiKE->SetName("rwiKE");
      TH1D* h1_rw = (TH1D*)rwoKE->Clone("h1_rw");
      h1_rw->Divide(rwiKE);
      //apply EPR correction
      applyR2(h1_rw, rwiKE, rwoKE, R2);
      TH1D* rw = (TH1D*)h1_rw->Clone("h1_rw_copy");
      histvec.push_back(rw);
      std::pair<TH1D*, TH1D*> weight = CalculateDeltaCut(nom, rw, rw, true);
      histvec_up.push_back(weight.first);
      histvec_down.push_back(weight.second);
  

      //=============systematic for ep angle======================
     
      //only saved for lower uncertainty 
      //retrieve 2D hist
      TH2D* h2_gluon = (TH2D*)fnom->Get(histname.c_str()); 
      //project in-plane
      TH1D* h1_glui = (TH1D*)h2_gluon->ProjectionX("h1_glui", 4, -1);
      TH1D* glui;
      if (sigma == true)  glui = (TH1D*)h1_glui->Rebin(merge.size()-1, "h1_glui_copy", merge.data());
      if (sigma == false) glui = (TH1D*)h1_glui->Clone("h1_glui_copy");
      //project out-of-plane
      TH1D* h1_gluo = (TH1D*)h2_gluon->ProjectionX("h1_gluo", 0, 1);
      TH1D* gluo;
      if (sigma == true)  gluo = (TH1D*)h1_gluo->Rebin(merge.size()-1, "h1_gluo_copy", merge.data());
      if (sigma == false) gluo = (TH1D*)h1_gluo->Clone("h1_gluo_copy");
      //Get Kinematic efficiencies
      TH1D *hKinO_glu = new TH1D("hKinO_glu", "hKinO_glu", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      TH1D *hKinI_glu = new TH1D("hKinI_glu", "hKinI_glu", kBinsUnfolded.size()-1, kBinsUnfolded.data());
      GetKinematicEfficiency(file_nom, hKinO_glu, hKinI_glu, "EP");
      TH1D* gluoKE = (TH1D*)ScaleSpec(gluo, hKinO_glu, hreceff);
      gluoKE->SetName("gluoKE");
      TH1D* gluiKE = (TH1D*)ScaleSpec(glui, hKinI_glu, hreceff);
      gluiKE->SetName("gluiKE");
      //Divide hists to get nominal
      TH1D* h1_gluons = (TH1D*)gluoKE->Clone("h1_gluons");
      h1_gluons->Divide(gluiKE);
      //apply EPR correction
      applyR2(h1_gluons, gluiKE, gluoKE, (R2-(0.02*R2)));
      TH1D* gluons = (TH1D*)h1_gluons->Clone("h1_gluons_copy");
      if (epup == false)  {
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
      if (sigma == false){
        if (epup == false) name << "preSystematics_R0" << R << "0_3050_" << (R/2)+1 << "0_hiinout_Nov14.root";
        if (epup == true)  name << "preSystematics_R0" << R << "0_3050_" << (R/2)+1 << "0_hiinoutup_Nov14.root";}
      if (sigma == true){
        if (epup == false) name << "preSystematics_R0" << R << "0_sigma_3050_" << (R/2)+1 << "0_hiinout_Nov14.root";
        if (epup == true)  name << "preSystematics_R0" << R << "0_sigma_3050_" << (R/2)+1 << "0_hiinoutup_Nov14.root";}

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
  TCanvas *c1 = new TCanvas("hists", "hists", 1000, 800);
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
     pad1->SetBottomMargin(0.15); // Upper and lower plot are joined
     pad1->Draw();
   pad1->cd();               // pad1 becomes the current pad   
   nom->GetXaxis()->SetRangeUser(30.0, 140.0);
   nom->SetMinimum(0.01);
   nom->SetMaximum(2.99);
   nom->Draw("same");
     nom->SetLineColor(kRed+2);
     nom->SetMarkerColor(kRed+2);
     nom->SetMarkerStyle(20);
   iterhigh->Draw("same");
     iterhigh->SetLineColor(kOrange);
     iterhigh->SetMarkerColor(kOrange);
     iterhigh->SetMarkerStyle(20);
   iterlow->Draw("same");
     iterlow->SetLineColor(kYellow);
     iterlow->SetMarkerColor(kYellow);
     iterlow->SetMarkerStyle(20);
   rw->Draw("same");
     rw->SetLineColor(kPink);
     rw->SetMarkerColor(kPink);
     rw->SetMarkerStyle(20);
   plus5->Draw("same");
     plus5->SetLineColor(kGreen+2);
     plus5->SetMarkerColor(kGreen+2);
     plus5->SetMarkerStyle(20);
   min5->Draw("same");
     min5->SetLineColor(kViolet-1);
     min5->SetMarkerColor(kViolet-1);
     min5->SetMarkerStyle(20);
   eff->Draw("same");
     eff->SetLineColor(kOrange+7);
     eff->SetMarkerColor(kOrange+7);
     eff->SetMarkerStyle(20);
  gluons->Draw("same");
     gluons->SetLineColor(kBlue);
     gluons->SetMarkerColor(kBlue);
     gluons->SetMarkerStyle(20);

  TLegend *leg = new TLegend(0.6, 0.6, 0.8, 0.8);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->AddEntry(nom, "nom", "pl");
  leg->AddEntry(plus5, "plus5", "pl");
  leg->AddEntry(min5, "min5", "pl");
  leg->AddEntry(gluons, "ep", "pl");
  leg->AddEntry(eff, "track", "pl");
  leg->AddEntry(rw, "rewe", "pl");
  leg->AddEntry(iterhigh, "iterhigh", "pl");
  leg->AddEntry(iterlow, "iterlow", "pl");
  leg->Draw("same"); 

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


//========================================================
//========================================================
//=========== EP Resolution Correction ===================
//========================================================
//========================================================

void applyR2(TH1D* Nratio, TH1D* NinPre, TH1D* NoutPre, double EPR) {
      //calculate the v2 bin-by-bin
        for (int i = 1; i <= Nratio->GetNbinsX(); i++)   {
          double Nin = NinPre->GetBinContent(i);
          double Nout = NoutPre->GetBinContent(i);
          double v2 = (pi/(3*sqrt(3)))*(1./EPR)*(Nin-Nout)/(Nin+Nout);
          double Nrat = ((2*pi)-(6*sqrt(3)*v2))/((2*pi)+(6*sqrt(3)*v2));
          Nratio->SetBinContent(i, Nrat);
        }
 }




//========================================================
//========================================================
//=========== ScaleSpec Function Definition ==============
//========================================================
//========================================================
TH1D* ScaleSpec(TH1D* h1_Orig,  TH1D* heffnum, TH1D* hrecEff){
  TH1D* h1 = (TH1D*)h1_Orig->Clone(Form("%s_copy", h1_Orig->GetName()));

  //correct the spectrum for the kinematic efficiency and rec eff at once
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


