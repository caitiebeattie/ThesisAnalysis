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


//========================================
//=========== main function ==============
//========================================

void formatSplines (string rootfilename) {

  const char* lerootfile = rootfilename.c_str();
  stringstream foutname;

   //==================================================
   //======== 1. get file tag for legend 
   //==================================================
   vector<const char*> letters(0);
   stringstream ss;
   const char* delim = ".";
   for (int i = 0; i < rootfilename.length(); i++)    {
       const char *a = &rootfilename[i];
       letters.push_back(a);
       if (strncmp(letters[i], delim, 1) == 0)    break;
       else ss << rootfilename[i];
       }
  
   foutname << ss.str() << "_formatted.root";

  

  //===========================================================
  //========= 3. retrieve info from data file
  //===========================================================
      //GetFile
      TFile *fin = new TFile(lerootfile);
      TFile *fout = new TFile(foutname.str().c_str(), "RECREATE");
 
      TDirectory *dir = (TDirectory*)fout->mkdir("SplineListq2V0C", "");
      dir->cd(); 

     

      //TList* l1 = new TList(); 
      //l1->SetName("mylist");
      //l1->Print();
      //l1->Write();
      

      //TList* l2 = new TList(); 
      

      TSpline3 *spline;
      for (int i = 0; i < 90; i++)  {
         stringstream splinename;
         splinename << "sp_q2V0C_" << i;
         spline = (TSpline3*)fin->Get(splinename.str().c_str());
         spline->Write();
         //l1->Add(spline);
      }

      fout->Add(dir);
      fout->Print();
      fout->Write();
      //fout->Close();
      //l1->Add(l2);
      //l2->Print();
      //l2->Write();

      //delete fout;

}


