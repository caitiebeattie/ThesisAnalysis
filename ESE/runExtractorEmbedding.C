class AliEmcalJetTask;
class AliAnalysisTaskEmcalJetPerformance;
class AliMultSelectionTask;
class AliAnalysisTaskRho;
class AliAnalysisTaskEmcalEmbeddingHelper;
class AliTrackContainer;
class AliParticleContainer;
class AliAnalysisTaskJetExtractor;
class AliAnalysisTaskJetQnVectors;

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "PWGJE/EMCALJetTasks/macros/AddTaskJetExtractor.C"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"

void runExtractorEmbedding() {

 //parameterslist
 double JetRadius = 0.4;
 double MinJpT = 5.0;
 double kTrackPtCut = 0.15;
 double kClusPtCut = 0.3;
 double kGhostArea = 0.005;

 //header location
 gInterpreter->ProcessLine(".include $ROOTSYS/include");
 gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

 //create analysis manager
 AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisMyTask");
 AliAODInputHandler *aodH = new AliAODInputHandler();
 mgr->SetInputEventHandler(aodH);

 //compile the class (locally) with debug symbols
 //gInterpreter->LoadMacro("AliAnalysisTaskMyTask.cxx++g");
 
 //===================================================
 //================ Add Task List ====================
 //===================================================
 //CDB Connect
 cout << "CDB Conncect\n";
 AliTaskCDBconnect *CDBConnectTask=reinterpret_cast<AliTaskCDBconnect*>(
        gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C()"));
 //TMacro CDBadd (gSystem->ExpandPathName("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C"));
 //AliTaskCDBconnect* CDBConnectTask = reinterpret_cast<AliTaskCDBconnect*>(CDBadd.Exec());
 CDBConnectTask->SetFallBackToRaw(kTRUE);

 //Embedding Helper
   AliAnalysisTaskEmcalEmbeddingHelper * EmbeddingTask = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
     EmbeddingTask->SetFileListFilename("aodFilesEmbedData.txt");
     EmbeddingTask->SetRandomFileAccess(kTRUE);
     EmbeddingTask->Initialize();


   //Jet Finder anti-kt spectrum
   cout << "Jet Finder: AKT hybrid\n";
   AliEmcalJetTask *AKT = AliEmcalJetTask::AddTaskEmcalJet("tracks", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, kTrackPtCut, kClusPtCut, kGhostArea, AliJetContainer::pt_scheme, "hybJet", 0.0 , kFALSE, kFALSE);
     AKT->SetUseNewCentralityEstimation(kTRUE);
     AKT->SetNCentBins(5);
     AKT->SelectCollisionCandidates(AliVEvent::kAnyINT);
     AKT->SetForceBeamType(AliAnalysisTaskEmcal::kAA);
     AKT->SetNeedEmcalGeom(kFALSE);
     AKT->SetMinJetArea(0.30);  // 0.6*pi*R^2

     AliTrackContainer *embTracks = new AliTrackContainer("tracks");
     embTracks->SetIsEmbedding(kTRUE);
     AKT->AdoptTrackContainer(embTracks);
     AKT->GetTrackContainer(0)->SetParticlePtCut(0.15);

   cout << "Jet Finder: AKT detector\n";
   AliEmcalJetTask *AKTD = AliEmcalJetTask::AddTaskEmcalJet("tracks", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, kTrackPtCut, kClusPtCut, kGhostArea, AliJetContainer::pt_scheme, "detJet", 0.0 , kFALSE, kFALSE);
     AKTD->SetUseNewCentralityEstimation(kTRUE);
     AKTD->SetNCentBins(5);
     AKTD->SelectCollisionCandidates(AliVEvent::kAnyINT);
     AKTD->SetForceBeamType(AliAnalysisTaskEmcal::kAA);
     AKTD->SetNeedEmcalGeom(kFALSE);
    AKTD->SetMinJetArea(0.30);  // 0.6*pi*R^2

     AliTrackContainer *tracksDetLevel = AKTD->GetTrackContainer(0);
     tracksDetLevel->SetIsEmbedding(kTRUE);
     tracksDetLevel->SetParticlePtCut(0.15);

   cout << "Jet Finder: AKT truth\n";
   AliEmcalJetTask *AKTP = AliEmcalJetTask::AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, kTrackPtCut, kClusPtCut, kGhostArea, AliJetContainer::pt_scheme, "parJet", 0.0 , kFALSE, kFALSE);
     AKTP->SetUseNewCentralityEstimation(kTRUE);
     AKTP->SetNCentBins(5);
     AKTP->SelectCollisionCandidates(AliVEvent::kAnyINT);
     AKTP->SetForceBeamType(AliAnalysisTaskEmcal::kAA);
     AKTP->SetNeedEmcalGeom(kFALSE);
    AKTP->SetMinJetArea(0.30);  // 0.6*pi*R^2

    AliMCParticleContainer *partCont = AKTP->GetMCParticleContainer(0);
    partCont->SetIsEmbedding(kTRUE);
    partCont->SetParticlePtCut(0.15);

   //Jet Finder kt
   cout << "Jet Finder: KT\n";
   AliEmcalJetTask *KT = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, 0.4, AliJetContainer::kChargedJet, kTrackPtCut, kClusPtCut, kGhostArea, AliJetContainer::pt_scheme, "Jet", 0., kFALSE, kFALSE);
     KT->SetUseNewCentralityEstimation(kTRUE);
     KT->SetNCentBins(5);
     KT->SelectCollisionCandidates(AliVEvent::kAnyINT);
     KT->SetForceBeamType(AliAnalysisTaskEmcal::kAA);

   //QnVector Getter
   cout << "Qn Vector\n";
   //handler
   //task
   //TMacro QN (gSystem->ExpandPathName("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJEQnVectors.C"));
   //AliAnalysisTaskJetQnVectors* QNtask = reinterpret_cast<AliAnalysisTaskJetQnVectors*>(QN.Exec());
    AliAnalysisTaskJetQnVectors *QNtask=reinterpret_cast<AliAnalysisTaskJetQnVectors*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJEQnVectors.C(\"AliAnalysisTaskJetQnVectors\", \"_V2_TPCEta0dot8\", 2, 1, 0, \"calibV0TrklNoEtaCutRun2.root\", \"\")")); /*calibq2V0CRun2Vtx10P118q.root, calibV0TrklNoEtaCutRun2.root*/
   //AliAnalysisTaskJetQnVectors *QN = AliAnalysisTaskJetQnVectors::AddTaskJEQnVectors("JEQnVectorsV2", "_V2_TPCEta0dot8", 2, 1, 0, "alien:////alice/cern.ch/user/f/fgrosa/QnVectorCalibrations/calibV0TrklTPCNoEtaCutRun218rVtx14MRP2New.root","alien:///alice/cern.ch/user/f/fgrosa/q2Splines/Splines_q2_010_3050_LHC18r.root",kFALSE,kFALSE);

   //Multiplicity Selection
   cout << "Multiplicity Selection\n";
   AliMultSelectionTask* MultTask = AliMultSelectionTask::AddTaskMultSelection(kFALSE);
     MultTask->SelectCollisionCandidates(AliVEvent::kAnyINT);

   //Rho Task
   cout << "Rho Task\n";
   AliAnalysisTaskRho* RhoTask = AliAnalysisTaskRho::AddTaskRhoNew("usedefault", "", "Rho", 0.4, AliEmcalJet::kTPCfid, AliJetContainer::kChargedJet, kTRUE);
     RhoTask->SetExcludeLeadJets(2);
     RhoTask->SetUseNewCentralityEstimation(kTRUE);
     RhoTask->SetNCentBins(5); 
     RhoTask->SelectCollisionCandidates(AliVEvent::kAnyINT);
     RhoTask->SetForceBeamType(AliAnalysisTaskEmcal::kAA);
 //load the addtask macro and create the task
 //AliAnalysisTaskMyTask *task = reinterpret_cast<AliAnalysisTaskMyTask*>(gInterpreter->ExecuteMacro("AddMyTask.C"));
 
   //Jet Extractor
   cout << "Jet Extractor\n";
   AliAnalysisTaskJetExtractor* JExt = AddTaskJetExtractor("tracks", "", "hybJet_AKTChargedR040_tracks_pT0150_pt_scheme", "Rho", 0.4, "hybJet");
     JExt->SelectCollisionCandidates(AliVEvent::kAnyINT);
     JExt->SetNeedEmcalGeom(kFALSE);
     JExt->SetVzRange(-10,10);
     JExt->SetSaveQVector(kTRUE);
     JExt->SetQ2Detector(2);
     JExt->SetEPDetector(0);
     JExt->SetUseNewCentralityEstimation(kTRUE);
     JExt->SetSaveMCInformation(kTRUE);
     JExt->SetSaveConstituents(kFALSE);
     JExt->SetSaveJetShapes(kFALSE);
     //JExt->SetLightTreeMode(kTRUE);
     JExt->GetJetTree()->AddExtractionPercentage(-20,10, 0.0);
     JExt->GetJetTree()->AddExtractionPercentage(10,20, 1.0);
     JExt->GetJetTree()->AddExtractionPercentage(20,30, 1.0);
     JExt->GetJetTree()->AddExtractionPercentage(30,40, 1.0);
     JExt->GetJetTree()->AddExtractionPercentage(40,60, 1.0);
     JExt->GetJetTree()->AddExtractionPercentage(60,80, 1.0);
     JExt->GetJetTree()->AddExtractionPercentage(80,100, 1.0);
     JExt->GetJetTree()->AddExtractionPercentage(100,200, 1.0);
     JExt->GetJetContainer(0)->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
     JExt->GetJetContainer(0)->SetPtBiasJetTrack(5.0);
     JExt->GetJetContainer(0)->SetMaxTrackPt(100.0);

     AliParticleContainer * mcpartCont =  new AliTrackContainer("tracks");
     mcpartCont->SetIsEmbedding(kTRUE);
     JExt->AdoptParticleContainer(mcpartCont);

     JExt->GetJetContainer(0)->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
     JExt->GetJetContainer(0)->SetPtBiasJetTrack(5.0);
     JExt->GetJetContainer(0)->SetMaxTrackPt(100.0);
     JExt->AddJetContainer("detJet_AKTChargedR040_tracks_pT0150_pt_scheme");
     JExt->GetJetContainer(1)->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
     JExt->GetJetContainer(1)->SetPtBiasJetTrack(5.0);
     JExt->GetJetContainer(1)->SetMaxTrackPt(100.0);
     JExt->AddJetContainer("parJet_AKTChargedR040_mcparticles_pT0150_pt_scheme");
     JExt->GetJetContainer(2)->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
     JExt->GetJetContainer(2)->SetPtBiasJetTrack(5.0);
     JExt->GetJetContainer(2)->SetMaxTrackPt(100.0);


 if(!mgr->InitAnalysis()) return;
 mgr->SetDebugLevel(0);


 //Alternative Method for Building a Chain
 TChain* chain = new TChain("aodTree");
 chain->Add("~/longlivealidock/data/2015/AliAOD.root");  //"~/alice/data/2015/LHC15o/000246757/pass1/AOD/002/AliAOD.root"
 //chain->Add("~/alice/data/2015/LHC15o/000246757/pass1/AOD/007/AliAOD.root");
 //chain->Add("~/alice/data/2017/LHC17q/000282365/pass1_FAST/AOD/003/AliAOD.root");

 //start the analysis locally
 mgr->StartAnalysis("local", chain);


}
