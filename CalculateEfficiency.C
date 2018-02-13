// 'CalculateEfficiency.C'
// Derek Anderson
// 02.07.2018
//
// Use this to calculate the tracking
// efficiency using output from the
// 'StThirdJetMaker'
//
// NOTE: in data, what I mean by
//       'efficiency' is the no.
//       of tracks before QA cuts
//       to the no. of tracks after
//       QA cuts.


#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t  NTrkMax(5000);
static const UInt_t  NTwrMax(5000);
static const UInt_t  NMatchMax(10);
static const UInt_t  NHotTwr(41);
static const UInt_t  NBadRuns(45);
static const UInt_t  NTrgs(2);
static const UInt_t  NCut(2);
static const UInt_t  NVal(2);
static const TString sTreeDefault("GfmtoDst_gnt");
static const TString sInputDefault("../../MuDstMatching/output/merged/pt5.match.root");
static const TString sOutputDefault("pp200r9pt35u.default.root");



void CalculateEfficiency(UInt_t &nTrgPi0, UInt_t &nTrgGam, const Bool_t isInBatchMode=false, const Bool_t isDetectorLevel=false, const Bool_t isTriggered=false, const TString sInput=sInputDefault, const TString sTree=sTreeDefault, const TString sOutput=sOutputDefault) {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning purity calculation..." << endl;


  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);

  // trigger paramters
  const Int_t    adcMax(6004);
  const Double_t eStrMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t hTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  const Double_t tspPi0[2] = {0., 0.08};
  const Double_t tspGam[2] = {0.2, 0.6};

  // track parameters
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(20.);


  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fInput) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input tree
  TTree *tInput;
  fInput -> GetObject(sTree.Data(), tInput);
  if (!tInput) {
    cerr << "PANIC: couldn't grab input tree!" << endl;
    return;
  }
  if (isDetectorLevel)
    cout << "    Grabbed detector tree." << endl;
  else
    cout << "    Grabbed particle tree." << endl;


  // declare input leaf addresses
  UInt_t   fUniqueID;
  UInt_t   fBits;
  Long64_t runNumber;
  Long64_t eventNumber;
  Int_t    trigID;
  Int_t    nGlobalTracks;
  Int_t    nPrimaryTracks;
  Int_t    refMult;
  Double_t vpdVz;
  Double_t xVertex;
  Double_t yVertex;
  Double_t zVertex;
  Double_t bbcZVertex;
  Double_t zdcCoincidenceRate;
  Double_t bbcCoincidenceRate;
  Double_t backgroundRate;
  Double_t bbcBlueBackgroundRate;
  Double_t bbcYellowBackgroundRate;
  Double_t refMultPos;
  Double_t refMultNeg;
  Double_t bTOFTrayMultiplicity;
  Int_t    nVerticies;
  Double_t MagF;
  Double_t VrtxRank;
  Int_t    FlagEvent_TrgTrkMisMtch;
  Float_t  Etsp;
  Int_t    ETwrdidT;
  Int_t    ETwradc11;
  Float_t  ETwreneT0;
  Float_t  ETwreT;
  Float_t  ETwrENET0;
  Float_t  ETwrphT;
  Float_t  ETwrPTower;
  Float_t  ETwrpidTower;
  Int_t    ETwrmoduleT;
  Float_t  EClustEneT0;
  Float_t  EClustetav1;
  Float_t  EClustphiv1;
  Float_t  EEstrpen01;
  Float_t  EEstrpen02;
  Float_t  EEstrpen03;
  Float_t  EEstrpen0;
  Float_t  EEstrpen1;
  Float_t  EEstrpen2;
  Float_t  EEstrpen3;
  Float_t  EEstrpen4;
  Float_t  EEstrpen5;
  Float_t  EEstrpen6;
  Float_t  EEstrpen7;
  Float_t  EEstrpen8;
  Float_t  EEstrpen9;
  Float_t  EEstrpen10;
  Float_t  EEstrpen11;
  Float_t  EEstrpen12;
  Float_t  EEstrpen13;
  Float_t  EEstrpen14;
  Float_t  EEstrpen15;
  Int_t    ETwrdidE;
  Float_t  EPstripenp01;
  Float_t  EPstripenp02;
  Float_t  EPstripenp03;
  Float_t  EPstripenp0;
  Float_t  EPstripenp1;
  Float_t  EPstripenp2;
  Float_t  EPstripenp3;
  Float_t  EPstripenp4;
  Float_t  EPstripenp5;
  Float_t  EPstripenp6;
  Float_t  EPstripenp7;
  Float_t  EPstripenp8;
  Float_t  EPstripenp9;
  Float_t  EPstripenp10;
  Float_t  EPstripenp11;
  Float_t  EPstripenp12;
  Float_t  EPstripenp13;
  Float_t  EPstripenp14;
  Float_t  EPstripenp15;
  Float_t  EclustEnnq1;
  Float_t  EclustEnnq20;
  Float_t  EclustEnnq19;
  Float_t  EclustEnpq1;
  Float_t  EclustEnpq20;
  Float_t  EclustEnpq19;
  Float_t  EclustEnpq21;
  Int_t    PrimaryTrackArray_;
  UInt_t   PrimaryTrackArray_fUniqueID[NTrkMax];
  UInt_t   PrimaryTrackArray_fBits[NTrkMax];
  Double_t PrimaryTrackArray_nHitsFit[NTrkMax];
  Double_t PrimaryTrackArray_nHitsPoss[NTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[NTrkMax];
  Double_t PrimaryTrackArray_pZ[NTrkMax];
  Double_t PrimaryTrackArray_pX[NTrkMax];
  Double_t PrimaryTrackArray_pY[NTrkMax];
  Double_t PrimaryTrackArray_pT[NTrkMax];
  Double_t PrimaryTrackArray_dEdx[NTrkMax];
  Double_t PrimaryTrackArray_charge[NTrkMax];
  Double_t PrimaryTrackArray_tofBeta[NTrkMax];
  Double_t PrimaryTrackArray_eta[NTrkMax];
  Double_t PrimaryTrackArray_phi[NTrkMax];
  Double_t PrimaryTrackArray_nSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_nSigPion[NTrkMax];
  Double_t PrimaryTrackArray_nSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_nSigProton[NTrkMax];
  Double_t PrimaryTrackArray_dcag[NTrkMax];
  Double_t PrimaryTrackArray_nHits[NTrkMax];
  Double_t PrimaryTrackArray_dEdxHits[NTrkMax];
  Double_t PrimaryTrackArray_firstZPoint[NTrkMax];
  Double_t PrimaryTrackArray_lastZPoint[NTrkMax];
  Double_t PrimaryTrackArray_tofSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_tofSigPion[NTrkMax];
  Double_t PrimaryTrackArray_tofSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_tofSigProton[NTrkMax];
  Double_t PrimaryTrackArray_timeOfflight[NTrkMax];
  Double_t PrimaryTrackArray_pathLength[NTrkMax];
  Int_t    PrimaryTrackArray_trkIndex[NTrkMax];
  Int_t    TowerArray_;
  UInt_t   TowerArray_fUniqueID[NTwrMax];
  UInt_t   TowerArray_fBits[NTwrMax];
  Int_t    TowerArray_TwrId[NTwrMax];
  Float_t  TowerArray_TwrEng[NTwrMax];
  Float_t  TowerArray_TwrEta[NTwrMax];
  Float_t  TowerArray_TwrPhi[NTwrMax];
  Float_t  TowerArray_TwrADC[NTwrMax];
  Float_t  TowerArray_TwrPed[NTwrMax];
  Float_t  TowerArray_TwrRMS[NTwrMax];
  Int_t    TowerArray_TwrMatchIdnex[NTwrMax];
  Int_t    TowerArray_NoOfmatchedTrk[NTwrMax];
  Float_t  TowerArray_TwrMatchP[NTwrMax];
  Float_t  TowerArray_TwrPx[NTwrMax];
  Float_t  TowerArray_TwrPy[NTwrMax];
  Float_t  TowerArray_TwrPz[NTwrMax];
  Int_t    TowerArray_fNAssocTracks[NTwrMax];
  Int_t    TowerArray_fMatchedTracksArray_[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_P[NTwrMax][NMatchMax];

  // declare input branches
  TBranch *bEventList_fUniqueID;
  TBranch *bEventList_fBits;
  TBranch *bEventList_runNumber;
  TBranch *bEventList_eventNumber;
  TBranch *bEventList_trigID;
  TBranch *bEventList_nGlobalTracks;
  TBranch *bEventList_nPrimaryTracks;
  TBranch *bEventList_refMult;
  TBranch *bEventList_vpdVz;
  TBranch *bEventList_xVertex;
  TBranch *bEventList_yVertex;
  TBranch *bEventList_zVertex;
  TBranch *bEventList_bbcZVertex;
  TBranch *bEventList_zdcCoincidenceRate;
  TBranch *bEventList_bbcCoincidenceRate;
  TBranch *bEventList_backgroundRate;
  TBranch *bEventList_bbcBlueBackgroundRate;
  TBranch *bEventList_bbcYellowBackgroundRate;
  TBranch *bEventList_refMultPos;
  TBranch *bEventList_refMultNeg;
  TBranch *bEventList_bTOFTrayMultiplicity;
  TBranch *bEventList_nVerticies;
  TBranch *bEventList_MagF;
  TBranch *bEventList_VrtxRank;
  TBranch *bEventList_FlagEvent_TrgTrkMisMtch;
  TBranch *bEventList_Etsp;
  TBranch *bEventList_ETwrdidT;
  TBranch *bEventList_ETwradc11;
  TBranch *bEventList_ETwreneT0;
  TBranch *bEventList_ETwreT;
  TBranch *bEventList_ETwrENET0;
  TBranch *bEventList_ETwrphT;
  TBranch *bEventList_ETwrPTower;
  TBranch *bEventList_ETwrpidTower;
  TBranch *bEventList_ETwrmoduleT;
  TBranch *bEventList_EClustEneT0;
  TBranch *bEventList_EClustetav1;
  TBranch *bEventList_EClustphiv1;
  TBranch *bEventList_EEstrpen01;
  TBranch *bEventList_EEstrpen02;
  TBranch *bEventList_EEstrpen03;
  TBranch *bEventList_EEstrpen0;
  TBranch *bEventList_EEstrpen1;
  TBranch *bEventList_EEstrpen2;
  TBranch *bEventList_EEstrpen3;
  TBranch *bEventList_EEstrpen4;
  TBranch *bEventList_EEstrpen5;
  TBranch *bEventList_EEstrpen6;
  TBranch *bEventList_EEstrpen7;
  TBranch *bEventList_EEstrpen8;
  TBranch *bEventList_EEstrpen9;
  TBranch *bEventList_EEstrpen10;
  TBranch *bEventList_EEstrpen11;
  TBranch *bEventList_EEstrpen12;
  TBranch *bEventList_EEstrpen13;
  TBranch *bEventList_EEstrpen14;
  TBranch *bEventList_EEstrpen15;
  TBranch *bEventList_ETwrdidE;
  TBranch *bEventList_EPstripenp01;
  TBranch *bEventList_EPstripenp02;
  TBranch *bEventList_EPstripenp03;
  TBranch *bEventList_EPstripenp0;
  TBranch *bEventList_EPstripenp1;
  TBranch *bEventList_EPstripenp2;
  TBranch *bEventList_EPstripenp3;
  TBranch *bEventList_EPstripenp4;
  TBranch *bEventList_EPstripenp5;
  TBranch *bEventList_EPstripenp6;
  TBranch *bEventList_EPstripenp7;
  TBranch *bEventList_EPstripenp8;
  TBranch *bEventList_EPstripenp9;
  TBranch *bEventList_EPstripenp10;
  TBranch *bEventList_EPstripenp11;
  TBranch *bEventList_EPstripenp12;
  TBranch *bEventList_EPstripenp13;
  TBranch *bEventList_EPstripenp14;
  TBranch *bEventList_EPstripenp15;
  TBranch *bEventList_EclustEnnq1;
  TBranch *bEventList_EclustEnnq20;
  TBranch *bEventList_EclustEnnq19;
  TBranch *bEventList_EclustEnpq1;
  TBranch *bEventList_EclustEnpq20;
  TBranch *bEventList_EclustEnpq19;
  TBranch *bEventList_EclustEnpq21;
  TBranch *bEventList_PrimaryTrackArray_;
  TBranch *bPrimaryTrackArray_fUniqueID;
  TBranch *bPrimaryTrackArray_fBits;
  TBranch *bPrimaryTrackArray_nHitsFit;
  TBranch *bPrimaryTrackArray_nHitsPoss;
  TBranch *bPrimaryTrackArray_trackFlag;
  TBranch *bPrimaryTrackArray_pZ;
  TBranch *bPrimaryTrackArray_pX;
  TBranch *bPrimaryTrackArray_pY;
  TBranch *bPrimaryTrackArray_pT;
  TBranch *bPrimaryTrackArray_dEdx;
  TBranch *bPrimaryTrackArray_charge;
  TBranch *bPrimaryTrackArray_tofBeta;
  TBranch *bPrimaryTrackArray_eta;
  TBranch *bPrimaryTrackArray_phi;
  TBranch *bPrimaryTrackArray_nSigElectron;
  TBranch *bPrimaryTrackArray_nSigPion;
  TBranch *bPrimaryTrackArray_nSigKaon;
  TBranch *bPrimaryTrackArray_nSigProton;
  TBranch *bPrimaryTrackArray_dcag;
  TBranch *bPrimaryTrackArray_nHits;
  TBranch *bPrimaryTrackArray_dEdxHits;
  TBranch *bPrimaryTrackArray_firstZPoint;
  TBranch *bPrimaryTrackArray_lastZPoint;
  TBranch *bPrimaryTrackArray_tofSigElectron;
  TBranch *bPrimaryTrackArray_tofSigPion;
  TBranch *bPrimaryTrackArray_tofSigKaon;
  TBranch *bPrimaryTrackArray_tofSigProton;
  TBranch *bPrimaryTrackArray_timeOfflight;
  TBranch *bPrimaryTrackArray_pathLength;
  TBranch *bPrimaryTrackArray_trkIndex;
  TBranch *bEventList_TowerArray_;
  TBranch *bTowerArray_fUniqueID;
  TBranch *bTowerArray_fBits;
  TBranch *bTowerArray_TwrId;
  TBranch *bTowerArray_TwrEng;
  TBranch *bTowerArray_TwrEta;
  TBranch *bTowerArray_TwrPhi;
  TBranch *bTowerArray_TwrADC;
  TBranch *bTowerArray_TwrPed;
  TBranch *bTowerArray_TwrRMS;
  TBranch *bTowerArray_TwrMatchIdnex;
  TBranch *bTowerArray_NoOfmatchedTrk;
  TBranch *bTowerArray_TwrMatchP;
  TBranch *bTowerArray_TwrPx;
  TBranch *bTowerArray_TwrPy;
  TBranch *bTowerArray_TwrPz;
  TBranch *bTowerArray_fNAssocTracks;
  TBranch *bTowerArray_fMatchedTracksArray_;
  TBranch *bTowerArray_fMatchedTracksArray_P;

  // set input branches
  tInput -> SetMakeClass(1);
  tInput -> SetBranchAddress("fUniqueID", &fUniqueID, &bEventList_fUniqueID);
  tInput -> SetBranchAddress("fBits", &fBits, &bEventList_fBits);
  tInput -> SetBranchAddress("runNumber", &runNumber, &bEventList_runNumber);
  tInput -> SetBranchAddress("eventNumber", &eventNumber, &bEventList_eventNumber);
  tInput -> SetBranchAddress("trigID", &trigID, &bEventList_trigID);
  tInput -> SetBranchAddress("nGlobalTracks", &nGlobalTracks, &bEventList_nGlobalTracks);
  tInput -> SetBranchAddress("nPrimaryTracks", &nPrimaryTracks, &bEventList_nPrimaryTracks);
  tInput -> SetBranchAddress("refMult", &refMult, &bEventList_refMult);
  tInput -> SetBranchAddress("vpdVz", &vpdVz, &bEventList_vpdVz);
  tInput -> SetBranchAddress("xVertex", &xVertex, &bEventList_xVertex);
  tInput -> SetBranchAddress("yVertex", &yVertex, &bEventList_yVertex);
  tInput -> SetBranchAddress("zVertex", &zVertex, &bEventList_zVertex);
  tInput -> SetBranchAddress("bbcZVertex", &bbcZVertex, &bEventList_bbcZVertex);
  tInput -> SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate, &bEventList_zdcCoincidenceRate);
  tInput -> SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate, &bEventList_bbcCoincidenceRate);
  tInput -> SetBranchAddress("backgroundRate", &backgroundRate, &bEventList_backgroundRate);
  tInput -> SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate, &bEventList_bbcBlueBackgroundRate);
  tInput -> SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate, &bEventList_bbcYellowBackgroundRate);
  tInput -> SetBranchAddress("refMultPos", &refMultPos, &bEventList_refMultPos);
  tInput -> SetBranchAddress("refMultNeg", &refMultNeg, &bEventList_refMultNeg);
  tInput -> SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity, &bEventList_bTOFTrayMultiplicity);
  tInput -> SetBranchAddress("nVerticies", &nVerticies, &bEventList_nVerticies);
  tInput -> SetBranchAddress("MagF", &MagF, &bEventList_MagF);
  tInput -> SetBranchAddress("VrtxRank", &VrtxRank, &bEventList_VrtxRank);
  tInput -> SetBranchAddress("FlagEvent_TrgTrkMisMtch", &FlagEvent_TrgTrkMisMtch, &bEventList_FlagEvent_TrgTrkMisMtch);
  tInput -> SetBranchAddress("Etsp", &Etsp, &bEventList_Etsp);
  tInput -> SetBranchAddress("ETwrdidT", &ETwrdidT, &bEventList_ETwrdidT);
  tInput -> SetBranchAddress("ETwradc11", &ETwradc11, &bEventList_ETwradc11);
  tInput -> SetBranchAddress("ETwreneT0", &ETwreneT0, &bEventList_ETwreneT0);
  tInput -> SetBranchAddress("ETwreT", &ETwreT, &bEventList_ETwreT);
  tInput -> SetBranchAddress("ETwrENET0", &ETwrENET0, &bEventList_ETwrENET0);
  tInput -> SetBranchAddress("ETwrphT", &ETwrphT, &bEventList_ETwrphT);
  tInput -> SetBranchAddress("ETwrPTower", &ETwrPTower, &bEventList_ETwrPTower);
  tInput -> SetBranchAddress("ETwrpidTower", &ETwrpidTower, &bEventList_ETwrpidTower);
  tInput -> SetBranchAddress("ETwrmoduleT", &ETwrmoduleT, &bEventList_ETwrmoduleT);
  tInput -> SetBranchAddress("EClustEneT0", &EClustEneT0, &bEventList_EClustEneT0);
  tInput -> SetBranchAddress("EClustetav1", &EClustetav1, &bEventList_EClustetav1);
  tInput -> SetBranchAddress("EClustphiv1", &EClustphiv1, &bEventList_EClustphiv1);
  tInput -> SetBranchAddress("EEstrpen01", &EEstrpen01, &bEventList_EEstrpen01);
  tInput -> SetBranchAddress("EEstrpen02", &EEstrpen02, &bEventList_EEstrpen02);
  tInput -> SetBranchAddress("EEstrpen03", &EEstrpen03, &bEventList_EEstrpen03);
  tInput -> SetBranchAddress("EEstrpen0", &EEstrpen0, &bEventList_EEstrpen0);
  tInput -> SetBranchAddress("EEstrpen1", &EEstrpen1, &bEventList_EEstrpen1);
  tInput -> SetBranchAddress("EEstrpen2", &EEstrpen2, &bEventList_EEstrpen2);
  tInput -> SetBranchAddress("EEstrpen3", &EEstrpen3, &bEventList_EEstrpen3);
  tInput -> SetBranchAddress("EEstrpen4", &EEstrpen4, &bEventList_EEstrpen4);
  tInput -> SetBranchAddress("EEstrpen5", &EEstrpen5, &bEventList_EEstrpen5);
  tInput -> SetBranchAddress("EEstrpen6", &EEstrpen6, &bEventList_EEstrpen6);
  tInput -> SetBranchAddress("EEstrpen7", &EEstrpen7, &bEventList_EEstrpen7);
  tInput -> SetBranchAddress("EEstrpen8", &EEstrpen8, &bEventList_EEstrpen8);
  tInput -> SetBranchAddress("EEstrpen9", &EEstrpen9, &bEventList_EEstrpen9);
  tInput -> SetBranchAddress("EEstrpen10", &EEstrpen10, &bEventList_EEstrpen10);
  tInput -> SetBranchAddress("EEstrpen11", &EEstrpen11, &bEventList_EEstrpen11);
  tInput -> SetBranchAddress("EEstrpen12", &EEstrpen12, &bEventList_EEstrpen12);
  tInput -> SetBranchAddress("EEstrpen13", &EEstrpen13, &bEventList_EEstrpen13);
  tInput -> SetBranchAddress("EEstrpen14", &EEstrpen14, &bEventList_EEstrpen14);
  tInput -> SetBranchAddress("EEstrpen15", &EEstrpen15, &bEventList_EEstrpen15);
  tInput -> SetBranchAddress("ETwrdidE", &ETwrdidE, &bEventList_ETwrdidE);
  tInput -> SetBranchAddress("EPstripenp01", &EPstripenp01, &bEventList_EPstripenp01);
  tInput -> SetBranchAddress("EPstripenp02", &EPstripenp02, &bEventList_EPstripenp02);
  tInput -> SetBranchAddress("EPstripenp03", &EPstripenp03, &bEventList_EPstripenp03);
  tInput -> SetBranchAddress("EPstripenp0", &EPstripenp0, &bEventList_EPstripenp0);
  tInput -> SetBranchAddress("EPstripenp1", &EPstripenp1, &bEventList_EPstripenp1);
  tInput -> SetBranchAddress("EPstripenp2", &EPstripenp2, &bEventList_EPstripenp2);
  tInput -> SetBranchAddress("EPstripenp3", &EPstripenp3, &bEventList_EPstripenp3);
  tInput -> SetBranchAddress("EPstripenp4", &EPstripenp4, &bEventList_EPstripenp4);
  tInput -> SetBranchAddress("EPstripenp5", &EPstripenp5, &bEventList_EPstripenp5);
  tInput -> SetBranchAddress("EPstripenp6", &EPstripenp6, &bEventList_EPstripenp6);
  tInput -> SetBranchAddress("EPstripenp7", &EPstripenp7, &bEventList_EPstripenp7);
  tInput -> SetBranchAddress("EPstripenp8", &EPstripenp8, &bEventList_EPstripenp8);
  tInput -> SetBranchAddress("EPstripenp9", &EPstripenp9, &bEventList_EPstripenp9);
  tInput -> SetBranchAddress("EPstripenp10", &EPstripenp10, &bEventList_EPstripenp10);
  tInput -> SetBranchAddress("EPstripenp11", &EPstripenp11, &bEventList_EPstripenp11);
  tInput -> SetBranchAddress("EPstripenp12", &EPstripenp12, &bEventList_EPstripenp12);
  tInput -> SetBranchAddress("EPstripenp13", &EPstripenp13, &bEventList_EPstripenp13);
  tInput -> SetBranchAddress("EPstripenp14", &EPstripenp14, &bEventList_EPstripenp14);
  tInput -> SetBranchAddress("EPstripenp15", &EPstripenp15, &bEventList_EPstripenp15);
  tInput -> SetBranchAddress("EclustEnnq1", &EclustEnnq1, &bEventList_EclustEnnq1);
  tInput -> SetBranchAddress("EclustEnnq20", &EclustEnnq20, &bEventList_EclustEnnq20);
  tInput -> SetBranchAddress("EclustEnnq19", &EclustEnnq19, &bEventList_EclustEnnq19);
  tInput -> SetBranchAddress("EclustEnpq1", &EclustEnpq1, &bEventList_EclustEnpq1);
  tInput -> SetBranchAddress("EclustEnpq20", &EclustEnpq20, &bEventList_EclustEnpq20);
  tInput -> SetBranchAddress("EclustEnpq19", &EclustEnpq19, &bEventList_EclustEnpq19);
  tInput -> SetBranchAddress("EclustEnpq21", &EclustEnpq21, &bEventList_EclustEnpq21);
  tInput -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_, &bEventList_PrimaryTrackArray_);
  tInput -> SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID, &bPrimaryTrackArray_fUniqueID);
  tInput -> SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits, &bPrimaryTrackArray_fBits);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit, &bPrimaryTrackArray_nHitsFit);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss, &bPrimaryTrackArray_nHitsPoss);
  tInput -> SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag, &bPrimaryTrackArray_trackFlag);
  tInput -> SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ, &bPrimaryTrackArray_pZ);
  tInput -> SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX, &bPrimaryTrackArray_pX);
  tInput -> SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY, &bPrimaryTrackArray_pY);
  tInput -> SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT, &bPrimaryTrackArray_pT);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx, &bPrimaryTrackArray_dEdx);
  tInput -> SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge, &bPrimaryTrackArray_charge);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta, &bPrimaryTrackArray_tofBeta);
  tInput -> SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta, &bPrimaryTrackArray_eta);
  tInput -> SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi, &bPrimaryTrackArray_phi);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron, &bPrimaryTrackArray_nSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion, &bPrimaryTrackArray_nSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon, &bPrimaryTrackArray_nSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton, &bPrimaryTrackArray_nSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag, &bPrimaryTrackArray_dcag);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits, &bPrimaryTrackArray_nHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits, &bPrimaryTrackArray_dEdxHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint, &bPrimaryTrackArray_firstZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint, &bPrimaryTrackArray_lastZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron, &bPrimaryTrackArray_tofSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion, &bPrimaryTrackArray_tofSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon, &bPrimaryTrackArray_tofSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton, &bPrimaryTrackArray_tofSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight, &bPrimaryTrackArray_timeOfflight);
  tInput -> SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength, &bPrimaryTrackArray_pathLength);
  tInput -> SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex, &bPrimaryTrackArray_trkIndex);
  tInput -> SetBranchAddress("TowerArray", &TowerArray_, &bEventList_TowerArray_);
  tInput -> SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID, &bTowerArray_fUniqueID);
  tInput -> SetBranchAddress("TowerArray.fBits", TowerArray_fBits, &bTowerArray_fBits);
  tInput -> SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId, &bTowerArray_TwrId);
  tInput -> SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng, &bTowerArray_TwrEng);
  tInput -> SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta, &bTowerArray_TwrEta);
  tInput -> SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi, &bTowerArray_TwrPhi);
  tInput -> SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC, &bTowerArray_TwrADC);
  tInput -> SetBranchAddress("TowerArray.TwrPed", TowerArray_TwrPed, &bTowerArray_TwrPed);
  tInput -> SetBranchAddress("TowerArray.TwrRMS", TowerArray_TwrRMS, &bTowerArray_TwrRMS);
  tInput -> SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex, &bTowerArray_TwrMatchIdnex);
  tInput -> SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk, &bTowerArray_NoOfmatchedTrk);
  tInput -> SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP, &bTowerArray_TwrMatchP);
  tInput -> SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx, &bTowerArray_TwrPx);
  tInput -> SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy, &bTowerArray_TwrPy);
  tInput -> SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz, &bTowerArray_TwrPz);
  tInput -> SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks, &bTowerArray_fNAssocTracks);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_, &bTowerArray_fMatchedTracksArray_);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P, &bTowerArray_fMatchedTracksArray_P);
  cout << "    Set branches." << endl;


  // define bad run and hot tower lists
  const UInt_t badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};
  const UInt_t hotTwrList[NHotTwr] = {1, 35, 141, 187, 224, 341, 424, 594, 814, 899, 900, 1046, 1128, 1132, 1244, 1382, 1388, 1405, 1588, 1766, 1773, 2066, 2160, 2253, 2281, 2284, 2301, 2303, 2306, 2590, 3007, 3495, 3840, 4043, 4047, 4053, 4057, 4121, 4442, 4569, 4617};
  cout << "    Bad run and hot tower lists defined:\n"
       << "      " << NBadRuns << " bad runs, " << NHotTwr << " hot towers."
       << endl;


  // define histograms
  TH1D *hPhiTrk[NTrgs][NCut];
  TH1D *hPhiRaw[NTrgs][NCut];
  TH1D *hEtaTrk[NTrgs][NCut];
  TH1D *hEtaRaw[NTrgs][NCut];
  TH1D *hPtTrk[NTrgs][NCut];
  TH1D *hPtRaw[NTrgs][NCut];
  TH1D *hPhiEff[NTrgs];
  TH1D *hEtaEff[NTrgs];
  TH1D *hPtEff[NTrgs];

  // binning
  const UInt_t   nF    = 60;
  const UInt_t   nH    = 40;
  const UInt_t   nPt   = 40;
  const Double_t f[2]  = {-3.15, 3.15};
  const Double_t h[2]  = {-1., 1.};
  const Double_t pT[2] = {0., 20.};
  // create histograms
  hPhiTrk[0][0] = new TH1D("hPhiBeforeQA_pi", "", nF, f[0], f[1]);
  hPhiTrk[1][0] = new TH1D("hPhiBeforeQA_ga", "", nF, f[0], f[1]);
  hPhiTrk[0][1] = new TH1D("hPhiAfterQA_pi", "", nF, f[0], f[1]);
  hPhiTrk[1][1] = new TH1D("hPhiAfterQA_ga", "", nF, f[0], f[1]);
  hPhiRaw[0][0] = new TH1D("hPhiRawBeforeQA_pi", "", nF, f[0], f[1]);
  hPhiRaw[1][0] = new TH1D("hPhiRawBeforeQA_ga", "", nF, f[0], f[1]);
  hPhiRaw[0][1] = new TH1D("hPhiRawAfterQA_pi", "", nF, f[0], f[1]);
  hPhiRaw[1][1] = new TH1D("hPhiRawAfterQA_ga", "", nF, f[0], f[1]);
  hEtaTrk[0][0] = new TH1D("hEtaBeforeQA_pi", "", nH, h[0], h[1]);
  hEtaTrk[1][0] = new TH1D("hEtaBeforeQA_ga", "", nH, h[0], h[1]);
  hEtaTrk[0][1] = new TH1D("hEtaAfterQA_pi", "", nH, h[0], h[1]);
  hEtaTrk[1][1] = new TH1D("hEtaAfterQA_ga", "", nH, h[0], h[1]);
  hEtaRaw[0][0] = new TH1D("hEtaRawBeforeQA_pi", "", nH, h[0], h[1]);
  hEtaRaw[1][0] = new TH1D("hEtaRawBeforeQA_ga", "", nH, h[0], h[1]);
  hEtaRaw[0][1] = new TH1D("hEtaRawAfterQA_pi", "", nH, h[0], h[1]);
  hEtaRaw[1][1] = new TH1D("hEtaRawAfterQA_ga", "", nH, h[0], h[1]);
  hPtTrk[0][0]  = new TH1D("hPtBeforeQA_pi", "", nPt, pT[0], pT[1]);
  hPtTrk[1][0]  = new TH1D("hPtBeforeQA_ga", "", nPt, pT[0], pT[1]);
  hPtTrk[0][1]  = new TH1D("hPtAfterQA_pi", "", nPt, pT[0], pT[1]);
  hPtTrk[1][1]  = new TH1D("hPtAfterQA_ga", "", nPt, pT[0], pT[1]);
  hPtRaw[0][0]  = new TH1D("hPtRawBeforeQA_pi", "", nPt, pT[0], pT[1]);
  hPtRaw[1][0]  = new TH1D("hPtRawBeforeQA_ga", "", nPt, pT[0], pT[1]);
  hPtRaw[0][1]  = new TH1D("hPtRawAfterQA_pi", "", nPt, pT[0], pT[1]);
  hPtRaw[1][1]  = new TH1D("hPtRawAfterQA_ga", "", nPt, pT[0], pT[1]);
  hPhiEff[0]    = new TH1D("hPhiEfficiency_pi", "", nF, f[0], f[1]);
  hPhiEff[1]    = new TH1D("hPhiEfficiency_ga", "", nF, f[0], f[1]);
  hEtaEff[0]    = new TH1D("hEtaEfficiency_pi", "", nH, h[0], h[1]);
  hEtaEff[1]    = new TH1D("hEtaEfficiency_ga", "", nH, h[0], h[1]);
  hPtEff[0]     = new TH1D("hPtEfficiency_pi", "", nPt, pT[0], pT[1]);
  hPtEff[1]     = new TH1D("hPtEfficiency_ga", "", nPt, pT[0], pT[1]);
  // errors
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    hPhiEff[iTrg] -> Sumw2();
    hEtaEff[iTrg] -> Sumw2();
    hPtEff[iTrg]  -> Sumw2();
    for (UInt_t iCut = 0; iCut < NCut; iCut++) {
      hPhiTrk[iTrg][iCut] -> Sumw2();
      hPhiRaw[iTrg][iCut] -> Sumw2();
      hEtaTrk[iTrg][iCut] -> Sumw2();
      hEtaRaw[iTrg][iCut] -> Sumw2();
      hPtTrk[iTrg][iCut]  -> Sumw2();
      hPtRaw[iTrg][iCut]  -> Sumw2();
    }
  }



  const UInt_t nEvts = tInput -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " events to process." << endl;

  // no. of triggers
  UInt_t   nTrgBin[NTrgs];
  Double_t nTrgScaled[NTrgs];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    nTrgBin[iTrg]    = 0;
    nTrgScaled[iTrg] = 0.;
  }

  // event loop
  UInt_t bytes(0);
  UInt_t nBytes(0);
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    // load entry
    bytes   = tInput -> GetEntry(iEvt);
    nBytes += bytes;
    if (bytes < 0) {
      cerr << "WARNING: issue with entry " << iEvt << "!" << endl;
      break;
    }
    else {
      if (isInBatchMode) {
        cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
      }
      else {
        cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
        if ((iEvt + 1) == nEvts) cout << endl;
      }
    }

    // event info
    const UInt_t   run   = runNumber;
    const Long64_t nTrks = nPrimaryTracks;
    const Double_t xVtx  = xVertex;
    const Double_t yVtx  = yVertex;
    const Double_t zVtx  = zVertex;
    const Double_t rVtx  = TMath::Sqrt((xVtx * xVtx) + (yVtx * yVtx));

    Bool_t isGoodRun = true;
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (run == badRunList[iRun]) {
        isGoodRun = false;
        break;
      }
    }

    // event cuts
    const Bool_t isInRcut = (TMath::Abs(rVtx) < rVtxMax);
    const Bool_t isInZcut = (TMath::Abs(zVtx) < zVtxMax);
    if (!isGoodRun || !isInRcut || !isInZcut) continue;


    // trigger info
    const Int_t    adc   = ETwradc11;
    const UInt_t   idTrg = ETwrdidT;
    const Double_t tsp   = Etsp;
    const Double_t eH4   = EEstrpen4;
    const Double_t eF4   = EPstripenp4;
    const Double_t pProj = ETwrPTower;
    const Double_t hDet  = ETwreT;
    const Double_t hPhys = EClustetav1;
    const Double_t fTrg  = EClustphiv1;
    const Double_t eTrg  = EClustEneT0;
    const Double_t tTrg  = 2. * TMath::ATan(TMath::Exp(-1. * hPhys));
    const Double_t eTtrg = eTrg * TMath::Sin(tTrg);

    Bool_t isGoodTwr = true;
    for (UInt_t iTwr = 0; iTwr < NHotTwr; iTwr++) {
      if (idTrg == hotTwrList[iTwr]) {
        isGoodTwr = false;
        break;
      }
    }

    // trigger cuts
    Bool_t isPi0     = true;
    Bool_t isGam     = true;
    Bool_t isGoodTrg = true;
    if (isTriggered) {
      const Bool_t isInAdcCut    = (adc <= adcMax);
      const Bool_t isInStrCut    = ((eH4 >= eStrMin) && (eF4 >= eStrMin));
      const Bool_t isInProjCut   = (pProj < pProjMax);
      const Bool_t isInEtaTrgCut = (TMath::Abs(hDet) < hTrgMax);
      const Bool_t isInEtCut     = ((eTtrg > eTtrgMin) && (eTtrg < eTtrgMax));
      const Bool_t isInPi0cut    = ((tsp > tspPi0[0]) && (tsp < tspPi0[1]));
      const Bool_t isInGamCut    = ((tsp > tspGam[0]) && (tsp < tspGam[1]));
      const Bool_t isInTspCut    = (isInPi0cut || isInGamCut);
      if (!isGoodTwr || !isInAdcCut || !isInStrCut || !isInProjCut || !isInEtaTrgCut || !isInEtCut || !isInTspCut) {
        isGoodTrg = false;
      }
      else {
        isPi0     = isInPi0cut;
        isGam     = isInGamCut;
        isGoodTrg = true;
      }
    }  // end trigger check

    if (isGoodTrg) {
      if (isPi0) nTrgBin[0]++;
      if (isGam) nTrgBin[1]++;
    }
    else {
      continue;
    }


    // track loop
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const UInt_t   nFit  = PrimaryTrackArray_nHitsFit[iTrk];
      const UInt_t   nPoss = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFit  = (Double_t) nFit / (Double_t) nPoss;
      const Double_t dca   = PrimaryTrackArray_dcag[iTrk];
      const Double_t chrg  = PrimaryTrackArray_charge[iTrk];
      const Double_t hTrk  = PrimaryTrackArray_eta[iTrk];
      const Double_t pXtrk = PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk = PrimaryTrackArray_pY[iTrk];
      const Double_t pTtrk = PrimaryTrackArray_pT[iTrk];

      // calculate phi (if need be)
      Double_t fCalc = TMath::ATan(pYtrk / pXtrk);
      if (!isDetectorLevel) {
        Bool_t isInQuad1  = ((pXtrk > 0.) && (pYtrk > 0.));
        Bool_t isInQuad2  = ((pXtrk < 0.) && (pYtrk > 0.));
        Bool_t isInQuad3  = ((pXtrk < 0.) && (pYtrk < 0.));
        Bool_t isInQuad4  = ((pXtrk > 0.) && (pYtrk < 0.));
        Bool_t isInBottom = (fCalc < 0.);
        if (isInQuad2 && isInBottom)  fCalc += TMath::Pi();
        if (isInQuad3 && !isInBottom) fCalc -= TMath::Pi();
      }

      Double_t fTrk = PrimaryTrackArray_phi[iTrk];
      if (!isDetectorLevel) fTrk = fCalc;


      // fill before histograms
      if (isPi0) {
        hPhiTrk[0][0] -> Fill(fTrk);
        hPhiRaw[0][0] -> Fill(fTrk);
        hEtaTrk[0][0] -> Fill(hTrk);
        hEtaRaw[0][0] -> Fill(hTrk);
        hPtTrk[0][0]  -> Fill(pTtrk);
        hPtRaw[0][0]  -> Fill(pTtrk);
      }
      if (isGam) {
        hPhiTrk[1][0] -> Fill(fTrk);
        hPhiRaw[1][0] -> Fill(fTrk);
        hEtaTrk[1][0] -> Fill(hTrk);
        hEtaRaw[1][0] -> Fill(hTrk);
        hPtTrk[1][0]  -> Fill(pTtrk);
        hPtRaw[1][0]  -> Fill(pTtrk);
      }


      // track cuts
      const Bool_t isInFitCut    = (nFit > nFitMin);
      const Bool_t isInRatioCut  = (rFit > rFitMin);
      const Bool_t isInDcaCut    = (dca < dcaMax);
      const Bool_t isInChrgCut   = (chrg != 0.);
      const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtCut     = ((pTtrk > pTtrkMin) && (pTtrk < pTtrkMax));

      Bool_t isGoodTrk = true;
      if (isDetectorLevel) {
        if (!isInFitCut || !isInRatioCut || !isInDcaCut || !isInEtaTrkCut || !isInPtCut)
          isGoodTrk = false;
      }
      else {
        if (!isInChrgCut || !isInEtaTrkCut || !isInPtCut)
          isGoodTrk = false;
      }
      if (!isGoodTrk) continue;


      // fill after histograms
      if (isPi0) {
        hPhiTrk[0][1] -> Fill(fTrk);
        hPhiRaw[0][1] -> Fill(fTrk);
        hEtaTrk[0][1] -> Fill(hTrk);
        hEtaRaw[0][1] -> Fill(hTrk);
        hPtTrk[0][1]  -> Fill(pTtrk);
        hPtRaw[0][1]  -> Fill(pTtrk);
      }
      if (isGam) {
        hPhiTrk[1][1] -> Fill(fTrk);
        hPhiRaw[1][1] -> Fill(fTrk);
        hEtaTrk[1][1] -> Fill(hTrk);
        hEtaRaw[1][1] -> Fill(hTrk);
        hPtTrk[1][1]  -> Fill(pTtrk);
        hPtRaw[1][1]  -> Fill(pTtrk);
      }

    }  // end track loop
  }  // end event loop

  cout << "    Event loop finished\n"
       << "      " << nTrgBin[0] << " pi0 triggers,\n"
       << "      " << nTrgBin[1] << " gamma triggers."
       << endl;

  // assign triggers
  nTrgPi0 = nTrgBin[0];
  nTrgGam = nTrgBin[1];

  // calculate ratios
  const Float_t weight(1.);
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    hPhiEff[iTrg] -> Divide(hPhiTrk[iTrg][1], hPhiTrk[iTrg][0], weight, weight);
    hEtaEff[iTrg] -> Divide(hEtaTrk[iTrg][1], hEtaTrk[iTrg][0], weight, weight);
    hPtEff[iTrg]  -> Divide(hPtTrk[iTrg][1], hPtTrk[iTrg][0], weight, weight);
  }
  cout << "    Ratios calculated." << endl;


  // normalize histograms
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    for (UInt_t iCut = 0; iCut < NCut; iCut++) {
      hPhiTrk[iTrg][iCut] -> Scale(1. / (Double_t) nTrgBin[iTrg]);
      hEtaTrk[iTrg][iCut] -> Scale(1. / (Double_t) nTrgBin[iTrg]);
      hPtTrk[iTrg][iCut]  -> Scale(1. / (Double_t) nTrgBin[iTrg]);
    }
  }
  cout << "    Normalized track distributions." << endl;


  // fit efficiencies
  const UInt_t   fLinF(1);
  const UInt_t   fSizF(2);
  const UInt_t   fColF[NTrgs] = {859, 899};
  const TString  sTrg[NTrgs]  = {"Pi", "Ga"};
  const Double_t fEffGuess(0.5);
  const Double_t hEffGuess(0.5);
  const Double_t pEffGuess(0.87);
  const Double_t pSigGuess(4.0);
  const Double_t fFitRange[2] = {f[0], f[1]};
  const Double_t hFitRange[2] = {-0.7, 0.7};
  const Double_t pFitRange[2] = {1., 7.};

  TF1      *fPhiEff[NTrgs];
  TF1      *fEtaEff[NTrgs];
  TF1      *fPtEff[NTrgs];
  Double_t phiEff[NTrgs][NVal];
  Double_t etaEff[NTrgs][NVal];
  Double_t ptEff[NTrgs][NVal];
  Double_t ptSig[NTrgs][NVal];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {

    // make names
    TString sPhi("fPhiEff");
    TString sEta("fEtaEff");
    TString sPt("fPtEff");
    sPhi += sTrg[iTrg];
    sEta += sTrg[iTrg];
    sPt  += sTrg[iTrg];

    // define fits
    fPhiEff[iTrg] = new TF1(sPhi.Data(), "[0]", f[0], f[1]);
    fEtaEff[iTrg] = new TF1(sEta.Data(), "[0]", h[0], h[1]);
    fPtEff[iTrg]  = new TF1(sPt.Data(), "[0]*(1-exp(-1.*[1]*x))", pT[0], pT[1]);
    fPhiEff[iTrg] -> SetParameter(0, fEffGuess);
    fPhiEff[iTrg] -> SetLineColor(fColF[iTrg]);
    fPhiEff[iTrg] -> SetLineStyle(fLinF);
    fPhiEff[iTrg] -> SetLineWidth(fSizF);
    fEtaEff[iTrg] -> SetParameter(0, hEffGuess);
    fEtaEff[iTrg] -> SetLineColor(fColF[iTrg]);
    fEtaEff[iTrg] -> SetLineStyle(fLinF);
    fEtaEff[iTrg] -> SetLineWidth(fSizF);
    fPtEff[iTrg]  -> SetParameters(0, pEffGuess);
    fPtEff[iTrg]  -> SetParameters(1, pSigGuess);
    fPtEff[iTrg]  -> SetLineColor(fColF[iTrg]);
    fPtEff[iTrg]  -> SetLineStyle(fLinF);
    fPtEff[iTrg]  -> SetLineWidth(fSizF);

    // fit histograms
    hPhiEff[iTrg] -> Fit(fPhiEff[iTrg], "BMQ0", "", fFitRange[0], fFitRange[1]);
    hEtaEff[iTrg] -> Fit(fEtaEff[iTrg], "BMQ0", "", hFitRange[0], hFitRange[1]);
    hPtEff[iTrg]  -> Fit(fPtEff[iTrg], "BMQ0", "", pFitRange[0], pFitRange[1]);
    phiEff[iTrg][0] = fPhiEff[iTrg] -> GetParameter(0);
    phiEff[iTrg][1] = fPhiEff[iTrg] -> GetParError(0);
    etaEff[iTrg][0] = fEtaEff[iTrg] -> GetParameter(0);
    etaEff[iTrg][1] = fEtaEff[iTrg] -> GetParError(0);
    ptEff[iTrg][0]  = fPtEff[iTrg]  -> GetParameter(0);
    ptEff[iTrg][1]  = fPtEff[iTrg]  -> GetParError(0);
    ptSig[iTrg][0]  = fPtEff[iTrg]  -> GetParameter(1);
    ptSig[iTrg][1]  = fPtEff[iTrg]  -> GetParError(1);

    // reset visibility
    const Int_t kNotDraw = 1<<9;
    if (hPhiEff[iTrg] -> GetFunction(sPhi.Data()))
      hPhiEff[iTrg] -> GetFunction(sPhi.Data()) -> ResetBit(kNotDraw);
    if (hEtaEff[iTrg] -> GetFunction(sEta.Data()))
      hEtaEff[iTrg] -> GetFunction(sEta.Data()) -> ResetBit(kNotDraw);
    if (hPtEff[iTrg] -> GetFunction(sPt.Data()))
      hPtEff[iTrg] -> GetFunction(sPt.Data()) -> ResetBit(kNotDraw);

  }
  cout << "    Fit efficiencies." << endl;


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColE(1);
  const UInt_t  fMarE(1);
  const UInt_t  fColD[NCut] = {1, 879};
  const UInt_t  fMarD[NCut] = {7, 4};
  const Float_t fLab(0.03);
  const Float_t fOffsetX(1.);
  const Float_t fPhiRangeY[2] = {0., 0.47};
  const Float_t fEtaRangeY[2] = {0., 0.73};
  const Float_t fPtRangeY[2]  = {0.00003, 13.};
  const Float_t fEffRangeY[2] = {0., 1.13};
  const TString sTitleXF("#varphi^{trk}");
  const TString sTitleXH("#eta^{trk}");
  const TString sTitleXP("p_{T}^{trk}");
  const TString sTitleYF("(1/N^{trg}) dN^{trk}/d#varphi^{trk}");
  const TString sTitleYFE("#epsilon(#varphi^{trk}) = #epsilon_{#varphi}");
  const TString sTitleYH("(1/N^{trg}) dN^{trk}/d#eta^{trk}");
  const TString sTitleYHE("#epsilon(#eta^{trk}) = #epsilon_{#eta}");
  const TString sTitleYP("(1/N^{trg}) dN^{trk}/dp_{T}^{trk}");
  const TString sTitleYPE("#epsilon(p_{T}^{trk}) = #epsilon_{p} #upoint (1 - exp(-#sigma #upoint p_{T}^{trk}))");
  const TString sTitleF[NTrgs] = {"Track #varphi: #pi^{0} trigger", "Track #varphi: #gamma^{rich}"};
  const TString sTitleH[NTrgs] = {"Track #eta: #pi^{0} trigger", "Track #eta: #gamma^{rich}"};
  const TString sTitleP[NTrgs] = {"Track p_{T}: #pi^{0}", "Track p_{T}: #gamma^{rich}"};
  const TString sTitleE[NTrgs] = {"Track efficiency, #epsilon: #pi^{0} trigger", "Track efficiency, #epsilon: #gamma^{rich}"};
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    // phi efficiency histogram
    hPhiEff[iTrg] -> SetLineColor(fColE);
    hPhiEff[iTrg] -> SetMarkerColor(fColE);
    hPhiEff[iTrg] -> SetMarkerStyle(fMarE);
    hPhiEff[iTrg] -> SetTitle(sTitleE[iTrg].Data());
    hPhiEff[iTrg] -> SetTitleFont(fTxt);
    hPhiEff[iTrg] -> GetXaxis() -> SetTitle(sTitleXF.Data());
    hPhiEff[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hPhiEff[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hPhiEff[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hPhiEff[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hPhiEff[iTrg] -> GetYaxis() -> SetTitle(sTitleYFE.Data());
    hPhiEff[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hPhiEff[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hPhiEff[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hPhiEff[iTrg] -> GetYaxis() -> SetRangeUser(fEffRangeY[0], fEffRangeY[1]);
    // eta efficiency histogram
    hEtaEff[iTrg] -> SetLineColor(fColE);
    hEtaEff[iTrg] -> SetMarkerColor(fColE);
    hEtaEff[iTrg] -> SetMarkerStyle(fMarE);
    hEtaEff[iTrg] -> SetTitle(sTitleE[iTrg].Data());
    hEtaEff[iTrg] -> SetTitleFont(fTxt);
    hEtaEff[iTrg] -> GetXaxis() -> SetTitle(sTitleXH.Data());
    hEtaEff[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hEtaEff[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hEtaEff[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hEtaEff[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hEtaEff[iTrg] -> GetYaxis() -> SetTitle(sTitleYHE.Data());
    hEtaEff[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hEtaEff[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hEtaEff[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hEtaEff[iTrg] -> GetYaxis() -> SetRangeUser(fEffRangeY[0], fEffRangeY[1]);
    // pT efficiency histogram
    hPtEff[iTrg] -> SetLineColor(fColE);
    hPtEff[iTrg] -> SetMarkerColor(fColE);
    hPtEff[iTrg] -> SetMarkerStyle(fMarE);
    hPtEff[iTrg] -> SetTitle(sTitleE[iTrg].Data());
    hPtEff[iTrg] -> SetTitleFont(fTxt);
    hPtEff[iTrg] -> GetXaxis() -> SetTitle(sTitleXP.Data());
    hPtEff[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtEff[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hPtEff[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hPtEff[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hPtEff[iTrg] -> GetYaxis() -> SetTitle(sTitleYPE.Data());
    hPtEff[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtEff[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hPtEff[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hPtEff[iTrg] -> GetYaxis() -> SetRangeUser(fEffRangeY[0], fEffRangeY[1]);
    for (UInt_t iCut = 0; iCut < NCut; iCut++) {
      // phi histograms
      hPhiTrk[iTrg][iCut] -> SetLineColor(fColD[iCut]);
      hPhiTrk[iTrg][iCut] -> SetMarkerColor(fColD[iCut]);
      hPhiTrk[iTrg][iCut] -> SetMarkerStyle(fMarD[iCut]);
      hPhiTrk[iTrg][iCut] -> SetTitle(sTitleF[iTrg].Data());
      hPhiTrk[iTrg][iCut] -> SetTitleFont(fTxt);
      hPhiTrk[iTrg][iCut] -> GetXaxis() -> SetTitle(sTitleXF.Data());
      hPhiTrk[iTrg][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
      hPhiTrk[iTrg][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPhiTrk[iTrg][iCut] -> GetXaxis() -> CenterTitle(fCnt);
      hPhiTrk[iTrg][iCut] -> GetXaxis() -> SetLabelSize(fLab);
      hPhiTrk[iTrg][iCut] -> GetYaxis() -> SetTitle(sTitleYF.Data());
      hPhiTrk[iTrg][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
      hPhiTrk[iTrg][iCut] -> GetYaxis() -> CenterTitle(fCnt);
      hPhiTrk[iTrg][iCut] -> GetYaxis() -> SetLabelSize(fLab);
      hPhiTrk[iTrg][iCut] -> GetYaxis() -> SetRangeUser(fPhiRangeY[0], fPhiRangeY[1]);
      // eta histograms
      hEtaTrk[iTrg][iCut] -> SetLineColor(fColD[iCut]);
      hEtaTrk[iTrg][iCut] -> SetMarkerColor(fColD[iCut]);
      hEtaTrk[iTrg][iCut] -> SetMarkerStyle(fMarD[iCut]);
      hEtaTrk[iTrg][iCut] -> SetTitle(sTitleH[iTrg].Data());
      hEtaTrk[iTrg][iCut] -> SetTitleFont(fTxt);
      hEtaTrk[iTrg][iCut] -> GetXaxis() -> SetTitle(sTitleXH.Data());
      hEtaTrk[iTrg][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
      hEtaTrk[iTrg][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hEtaTrk[iTrg][iCut] -> GetXaxis() -> CenterTitle(fCnt);
      hEtaTrk[iTrg][iCut] -> GetXaxis() -> SetLabelSize(fLab);
      hEtaTrk[iTrg][iCut] -> GetYaxis() -> SetTitle(sTitleYH.Data());
      hEtaTrk[iTrg][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
      hEtaTrk[iTrg][iCut] -> GetYaxis() -> CenterTitle(fCnt);
      hEtaTrk[iTrg][iCut] -> GetYaxis() -> SetLabelSize(fLab);
      hEtaTrk[iTrg][iCut] -> GetYaxis() -> SetRangeUser(fEtaRangeY[0], fEtaRangeY[1]);
      // pT histograms
      hPtTrk[iTrg][iCut] -> SetLineColor(fColD[iCut]);
      hPtTrk[iTrg][iCut] -> SetMarkerColor(fColD[iCut]);
      hPtTrk[iTrg][iCut] -> SetMarkerStyle(fMarD[iCut]);
      hPtTrk[iTrg][iCut] -> SetTitle(sTitleP[iTrg].Data());
      hPtTrk[iTrg][iCut] -> SetTitleFont(fTxt);
      hPtTrk[iTrg][iCut] -> GetXaxis() -> SetTitle(sTitleXP.Data());
      hPtTrk[iTrg][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
      hPtTrk[iTrg][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPtTrk[iTrg][iCut] -> GetXaxis() -> CenterTitle(fCnt);
      hPtTrk[iTrg][iCut] -> GetXaxis() -> SetLabelSize(fLab);
      hPtTrk[iTrg][iCut] -> GetYaxis() -> SetTitle(sTitleYP.Data());
      hPtTrk[iTrg][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
      hPtTrk[iTrg][iCut] -> GetYaxis() -> CenterTitle(fCnt);
      hPtTrk[iTrg][iCut] -> GetYaxis() -> SetLabelSize(fLab);
      hPtTrk[iTrg][iCut] -> GetYaxis() -> SetRangeUser(fPtRangeY[0], fPtRangeY[1]);
    }  // end cut loop
  }  // end trigger loop
  cout << "    Set styles." << endl;


  // make labels
  const UInt_t  nDec(3);
  const UInt_t  fAlign(12);
  const UInt_t  fColP(0);
  const UInt_t  fColT[NTrgs] = {859, 899};
  const Float_t xPav[2]      = {0.7, 0.9};
  const Float_t xLeg[2]      = {0.7, 0.9};
  const Float_t yPav[2]      = {0.7, 0.9};
  const Float_t yLeg[2]      = {0.5, 0.7};
  const TString sSystem("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrgKin("E_{T}^{trg} #in (9, 20) GeV, |#eta^{trg}| < 0.9");
  const TString sPhiEff("#epsilon_{#varphi} = ");
  const TString sEtaEff("#epsilon_{#eta} = ");
  const TString sPtEff("#epsilon_{p} = ");
  const TString sPtSig("#sigma_{p} = ");
  const TString sLegend[2] = {"before QA cuts", "after QA cuts"};

  TLegend   *lPhi[NTrgs];
  TLegend   *lEta[NTrgs];
  TLegend   *lPt[NTrgs];
  TPaveText *pPhi[NTrgs];
  TPaveText *pEta[NTrgs];
  TPaveText *pPt[NTrgs];
  TString sFraw[NVal];
  TString sHraw[NVal];
  TString sPEraw[NVal];
  TString sPSraw[NVal];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {

    TString sFtxt(sPhiEff.Data());
    TString sHtxt(sEtaEff.Data());
    TString sPEtxt(sPtEff.Data());
    TString sPStxt(sPtSig.Data());
    for (UInt_t iVal = 0; iVal < NVal; iVal++) {
      sFraw[iVal]   = "";
      sHraw[iVal]   = "";
      sPEraw[iVal]  = "";
      sPSraw[iVal]  = "";
      sFraw[iVal]  += phiEff[iTrg][iVal];
      sHraw[iVal]  += etaEff[iTrg][iVal];
      sPEraw[iVal] += ptEff[iTrg][iVal];
      sPSraw[iVal] += ptSig[iTrg][iVal];

      const UInt_t nFraw  = sFraw[iVal].First(".");
      const UInt_t nHraw  = sHraw[iVal].First(".");
      const UInt_t nPEraw = sPEraw[iVal].First(".");
      const UInt_t nPSraw = sPSraw[iVal].First(".");
      const UInt_t nFtxt  = (nFraw + nDec) + 1;
      const UInt_t nHtxt  = (nHraw + nDec) + 1;
      const UInt_t nPEtxt = (nPEraw + nDec) + 1;
      const UInt_t nPStxt = (nPSraw + nDec) + 1;
      if (iVal == 0) {
        sFtxt.Append(sFraw[iVal].Data(), nFtxt);
        sHtxt.Append(sHraw[iVal].Data(), nHtxt);
        sPEtxt.Append(sPEraw[iVal].Data(), nPEtxt);
        sPStxt.Append(sPSraw[iVal].Data(), nPStxt);
      } 
      else {
        sFtxt  += " #pm ";
        sHtxt  += " #pm ";
        sPEtxt += " #pm ";
        sPStxt += " #pm ";
        sFtxt.Append(sFraw[iVal].Data(), nFtxt);
        sHtxt.Append(sHraw[iVal].Data(), nHtxt);
        sPEtxt.Append(sPEraw[iVal].Data(), nPEtxt);
        sPStxt.Append(sPSraw[iVal].Data(), nPStxt);
      }
    }  // end value loop

    // phi labels
    pPhi[iTrg] = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
    pPhi[iTrg] -> SetFillColor(fColP);
    pPhi[iTrg] -> SetLineColor(fColP);
    pPhi[iTrg] -> SetTextColor(fColT[iTrg]);
    pPhi[iTrg] -> SetTextFont(fTxt);
    pPhi[iTrg] -> SetTextAlign(fAlign);
    pPhi[iTrg] -> AddText(sSystem.Data());
    pPhi[iTrg] -> AddText(sTrgKin.Data());
    pPhi[iTrg] -> AddText(sFtxt.Data());
    // eta labels
    pEta[iTrg] = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
    pEta[iTrg] -> SetFillColor(fColP);
    pEta[iTrg] -> SetLineColor(fColP);
    pEta[iTrg] -> SetTextColor(fColT[iTrg]);
    pEta[iTrg] -> SetTextFont(fTxt);
    pEta[iTrg] -> SetTextAlign(fAlign);
    pEta[iTrg] -> AddText(sSystem.Data());
    pEta[iTrg] -> AddText(sTrgKin.Data());
    pEta[iTrg] -> AddText(sHtxt.Data());
    // pT labels
    pPt[iTrg]  = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
    pPt[iTrg] -> SetFillColor(fColP);
    pPt[iTrg] -> SetLineColor(fColP);
    pPt[iTrg] -> SetTextColor(fColT[iTrg]);
    pPt[iTrg] -> SetTextFont(fTxt);
    pPt[iTrg] -> SetTextAlign(fAlign);
    pPt[iTrg] -> AddText(sSystem.Data());
    pPt[iTrg] -> AddText(sTrgKin.Data());
    pPt[iTrg] -> AddText(sPEtxt.Data());
    pPt[iTrg] -> AddText(sPStxt.Data());


    // phi legend
    lPhi[iTrg] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lPhi[iTrg] -> SetFillColor(fColP);
    lPhi[iTrg] -> SetLineColor(fColP);
    lPhi[iTrg] -> SetTextFont(fTxt);
    lPhi[iTrg] -> AddEntry(hPhiTrk[iTrg][0], sLegend[0].Data());
    lPhi[iTrg] -> AddEntry(hPhiTrk[iTrg][1], sLegend[1].Data());
    // eta legend
    lEta[iTrg] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lEta[iTrg] -> SetFillColor(fColP);
    lEta[iTrg] -> SetLineColor(fColP);
    lEta[iTrg] -> SetTextFont(fTxt);
    lEta[iTrg] -> AddEntry(hEtaTrk[iTrg][0], sLegend[0].Data());
    lEta[iTrg] -> AddEntry(hEtaTrk[iTrg][1], sLegend[1].Data());
    // pT legend
    lPt[iTrg] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lPt[iTrg] -> SetFillColor(fColP);
    lPt[iTrg] -> SetLineColor(fColP);
    lPt[iTrg] -> SetTextFont(fTxt);
    lPt[iTrg] -> AddEntry(hPtTrk[iTrg][0], sLegend[0].Data());
    lPt[iTrg] -> AddEntry(hPtTrk[iTrg][1], sLegend[1].Data());

  }  // end trigger loop
  cout << "    Made labels." << endl;



  // make directories and save histograms
  const TString sDir[NTrgs] = {"pi0", "gamma"};

  TDirectory *dirs[NTrgs];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    dirs[iTrg] = (TDirectory*) fOutput -> mkdir(sDir[iTrg].Data());
    dirs[iTrg] -> cd();
    for (UInt_t iCut = 0; iCut < NCut; iCut++) {
      hPhiTrk[iTrg][iCut] -> Write();
      hPhiRaw[iTrg][iCut] -> Write();
      hEtaTrk[iTrg][iCut] -> Write();
      hEtaRaw[iTrg][iCut] -> Write();
      hPtTrk[iTrg][iCut]  -> Write();
      hPtRaw[iTrg][iCut]  -> Write();
    }
    hPhiEff[iTrg] -> Write();
    hEtaEff[iTrg] -> Write();
    hPtEff[iTrg]  -> Write();
  }
  cout << "    Made directories." << endl;


  // make plots
  const UInt_t  width(1500);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  ticks(1);
  const UInt_t  log(1);
  const Float_t margin(0.);
  const Float_t xPad[NTrgs + 2] = {0., 0.5, 0.5, 1.};
  const TString sPads[NTrgs]    = {"pPi0", "pGamma"};

  TPad    *pPhiTrk[NTrgs];
  TPad    *pEtaTrk[NTrgs];
  TPad    *pPtTrk[NTrgs];
  TPad    *pPhiEff[NTrgs];
  TPad    *pEtaEff[NTrgs];
  TPad    *pPtEff[NTrgs];
  TCanvas *cPhiTrk;
  TCanvas *cEtaTrk;
  TCanvas *cPtTrk;
  TCanvas *cPhiEff;
  TCanvas *cEtaEff;
  TCanvas *cPtEff;
  fOutput -> cd();
  // phi plots
  cPhiTrk    = new TCanvas("cPhiTrk", "", width, height);
  pPhiTrk[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPhiTrk[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPhiTrk[0]    -> SetGrid(grid, grid);
  pPhiTrk[0]    -> SetTicks(ticks, ticks);
  pPhiTrk[0]    -> SetRightMargin(margin);
  pPhiTrk[1]    -> SetGrid(grid, grid);
  pPhiTrk[1]    -> SetTicks(ticks, ticks);
  pPhiTrk[1]    -> SetLeftMargin(margin);
  cPhiTrk       -> cd();
  pPhiTrk[0]    -> Draw();
  pPhiTrk[1]    -> Draw();
  pPhiTrk[0]    -> cd();
  hPhiTrk[0][0] -> Draw();
  hPhiTrk[0][1] -> Draw("same");
  lPhi[0]       -> Draw();
  pPhi[0]       -> Draw();
  pPhiTrk[1]    -> cd();
  hPhiTrk[1][0] -> Draw();
  hPhiTrk[1][1] -> Draw("same");
  lPhi[1]       -> Draw();
  pPhi[1]       -> Draw();
  cPhiTrk       -> Write();
  cPhiTrk       -> Close();
  // eta plots
  cEtaTrk    = new TCanvas("cEtaTrk", "", width, height);
  pEtaTrk[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pEtaTrk[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pEtaTrk[0]    -> SetGrid(grid, grid);
  pEtaTrk[0]    -> SetTicks(ticks, ticks);
  pEtaTrk[0]    -> SetRightMargin(margin);
  pEtaTrk[1]    -> SetGrid(grid, grid);
  pEtaTrk[1]    -> SetTicks(ticks, ticks);
  pEtaTrk[1]    -> SetLeftMargin(margin);
  cEtaTrk       -> cd();
  pEtaTrk[0]    -> Draw();
  pEtaTrk[1]    -> Draw();
  pEtaTrk[0]    -> cd();
  hEtaTrk[0][0] -> Draw();
  hEtaTrk[0][1] -> Draw("same");
  lEta[0]       -> Draw();
  pEta[0]       -> Draw();
  pEtaTrk[1]    -> cd();
  hEtaTrk[1][0] -> Draw();
  hEtaTrk[1][1] -> Draw("same");
  lEta[1]       -> Draw();
  pEta[1]       -> Draw();
  cEtaTrk       -> Write();
  cEtaTrk       -> Close();
  // pT plots
  cPtTrk    = new TCanvas("cPtTrk", "", width, height);
  pPtTrk[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPtTrk[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPtTrk[0]    -> SetGrid(grid, grid);
  pPtTrk[0]    -> SetTicks(ticks, ticks);
  pPtTrk[0]    -> SetRightMargin(margin);
  pPtTrk[0]    -> SetLogy(log);
  pPtTrk[1]    -> SetGrid(grid, grid);
  pPtTrk[1]    -> SetTicks(ticks, ticks);
  pPtTrk[1]    -> SetLeftMargin(margin);
  pPtTrk[1]    -> SetLogy(log);
  cPtTrk       -> cd();
  pPtTrk[0]    -> Draw();
  pPtTrk[1]    -> Draw();
  pPtTrk[0]    -> cd();
  hPtTrk[0][0] -> Draw();
  hPtTrk[0][1] -> Draw("same");
  lPt[0]       -> Draw();
  pPt[0]       -> Draw();
  pPtTrk[1]    -> cd();
  hPtTrk[1][0] -> Draw();
  hPtTrk[1][1] -> Draw("same");
  lPt[1]       -> Draw();
  pPt[1]       -> Draw();
  cPtTrk       -> Write();
  cPtTrk       -> Close();
  // phi efficiency plots
  cPhiEff    = new TCanvas("cPhiEff", "", width, height);
  pPhiEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPhiEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPhiEff[0] -> SetGrid(grid, grid);
  pPhiEff[0] -> SetTicks(ticks, ticks);
  pPhiEff[0] -> SetRightMargin(margin);
  pPhiEff[1] -> SetGrid(grid, grid);
  pPhiEff[1] -> SetTicks(ticks, ticks);
  pPhiEff[1] -> SetLeftMargin(margin);
  cPhiEff    -> cd();
  pPhiEff[0] -> Draw();
  pPhiEff[1] -> Draw();
  pPhiEff[0] -> cd();
  hPhiEff[0] -> Draw();
  pPhi[0]    -> Draw();
  pPhiEff[1] -> cd();
  hPhiEff[1] -> Draw();
  pPhi[1]    -> Draw();
  cPhiEff    -> Write();
  cPhiEff    -> Close();
  // eta efficiency plots
  cEtaEff    = new TCanvas("cEtaEff", "", width, height);
  pEtaEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pEtaEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pEtaEff[0] -> SetGrid(grid, grid);
  pEtaEff[0] -> SetTicks(ticks, ticks);
  pEtaEff[0] -> SetRightMargin(margin);
  pEtaEff[1] -> SetGrid(grid, grid);
  pEtaEff[1] -> SetTicks(ticks, ticks);
  pEtaEff[1] -> SetLeftMargin(margin);
  cEtaEff    -> cd();
  pEtaEff[0] -> Draw();
  pEtaEff[1] -> Draw();
  pEtaEff[0] -> cd();
  hEtaEff[0] -> Draw();
  pEta[0]    -> Draw();
  pEtaEff[1] -> cd();
  hEtaEff[1] -> Draw();
  pEta[1]    -> Draw();
  cEtaEff    -> Write();
  cEtaEff    -> Close();
  // pT efficiency plots
  cPtEff    = new TCanvas("cPtEff", "", width, height);
  pPtEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPtEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPtEff[0] -> SetGrid(grid, grid);
  pPtEff[0] -> SetTicks(ticks, ticks);
  pPtEff[0] -> SetRightMargin(margin);
  pPtEff[1] -> SetGrid(grid, grid);
  pPtEff[1] -> SetTicks(ticks, ticks);
  pPtEff[1] -> SetLeftMargin(margin);
  cPtEff    -> cd();
  pPtEff[0] -> Draw();
  pPtEff[1] -> Draw();
  pPtEff[0] -> cd();
  hPtEff[0] -> Draw();
  pPt[0]    -> Draw();
  pPtEff[1] -> cd();
  hPtEff[1] -> Draw();
  pPt[1]    -> Draw();
  cPtEff    -> Write();
  cPtEff    -> Close();
  cout << "    Drew plots." << endl;


  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Calculation finished!\n" << endl;

}

// End ------------------------------------------------------------------------
