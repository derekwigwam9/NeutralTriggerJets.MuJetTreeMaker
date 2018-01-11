// 'MakeMuDstJetTree.C'
// Derek Anderson
// 11.13.2017
//
// This produces a tree of jets from the specified
// input with the specified parameters.
//
// NOTE: type = 0, makes charged jets
//       type = 1, makes full jets.


#include "TString.h"
#include "TSystem.h"
#include "TDatime.h"

using namespace std;

class StMuDstJetTreeMaker;


// i/o parameters
static const TString sInput("../../ThirdJetMaker/output/merged/pt35.thirdJetMaker.root");
static const TString sOuput("pp200r12pt35u.r04rm1hc100full.d6m1y2018.root");
// jet parameters
static const UInt_t   type(1);
static const UInt_t   nRepeat(1);
static const UInt_t   nRemove(1);
static const Double_t rJet(0.4);
static const Double_t aGhost(0.01);
static const Double_t pTjetMin(0.2);
static const Double_t etaGhostMax(1.0 + rJet);
static const Double_t etaJetMax(1.0 - rJet);
static const Double_t etaBkgdMax(1.0);


void MakeMuDstJetTree(const Bool_t isInBatchMode=false) {

  gSystem -> Load("/opt/star/Xsl64_gcc482/lib/libfastjet.so");
  gSystem -> Load("/opt/star/Xsl64_gcc482/lib/libfastjettools.so");
  gSystem -> Load("StMuDstJetTreeMaker");

  // event/trigger parameters
  const Int_t    adcMax(6004);
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);
  const Double_t eEtaMin(0.5);
  const Double_t ePhiMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t etaTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  const Double_t tspPi0Min(0.);
  const Double_t tspPi0Max(0.08);
  const Double_t tspGamMin(0.2);
  const Double_t tspGamMax(0.6);
  // track parameters
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.0);
  const Double_t etaTrkMax(1.0);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(20.);
  // tower parameters
  const Double_t etaTwrMax(0.9);
  const Double_t eTwrMin(0.2);
  const Double_t eTwrMax(100.);
  const Double_t eCorrMin(0.2);
  const Double_t eCorrMax(20.);


  StMuDstJetTreeMaker *jetMaker = new StMuDstJetTreeMaker(isInBatchMode);
  // set parameters
  jetMaker -> SetInputAndOutputFiles(sInput.Data(), sOuput.Data());
  jetMaker -> SetEventParameters(rVtxMax, zVtxMax);
  jetMaker -> SetTriggerParameters(adcMax, eEtaMin, ePhiMin, pProjMax, etaTrgMax, eTtrgMin, eTtrgMax, tspPi0Min, tspPi0Max, tspGamMin, tspGamMax);
  jetMaker -> SetTrackParameters(nFitMin, rFitMin, dcaMax, etaTrkMax, pTtrkMin, pTtrkMax);
  jetMaker -> SetTowerParameters(etaTwrMax, eTwrMin, eTwrMax, eCorrMin, eCorrMax);
  jetMaker -> SetJetParameters(type, nRepeat, nRemove, rJet, aGhost, pTjetMin, etaGhostMax, etaJetMax, etaBkgdMax);
  // find jets
  jetMaker -> Init();
  jetMaker -> Make();
  jetMaker -> Finish();

}

// End ------------------------------------------------------------------------
