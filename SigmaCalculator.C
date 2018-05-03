// 'SigmaCalculator.C'
// Derek Anderson
// 05.02.2018
//
// Use this calculate (and plot) the
// resolution width.


#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TString.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

using namespace std;


// global constatns
static const UInt_t NBins(9);
static const UInt_t NFiles(10);
static const UInt_t NEmbed(10);



void SigmaCalculator() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sigma calculation..." << endl;


  // io parameters
  const TString sOutput("test.root");
  const TString sInput[NFiles] = {"pp200r9pt2.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt3.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt4.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt5.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt7.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt9.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt11.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt15.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt25.resNoPQaCuts.noQaTruth50.d30m4y2018.root", "pp200r9pt35.resNoPQaCuts.noQaTruth50.d30m4y2018.root"};

  // constants
  const TString  sHist("hPtTrk");
  const TString  sTuple("nTrkMatch");
  const Float_t  loEnds[NBins]   = {0, 2, 4, 6, 8, 10, 12, 14, 16};
  const Float_t  hiEnds[NBins]   = {2, 4, 6, 8, 10, 12, 14, 16, 20};
  const Double_t weights[NEmbed] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};


  // open output
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!" << endl;
    return;
  }
  cout << "    Opened output." << endl;

  // create histograms
  const UInt_t   nPtBins(300);
  const Double_t pTbins[2] = {-5., 25.};

  TH1D *hPtTrk[NBins];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    TString sName(sHist.Data());
    sName += loEnds[iBin];
    hPtTrk[iBin] = new TH1D(sName.Data(), "", nPtBins, pTbins[0], pTbins[1]);
    hPtTrk[iBin] -> Sumw2();
  }
  cout << "    Declared histograms." << endl;


  // leaf addresses
  Float_t mcIdTruth;
  Float_t mcIdVtx;
  Float_t mcVtxEnd;
  Float_t mcIdGeant;
  Float_t mcEta;
  Float_t mcChrg;
  Float_t mcPt;
  Float_t trkIdTruth;
  Float_t trkQaTruth;
  Float_t trkNumFit;
  Float_t trkNumPoss;
  Float_t trkDca;
  Float_t trkEta;
  Float_t trkPt;

  // leaf pointers
  TLeaf *lMcIdTruth;
  TLeaf *lMcIdVtx;
  TLeaf *lMcVtxEnd;
  TLeaf *lMcIdGeant;
  TLeaf *lMcEta;
  TLeaf *lMcChrg;
  TLeaf *lMcPt;
  TLeaf *lTrkIdTruth;
  TLeaf *lTrkQaTruth;
  TLeaf *lTrkNumFit;
  TLeaf *lTrkNumPoss;
  TLeaf *lTrkDca;
  TLeaf *lTrkEta;
  TLeaf *lTrkPt;


  // file loop
  Int_t error(0);
  cout << "    Beginning file loop: " << NFiles << " files." << endl;

  TFile   *fInput[NFiles];
  TNtuple *tMatch[NFiles];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {

    // open input
    fInput[iFile] = new TFile(sInput[iFile].Data(), "read");
    if (!fInput[iFile]) {
      cerr << "PANIC: couldn't open file #" << iFile << "!" << endl;
      error = 1;
      break;
    }

    // grab tuple
    tMatch[iFile] = (TNtuple*) fInput[iFile] -> Get(sTuple.Data());
    if (!tMatch) {
      cerr << "PANIC: couldn't grab tuple #" << iFile << "!" << endl;
      error = 1;
      break;
    }

    // get leaves
    lMcIdTruth  = tMatch[iFile] -> GetLeaf("mcIdTruth");
    lMcIdVtx    = tMatch[iFile] -> GetLeaf("mcIdVtx");
    lMcVtxEnd   = tMatch[iFile] -> GetLeaf("mcVtxEnd");
    lMcIdGeant  = tMatch[iFile] -> GetLeaf("mcIdGeant");
    lMcEta      = tMatch[iFile] -> GetLeaf("mcEta");
    lMcChrg     = tMatch[iFile] -> GetLeaf("mcChrg");
    lMcPt       = tMatch[iFile] -> GetLeaf("mcPt");
    lTrkIdTruth = tMatch[iFile] -> GetLeaf("trkIdTruth");
    lTrkQaTruth = tMatch[iFile] -> GetLeaf("trkQaTruth");
    lTrkNumFit  = tMatch[iFile] -> GetLeaf("trkNumFit");
    lTrkNumPoss = tMatch[iFile] -> GetLeaf("trkNumPoss");
    lTrkDca     = tMatch[iFile] -> GetLeaf("trkDca");
    lTrkEta     = tMatch[iFile] -> GetLeaf("trkEta");
    lTrkPt      = tMatch[iFile] -> GetLeaf("trkPt");

    // set leaf addresses
    lMcIdTruth  -> SetAddress(&mcIdTruth);
    lMcIdVtx    -> SetAddress(&mcIdVtx);
    lMcVtxEnd   -> SetAddress(&mcVtxEnd);
    lMcIdGeant  -> SetAddress(&mcIdGeant);
    lMcEta      -> SetAddress(&mcEta);
    lMcChrg     -> SetAddress(&mcChrg);
    lMcPt       -> SetAddress(&mcPt);
    lTrkIdTruth -> SetAddress(&trkIdTruth);
    lTrkQaTruth -> SetAddress(&trkQaTruth);
    lTrkNumFit  -> SetAddress(&trkNumFit);
    lTrkNumPoss -> SetAddress(&trkNumPoss);
    lTrkDca     -> SetAddress(&trkDca);
    lTrkEta     -> SetAddress(&trkEta);
    lTrkPt      -> SetAddress(&trkPt);


    // track loop
    const UInt_t nTrks   = tMatch[iFile] -> GetEntries();
    const UInt_t iWeight = (NEmbed - NFiles) + iFile;
    cout << "      File " << iFile << ": " << nTrks << " tracks." << endl;

    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      const Int_t bytes = tMatch[iFile] -> GetEntry(iTrk);
      if (bytes < 0) {
        cerr << "PANIC: something weird in file #" << iFile << " at entry #" << iTrk << "!" << endl;
        error = 2;
        break;
      }

      // fill histogram
      UInt_t bin(0);
      Bool_t isInBin(false);
      for (UInt_t iBin = 0; iBin < NBins; iBin++) {
        const Bool_t isInMcPtCut = ((mcPt > loEnds[iBin]) && (mcPt <= hiEnds[iBin]));
        if (isInMcPtCut) {
          bin     = iBin;
          isInBin = true;
          break;
        }
      }  // end bin loop
      if (isInBin) hPtTrk[bin] -> Fill(trkPt, weights[iWeight]);

    }  // end track loop
  }  // end file loop

  if (error > 0) {
    cerr << "PANIC: error! (code = " << error << ")" << endl;
    return;
  }
  cout << "    File loop finished!" << endl;


  // fit histograms
  Double_t muVal[NBins];
  Double_t muErr[NBins];
  Double_t sigVal[NBins];
  Double_t sigErr[NBins];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    hPtTrk[iBin] -> Fit("gaus", "Q0");
    muVal[iBin]  = hPtTrk[iBin] -> GetFunction("gaus") -> GetParameter(1);
    muErr[iBin]  = hPtTrk[iBin] -> GetFunction("gaus") -> GetParError(1);
    sigVal[iBin] = hPtTrk[iBin] -> GetFunction("gaus") -> GetParameter(2);
    sigErr[iBin] = hPtTrk[iBin] -> GetFunction("gaus") -> GetParError(2);
    hPtTrk[iBin] -> GetFunction("gaus") -> ResetBit((Int_t) 1<<9);
  }
  fOutput -> cd();
  cout << "    Fit histograms." << endl;

  // create graph
  const UInt_t  cSig(898);
  const UInt_t  cFit(858);
  const UInt_t  lSig(2);
  const UInt_t  lFit(1);
  const UInt_t  mSig(8);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.02);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  const TString sNamSig("gSigma");
  const TString sTtlSig("Resolution #sigma(#LTp_{T}^{reco}#GT)");
  const TString sTtlSigX("#LTp_{T}^{reco}#GT");
  const TString sTtlSigY("#sigma");

  TGraphErrors *gSigma = new TGraphErrors(NBins, muVal, sigVal, muErr, sigErr);
  gSigma -> SetLineColor(cSig);
  gSigma -> SetLineStyle(lSig);
  gSigma -> SetMarkerColor(cSig);
  gSigma -> SetMarkerStyle(mSig);
  gSigma -> SetName(sNamSig.Data());
  gSigma -> SetTitle(sTtlSig.Data());
  gSigma -> GetXaxis() -> SetTitle(sTtlSigX.Data());
  gSigma -> GetXaxis() -> SetTitleFont(fTxt);
  gSigma -> GetXaxis() -> SetTitleOffset(fOffX);
  gSigma -> GetXaxis() -> CenterTitle(fCnt);
  gSigma -> GetXaxis() -> SetLabelSize(fLab);
  gSigma -> GetYaxis() -> SetTitle(sTtlSigY.Data());
  gSigma -> GetYaxis() -> SetTitleFont(fTxt);
  gSigma -> GetYaxis() -> SetTitleOffset(fOffY);
  gSigma -> GetYaxis() -> CenterTitle(fCnt);
  gSigma -> GetYaxis() -> SetLabelSize(fLab);
  cout << "    Made graph." << endl;

  // fit graph
  const Float_t fitRange[2] = {muVal[0], muVal[NBins - 2]};
  gSigma -> Fit("pol2", "0", "", fitRange[0], fitRange[1]);
  gSigma -> GetFunction("pol2") -> SetLineColor(cFit);
  gSigma -> GetFunction("pol2") -> SetLineStyle(lFit);
  gSigma -> GetFunction("pol2") -> ResetBit((Int_t) 1<<9);


  // make plot
  const UInt_t  wSig(750);
  const UInt_t  hSig(750);
  const UInt_t  grid(0);
  const TString sCanSig("cSigma");

  TCanvas *cSigma = new TCanvas(sCanSig.Data(), "", wSig, hSig);
  cSigma -> SetGrid(grid, grid);
  gSigma -> Draw("ALP");
  cSigma -> Write();
  cSigma -> Close();
  cout << "    Made plot." << endl;


  // save histograms
  fOutput -> cd();
  gSigma  -> Write();
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    hPtTrk[iBin] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  cout << "  Calculation finished!\n" << endl;

}
