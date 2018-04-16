// 'CalculateResolution.C'
// Derek Anderson
// 04.16.2018
//
// Use this to calculate the resolution for a set of embedding
// files.  Uses (0) 'CalculateDataEfficiency.C' or (1)
// 'CalculateEmbeddingEfficiency.C' to accumulate the track
// distributions for each file.
//
// NOTE: (0) corresponds to the output
//       of 'StJetTreeThirdMaker' and
//       (1) to the output of
//       'StJetTreeMcMaker'.
//
// NOTE: this assumes the set of files
//       you're looking at is the LAST
//       NFiles out of NEmbed possible
//       files.


#include <vector>
#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TROOT.h"
#include "TFile.h"
#include "TError.h"
#include "TSystem.h"
#include "TString.h"
#include "TCanvas.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t NPhiBins(61);
static const UInt_t NEtaBins(41);
static const UInt_t NPtBins(48);
static const UInt_t NEmbed(10);
static const UInt_t NFiles(10);
static const UInt_t NLevel(2);
static const UInt_t NHist(3);
static const UInt_t NCut(2);
static const UInt_t NVal(2);
static const UInt_t NTrg(2);



void CalculateResolution() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning resolution calculation..." << endl;


  // io parameters
  const UInt_t  maker(1);
  const TString sOutput("pp200r9embed.resolution.d16m4y2018.root");
  const TString sInput[NFiles] = {"../../MuDstMatching/output/merged/pt2.matchWithMc.root", "../../MuDstMatching/output/merged/pt3.matchWithMc.root", "../../MuDstMatching/output/merged/pt4.matchWithMc.root", "../../MuDstMatching/output/merged/pt5.matchWithMc.root", "../../MuDstMatching/output/merged/pt7.matchWithMc.root", "../../MuDstMatching/output/merged/pt9.matchWithMc.root", "../../MuDstMatching/output/merged/pt11.matchWithMc.root", "../../MuDstMatching/output/merged/pt15.matchWithMc.root", "../../MuDstMatching/output/merged/pt25.matchWithMc.root", "../../MuDstMatching/output/merged/pt35.matchWithMc.root"};

  // parameters for individual files
  const TString sDate("d16m4y2018");
  const TString sSystem("pp200r9");
  const TString sOutLabel("resolution");
  const TString sTree[NLevel]     = {"McTracks", "GfmtoDst_mu"};
  const TString sBinLabel[NFiles] = {"pt2", "pt3", "pt4", "pt5", "pt7", "pt9", "pt11", "pt15", "pt25", "pt35"};

  // other parameters
  const Bool_t   batch(false);
  const Bool_t   trigger(false);
  const Bool_t   level[NLevel]   = {false, true};
  const Double_t weights[NEmbed] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};


  // load macro
  UInt_t err(0);
  switch (maker) {
    case 0:
      gROOT -> ProcessLine(".L CalculateDataEfficiency.C");
      break;
    case 1:
      gROOT -> ProcessLine(".L CalculateEmbeddingEfficiency.C");
      break;
    default:
      cerr << "PANIC: check the value of 'maker.' It should be 0 or 1..." << endl;
      err = 1;
      break;
  }
  if (err == 1)
    assert(0);
  else
    cout << "    Beginning file loop: " << NFiles << " to process." << endl;


  // file loop
  UInt_t   nTrgBin[NTrg] = {0, 0};
  Double_t nTrgTot[NTrg] = {0., 0.};
  TString  sName("");
  TString  sOutBin[NFiles];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {

    // generate output names
    sName  = sSystem.Data();
    sName += sBinLabel[iFile].Data();
    sName += ".";
    sName += sOutLabel.Data();
    sName += ".";
    sName += sDate.Data();
    sName += ".root";

    // assign name and do calculation
    sOutBin[iFile] = sName;
    switch (maker) {
      case 0:
        CalculateDataEfficiency(nTrgBin[0], nTrgBin[1], batch, level[0], trigger, sInput[iFile], sTree[0], sOutBin[iFile]);
        break;
      case 1:
        CalculateEmbeddingEfficiency(nTrgBin[0], nTrgBin[1], batch, trigger, sInput[iFile], sOutBin[iFile]);
        break;
    }

    // scale triggers
    const UInt_t iWeight = (NEmbed - NFiles) + iFile;
    nTrgTot[0] += nTrgBin[0] * weights[iWeight];
    nTrgTot[1] += nTrgBin[1] * weights[iWeight];
    cout << "      Finished file " << iFile << ":\n"
         << "        nPi0 = " << nTrgBin[0] << ", nGam = " << nTrgBin[1] << " triggers."
         << endl;

  }  // end file loop

  cout << "    File loop finished:\n"
       << "      nPi0 = " << nTrgTot[0] << ", nGam = " << nTrgTot[1] << " scaled triggers."
       << endl;


   // open output file
   TFile *fOutput = new TFile(sOutput.Data(), "recreate");
   if (!fOutput) {
     cerr << "PANIC: couldn't open output file." << endl;
     assert(0);
   }

  // open input files
  TFile *fInput[NFiles];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    fInput[iFile] = new TFile(sOutBin[iFile], "read");
    if (!fInput[iFile]) {
      cerr << "PANIC: couldn't open input file no. " << iFile << endl;
      err = 2;
      break;
    }
  }  // end file loop
  if (err == 2)
    assert(0);
  else
    cout << "    Opened files." << endl;

  // grab histograms
  const TString sRes1D[NHist] = {"pi0/hPhiRes_pi", "pi0/hEtaRes_pi", "pi0/hPtRes_pi"};
  const TString sRes2D[NHist] = {"pi0/hPhiResVsPhi_pi", "pi0/hEtaResVsEta_pi", "pi0/hPtResVsPt_pi"};
  const TString sResPR[NHist] = {"pi0/pPhiResVsPhi_pi", "pi0/pEtaResVsEta_pi", "pi0/pPtResVsPt_pi"};

  TH1D     *hPhiRes[NFiles];
  TH1D     *hEtaRes[NFiles];
  TH1D     *hPtRes[NFiles];
  TH2D     *hResVsPhi[NFiles];
  TH2D     *hResVsEta[NFiles];
  TH2D     *hResVsPt[NFiles];
  TProfile *pResVsPhi[NFiles];
  TProfile *pResVsEta[NFiles];
  TProfile *pResVsPt[NFiles];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    // grab 1d histograms
    hPhiRes[iFile] = (TH1D*) fInput[iFile] -> Get(sRes1D[0].Data());
    hEtaRes[iFile] = (TH1D*) fInput[iFile] -> Get(sRes1D[1].Data());
    hPtRes[iFile]  = (TH1D*) fInput[iFile] -> Get(sRes1D[2].Data());
    if (!hPhiRes[iFile] || !hEtaRes[iFile] || !hPtRes[iFile]) {
      cerr << "PANIC: couldn't grab 1D resolution histogram!\n"
           << "       check file " << iFile << "..."
           << endl;
      err = 4;
      break;
    }

    // grab 2d histograms
    hResVsPhi[iFile] = (TH2D*) fInput[iFile] -> Get(sRes2D[0].Data());
    hResVsEta[iFile] = (TH2D*) fInput[iFile] -> Get(sRes2D[1].Data());
    hResVsPt[iFile]  = (TH2D*) fInput[iFile] -> Get(sRes2D[2].Data());
    if (!hResVsPhi[iFile] || !hResVsEta[iFile] || !hResVsPt[iFile]) {
      cerr << "PANIC: couldn't grab 2D resolution histogram!\n"
           << "       check file " << iFile << "..."
           << endl;
      err = 4;
      break;
    }

    // grab profiles
    pResVsPhi[iFile] = (TProfile*) fInput[iFile] -> Get(sResPR[0].Data());
    pResVsEta[iFile] = (TProfile*) fInput[iFile] -> Get(sResPR[1].Data());
    pResVsPt[iFile]  = (TProfile*) fInput[iFile] -> Get(sResPR[2].Data());
    if (!pResVsPhi[iFile] || !pResVsEta[iFile] || !pResVsPt[iFile]) {
      cerr << "PANIC: couldn't grab resolution profile!\n"
           << "       check file " << iFile << "..."
           << endl;
      err = 4;
      break;
    }
  }  // end file loop
  if (err == 4)
    assert(0);
  else
    cout << "    Grabbed histograms." << endl;


  // binning
  const UInt_t   nPhiResBins = hPhiRes[0]   -> GetNbinsX();
  const UInt_t   nPhiTrkBins = hResVsPhi[0] -> GetNbinsX();
  const UInt_t   nEtaResBins = hEtaRes[0]   -> GetNbinsX();
  const UInt_t   nEtaTrkBins = hResVsEta[0] -> GetNbinsX();
  const UInt_t   nPtResBins  = hPtRes[0]    -> GetNbinsX();
  const UInt_t   nPtTrkBins  = hResVsPt[0]  -> GetNbinsX();
  const Double_t phiRes[2]   = {hPhiRes[0] -> GetBinLowEdge(0), hPhiRes[0] -> GetBinLowEdge(nPhiResBins + 1)};
  const Double_t phiTrk[2]   = {hResVsPhi[0] -> GetXaxis() -> GetBinLowEdge(1), hResVsPhi[0] -> GetXaxis() -> GetBinLowEdge(nPhiTrkBins + 1)};
  const Double_t etaRes[2]   = {hEtaRes[0] -> GetBinLowEdge(0), hEtaRes[0] -> GetBinLowEdge(nEtaResBins + 1)};
  const Double_t etaTrk[2]   = {hResVsEta[0] -> GetXaxis() -> GetBinLowEdge(1), hResVsEta[0] -> GetXaxis() -> GetBinLowEdge(nEtaTrkBins + 1)};
  const Double_t ptRes[2]    = {hPtRes[0] -> GetBinLowEdge(0), hPtRes[0] -> GetBinLowEdge(nPtResBins + 1)};
  const Double_t ptTrk[2]    = {hResVsPt[0] -> GetXaxis() -> GetBinLowEdge(1), hResVsPt[0] -> GetXaxis() -> GetBinLowEdge(nPtTrkBins + 1)};
  cout << "CHECK: pTx = (" << ptTrk[0] << ", " << ptTrk[1] << ")" << endl;

  // names
  const TString sResBase1D[NHist] = {"hPhiRes", "hEtaRes", "hPtRes"};
  const TString sResBase2D[NHist] = {"hResVsPhi", "hResVsEta", "hResVsPt"};
  const TString sResBasePR[NHist] = {"pResVsPhi", "pResVsEta", "pResVsPt"};

  // sum histograms
  TH1D     *hRes1D[NHist];
  TH2D     *hRes2D[NHist];
  TProfile *pRes2D[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    // declare histograms
    switch (iHist) {
      case 0:
        hRes1D[iHist] = new TH1D(sResBase1D[iHist].Data(), "", nPhiResBins, phiRes[0], phiRes[1]);
        hRes2D[iHist] = new TH2D(sResBase2D[iHist].Data(), "", nPhiTrkBins, phiTrk[0], phiTrk[1], nPhiResBins, phiRes[0], phiRes[1]);
        pRes2D[iHist] = new TProfile(sResBasePR[iHist].Data(), "", nPhiTrkBins, phiTrk[0], phiTrk[1], "S");
        break;
      case 1:
        hRes1D[iHist] = new TH1D(sResBase1D[iHist].Data(), "", nEtaResBins, etaRes[0], etaRes[1]);
        hRes2D[iHist] = new TH2D(sResBase2D[iHist].Data(), "", nEtaTrkBins, etaTrk[0], etaTrk[1], nEtaResBins, etaRes[0], etaRes[1]);
        pRes2D[iHist] = new TProfile(sResBasePR[iHist].Data(), "", nEtaTrkBins, etaTrk[0], etaTrk[1], "S");
        break;
      case 2:
        hRes1D[iHist] = new TH1D(sResBase1D[iHist].Data(), "", nPtResBins, ptRes[0], ptRes[1]);
        hRes2D[iHist] = new TH2D(sResBase2D[iHist].Data(), "", nPtTrkBins, ptTrk[0], ptTrk[1], nPtResBins, ptRes[0], ptRes[1]);
        pRes2D[iHist] = new TProfile(sResBasePR[iHist].Data(), "", nPtTrkBins, ptTrk[0], ptTrk[1], "S");
        break;
    }
    for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
      const UInt_t iWeight = (NEmbed - NFiles) + iFile;
      switch (iHist) {
        case 0:
          hRes1D[iHist] -> Add(hPhiRes[iFile], weights[iWeight]);
          hRes2D[iHist] -> Add(hResVsPhi[iFile], weights[iWeight]);
          pRes2D[iHist] -> Add(pResVsPhi[iFile], weights[iWeight]);
          break;
        case 1:
          hRes1D[iHist] -> Add(hEtaRes[iFile], weights[iWeight]);
          hRes2D[iHist] -> Add(hResVsEta[iFile], weights[iWeight]);
          pRes2D[iHist] -> Add(pResVsEta[iFile], weights[iWeight]);
          break;
        case 2:
          hRes1D[iHist] -> Add(hPtRes[iFile], weights[iWeight]);
          hRes2D[iHist] -> Add(hResVsPt[iFile], weights[iWeight]);
          pRes2D[iHist] -> Add(pResVsPt[iFile], weights[iWeight]);
          break;
      }
    }  // end file loop
  }  // end variable loop
  cout << "    Summed histograms." << endl;


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColR(1);
  const UInt_t  fMarR(8);
  const Float_t fLab(0.025);
  const Float_t fOffsetX(1.);
  const TString sTitleR1D[NHist] = {"#varphi^{trk} resolution, #sigma(#varphi^{trk})", "#eta^{trk} resolution, #sigma(#eta^{trk})", "p_{T}^{trk} resolution, #sigma(p_{T}^{trk})"};
  const TString sTitleR2D[NHist] = {"#sigma(#varphi^{trk}) vs. #varphi^{MC}", "#sigma(#eta^{trk}) vs. #eta^{MC}", "#sigma(p_{T}^{trk}) vs. p_{T}^{MC}"};
  const TString sTitleR1X[NHist] = {"#sigma(#varphi^{trk}) = (#varphi^{MC} - #varphi^{trk}) / #varphi^{MC}", "#sigma(#eta^{trk}) = (#eta^{MC} - #eta^{trk}) / #eta^{MC}", "#sigma(p_{T}^{trk}) = (p_{T}^{MC} - p_{T}^{trk}) / p_{T}^{MC}"};
  const TString sTitleR1Y[NHist] = {"dN^{trk}/d#sigma(#varphi^{trk})", "dN^{trk}/d#sigma(#eta^{trk})", "dN^{trk}/d#sigma(p_{T}^{trk})"};
  const TString sTitleR2X[NHist] = {"#varphi^{MC}", "#eta^{MC}", "p_{T}^{MC}"};
  const TString sTitleR2Y[NHist] = {"#sigma(#varphi^{trk})", "#sigma(#eta^{trk})", "#sigma(p_{T}^{trk})"};
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hRes1D[iHist] -> SetLineColor(fColR);
    hRes1D[iHist] -> SetMarkerColor(fColR);
    hRes1D[iHist] -> SetMarkerStyle(fMarR);
    hRes1D[iHist] -> SetTitle(sTitleR1D[iHist]);
    hRes1D[iHist] -> SetTitleFont(fTxt);
    hRes1D[iHist] -> GetXaxis() -> SetTitle(sTitleR1X[iHist].Data());
    hRes1D[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hRes1D[iHist] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hRes1D[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hRes1D[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hRes1D[iHist] -> GetYaxis() -> SetTitle(sTitleR1Y[iHist].Data());
    hRes1D[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hRes1D[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hRes1D[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hRes2D[iHist] -> SetTitle(sTitleR2D[iHist]);
    hRes2D[iHist] -> SetTitleFont(fTxt);
    hRes2D[iHist] -> GetXaxis() -> SetTitle(sTitleR2X[iHist].Data());
    hRes2D[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hRes2D[iHist] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hRes2D[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hRes2D[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hRes2D[iHist] -> GetYaxis() -> SetTitle(sTitleR2Y[iHist].Data());
    hRes2D[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hRes2D[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hRes2D[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    pRes2D[iHist] -> SetTitle(sTitleR2D[iHist]);
    pRes2D[iHist] -> SetTitleFont(fTxt);
    pRes2D[iHist] -> GetXaxis() -> SetTitle(sTitleR2X[iHist].Data());
    pRes2D[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    pRes2D[iHist] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    pRes2D[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    pRes2D[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    pRes2D[iHist] -> GetYaxis() -> SetTitle(sTitleR2Y[iHist].Data());
    pRes2D[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    pRes2D[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    pRes2D[iHist] -> GetYaxis() -> SetLabelSize(fLab);
  }
  cout << "    Set styles." << endl;


  // make plots
  const UInt_t  width(1500);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  ticks(0);
  const UInt_t  logRes(1);
  const Float_t xPad[4]         = {0., 0.5, 0.5, 1.};
  const Float_t yPad[4]         = {0., 1., 0., 1.};
  const TString sResPads[2]     = {"pRes1D", "pRes2D"};
  const TString sResPlot[NHist] = {"cPhiResolution", "cEtaResolution", "cPtResolution"};

  TPad    *pRes[NHist][2];
  TCanvas *cRes[NHist];
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    cRes[iHist]    = new TCanvas(sResPlot[iHist].Data(), "", width, height);
    pRes[iHist][0] = new TPad(sResPads[0].Data(), "", xPad[0], yPad[0], xPad[1], yPad[1]);
    pRes[iHist][1] = new TPad(sResPads[1].Data(), "", xPad[2], yPad[2], xPad[3], yPad[3]);
    pRes[iHist][0] -> SetGrid(grid, grid);
    pRes[iHist][0] -> SetLogy(logRes);
    pRes[iHist][1] -> SetGrid(grid, grid);
    pRes[iHist][1] -> SetLogz(logRes);
    cRes[iHist]    -> cd();
    pRes[iHist][0] -> Draw();
    pRes[iHist][1] -> Draw();
    pRes[iHist][0] -> cd();
    hRes1D[iHist]  -> Draw();
    pRes[iHist][1] -> cd();
    hRes2D[iHist]  -> Draw("colz");
    pRes2D[iHist]  -> Draw("same");
    cRes[iHist]    -> Write();
    cRes[iHist]    -> Close();
  }
  cout << "    Made plots." << endl;


  // save histograms
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hRes1D[iHist] -> Write();
    hRes2D[iHist] -> Write();
    pRes2D[iHist] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      fInput[iFile] -> cd();
      fInput[iFile] -> Close();
    }
  }
  cout << "  Calculation finished!\n" << endl;

}

// End ------------------------------------------------------------------------
