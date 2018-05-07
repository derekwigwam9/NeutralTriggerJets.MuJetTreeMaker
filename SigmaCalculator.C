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
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

using namespace std;


// global constants
static const UInt_t  NBins(81);
static const UInt_t  NDraw(6);
static const UInt_t  NFiles(10);
static const UInt_t  NEmbed(10);
static const TString SFuncComp("-0.026+0.020*x+0.003*x*x");



void SigmaCalculator() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sigma calculation..." << endl;


  // io parameters
  const TString sOutput("pp200r9rff.sigmaPtDiff.withPtPlots.d7m5y2018.root");
  const TString sInput[NFiles] = {"output/ResolutionStudy/pp200r9pt2.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt3.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt4.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt5.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt7.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt9.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt11.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt15.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt25.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root", "output/ResolutionStudy/pp200r9pt35.resQaCutsButNoPtRecoOrQaTruthWithRFF.d3m5y2018.root"};

  // constants
  const TString  sHist("hPtTrk");
  const TString  sTuple("nTrkMatch");
  const Float_t  loEnds[NBins]   = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6, 12.8, 13.0, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.4, 14.6, 14.8, 15.0, 15.2, 15.4, 15.6, 15.8, 16.0};
  const Float_t  hiEnds[NBins]   = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6, 11.8, 12.0, 12.2, 12.4, 12.6, 12.8, 13.0, 13.2, 13.4, 13.6, 13.8, 14.0, 14.2, 14.4, 14.6, 14.8, 15.0, 15.2, 15.4, 15.6, 15.8, 16.0, 16.2};
  const Double_t weights[NEmbed] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};


  // open output
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!" << endl;
    return;
  }
  cout << "    Opened output." << endl;

  // create histograms
  //const UInt_t   nPtBins(300);
  const UInt_t   nPtBins(500);
  //const Double_t pTbins[2] = {-5., 25.};
  const Double_t pTbins[2] = {-25., 25.};

  TH1D *hPtTrk[NBins];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    TString sName(sHist.Data());
    sName += (Int_t) (loEnds[iBin] * 10.);
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
      //if (isInBin) hPtTrk[bin] -> Fill(trkPt, weights[iWeight]);
      if (isInBin) hPtTrk[bin] -> Fill(mcPt - trkPt, weights[iWeight]);

    }  // end track loop
  }  // end file loop

  if (error > 0) {
    cerr << "PANIC: error! (code = " << error << ")" << endl;
    return;
  }
  cout << "    File loop finished!" << endl;


  // fit histograms
  Double_t mcVal[NBins];
  Double_t mcErr[NBins];
  Double_t sigVal[NBins];
  Double_t sigErr[NBins];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    // normalize histogram
    const Double_t nEntries = (Double_t) hPtTrk[iBin] -> GetEntries();
    if (nEntries > 0.) hPtTrk[iBin] -> Scale(1. / nEntries);

    // do fitting
    hPtTrk[iBin] -> Fit("gaus", "Q0");
    if (hPtTrk[iBin] -> GetFunction("gaus")) {
      mcVal[iBin]  = loEnds[iBin] + ((hiEnds[iBin] - loEnds[iBin]) / 2.);
      mcErr[iBin]  = 0.;
      sigVal[iBin] = hPtTrk[iBin] -> GetFunction("gaus") -> GetParameter(2);
      sigErr[iBin] = hPtTrk[iBin] -> GetFunction("gaus") -> GetParError(2);
      hPtTrk[iBin] -> GetFunction("gaus") -> ResetBit((Int_t) 1<<9);
    }
    else {
      mcVal[iBin]  = loEnds[iBin] + ((hiEnds[iBin] - loEnds[iBin]) / 2.);
      mcErr[iBin]  = 0.;
      sigVal[iBin] = 0.;
      sigErr[iBin] = 0.;
    }
  }
  fOutput -> cd();
  cout << "    Fit (and normalized) histograms." << endl;

  // create graph
  const UInt_t  cSig(898);
  const UInt_t  cFit(858);
  const UInt_t  cComp(818);
  const UInt_t  lSig(2);
  const UInt_t  lFit(1);
  const UInt_t  mSig(8);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.02);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  const TString sNamSig("gSigma");
  //const TString sTtlSig("Resolution #sigma(p_{T}^{reco} | p_{T}^{MC})");
  const TString sTtlSig("Resolution #sigma(#Deltap_{T} | p_{T}^{MC}), #Deltap_{T} = p_{T}^{MC} - p_{T}^{reco}");
  const TString sTtlSigX("#LTp_{T}^{MC}#GT");
  const TString sTtlSigY("#sigma");

  TGraphErrors *gSigma = new TGraphErrors(NBins, mcVal, sigVal, mcErr, sigErr);
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


  // select histograms to draw
  const UInt_t  cDraw[NDraw]   = {794, 814, 834, 854, 874, 894};
  const UInt_t  mDraw[NDraw]   = {4, 4, 4, 4, 4, 4};
  const Float_t pTdraw[NDraw]  = {0.9, 3.9, 6.9, 9.9, 12.9, 15.9};
  const TString sPtDraw[NDraw] = {"0.9 GeV/c", "3.9 GeV/c", "6.9 GeV/c", "9.9 GeV/c", "12.9 GeV/c", "15.9 GeV/c"};

  UInt_t iSelect(0);
  UInt_t iDraw[NDraw];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    const Bool_t isInBin = ((pTdraw[iSelect] >= loEnds[iBin]) && (pTdraw[iSelect] <= hiEnds[iBin]));
    if (isInBin) {
      iDraw[iSelect] = iBin;
      iSelect++;
    }
    if (iSelect == NDraw) break;
  }
  cout << "    Selected histograms." << endl;

  // set styles
  const UInt_t  cLegPt(0);
  const UInt_t  lWeight(2);
  //const TString sTtlPt("Matched p_{T}^{reco}");
  const TString sTtlPt("Matched #Deltap_{T}");
  //const TString sTtlPtX("p_{T}^{reco}");
  const TString sTtlPtX("#Deltap_{T} = p_{T}^{MC} - p_{T}^{reco}");
  const TString sTtlPtY("arb.");
  const Float_t xyLegPt[4] = {0.3, 0.1, 0.5, 0.3};

  TLegend *lPtTrk = new TLegend(xyLegPt[0], xyLegPt[1], xyLegPt[2], xyLegPt[3], "p_{T}^{MC} #pm 0.1 GeV/c");
  lPtTrk -> SetFillColor(cLegPt);
  lPtTrk -> SetLineColor(cLegPt);
  lPtTrk -> SetTextFont(fTxt);
  for (UInt_t iHist = 0; iHist < NDraw; iHist++) {
    // set histogram styles
    hPtTrk[iDraw[iHist]] -> SetLineColor(cDraw[iHist]);
    hPtTrk[iDraw[iHist]] -> SetMarkerColor(cDraw[iHist]);
    hPtTrk[iDraw[iHist]] -> SetMarkerStyle(mDraw[iHist]);
    hPtTrk[iDraw[iHist]] -> SetTitle(sTtlPt.Data());
    hPtTrk[iDraw[iHist]] -> SetTitleFont(fTxt);
    hPtTrk[iDraw[iHist]] -> GetXaxis() -> SetTitle(sTtlPtX.Data());
    hPtTrk[iDraw[iHist]] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtTrk[iDraw[iHist]] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPtTrk[iDraw[iHist]] -> GetXaxis() -> CenterTitle(fCnt);
    hPtTrk[iDraw[iHist]] -> GetXaxis() -> SetLabelSize(fLab);
    hPtTrk[iDraw[iHist]] -> GetYaxis() -> SetTitle(sTtlPtY.Data());
    hPtTrk[iDraw[iHist]] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtTrk[iDraw[iHist]] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPtTrk[iDraw[iHist]] -> GetYaxis() -> CenterTitle(fCnt);
    hPtTrk[iDraw[iHist]] -> GetYaxis() -> SetLabelSize(fLab);

    // set function color
    if (hPtTrk[iDraw[iHist]] -> GetFunction("gaus")) {
      hPtTrk[iDraw[iHist]] -> GetFunction("gaus") -> SetLineColor(cDraw[iHist]);
      hPtTrk[iDraw[iHist]] -> GetFunction("gaus") -> SetLineWidth(lWeight);
    }

    // add to legend
    lPtTrk -> AddEntry(hPtTrk[iDraw[iHist]], sPtDraw[iHist].Data());
  }
  cout << "    Set styles." << endl;


  // fit graph
  const UInt_t  cLegSig(0);
  const TString sLegFit("Fit");
  const TString sLegComp("Comparison");
  const Float_t xComp[2]    = {0., 20.};
  const Float_t fitRange[2] = {mcVal[0], mcVal[NBins - 2]};
  const Float_t xyLegSig[4] = {0.3, 0.1, 0.5, 0.3};
  gSigma -> Fit("pol2", "0", "", fitRange[0], fitRange[1]);
  gSigma -> GetFunction("pol2") -> SetLineColor(cFit);
  gSigma -> GetFunction("pol2") -> SetLineStyle(lFit);
  gSigma -> GetFunction("pol2") -> SetLineWidth(lWeight);
  gSigma -> GetFunction("pol2") -> ResetBit((Int_t) 1<<9);

  // make comparison
  TF1 *fComp = new TF1("fComp", SFuncComp.Data(), xComp[0], xComp[1]);
  fComp -> SetLineColor(cComp);
  fComp -> SetLineStyle(lFit);
  fComp -> SetLineWidth(lWeight);

  // make legend
  TLegend *lSigma = new TLegend(xyLegSig[0], xyLegSig[1], xyLegSig[2], xyLegSig[3]);
  lSigma -> SetFillColor(cLegSig);
  lSigma -> SetLineColor(cLegSig);
  lSigma -> SetTextFont(fTxt);
  lSigma -> AddEntry(gSigma -> GetFunction("pol2"), sLegFit.Data());
  lSigma -> AddEntry(fComp, sLegComp.Data());


  // make plots
  const UInt_t  wSig(750);
  const UInt_t  wPt(750);
  const UInt_t  hSig(750);
  const UInt_t  hPt(750);
  const UInt_t  grid(0);
  const UInt_t  log(1);
  const TString sCanSig("cSigma");
  const TString sCanPt("cPtTrk");

  TCanvas *cSigma = new TCanvas(sCanSig.Data(), "", wSig, hSig);
  cSigma -> SetGrid(grid, grid);
  cSigma -> SetLogy(log);
  gSigma -> Draw("ALP");
  fComp  -> Draw("same");
  lSigma -> Draw();
  cSigma -> Write();
  cSigma -> Close();

  TCanvas *cPtTrk = new TCanvas(sCanPt.Data(), "", wPt, hPt);
  cPtTrk -> SetGrid(grid, grid);
  cPtTrk -> SetLogy(log);
  for (UInt_t iHist = 0; iHist < NDraw; iHist++) {
    //if (iHist == 0)
      //hPtTrk[iDraw[iHist]] -> Draw();
    //else
      //hPtTrk[iDraw[iHist]] -> Draw("same");
    if (iHist == 0)
      hPtTrk[iDraw[iHist]] -> Draw("hist");
    else
      hPtTrk[iDraw[iHist]] -> Draw("same hist");
  }
  lPtTrk -> Draw();
  cPtTrk -> Write();
  cPtTrk -> Close();
  cout << "    Made plots." << endl;


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
