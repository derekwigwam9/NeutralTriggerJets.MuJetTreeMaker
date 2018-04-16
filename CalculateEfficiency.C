// 'CalculateEfficiency.C'
// Derek Anderson
// 02.13.2018
//
// Use this to calculate the efficiency for a set of embedding
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



void CalculateEfficiency() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning (embedding) efficiency calculation..." << endl;


  // io parameters
  const UInt_t  maker(1);
  const TString sOutput("pp200r9embed.efficiency.withResolution.d13m4y2018.root");
  const TString sInput[NFiles] = {"../../MuDstMatching/output/merged/pt2.matchWithMc.root", "../../MuDstMatching/output/merged/pt3.matchWithMc.root", "../../MuDstMatching/output/merged/pt4.matchWithMc.root", "../../MuDstMatching/output/merged/pt5.matchWithMc.root", "../../MuDstMatching/output/merged/pt7.matchWithMc.root", "../../MuDstMatching/output/merged/pt9.matchWithMc.root", "../../MuDstMatching/output/merged/pt11.matchWithMc.root", "../../MuDstMatching/output/merged/pt15.matchWithMc.root", "../../MuDstMatching/output/merged/pt25.matchWithMc.root", "../../MuDstMatching/output/merged/pt35.matchWithMc.root"};

  // parameters for individual files
  const TString sDate("d13m4y2018");
  const TString sSystem("pp200r9");
  const TString sOutLabel("withResolution");
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

  // grab efficiency histograms
  const TString sEff[NHist]  = {"hPhiForEff_pi", "hEtaForEff_pi", "hPtForEff_pi"};
  const TString sDir[NLevel] = {"pi0/particle/", "pi0/detector/"};
  const TString sLvl[NLevel] = {"Par", "Det"};

  TH1D    *hPhi[NFiles][NLevel];
  TH1D    *hEta[NFiles][NLevel];
  TH1D    *hPt[NFiles][NLevel];
  TString sHistName[NHist];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      // make name
      sHistName[0] = sDir[iLevel];
      sHistName[1] = sDir[iLevel];
      sHistName[2] = sDir[iLevel];
      sHistName[0].Append(sEff[0].Data());
      sHistName[1].Append(sEff[1].Data());
      sHistName[2].Append(sEff[2].Data());
      sHistName[0].Append(sLvl[iLevel].Data());
      sHistName[1].Append(sLvl[iLevel].Data());
      sHistName[2].Append(sLvl[iLevel].Data());
      // grab histogram
      hPhi[iFile][iLevel] = (TH1D*) fInput[iFile] -> Get(sHistName[0].Data());
      hEta[iFile][iLevel] = (TH1D*) fInput[iFile] -> Get(sHistName[1].Data());
      hPt[iFile][iLevel]  = (TH1D*) fInput[iFile] -> Get(sHistName[2].Data());
      if (!hPhi[iFile][iLevel] || !hEta[iFile][iLevel] || !hPt[iFile][iLevel]) {
        cerr << "PANIC: couldn't grab input histogram!\n"
             << "       code = (" << iFile << "," << iLevel << ")"
             << endl;
        err = 3;
        break;
      }
      if (err == 3) break;
    }  // end level loop
    if (err == 3) break;
  }  // end file loop
  if (err == 3)
    assert(0);
  else
    cout << "    Grabbed efficiency histograms." << endl;


  // binning
  Double_t phiBins[NPhiBins];
  Double_t etaBins[NEtaBins];
  Double_t pTbins[NPtBins];

  // get efficiency bin edges
  const Double_t nPhiHistBins = hPhi[0][0] -> GetNbinsX();
  const Double_t nEtaHistBins = hEta[0][0] -> GetNbinsX();;
  const Double_t nPtHistBins  = hPt[0][0]  -> GetNbinsX();
  const Double_t phiMaxEdge   = hPhi[0][0] -> GetBinLowEdge(NPhiBins);
  const Double_t etaMaxEdge   = hEta[0][0] -> GetBinLowEdge(NEtaBins);
  const Double_t pTmaxEdge    = hPt[0][0]  -> GetBinLowEdge(NPtBins);
  hPhi[0][0] -> GetLowEdge(phiBins);
  hEta[0][0] -> GetLowEdge(etaBins);
  hPt[0][0]  -> GetLowEdge(pTbins);
  phiBins[nPhiHistBins] = phiMaxEdge;
  etaBins[nEtaHistBins] = etaMaxEdge;
  pTbins[nPtHistBins]   = pTmaxEdge;

  // names
  const TString sSumBase[NHist]   = {"hPhiSum", "hEtaSum", "hPtSum"};
  const TString sSumLvl[NLevel]   = {"_par", "_det"};

  // sum histograms
  TH1D    *hSum[NHist][NLevel];
  TString sSumName[NHist];
  for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
    // make names
    sSumName[0] = sSumBase[0].Data();
    sSumName[1] = sSumBase[1].Data();
    sSumName[2] = sSumBase[2].Data();
    sSumName[0].Append(sSumLvl[iLevel].Data());
    sSumName[1].Append(sSumLvl[iLevel].Data());
    sSumName[2].Append(sSumLvl[iLevel].Data());
    // declare histograms
    hSum[0][iLevel] = new TH1D(sSumName[0].Data(), "", nPhiHistBins, phiBins);
    hSum[1][iLevel] = new TH1D(sSumName[1].Data(), "", nEtaHistBins, etaBins);
    hSum[2][iLevel] = new TH1D(sSumName[2].Data(), "", nPtHistBins, pTbins);
    hSum[0][iLevel] -> Sumw2();
    hSum[1][iLevel] -> Sumw2();
    hSum[2][iLevel] -> Sumw2();
    // add files
    for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
      const UInt_t iWeight = (NEmbed - NFiles) + iFile;
      hSum[0][iLevel] -> Add(hPhi[iFile][iLevel], weights[iWeight]);
      hSum[1][iLevel] -> Add(hEta[iFile][iLevel], weights[iWeight]);
      hSum[2][iLevel] -> Add(hPt[iFile][iLevel], weights[iWeight]);
    }  // end file loop
  }  // end level loop
  cout << "    Summed histograms." << endl;


  // normalize histograms
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      const UInt_t   bins = hSum[iHist][iLevel] -> GetNbinsX();
      const Double_t norm = nTrgTot[0];
      hSum[iHist][iLevel] -> Scale(1. / norm);
      for (UInt_t iBin = 1; iBin <= bins; iBin++) {
        const Double_t binW = hSum[iHist][iLevel] -> GetBinWidth(iBin);
        const Double_t binV = hSum[iHist][iLevel] -> GetBinContent(iBin);
        const Double_t binE = hSum[iHist][iLevel] -> GetBinError(iBin);
        const Double_t newV = binV / binW;
        const Double_t newE = binE / binW;
        hSum[iHist][iLevel] -> SetBinContent(iBin, newV);
        hSum[iHist][iLevel] -> SetBinError(iBin, newE);
      }
    }  // end level loop
  }  // end variable loop
  cout << "    Normalized histograms." << endl;


  // calculate efficiencies
  const TString sEffName[NHist] = {"hPhiEfficiency", "hEtaEfficiency", "hPtEfficiency"};

  TH1D    *hEff[NHist];
  Float_t weight(1.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (iHist) {
      case 0:
        hEff[iHist] = new TH1D(sEffName[iHist].Data(), "", nPhiHistBins, phiBins);
        break;
      case 1:
        hEff[iHist] = new TH1D(sEffName[iHist].Data(), "", nEtaHistBins, etaBins);
        break;
      case 2:
        hEff[iHist] = new TH1D(sEffName[iHist].Data(), "", nPtHistBins, pTbins);
        break;
    }
    hEff[iHist] -> Sumw2();
    hEff[iHist] -> Divide(hSum[iHist][1], hSum[iHist][0], weight, weight);
  }
  cout << "    Calculated efficiencies." << endl;


  // fit efficiencies
  const UInt_t   fLinF(1);
  const UInt_t   fSizF(2);
  const UInt_t   fColF(859);
  const TString  sDraw("BMQ0");
  const TString  sFunc("Eff");
  const TString  sVar[NHist]  = {"fPhi", "fEta", "fPt"};
  const TString  sForm[NHist] = {"[0]", "[0]", "[0]*(1-exp(-1.*[1]*x))"};
  const Double_t guess[NHist] = {0.5, 0.5, 0.87};
  const Double_t sigGuess(4.0);
  const Double_t fFitRange[2] = {-3.15, 3.15};
  const Double_t hFitRange[2] = {-0.7, 0.7};
  const Double_t pFitRange[2] = {1., 7.};
  const Double_t start[NHist] = {fFitRange[0], hFitRange[0], pFitRange[0]};
  const Double_t stop[NHist]  = {fFitRange[1], hFitRange[1], pFitRange[1]};

  TF1      *fFit[NHist];
  Double_t amplitude[NHist][NVal];
  Double_t sigma[NVal];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {

    // make names
    TString sFuncName(sVar[iHist].Data());
    sFuncName += sFunc.Data();

    // define fits
    fFit[iHist] = new TF1(sFuncName.Data(), sForm[iHist].Data(), start[iHist], stop[iHist]);
    fFit[iHist] -> SetLineColor(fColF);
    fFit[iHist] -> SetLineStyle(fLinF);
    fFit[iHist] -> SetLineWidth(fSizF);
    fFit[iHist] -> SetParameter(0, guess[iHist]);
    if (iHist == 2)
      fFit[iHist] -> SetParameter(1, sigGuess);

    // fit histograms
    const Int_t kNotDraw = 1<<9;
    hEff[iHist] -> Fit(fFit[iHist], sDraw.Data(), "", start[iHist], stop[iHist]);
    amplitude[iHist][0] = fFit[iHist] -> GetParameter(0);
    amplitude[iHist][1] = fFit[iHist] -> GetParError(0);
    if (iHist == 2) {
      sigma[0] = fFit[iHist] -> GetParameter(1);
      sigma[1] = fFit[iHist] -> GetParError(1);
    }
    if (hEff[iHist] -> GetFunction(sFuncName.Data()))
      hEff[iHist] -> GetFunction(sFuncName.Data()) -> ResetBit(kNotDraw);

  }  // end variable loop
  cout << "    Fit efficiencies." << endl;


  // set efficiency styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColE(1);
  const UInt_t  fMarE(1);
  const UInt_t  fColL[NLevel] = {859, 899};
  const UInt_t  fColC[NCut]   = {1, 879};
  const UInt_t  fMarL[NLevel] = {7, 4};
  const UInt_t  fMarC[NCut]   = {7, 4};
  const Float_t fLab(0.025);
  const Float_t fOffsetX(1.);
  const TString sTitleX[NHist]   = {"#varphi^{trk}", "#eta^{trk}", "p_{T}^{trk}"};
  const TString sTitleY[NHist]   = {"(1/N^{trg}_{eff}) dN^{trk}/d#varphi^{trk}", "(1/N^{trg}_{eff}) dN^{trk}/d#eta^{trk}", "(1/N^{trg}_{eff}) dN^{trk}/dp_{T}^{trk}"};
  const TString sTitleYR[NHist]  = {"#tilde{#epsilon}(#varphi^{trk}) = #epsilon_{#varphi}", "#tilde{#epsilon}(#eta^{trk}) = #epsilon_{#eta}", "#tilde{#epsilon}(p_{T}^{trk}) = #epsilon_{p} (1 - exp(-#sigma #upoint p_{T}^{trk}))"};
  const TString sTitleYE[NHist]  = {"#epsilon(#varphi^{trk}) = #epsilon_{#varphi}", "#epsilon(#eta^{trk}) = #epsilon_{#eta}", "#epsilon(p_{T}^{trk}) = #epsilon_{p} (1 - e^{-#sigma #upoint p_{T}^{trk}})"};
  const TString sSumTitle[NHist] = {"Track #varphi", "Track #eta", "Track p_{T}"};
  const TString sResTitle[NHist] = {"Cut response, #tilde{#epsilon}(#varphi^{trk})", "Cut response, #tilde{#epsilon}(#eta^{trk})", "Cut response, #tilde{#epsilon}(p_{T}^{trk})"};
  const TString sEffTitle[NHist] = {"Track efficiency, #epsilon(#varphi^{trk})", "Track efficiency, #epsilon(#eta^{trk})", "Track efficiency, #epsilon(p_{T}^{trk})"};
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hEff[iHist] -> SetLineColor(fColE);
    hEff[iHist] -> SetMarkerColor(fColE);
    hEff[iHist] -> SetMarkerStyle(fMarE);
    hEff[iHist] -> SetTitle(sEffTitle[iHist].Data());
    hEff[iHist] -> SetTitleFont(fTxt);
    hEff[iHist] -> GetXaxis() -> SetTitle(sTitleX[iHist].Data());
    hEff[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hEff[iHist] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hEff[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hEff[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hEff[iHist] -> GetYaxis() -> SetTitle(sTitleYE[iHist].Data());
    hEff[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hEff[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hEff[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      hSum[iHist][iLevel] -> SetLineColor(fColL[iLevel]);
      hSum[iHist][iLevel] -> SetMarkerColor(fColL[iLevel]);
      hSum[iHist][iLevel] -> SetMarkerStyle(fMarL[iLevel]);
      hSum[iHist][iLevel] -> SetTitle(sSumTitle[iHist].Data());
      hSum[iHist][iLevel] -> SetTitleFont(fTxt);
      hSum[iHist][iLevel] -> GetXaxis() -> SetTitle(sTitleX[iHist].Data());
      hSum[iHist][iLevel] -> GetXaxis() -> SetTitleFont(fTxt);
      hSum[iHist][iLevel] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hSum[iHist][iLevel] -> GetXaxis() -> CenterTitle(fCnt);
      hSum[iHist][iLevel] -> GetXaxis() -> SetLabelSize(fLab);
      hSum[iHist][iLevel] -> GetYaxis() -> SetTitle(sTitleY[iHist].Data());
      hSum[iHist][iLevel] -> GetYaxis() -> SetTitleFont(fTxt);
      hSum[iHist][iLevel] -> GetYaxis() -> CenterTitle(fCnt);
      hSum[iHist][iLevel] -> GetYaxis() -> SetLabelSize(fLab);
    }  // end level loop
  }  // end variable loop
  cout << "    Set styles." << endl;


  // make labels
  const UInt_t  nDec(3);
  const UInt_t  fAlign(12);
  const UInt_t  fColP(0);
  const UInt_t  fColT(879);
  const Float_t xPav[2] = {0.7, 0.9};
  const Float_t xLeg[2] = {0.7, 0.9};
  const Float_t yPav[2] = {0.7, 0.9};
  const Float_t yLeg[2] = {0.5, 0.7};
  const TString sSystem("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrgKin("Embedding (hard 2-to-2 QCD)");
  const TString sSig("#sigma = ");
  const TString sAmp[NHist]     = {"#epsilon_{#varphi} = ", "#epsilon_{#eta} = ", "#epsilon_{p} = "};
  const TString sEffLeg[NLevel] = {"particle level", "detector level"};

  TString   sAmpRaw[NVal];
  TString   sSigRaw[NVal];
  TLegend   *lTrks[NHist];
  TPaveText *pInfo[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {

    TString sAmpTxt(sAmp[iHist].Data());
    TString sSigTxt(sSig.Data());
    for (UInt_t iVal = 0; iVal < NVal; iVal++) {
      sAmpRaw[iVal]  = "";
      sSigRaw[iVal]  = "";
      sAmpRaw[iVal] += amplitude[iHist][iVal];
      sSigRaw[iVal] += sigma[iVal];

      const UInt_t nAmpRaw = sAmpRaw[iVal].First(".");
      const UInt_t nSigRaw = sSigRaw[iVal].First(".");
      const UInt_t nAmpTxt = (nAmpRaw + nDec) + 1;
      const UInt_t nSigTxt = (nSigRaw + nDec) + 1;
      if (iVal == 0) {
        sAmpTxt.Append(sAmpRaw[iVal].Data(), nAmpTxt);
        sSigTxt.Append(sSigRaw[iVal].Data(), nSigTxt);
      } 
      else {
        sAmpTxt += " #pm ";
        sSigTxt += " #pm ";
        sAmpTxt.Append(sAmpRaw[iVal].Data(), nAmpTxt);
        sSigTxt.Append(sSigRaw[iVal].Data(), nSigTxt);
      }
    }  // end value loop

    // information
    pInfo[iHist] = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
    pInfo[iHist] -> SetFillColor(fColP);
    pInfo[iHist] -> SetLineColor(fColP);
    pInfo[iHist] -> SetTextColor(fColT);
    pInfo[iHist] -> SetTextFont(fTxt);
    pInfo[iHist] -> SetTextAlign(fAlign);
    pInfo[iHist] -> AddText(sSystem.Data());
    pInfo[iHist] -> AddText(sTrgKin.Data());
    pInfo[iHist] -> AddText(sAmpTxt.Data());
    if (iHist == 2)
      pInfo[iHist] -> AddText(sSigTxt.Data());

    // legends
    lTrks[iHist] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lTrks[iHist] -> SetFillColor(fColP);
    lTrks[iHist] -> SetLineColor(fColP);
    lTrks[iHist] -> SetTextFont(fTxt);
    lTrks[iHist] -> AddEntry(hSum[iHist][0], sEffLeg[0].Data());
    lTrks[iHist] -> AddEntry(hSum[iHist][1], sEffLeg[1].Data());

  }  // end variable loop
  cout << "    Made labels." << endl;


  // make plots
  const UInt_t  width(1500);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  ticks(0);
  const UInt_t  log[NHist]      = {0, 0, 1};
  const Float_t xPad[4]         = {0., 0.5, 0.5, 1.};
  const Float_t yPad[4]         = {0., 1., 0., 1.};
  const TString sEffPads[2]     = {"pRatio", "pTracks"};
  const TString sEffPlot[NHist] = {"cPhiEfficiency", "cEtaEfficiency", "cPtEfficiency"};

  TPad    *pEff[NHist][2];
  TCanvas *cEff[NHist];
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    cEff[iHist]    = new TCanvas(sEffPlot[iHist].Data(), "", width, height);
    pEff[iHist][0] = new TPad(sEffPads[0].Data(), "", xPad[0], yPad[0], xPad[1], yPad[1]);
    pEff[iHist][1] = new TPad(sEffPads[1].Data(), "", xPad[2], yPad[2], xPad[3], yPad[3]);
    pEff[iHist][0] -> SetGrid(grid, grid);
    pEff[iHist][0] -> SetTicks(ticks, ticks);
    pEff[iHist][1] -> SetGrid(grid, grid);
    pEff[iHist][1] -> SetTicks(ticks, ticks);
    pEff[iHist][1] -> SetLogy(log[iHist]);
    cEff[iHist]    -> cd();
    pEff[iHist][0] -> Draw();
    pEff[iHist][1] -> Draw();
    pEff[iHist][0] -> cd();
    hEff[iHist]    -> Draw();
    pInfo[iHist]   -> Draw();
    pEff[iHist][1] -> cd();
    hSum[iHist][0] -> Draw();
    hSum[iHist][1] -> Draw("same");
    lTrks[iHist]   -> Draw();
    cEff[iHist]    -> Write();
    cEff[iHist]    -> Close();
  }
  cout << "    Made plots." << endl;


  // make directories, save histograms
  const TString sDirPar("particle");
  const TString sDirDet("detector");

  TDirectory *dLvl[NLevel];
  for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
    if (iLevel == 0)
      dLvl[iLevel] = (TDirectory*) fOutput -> mkdir(sDirPar.Data());
    else
      dLvl[iLevel] = (TDirectory*) fOutput -> mkdir(sDirDet.Data());
    dLvl[iLevel] -> cd();
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      hSum[iHist][iLevel] -> Write();
    }  // end variable loop
  }  // end level loop
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hEff[iHist] -> Write();
  }
  cout << "    Made directories and saved histograms." << endl;


  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    fInput[iFile] -> cd();
    fInput[iFile] -> Close();
  }
  cout << "  Calculation finished!\n" << endl;

}

// End ------------------------------------------------------------------------
