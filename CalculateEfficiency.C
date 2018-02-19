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


#include <vector>
#include <cassert>
#include "TH1.h"
#include "TPad.h"
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TCanvas.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t NEmbed(10);
static const UInt_t NFiles(7);
static const UInt_t NLevel(2);
static const UInt_t NHist(3);
static const UInt_t NFunc(2);
static const UInt_t NCut(2);
static const UInt_t NVal(2);



void CalculateEfficiency() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning (embedding) efficiency calculation..." << endl;


  // io parameters
  const UInt_t  maker(1);
  const TString sOutput("pp200r9embed.efficiencyWeightCheck.d172y2018.root");
  const TString sInput[NFiles] = {"../../MuDstMatching/output/merged/pt5.match.root", "../../MuDstMatching/output/merged/pt7.match.root", "../../MuDstMatching/output/merged/pt9.match.root", "../../MuDstMatching/output/merged/pt11.match.root", "../../MuDstMatching/output/merged/pt15.match.root", "../../MuDstMatching/output/merged/pt25.match.root", "../../MuDstMatching/output/merged/pt35.match.root"};

  // parameters for individual files
  const TString sDate("d17m2y2018");
  const TString sSystem("pp200r9");
  const TString sOutLabel("cutResponseWeightCheck");
  const TString sTree[NLevel]     = {"GfmtoDst_gnt", "GfmtoDst_mu"};
  const TString sLvlLabel[NLevel] = {"g", "u"};
  const TString sBinLabel[NFiles] = {"pt5", "pt7", "pt9", "pt11", "pt15", "pt25", "pt35"};

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

  // declare list of bad events
  vector<UInt_t> badTreeIndices;
  vector<UInt_t> parTreeIndices;
  badTreeIndices.clear();
  parTreeIndices.clear();


  // file loop
  UInt_t   nPi0Bin[NLevel] = {0, 0};
  UInt_t   nGamBin[NLevel] = {0, 0};
  Double_t nPi0Tot[NLevel] = {0., 0.};
  Double_t nGamTot[NLevel] = {0., 0.};
  TString  sName[NLevel];
  TString  sOutBin[NFiles][NLevel];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {

    badTreeIndices.clear();
    parTreeIndices.clear();
    cout << "      Processing file " << iFile << "..." << endl;

    // generate output names
    sName[0]  = sSystem.Data();
    sName[0] += sBinLabel[iFile].Data();
    sName[0] += sLvlLabel[0].Data();
    sName[0] += ".";
    sName[0] += sOutLabel.Data();
    sName[0] += ".";
    sName[0] += sDate.Data();
    sName[0] += ".root";
    sName[1]  = sSystem.Data();
    sName[1] += sBinLabel[iFile].Data();
    sName[1] += sLvlLabel[1].Data();
    sName[1] += ".";
    sName[1] += sOutLabel.Data();
    sName[1] += ".";
    sName[1] += sDate.Data();
    sName[1] += ".root";

    // assign name and do calculation
    sOutBin[iFile][0] = sName[0];
    sOutBin[iFile][1] = sName[1];
    switch (maker) {
      case 0:
        CalculateDataEfficiency(nPi0Bin[1], nGamBin[1], batch, level[1], trigger, sInput[iFile], sTree[1], sOutBin[iFile][1]);
        CalculateDataEfficiency(nPi0Bin[0], nGamBin[0], batch, level[0], trigger, sInput[iFile], sTree[0], sOutBin[iFile][0]);
        break;
      case 1:
        badTreeIndices = CalculateEmbeddingEfficiency(nPi0Bin[1], nGamBin[1], batch, level[1], trigger, sInput[iFile], sTree[1], sOutBin[iFile][1]);
        parTreeIndices = CalculateEmbeddingEfficiency(nPi0Bin[0], nGamBin[0], batch, level[0], trigger, sInput[iFile], sTree[0], sOutBin[iFile][0], badTreeIndices);
        break;
    }

    // scale triggers
    const UInt_t iWeight = (NEmbed - NFiles) + iFile;
    //nPi0Tot[0] += nPi0Bin[0] * weights[iWeight];
    nPi0Tot[0] += nPi0Bin[0] * 1.;
    nPi0Tot[1] += nPi0Bin[1] * weights[iWeight];
    //nGamTot[0] += nGamBin[0] * weights[iWeight];
    nGamTot[0] += nGamBin[0] * 1.;
    nGamTot[1] += nGamBin[1] * weights[iWeight];
    cout << "      Finished file " << iFile << ":\n"
         << "        [particle] nPi0 = " << nPi0Bin[0] << ", nGam = " << nGamBin[0] << " triggers\n"
         << "        [detector] nPi0 = " << nPi0Bin[1] << ", nGam = " << nGamBin[1] << " triggers."
         << endl;

  }  // end file loop

  cout << "    File loop finished:\n"
       << "      [particle] nPi0 = " << nPi0Tot[0] << ", nGam = " << nGamTot[0] << " scaled triggers\n"
       << "      [detector] nPi0 = " << nPi0Tot[1] << ", nGam = " << nGamTot[1] << " scaled triggers."
       << endl;


   // open output file
   TFile *fOutput = new TFile(sOutput.Data(), "recreate");
   if (!fOutput) {
     cerr << "PANIC: couldn't open output file." << endl;
     assert(0);
   }

  // open input files
  TFile *fInput[NFiles][NLevel];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      fInput[iFile][iLevel] = new TFile(sOutBin[iFile][iLevel], "read");
      if (!fInput[iFile][iLevel]) {
        cerr << "PANIC: couldn't open input file no. " << iFile << ", level " << iLevel << endl;
        err = 2;
        break;
      }
    }  // end level loop
  }  // end file loop
  if (err == 2)
    assert(0);
  else
    cout << "    Opened files." << endl;

  // grab histograms
  const TString sBase[NHist] = {"pi0/hPhiRaw", "pi0/hEtaRaw", "pi0/hPtRaw"};
  const TString sCut[NCut]   = {"BeforeQA_pi", "AfterQA_pi"};

  TH1D    *hPhi[NFiles][NLevel][NCut];
  TH1D    *hEta[NFiles][NLevel][NCut];
  TH1D    *hPt[NFiles][NLevel][NCut];
  TString sHistName[NHist];
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        // make name
        sHistName[0] = sBase[0];
        sHistName[1] = sBase[1];
        sHistName[2] = sBase[2];
        sHistName[0].Append(sCut[iCut].Data());
        sHistName[1].Append(sCut[iCut].Data());
        sHistName[2].Append(sCut[iCut].Data());
        // grab histogram
        hPhi[iFile][iLevel][iCut] = (TH1D*) fInput[iFile][iLevel] -> Get(sHistName[0].Data());
        hEta[iFile][iLevel][iCut] = (TH1D*) fInput[iFile][iLevel] -> Get(sHistName[1].Data());
        hPt[iFile][iLevel][iCut]  = (TH1D*) fInput[iFile][iLevel] -> Get(sHistName[2].Data());
        if (!hPhi[iFile][iLevel][iCut] || !hEta[iFile][iLevel][iCut] || !hPt[iFile][iLevel][iCut]) {
          cerr << "PANIC: couldn't grab input phi histogram!\n"
               << "       code = (" << iFile << "," << iLevel << "," << iCut << ")"
               << endl;
          err = 3;
          break;
        }
      }  // end cut loop
      if (err == 3) break;
    }  // end level loop
    if (err == 3) break;
  }  // end file loop
  if (err == 3)
    assert(0);
  else
    cout << "    Grabbed histograms." << endl;


  // sum histograms
  const UInt_t  nF = hPhi[0][0][0] -> GetNbinsX();
  const UInt_t  nH = hEta[0][0][0] -> GetNbinsX();
  const UInt_t  nP = hEta[0][0][0] -> GetNbinsX();
  const Float_t f1 = hPhi[0][0][0] -> GetBinLowEdge(1);
  const Float_t f2 = hPhi[0][0][0] -> GetBinLowEdge(nF + 1);
  const Float_t h1 = hEta[0][0][0] -> GetBinLowEdge(1);
  const Float_t h2 = hEta[0][0][0] -> GetBinLowEdge(nH + 1);
  const Float_t p1 = hPt[0][0][0]  -> GetBinLowEdge(1);
  const Float_t p2 = hPt[0][0][0]  -> GetBinLowEdge(nP + 1);

  const TString sSumBase[NHist] = {"hPhiSum", "hEtaSum", "hPtSum"};
  const TString sSumCut[NCut]   = {"BeforeQA", "AfterQA"};
  const TString sSumLvl[NLevel] = {"_par", "_det"};

  TH1D    *hSum[NHist][NLevel][NCut];
  TString sSumName[NHist];
  for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
    for (UInt_t iCut = 0; iCut < NCut; iCut++) {
      // make names
      sSumName[0] = sSumBase[0].Data();
      sSumName[1] = sSumBase[1].Data();
      sSumName[2] = sSumBase[2].Data();
      sSumName[0].Append(sSumCut[iCut].Data());
      sSumName[0].Append(sSumLvl[iLevel].Data());
      sSumName[1].Append(sSumCut[iCut].Data());
      sSumName[1].Append(sSumLvl[iLevel].Data());
      sSumName[2].Append(sSumCut[iCut].Data());
      sSumName[2].Append(sSumLvl[iLevel].Data());
      // declare histograms
      hSum[0][iLevel][iCut] = new TH1D(sSumName[0].Data(), "", nF, f1, f2);
      hSum[1][iLevel][iCut] = new TH1D(sSumName[1].Data(), "", nH, h1, h2);
      hSum[2][iLevel][iCut] = new TH1D(sSumName[2].Data(), "", nP, p1, p2);
      hSum[0][iLevel][iCut] -> Sumw2();
      hSum[1][iLevel][iCut] -> Sumw2();
      hSum[2][iLevel][iCut] -> Sumw2();
      // add files
      for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
        const UInt_t iWeight = (NEmbed - NFiles) + iFile;
        //hSum[0][iLevel][iCut] -> Add(hPhi[iFile][iLevel][iCut], weights[iWeight]);
        //hSum[1][iLevel][iCut] -> Add(hEta[iFile][iLevel][iCut], weights[iWeight]);
        //hSum[2][iLevel][iCut] -> Add(hPt[iFile][iLevel][iCut], weights[iWeight]);
        if (iLevel == 0) {
          hSum[0][iLevel][iCut] -> Add(hPhi[iFile][iLevel][iCut]);
          hSum[1][iLevel][iCut] -> Add(hEta[iFile][iLevel][iCut]);
          hSum[2][iLevel][iCut] -> Add(hPt[iFile][iLevel][iCut]);
        }
        else {
          hSum[0][iLevel][iCut] -> Add(hPhi[iFile][iLevel][iCut], weights[iWeight]);
          hSum[1][iLevel][iCut] -> Add(hEta[iFile][iLevel][iCut], weights[iWeight]);
          hSum[2][iLevel][iCut] -> Add(hPt[iFile][iLevel][iCut], weights[iWeight]);
        }
      }  // end file loop
    }  // end cut loop
  }  // end level loop
  cout << "    Summed histograms." << endl;


  // normalize histograms
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        const Double_t nTrgTot = nPi0Tot[iLevel];
        hSum[iHist][iLevel][iCut] -> Scale(1. / nTrgTot);
      }  // end cut loop
    }  // end level loop
  }  // end hist loop
  cout << "    Normalized histograms." << endl;


  // calculate efficiencies
  const UInt_t  nBins[NHist]    = {nF, nH, nP};
  const Float_t bin1[NHist]     = {f1, h1, p1};
  const Float_t bin2[NHist]     = {f2, h2, p2};
  const TString sResName[NHist] = {"hPhiResponse", "hEtaResponse", "hPtResponse"};
  const TString sEffName[NHist] = {"hPhiEfficiency", "hEtaEfficiency", "hPtEfficiency"};
  const TString sCopyBe[NHist]  = {"hPhiBeforeCopy", "hEtaBeforeCopy", "hPtBeforeCopy"};
  const TString sCopyAf[NHist]  = {"hPhiAfterCopy", "hEtaAfterCopy", "hPtAfterCopy"};

  TH1D    *hRes[NHist];
  TH1D    *hEff[NHist];
  TH1D    *hCopy[NHist][NLevel];
  Float_t weight(1.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hRes[iHist]     = new TH1D(sResName[iHist].Data(), "", nBins[iHist], bin1[iHist], bin2[iHist]);
    hEff[iHist]     = new TH1D(sEffName[iHist].Data(), "", nBins[iHist], bin1[iHist], bin2[iHist]);
    hRes[iHist]     -> Sumw2();
    hEff[iHist]     -> Sumw2();
    hRes[iHist]     -> Divide(hSum[iHist][1][1], hSum[iHist][1][0], weight, weight);
    hEff[iHist]     -> Divide(hSum[iHist][1][1], hSum[iHist][0][1], weight, weight);
    hCopy[iHist][0] = (TH1D*) hSum[iHist][1][0] -> Clone();
    hCopy[iHist][1] = (TH1D*) hSum[iHist][1][1] -> Clone();
    hCopy[iHist][0] -> SetName(sCopyBe[iHist].Data());
    hCopy[iHist][1] -> SetName(sCopyAf[iHist].Data());
  }
  cout << "    Calculated efficiencies." << endl;


  // fit efficiencies
  const UInt_t   fLinF(1);
  const UInt_t   fSizF(2);
  const UInt_t   fColF(859);
  const TString  sDraw("BMQ0");
  const TString  sVar[NHist]  = {"fPhi", "fEta", "fPt"};
  const TString  sFunc[NFunc] = {"Res", "Eff"};
  const TString  sForm[NHist] = {"[0]", "[0]", "[0]*(1-exp(-1.*[1]*x))"};
  const Double_t guess[NHist] = {0.5, 0.5, 0.87};
  const Double_t sigGuess(4.0);
  const Double_t fFitRange[2] = {f1, f2};
  const Double_t hFitRange[2] = {-0.7, 0.7};
  const Double_t pFitRange[2] = {1., 7.};
  const Double_t start[NHist] = {fFitRange[0], hFitRange[0], pFitRange[0]};
  const Double_t stop[NHist]  = {fFitRange[1], hFitRange[1], pFitRange[1]};

  TF1      *fFit[NHist][NFunc];
  Double_t amplitude[NHist][NFunc][NVal];
  Double_t sigma[NFunc][NVal];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {

      // make names
      TString sFuncName(sVar[iHist].Data());
      sFuncName += sFunc[iFunc].Data();

      // define fits
      fFit[iHist][iFunc] = new TF1(sFuncName.Data(), sForm[iHist].Data(), bin1[iHist], bin2[iHist]);
      fFit[iHist][iFunc] -> SetLineColor(fColF);
      fFit[iHist][iFunc] -> SetLineStyle(fLinF);
      fFit[iHist][iFunc] -> SetLineWidth(fSizF);
      fFit[iHist][iFunc] -> SetParameter(0, guess[iHist]);
      if (iHist == 2)
        fFit[iHist][iFunc] -> SetParameter(1, sigGuess);

      // fit histograms
      const Int_t kNotDraw = 1<<9;
      if (iFunc == 0) {
        hRes[iHist] -> Fit(fFit[iHist][iFunc], sDraw.Data(), "", start[iHist], stop[iHist]);
        amplitude[iHist][iFunc][0] = fFit[iHist][iFunc] -> GetParameter(0);
        amplitude[iHist][iFunc][1] = fFit[iHist][iFunc] -> GetParError(0);
        if (iHist == 2) {
          sigma[iFunc][0] = fFit[iHist][iFunc] -> GetParameter(1);
          sigma[iFunc][1] = fFit[iHist][iFunc] -> GetParError(1);
        }
        if (hRes[iHist] -> GetFunction(sFuncName.Data()))
          hRes[iHist] -> GetFunction(sFuncName.Data()) -> ResetBit(kNotDraw);
      }
      else {
        hEff[iHist] -> Fit(fFit[iHist][iFunc], sDraw.Data(), "", start[iHist], stop[iHist]);
        amplitude[iHist][iFunc][0] = fFit[iHist][iFunc] -> GetParameter(0);
        amplitude[iHist][iFunc][1] = fFit[iHist][iFunc] -> GetParError(0);
        if (iHist == 2) {
          sigma[iFunc][0] = fFit[iHist][iFunc] -> GetParameter(1);
          sigma[iFunc][1] = fFit[iHist][iFunc] -> GetParError(1);
        }
        if (hEff[iHist] -> GetFunction(sFuncName.Data()))
          hEff[iHist] -> GetFunction(sFuncName.Data()) -> ResetBit(kNotDraw);
      }

    }  // end function loop
  }  // end variable loop
  cout << "    Fit efficiencies." << endl;


  // set styles
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
    // cut responses
    hRes[iHist] -> SetLineColor(fColE);
    hRes[iHist] -> SetMarkerColor(fColE);
    hRes[iHist] -> SetMarkerStyle(fMarE);
    hRes[iHist] -> SetTitle(sResTitle[iHist].Data());
    hRes[iHist] -> SetTitleFont(fTxt);
    hRes[iHist] -> GetXaxis() -> SetTitle(sTitleX[iHist].Data());
    hRes[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hRes[iHist] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hRes[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hRes[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hRes[iHist] -> GetYaxis() -> SetTitle(sTitleYR[iHist].Data());
    hRes[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hRes[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hRes[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    // efficiencies
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
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        hSum[iHist][iLevel][iCut] -> SetLineColor(fColL[iLevel]);
        hSum[iHist][iLevel][iCut] -> SetMarkerColor(fColL[iLevel]);
        hSum[iHist][iLevel][iCut] -> SetMarkerStyle(fMarL[iLevel]);
        hSum[iHist][iLevel][iCut] -> SetTitle(sSumTitle[iHist].Data());
        hSum[iHist][iLevel][iCut] -> SetTitleFont(fTxt);
        hSum[iHist][iLevel][iCut] -> GetXaxis() -> SetTitle(sTitleX[iHist].Data());
        hSum[iHist][iLevel][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
        hSum[iHist][iLevel][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hSum[iHist][iLevel][iCut] -> GetXaxis() -> CenterTitle(fCnt);
        hSum[iHist][iLevel][iCut] -> GetXaxis() -> SetLabelSize(fLab);
        hSum[iHist][iLevel][iCut] -> GetYaxis() -> SetTitle(sTitleY[iHist].Data());
        hSum[iHist][iLevel][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
        hSum[iHist][iLevel][iCut] -> GetYaxis() -> CenterTitle(fCnt);
        hSum[iHist][iLevel][iCut] -> GetYaxis() -> SetLabelSize(fLab);
      }  // end cut loop 1
    }  // end level loop
    for (UInt_t iCut = 0; iCut < NCut; iCut++) {
      hCopy[iHist][iCut] -> SetLineColor(fColC[iCut]);
      hCopy[iHist][iCut] -> SetMarkerColor(fColC[iCut]);
      hCopy[iHist][iCut] -> SetMarkerStyle(fMarC[iCut]);
      hCopy[iHist][iCut] -> SetTitle(sSumTitle[iHist].Data());
      hCopy[iHist][iCut] -> SetTitleFont(fTxt);
      hCopy[iHist][iCut] -> GetXaxis() -> SetTitle(sTitleX[iHist].Data());
      hCopy[iHist][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
      hCopy[iHist][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hCopy[iHist][iCut] -> GetXaxis() -> CenterTitle(fCnt);
      hCopy[iHist][iCut] -> GetXaxis() -> SetLabelSize(fLab);
      hCopy[iHist][iCut] -> GetYaxis() -> SetTitle(sTitleY[iHist].Data());
      hCopy[iHist][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
      hCopy[iHist][iCut] -> GetYaxis() -> CenterTitle(fCnt);
      hCopy[iHist][iCut] -> GetYaxis() -> SetLabelSize(fLab);
    }  // end cut loop 2
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
  const TString sSig("#sigma_{p} = ");
  const TString sAmp[NHist]     = {"#epsilon_{#varphi} = ", "#epsilon_{#eta} = ", "#epsilon_{p} = "};
  const TString sResLeg[NCut]   = {"detector, before QA cuts", "detector, after QA cuts"};
  const TString sEffLeg[NLevel] = {"particle level", "detector level"};

  TString   sAmpRaw[NVal];
  TString   sSigRaw[NVal];
  TLegend   *lTrks[NHist][NFunc];
  TPaveText *pInfo[NHist][NFunc];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {

      TString sAmpTxt(sAmp[iHist].Data());
      TString sSigTxt(sSig.Data());
      for (UInt_t iVal = 0; iVal < NVal; iVal++) {
        sAmpRaw[iVal]  = "";
        sSigRaw[iVal]  = "";
        sAmpRaw[iVal] += amplitude[iHist][iFunc][iVal];
        sSigRaw[iVal] += sigma[iFunc][iVal];

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
      pInfo[iHist][iFunc] = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
      pInfo[iHist][iFunc] -> SetFillColor(fColP);
      pInfo[iHist][iFunc] -> SetLineColor(fColP);
      pInfo[iHist][iFunc] -> SetTextColor(fColT);
      pInfo[iHist][iFunc] -> SetTextFont(fTxt);
      pInfo[iHist][iFunc] -> SetTextAlign(fAlign);
      pInfo[iHist][iFunc] -> AddText(sSystem.Data());
      pInfo[iHist][iFunc] -> AddText(sTrgKin.Data());
      pInfo[iHist][iFunc] -> AddText(sAmpTxt.Data());
      if (iHist == 2)
        pInfo[iHist][iFunc] -> AddText(sSigTxt.Data());

      // legends
      if (iFunc == 0) {
        lTrks[iHist][iFunc] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
        lTrks[iHist][iFunc] -> SetFillColor(fColP);
        lTrks[iHist][iFunc] -> SetLineColor(fColP);
        lTrks[iHist][iFunc] -> SetTextFont(fTxt);
        lTrks[iHist][iFunc] -> AddEntry(hCopy[iHist][0], sResLeg[0].Data());
        lTrks[iHist][iFunc] -> AddEntry(hCopy[iHist][1], sResLeg[1].Data());
      }
      else {
        lTrks[iHist][iFunc] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
        lTrks[iHist][iFunc] -> SetFillColor(fColP);
        lTrks[iHist][iFunc] -> SetLineColor(fColP);
        lTrks[iHist][iFunc] -> SetTextFont(fTxt);
        lTrks[iHist][iFunc] -> AddEntry(hSum[iHist][0][1], sEffLeg[0].Data());
        lTrks[iHist][iFunc] -> AddEntry(hSum[iHist][1][1], sEffLeg[1].Data());
      }

    }  // end function loop
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
  const TString sPads[2]        = {"pRatio", "pTracks"};
  const TString sResPlot[NHist] = {"cPhiResponse", "cEtaResponse", "cPtResponse"};
  const TString sEffPlot[NHist] = {"cPhiEfficiency", "cEtaEfficiency", "cPtEfficiency"};

  TPad    *pPad[NHist][2];
  TCanvas *cRes[NHist];
  TCanvas *cEff[NHist];
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    // response
    cRes[iHist]    = new TCanvas(sResPlot[iHist].Data(), "", width, height);
    pPad[iHist][0] = new TPad(sPads[0].Data(), "", xPad[0], yPad[0], xPad[1], yPad[1]);
    pPad[iHist][1] = new TPad(sPads[1].Data(), "", xPad[2], yPad[2], xPad[3], yPad[3]);
    pPad[iHist][0]  -> SetGrid(grid, grid);
    pPad[iHist][0]  -> SetTicks(ticks, ticks);
    pPad[iHist][1]  -> SetGrid(grid, grid);
    pPad[iHist][1]  -> SetTicks(ticks, ticks);
    pPad[iHist][1]  -> SetLogy(log[iHist]);
    cRes[iHist]     -> cd();
    pPad[iHist][0]  -> Draw();
    pPad[iHist][1]  -> Draw();
    pPad[iHist][0]  -> cd();
    hRes[iHist]     -> Draw();
    pInfo[iHist][0] -> Draw();
    pPad[iHist][1]  -> cd();
    hCopy[iHist][0] -> Draw();
    hCopy[iHist][1] -> Draw("same");
    lTrks[iHist][0] -> Draw();
    cRes[iHist]     -> Write();
    cRes[iHist]     -> Close();

    // efficiency
    cEff[iHist]    = new TCanvas(sEffPlot[iHist].Data(), "", width, height);
    pPad[iHist][0] = new TPad(sPads[0].Data(), "", xPad[0], yPad[0], xPad[1], yPad[1]);
    pPad[iHist][1] = new TPad(sPads[1].Data(), "", xPad[2], yPad[2], xPad[3], yPad[3]);
    pPad[iHist][0]    -> SetGrid(grid, grid);
    pPad[iHist][0]    -> SetTicks(ticks, ticks);
    pPad[iHist][1]    -> SetGrid(grid, grid);
    pPad[iHist][1]    -> SetTicks(ticks, ticks);
    pPad[iHist][1]    -> SetLogy(log[iHist]);
    cEff[iHist]       -> cd();
    pPad[iHist][0]    -> Draw();
    pPad[iHist][1]    -> Draw();
    pPad[iHist][0]    -> cd();
    hEff[iHist]       -> Draw();
    pInfo[iHist][1]   -> Draw();
    pPad[iHist][1]    -> cd();
    hSum[iHist][0][1] -> Draw();
    hSum[iHist][1][1] -> Draw("same");
    lTrks[iHist][1]   -> Draw();
    cEff[iHist]       -> Write();
    cEff[iHist]       -> Close();
  }


  // make directories, save histograms
  const TString sDir[NLevel] = {"particle", "detector"};

  TDirectory *dLvl[NLevel];
  for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
    dLvl[iLevel] = (TDirectory*) fOutput -> mkdir(sDir[iLevel].Data());
    dLvl[iLevel] -> cd();
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      hRes[iHist] -> Write();
      hEff[iHist] -> Write();
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        hSum[iHist][iLevel][iCut] -> Write();
      }  // end cut loop
    }  // end hist loop
  }  // end level loop
  cout << "    Made directories and saved histograms." << endl;


  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iFile = 0; iFile < NFiles; iFile++) {
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      fInput[iFile][iLevel] -> cd();
      fInput[iFile][iLevel] -> Close();
    }
  }
  cout << "  Calculation finished!\n" << endl;

}

// End ------------------------------------------------------------------------
