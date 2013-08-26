#ifndef _DRAWDATAMC
#define _DRAWDATAMC
#include <TLegend.h>
#include <TStyle.h>
#include <TFormula.h>
#include <TF1.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TObject.h>
#include <string.h>
#include <algorithm>
#include <TBranch.h>
#include <TString.h>
#include <TBits.h>
#include <TMath.h>
#include "TROOT.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <TLatex.h>


/*
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootBardak.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootBeamSpot.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootCluster.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootDummyEvent.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootEcalRecHit.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootElectron.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootEvent.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootJet.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMCParticle.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMCPhoton.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMET.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootMuon.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootParticle.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootPhoton.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootRun.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootSignalEvent.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootSuperCluster.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootTopTop.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootTrack.h"
#include "/sps/cms/obondu/CMSSW_4_1_2/src/UserCode/IpnTreeProducer/interface/TRootVertex.h"
*/

#pragma optimize 0

using namespace std;

//void DrawDataMCplot(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *QCDMu_miniTree, TTree* TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, bool doFit = false);
void DrawDataMCplot(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree* TTJets_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, double integratedLuminosity, bool doFit = false);
//void DrawDataMCplot(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, bool doFit = false);

//void DrawDataMCplot_TH1I(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var, string pic, string limits, string cut, string name, string Title, bool inlog, bool drawUnderOverFsubleading, TCanvas *c1, bool doFit = false);

//void DrawDataMCplot_TH2F(TTree *Data_miniTree, TTree *FSR_DYToMuMu_miniTree, TTree *nonFSR_DYToMuMu_miniTree, TTree *TTJets_miniTree, TTree *WJetsToLNu_miniTree, string var1, string var2, string limits, string cut, string name, string Title_var1, string Title_var2, string pic, bool inlog_var1, bool inlog_var2, TCanvas *c1, bool doFit = false);

#endif
