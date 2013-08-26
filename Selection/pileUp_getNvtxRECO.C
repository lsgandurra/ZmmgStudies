// ROOT HEADERS
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TBits.h"
#include "TMath.h"
#include "TSystem.h"

// C++ HEADERS
#include <sstream>
#include <iostream>
#include <fstream>
#include <utility>
//#pragma optimize 0

// TOTOANA HEADERS
#include "interface/TRootBardak.h"
#include "interface/TRootBeamSpot.h"
#include "interface/TRootCluster.h"
#include "interface/TRootDummyEvent.h"
#include "interface/TRootEcalRecHit.h"
#include "interface/TRootElectron.h"
#include "interface/TRootEvent.h"
#include "interface/TRootJet.h"
#include "interface/TRootMCParticle.h"
#include "interface/TRootMCPhoton.h"
#include "interface/TRootMET.h"
#include "interface/TRootMuon.h"
#include "interface/TRootParticle.h"
#include "interface/TRootPhoton.h"
#include "interface/TRootRun.h"
#include "interface/TRootSignalEvent.h"
#include "interface/TRootSuperCluster.h"
#include "interface/TRootTopTop.h"
#include "interface/TRootTrack.h"
#include "interface/TRootVertex.h"


// *****************************************************************************************************
// ******************* Execute bash command line and get output
// *****************************************************************************************************
std::string exec(char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
                result += buffer;
    }
    pclose(pipe);
    return result;
}



Int_t nVertices;

//int Selection_miniTree()
int main(int argc, char *argv[])
//int main()
{

	cout << "argc= " << argc << endl;
	for(int iarg = 0 ; iarg < argc; iarg++)
	{
		cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
	}
  if( argc == 1 )
  {
    cerr << "arguments should be passed !! sample (ntotjob) (outputname) (ijob) (isZgammaMC)" << endl;
    return 1;
  }



	
 // ******************************************
  // First argument is sample
  // ******************************************
  char* sample_char = argv[1];
//  char* sample_char2 = argv[1];

  // ******************************************
  // Optional argument : output root file
  // ******************************************
  string sample = argv[1];
  if( argc > 2 )
  {
    sample = argv[2];
  }

  // ******************************************
  // Optional argument : ntotjob
  // ******************************************
  int ntotjob = -1;
  if( argc > 3 )
  {
    std::stringstream ss ( argv[3] );
    ss >> ntotjob;
  }
  


  // ******************************************
  // Optional argument : ijob
  // ******************************************
  int ijob = -1;
  if( argc > 4 )
  {
    std::stringstream ss ( argv[4] );
    ss >> ijob;
  }

  // ******************************************
  // Optional argument : isZgammaMC
  // ******************************************
  int isZgammaMC = 0;
  if( argc > 5 )
  {
    std::stringstream ss ( argv[5] );
    ss >> isZgammaMC;
  }



	gSystem->Load("libToto.so");
	bool doHLT										= false;
	bool doMC										 = true;
	bool doJetMC									= false;
	bool doMETMC									= false;
	bool doPDFInfo								= false;
	bool doSignalMuMuGamma				= true;
	bool doSignalTopTop					 = false;
	bool doPhotonConversionMC		 = false;
	bool doBeamSpot							 = false;
	bool doPrimaryVertex					= true;
	bool doZeePrimaryVertex			 = false;
	bool doTrack									= false;
	bool doJet										= false;
	bool doMuon									 = true;
	bool doElectron							 = false;
	bool doPhoton								 = true;
	bool doCluster								= false;
	bool doPhotonConversion			 = false;
	bool doMET										= false;
	bool doBardak								 = false;
	bool doPhotonVertexCorrection = false;
	bool doPhotonIsolation				= false;

	// DATASET	
	TChain *inputEventTree = new TChain("eventTree");
	TChain *inputRunTree = new TChain("runTree");

  string line;
  string filename = Form("listFiles_%s", sample_char);
  ifstream myfile(filename.c_str());
  string nlines_ = exec(Form("wc -l %s | awk '{print $1}'", filename.c_str() ));
  int nlines = atoi(nlines_.c_str());
  string protocol = "dcap://ccdcapcms.in2p3.fr:22125";
//  int nfilesPerJob = ceil( nlines / ntotjob ); 
//  if( nlines % ntotjob != 0 ) nfilesPerJob++;
//  int ilineBegin = nfilesPerJob * ijob + 1;
//  int ilineEnd = min( nlines + 1  ,   nfilesPerJob * (ijob + 1) + 1);
  int nJobsWithExtraFile = nlines % ntotjob;
  int nfilesPerJob = ceil( nlines / ntotjob );
  int ilineBegin = nfilesPerJob * ijob + 1;
  int ilineEnd = min( nlines + 1  ,   nfilesPerJob * (ijob + 1) + 1);
  if( nJobsWithExtraFile != 0 )
  {
    if( ijob < nJobsWithExtraFile )
    {
      ilineBegin = (nfilesPerJob + 1 ) * ijob + 1;
      ilineEnd = min( nlines + 1  ,   (nfilesPerJob + 1) * (ijob + 1) + 1);
    } else {
      ilineBegin = nJobsWithExtraFile * (nfilesPerJob + 1) + nfilesPerJob * ( ijob - nJobsWithExtraFile ) + 1;
      ilineEnd = min( nlines + 1  ,  nJobsWithExtraFile * (nfilesPerJob + 1) + nfilesPerJob * (ijob - nJobsWithExtraFile + 1) + 1);
    }
  }

  cout << "ijob= " << ijob << endl;
  cout << "ntotjob= " << ntotjob << endl;
  cout << "nlines= " << nlines << endl;
  cout << "nJobsWithExtraFile= " << nJobsWithExtraFile << endl;
  cout << "ilineBegin= " << ilineBegin << endl;
  cout << "ilineEnd= " << ilineEnd << endl;

int iline = 0;
if( ntotjob == 9999 )
{
  //inputEventTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
  //inputRunTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
	inputEventTree->Add("DYToMuMu_test_2011.root");
	inputRunTree->Add("DYToMuMu_test_2011.root");
	cout << "addition ok "<< endl;
} else {
  if (myfile.is_open())
  {
    while( myfile.good() )
    {
      iline++;
      getline(myfile, line);
      if( line == "") continue; // EOF !
      if( iline >= ilineBegin && iline < ilineEnd )
      {
        cout << "Adding file #" << iline << " ( / " << nlines << ") : " << line << endl;
        inputEventTree->Add(Form("%s%s", protocol.c_str(), line.c_str()));
        inputRunTree->Add(Form("%s%s", protocol.c_str(), line.c_str()));
      }
    }
    myfile.close();
  } else {
    cout << "Unable to open file" << endl;
    return 987;
  }
}


//	inputEventTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
//	inputRunTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));

// INSERTFILES

	TFile* OutputRootFile = new TFile(Form("miniTree_pileUp_%s_part%i.root", sample.c_str(), ijob), "RECREATE");
	
	TBranch* event_br = 0;
	TRootEvent* event = 0;
	inputEventTree->SetBranchAddress("Event", &event, &event_br);
	inputEventTree->SetBranchStatus("Event", 1);

  TBranch* run_br = 0;
  TRootRun* runInfos = 0;
  inputRunTree->SetBranchAddress("runInfos", &runInfos, &run_br);
  inputRunTree->SetBranchStatus("runInfos", 1);
	
	TBranch* mcParticles_br = 0;
	TClonesArray* mcParticles = new TClonesArray("TRootMCParticle", 0);
	if(doMC)
	{
		inputEventTree->SetBranchAddress("MCParticles", &mcParticles, &mcParticles_br);
		inputEventTree->SetBranchStatus("MCParticles", 1);
	}
/*	
	TBranch* genJets_br = 0;
	TClonesArray* genJets = new TClonesArray("TRootParticle", 0);
	if(doJetMC)
	{
		inputEventTree->SetBranchAddress("genJets", &genJets, &genJets_br);
		inputEventTree->SetBranchStatus("genJets", 1);
	}
	
		TBranch* genMETs_br = 0;
		TClonesArray* genMETs = new TClonesArray("TRootParticle", 0);
	if(doMETMC)
	{
		inputEventTree->SetBranchAddress("genMETs", &genMETs, &genMETs_br);
		inputEventTree->SetBranchStatus("genMETs", 1);
	}
*/	
		TBranch* mcSignalMuMuGamma_br = 0;
		TRootSignalEvent* mcMuMuGammaEvent = 0;
	if(doSignalMuMuGamma)
	{
		inputEventTree->SetBranchAddress("MuMuGamma", &mcMuMuGammaEvent, &mcSignalMuMuGamma_br);
		inputEventTree->SetBranchStatus("MuMuGamma", 1);
	}
/*	
		TBranch* mcTopTopEvent_br = 0;
		TRootSignalEvent* mcTopTopEvent = 0;
	if(doSignalTopTop)
	{
		inputEventTree->SetBranchAddress("rootMCTopTop", &mcTopTopEvent, &mcTopTopEvent_br);
		inputEventTree->SetBranchStatus("rootMCTopTop", 1);
	}
*/	
		TBranch* mcPhotons_br = 0;
		TClonesArray* mcPhotons = new TClonesArray("TRootMCPhoton", 0);
	if(doPhotonConversionMC)
	{
		inputEventTree->SetBranchAddress("MCPhotons", &mcPhotons, &mcPhotons_br);
		inputEventTree->SetBranchStatus("MCPhotons", 1);
	}
/*	
		TBranch* beamSpot_br = 0;
		TRootBeamSpot* beamSpot = 0;
	if(doBeamSpot)
	{
		inputEventTree->SetBranchAddress("BeamSpot", &beamSpot, &beamSpot_br);
		inputEventTree->SetBranchStatus("BeamSpot", 1);
	}
*/
		TBranch* vertices_br = 0;
		TClonesArray* vertices = new TClonesArray("TRootVertex", 0);
	if(doPrimaryVertex)
	{
		inputEventTree->SetBranchAddress("Vertices", &vertices, &vertices_br);
		inputEventTree->SetBranchStatus("Vertices", 1);
	}
/*	
		TBranch* zeeVertices_br = 0;
		TClonesArray* zeeVertices = new TClonesArray("TRootVertex", 0);
	if(doZeePrimaryVertex)
	{
		inputEventTree->SetBranchAddress("ZeeVertices", &zeeVertices, &zeeVertices_br);
		inputEventTree->SetBranchStatus("ZeeVertices", 1);
	}
	
		TBranch* tracks_br = 0;
		TClonesArray* tracks = new TClonesArray("TRootTrack", 0);
	if(doTrack)
	{
		inputEventTree->SetBranchAddress("Tracks", &tracks, &tracks_br);
		inputEventTree->SetBranchStatus("Tracks", 1);
	}
	
		TBranch* jets_br = 0;
		TClonesArray* jets = new TClonesArray("TRootJet", 0);
	if(doJet)	{
		inputEventTree->SetBranchAddress("ak5CaloJets", &jets, &jets_br);
		inputEventTree->SetBranchStatus("ak5CaloJets", 1);
	}
*/	
		TBranch* muons_br = 0;
		TClonesArray* muons = new TClonesArray("TRootMuon", 0);
	if(doMuon)
	{
		inputEventTree->SetBranchAddress("muons", &muons, &muons_br);
		inputEventTree->SetBranchStatus("muons", 1);
	}
/*	
		TBranch* electrons_br = 0;
		TClonesArray* electrons = new TClonesArray("TRootElectron", 0);
	if(doElectron)
	{
		inputEventTree->SetBranchAddress("gsfElectrons", &electrons, &electrons_br);
		inputEventTree->SetBranchStatus("gsfElectrons", 1);
	}
*/	
		TBranch* photons_br = 0;
		TClonesArray* photons = new TClonesArray("TRootPhoton", 0);
	if(doPhoton)
	{
		inputEventTree->SetBranchAddress("photons", &photons, &photons_br);
		inputEventTree->SetBranchStatus("photons", 1);
	}
/*	
		TBranch* clusters_br = 0;
		TClonesArray* clusters = new TClonesArray("TRootCluster", 0);
	if(doCluster)
	{
		inputEventTree->SetBranchAddress("BasicClusters", &clusters, &clusters_br);
		inputEventTree->SetBranchStatus("BasicClusters", 1);
	}
	
		TBranch* superClusters_br = 0;
		TClonesArray* superClusters = new TClonesArray("TRootSuperCluster", 0);
	if(doCluster)
	{
		inputEventTree->SetBranchAddress("SuperClusters", &superClusters, &superClusters_br);
		inputEventTree->SetBranchStatus("SuperClusters", 1);
	}
	
		TBranch* conversions_br = 0;
		TClonesArray* conversionTracks = new TClonesArray("TRootTrack", 0);
	if(doPhotonConversion)
	{
		inputEventTree->SetBranchAddress("ConversionTracks", &conversionTracks, &conversions_br);
		inputEventTree->SetBranchStatus("ConversionTracks", 1);
	}
	
		TBranch* met_br = 0;
		TClonesArray* met = new TClonesArray("TRootMET", 0);
	if(doMET)
	{
		inputEventTree->SetBranchAddress("met", &met, &met_br);
		inputEventTree->SetBranchStatus("met", 1);
	}
	
		TBranch* bardak_br = 0;
		TRootBardak* bardak = 0;
	if(doBardak)
	{
		inputEventTree->SetBranchAddress("Bardak", &bardak, &bardak_br);
		inputEventTree->SetBranchStatus("Bardak", 1);
	}
*/
	TTree* miniTree = new TTree("miniTree","NVertex information");

	// ____________________________________________
	// Event information
	// ____________________________________________
	miniTree->Branch("nVertices", &nVertices, "nVertices/I");
	
	// SETUP PARAMETERS	
	unsigned int NbEvents = (int)inputEventTree->GetEntries();
//	unsigned int NbEvents = 10000;
//	bool powheg = false;
	bool powheg = true;
	bool signal = false;
	bool stew = false;
	bool zjet_veto = (bool)(isZgammaMC >= 1);
	cout << "Nb of events : " << NbEvents << endl;
	cout << "Signal is: " << signal <<endl;
	cout << "Stew is: " << stew << endl;
	cout << "ZJet veto is: " << zjet_veto << endl;
	int nBeforeAllCuts, nAfterCutPthatFilter, nAfterCutCSA07ID, nAfterCutZJETVETO, nAfterVeryLooseMMG, nAfterLooseMMG, nAfterCut1c, nAfterTightMMG, nAfterCut1e, nAfterCut2a, nAfterCut2b, nAfterCut2c, nAfterCut3, nAfterCut4, nAfterCut5, nAfterCut6, nAfterCut7, nAfterCut8, nAfterCut9, nAfterCut10, nSelected;
	nBeforeAllCuts = nAfterCutPthatFilter = nAfterCutCSA07ID = nAfterCutZJETVETO = nAfterVeryLooseMMG = nAfterLooseMMG = nAfterCut1c = nAfterTightMMG = nAfterCut1e = nAfterCut2a = nAfterCut2b = nAfterCut2c = nAfterCut3 = nAfterCut4 = nAfterCut5 = nAfterCut6 = nAfterCut7 = nAfterCut8 = nAfterCut9 = nAfterCut10 = nSelected = 0;

	double minPtHat = -100;
  double maxPtHat = 1000000;
  int verbosity = 0;
	int TOTALnbMuonsAfterID[12] = {0};
	int TOTALnbEventsAfterMuonID[12] = {0};
	int TOTALnbDimuonsAfterID[3] = {0};
	int TOTALnbEventsAfterDimuonID[3] = {0};
	int TOTALnbPhotonsAfterID[6] = {0};
	int TOTALnbEventsAfterPhotonID[6] = {0};
	int TOTALnbMuMuGammaAfterID[8] = {0};
	int TOTALnbEventsAfterMuMuGammaID[8] = {0};

//	TH1D* pileup = new TH1D("pileup", "pileup", 51, -0.5, 50.5);
//	TH1D* pileup = new TH1D("pileup", "pileup", 36, -0.5, 35.5);
	TH1D* pileup = new TH1D("pileup", "pileup", 1000, 0, 25);


  int NbEventsPerJob = NbEvents;
  int NbEventsBegin = 0;
  int NbEventsEnd = NbEvents;
  if( ntotjob == 9999 && ijob != -1 )
  {
    NbEventsPerJob = 200000;
//    NbEventsPerJob = 100;
    NbEventsBegin = ijob * NbEventsPerJob;
    NbEventsEnd = min( (ijob + 1)* NbEventsPerJob , (int)NbEvents);
    NbEvents = NbEventsEnd - NbEventsBegin ;
  }
  cout << "NbEventsBegin= " << NbEventsBegin << "\tNbEventsEnd= " << NbEventsEnd << "\tNbEventsPerJob= " << NbEventsPerJob << endl;


// ***************************************************************************************************
// ***************************************************************************************************
// ********************************** LOOP over events ***********************************************
// ***************************************************************************************************
// ***************************************************************************************************
  for(unsigned int ievt=NbEventsBegin; ievt<NbEventsEnd; ievt++)
//	for(unsigned int ievt=0; ievt<NbEvents; ievt++)
//	for(unsigned int ievt=0; ievt<20000; ievt++)
	{
		if(verbosity>4) cout << "analysing event ievt= " << ievt << endl;
		nBeforeAllCuts++;
		int nprint = (int)((double)NbEvents/(double)100.0);
		if( (ievt % nprint)==0 ){ cout<< ievt <<" events done over "<<NbEvents<<" ( "<<ceil((double)ievt/(double)NbEvents*100)<<" \% )"<<endl; }

		inputEventTree->GetEvent(ievt);

		// ____________________________________________
		// Event information
		// ____________________________________________
//		nVertices = vertices->GetEntries();
//		if( (isZgammaMC >= 1) )
//    {
			//nVertices = event->nInTimePUVertices();
			nVertices = event->pu_TrueNumInteractions();
//    }
		pileup->Fill(nVertices);
		miniTree->Fill();

	} // fin boucle sur evts LOOP
// Writing stuff out
	pileup->Write();
	OutputRootFile->Write();
	OutputRootFile->Close();

	delete OutputRootFile;
	delete inputEventTree;
	delete inputRunTree;



	return 0;

}
