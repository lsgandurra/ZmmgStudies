#include "Muons_miniTree.h"
#include "rochcor_v4_new.C"

	// ____________________________________________
	// Event information
	// ____________________________________________
	ULong64_t iEvent, iEventID, iLumiID, iRunID;

	Int_t isMM;

	Int_t nVertices;
	Int_t nGenVertices;
	Float_t weight_pileUp, weight_Xsection;

	// ____________________________________________
	// Muon variables
	// ____________________________________________
	Int_t NbMuons;

	Float_t Pt_allMuons, Eta_allMuons, Phi_allMuons, Charge_allMuons;
// (M minus charge, P plus charge), (F far, N near), (L leading, S subleading)
	Float_t MuonM_Pt, MuonP_Pt, MuonL_Pt, MuonS_Pt;
	Float_t MuonM_Eta, MuonP_Eta, MuonL_Eta, MuonS_Eta;
	Float_t MuonM_Phi, MuonP_Phi, MuonL_Phi, MuonS_Phi;
	Int_t MuonL_Charge, MuonS_Charge;

	Float_t MuonM_isoR03_emEt, MuonP_isoR03_emEt, MuonL_isoR03_emEt, MuonS_isoR03_emEt;
	Float_t MuonM_isoR03_hadEt, MuonP_isoR03_hadEt, MuonL_isoR03_hadEt, MuonS_isoR03_hadEt;
	Float_t MuonM_isoR03_hoEt, MuonP_isoR03_hoEt, MuonL_isoR03_hoEt, MuonS_isoR03_hoEt;
	Float_t MuonM_isoR03_nJets, MuonP_isoR03_nJets, MuonL_isoR03_nJets, MuonS_isoR03_nJets;
	Float_t MuonM_isoR03_nTracks, MuonP_isoR03_nTracks, MuonL_isoR03_nTracks, MuonS_isoR03_nTracks;
	Float_t MuonM_isoR03_sumPt, MuonP_isoR03_sumPt, MuonL_isoR03_sumPt, MuonS_isoR03_sumPt;

	Float_t MuonM_isoR05_emEt, MuonP_isoR05_emEt, MuonL_isoR05_emEt, MuonS_isoR05_emEt;
	Float_t MuonM_isoR05_hadEt, MuonP_isoR05_hadEt, MuonL_isoR05_hadEt, MuonS_isoR05_hadEt;
	Float_t MuonM_isoR05_hoEt, MuonP_isoR05_hoEt, MuonL_isoR05_hoEt, MuonS_isoR05_hoEt;
	Float_t MuonM_isoR05_nJets, MuonP_isoR05_nJets, MuonL_isoR05_nJets, MuonS_isoR05_nJets;
	Float_t MuonM_isoR05_nTracks, MuonP_isoR05_nTracks, MuonL_isoR05_nTracks, MuonS_isoR05_nTracks;
	Float_t MuonM_isoR05_sumPt, MuonP_isoR05_sumPt, MuonL_isoR05_sumPt, MuonS_isoR05_sumPt;

	Float_t MuonM_E, MuonP_E, MuonL_E, MuonS_E;
	Float_t MuonM_Px, MuonP_Px, MuonL_Px, MuonS_Px;
	Float_t MuonM_Py, MuonP_Py, MuonL_Py, MuonS_Py;
	Float_t MuonM_Pz, MuonP_Pz, MuonL_Pz, MuonS_Pz;

	// ____________________________________________
	// mumu information
	// ____________________________________________
	
	Float_t Mmumu;
	Float_t Ptmumu;

	// ____________________________________________
	// MC Truth
	// ___________________________________________

	Float_t MuonM_MC_E, MuonM_MC_Px, MuonM_MC_Py, MuonM_MC_Pz, MuonM_MC_Phi, MuonM_MC_Eta, MuonM_MC_Pt;
	Float_t MuonP_MC_E, MuonP_MC_Px, MuonP_MC_Py, MuonP_MC_Pz, MuonP_MC_Phi, MuonP_MC_Eta, MuonP_MC_Pt;
	Float_t MuonL_MC_E, MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_Phi, MuonL_MC_Eta, MuonL_MC_Pt;
	Float_t MuonS_MC_E, MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_Phi, MuonS_MC_Eta, MuonS_MC_Pt;
	Float_t Mmumu_Muons_MC;

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
		cerr << "arguments should be passed !! sample (outputname) (ntotjob) (ijob) (isZgammaMC) (lumi_set) (pu_set) (low m_mumu cut) (high m_mumu cut) (applyMuonScaleCorrection) (muon correction sys)" << endl;
		return 1;
	}

	// ******************************************
	// First argument is sample
	// ******************************************
	char* sample_char = argv[1];
//	char* sample_char2 = argv[1];

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
	// Optional argument : isZgammaMC (1: FSR -- 2: nonFSR -- 3: MC info)
	// ******************************************
	int isZgammaMC = 0;
	if( argc > 5 )
	{
		std::stringstream ss ( argv[5] );
		ss >> isZgammaMC;
	}

	// ******************************************
	// Optional argument : lumi set used
	// ******************************************
	string lumi_set = "";
	double integratedLuminosity = 1.0;
	int runopt = 0;
	if( argc > 6 )
	{
		lumi_set = argv[6];
// LAST UPDATE: May 02nd 2012 with pixel lumi from https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc: cvs co	-r V03-05-05 RecoLuminosity/LumiDB
		if( lumi_set == "May10" ) integratedLuminosity = 215.552;
		if( lumi_set == "Promptv4" ) integratedLuminosity = 951.716;
		if( lumi_set == "July05" ) integratedLuminosity = 1.157*1000.0;
		if( lumi_set == "Aug05" ) integratedLuminosity = 389.876;
		if( lumi_set == "Oct03" ) integratedLuminosity = 706.719;
		if( lumi_set == "2011A" ){ integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719; runopt = 0;}
		if( lumi_set == "2011A_rereco" ){ integratedLuminosity = 2.221*1000.0; runopt = 0;}
		if( lumi_set == "2011B" ){ integratedLuminosity = 2.714*1000.0; runopt = 1;}
		if( lumi_set == "2011" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0;
		if( lumi_set == "2011_rereco" ) integratedLuminosity = 2.221*1000.0 +	2.714*1000.0;
	}

	// ******************************************
	// Optional argument : pile-up scenario used
	// ******************************************
	string pu_set = "";
	if( argc > 7 )
	{
		pu_set = argv[7];
	}

	// ******************************************
	// Optional argument : low_m_mumu
	// ******************************************
	double low_m_mumu = 60.0;
	if( argc > 8 )
	{
		std::stringstream ss ( argv[8] );
		ss >> low_m_mumu;
	}

	// ******************************************
	// Optional argument : high_m_mumu
	// ******************************************
	double high_m_mumu = 120.0;
	if( argc > 9 )
	{
		std::stringstream ss ( argv[9] );
		ss >> high_m_mumu;
	}
	
	// ******************************************
	// Optional argument : applyMuonScaleCorrection: 0)nothing 1) MuScleFit 2)SIDRA 3)Rochester(to be implemented)
	// ******************************************
	int applyMuonScaleCorrection = 0;
	if( argc > 10 )
	{
		std::stringstream ss ( argv[10] );
		ss >> applyMuonScaleCorrection;
				if( applyMuonScaleCorrection == 0 ) cout << "No muon correction will be applied" << endl;
				if( applyMuonScaleCorrection == 1 ) cout << "MuScleFit muon corrections will be applied upstream" << endl;
				if( applyMuonScaleCorrection == 2 ) cout << "SIDRA muon corrections will be applied upstream" << endl;
				if( applyMuonScaleCorrection == 3 ) cout << "Rochester muon corrections will be applied upstream" << endl;
	}

	// ******************************************
	// Optional argument is extra muon correction smearing (Rochester only)
	// ******************************************
	double sysdev = 0.0;
	if( argc > 11 )
	{
		std::stringstream ss ( argv[11] );
		ss >> sysdev;
	}




	TTimeStamp *time = new TTimeStamp();
	TRandom3* generator = new TRandom3(time->GetNanoSec());
		// TODO

//	char* inputfile = argv[6];

//	TProof * p = TProof::Open("ccaplmaster.in2p3.fr");
	gSystem->Load("libToto.so");
	gROOT->ProcessLine(".L rochcor_v4_new.h+");
//	gSystem->Load("libFWCoreFWLite.so");
//	gSystem->Load("libDataFormatsFWLite.so");
//	AutoLibraryLoader::enable();
//	#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
//	#include "AutoLibraryLoader.h"
//	TruncatedPyramid a;

	bool doHLT										= false;
	bool doMC										 = (bool)(isZgammaMC >= 1);
	bool doJetMC									= false;
	bool doMETMC									= false;
	bool doPDFInfo								= false;
	bool doSignalMuMuGamma				= (bool)(isZgammaMC == 1 || isZgammaMC == 2);
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
	bool doCluster								= true;
	bool doPhotonConversion			 = true;
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
//	int nfilesPerJob = ceil( nlines / ntotjob ); 
//	if( nlines % ntotjob != 0 ) nfilesPerJob++;
//	int ilineBegin = nfilesPerJob * ijob + 1;
//	int ilineEnd = min( nlines + 1	,	 nfilesPerJob * (ijob + 1) + 1);
	int nJobsWithExtraFile = nlines % ntotjob;
	int nfilesPerJob = ceil( nlines / ntotjob );
	int ilineBegin = nfilesPerJob * ijob + 1;
	int ilineEnd = min( nlines + 1	,	 nfilesPerJob * (ijob + 1) + 1);
	if( nJobsWithExtraFile != 0 )
	{
		if( ijob < nJobsWithExtraFile )
		{
			ilineBegin = (nfilesPerJob + 1 ) * ijob + 1;
			ilineEnd = min( nlines + 1	,	 (nfilesPerJob + 1) * (ijob + 1) + 1);
		} else {
			ilineBegin = nJobsWithExtraFile * (nfilesPerJob + 1) + nfilesPerJob * ( ijob - nJobsWithExtraFile ) + 1;
			ilineEnd = min( nlines + 1	,	nJobsWithExtraFile * (nfilesPerJob + 1) + nfilesPerJob * (ijob - nJobsWithExtraFile + 1) + 1);
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
	inputEventTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
	inputRunTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
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
//	inputEventTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/SkimmedTotoSamples/%s/%s*root", sample_char, sample_char));
//	inputRunTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/SkimmedTotoSamples/%s/%s*root", sample_char, sample_char));
// INSERTFILES

	TFile* OutputRootFile = new TFile(Form("miniTreeMuons_%s_part%i.root", sample.c_str(), ijob), "RECREATE");
	
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
/*	
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

	int* NumWantedHLTnames;

	

	// ____________________________________________
	// preparing the tree
	// ____________________________________________

	TTree* miniTree = new TTree("miniTree","Mu Mu informations");

	// ____________________________________________
	// Event information
	// ____________________________________________
	
	miniTree->Branch("iEvent", &iEvent, "iEvent/l");
	miniTree->Branch("iEventID", &iEventID, "iEventID/l");
	miniTree->Branch("iLumiID", &iLumiID, "iLumiID/l");
	miniTree->Branch("iRunID", &iRunID, "iRunID/l");

	miniTree->Branch("isMM", &isMM, "isMM/I");
	miniTree->Branch("nVertices", &nVertices, "nVertices/I");
	miniTree->Branch("nGenVertices", &nGenVertices, "nGenVertices/I");
	miniTree->Branch("weight_pileUp", &weight_pileUp, "weight_pileUp/F");
	miniTree->Branch("weight_Xsection", &weight_Xsection, "weight_Xsection/F");

	// ____________________________________________
	// Muon variables
	// ____________________________________________

	miniTree->Branch("NbMuons", &NbMuons, "NbMuons/I");

	miniTree->Branch("MuonM_Pt", &MuonM_Pt, "MuonM_Pt/F");
	miniTree->Branch("MuonP_Pt", &MuonP_Pt, "MuonP_Pt/F");
	miniTree->Branch("MuonL_Pt", &MuonL_Pt, "MuonL_Pt/F");
	miniTree->Branch("MuonS_Pt", &MuonS_Pt, "MuonS_Pt/F");

	miniTree->Branch("MuonM_Eta", &MuonM_Eta, "MuonM_Eta/F");
	miniTree->Branch("MuonP_Eta", &MuonP_Eta, "MuonP_Eta/F");
	miniTree->Branch("MuonL_Eta", &MuonL_Eta, "MuonL_Eta/F");
	miniTree->Branch("MuonS_Eta", &MuonS_Eta, "MuonS_Eta/F");

	miniTree->Branch("MuonM_Phi", &MuonM_Phi, "MuonM_Phi/F");
	miniTree->Branch("MuonP_Phi", &MuonP_Phi, "MuonP_Phi/F");
	miniTree->Branch("MuonL_Phi", &MuonL_Phi, "MuonL_Phi/F");
	miniTree->Branch("MuonS_Phi", &MuonS_Phi, "MuonS_Phi/F");

	miniTree->Branch("MuonL_Charge", &MuonL_Charge, "MuonL_Charge/I");
	miniTree->Branch("MuonS_Charge", &MuonS_Charge, "MuonS_Charge/I");

	miniTree->Branch("MuonM_isoR03_emEt", &MuonM_isoR03_emEt, "MuonM_isoR03_emEt/F");
	miniTree->Branch("MuonP_isoR03_emEt", &MuonP_isoR03_emEt, "MuonP_isoR03_emEt/F");
	miniTree->Branch("MuonL_isoR03_emEt", &MuonL_isoR03_emEt, "MuonL_isoR03_emEt/F");
	miniTree->Branch("MuonS_isoR03_emEt", &MuonS_isoR03_emEt, "MuonS_isoR03_emEt/F");

	miniTree->Branch("MuonM_isoR03_hadEt", &MuonM_isoR03_hadEt, "MuonM_isoR03_hadEt/F");
	miniTree->Branch("MuonP_isoR03_hadEt", &MuonP_isoR03_hadEt, "MuonP_isoR03_hadEt/F");
	miniTree->Branch("MuonL_isoR03_hadEt", &MuonL_isoR03_hadEt, "MuonL_isoR03_hadEt/F");
	miniTree->Branch("MuonS_isoR03_hadEt", &MuonS_isoR03_hadEt, "MuonS_isoR03_hadEt/F");

	miniTree->Branch("MuonM_isoR03_hoEt", &MuonM_isoR03_hoEt, "MuonM_isoR03_hoEt/F");
	miniTree->Branch("MuonP_isoR03_hoEt", &MuonP_isoR03_hoEt, "MuonP_isoR03_hoEt/F");
	miniTree->Branch("MuonL_isoR03_hoEt", &MuonL_isoR03_hoEt, "MuonL_isoR03_hoEt/F");
	miniTree->Branch("MuonS_isoR03_hoEt", &MuonS_isoR03_hoEt, "MuonS_isoR03_hoEt/F");

	miniTree->Branch("MuonM_isoR03_nJets", &MuonM_isoR03_nJets, "MuonM_isoR03_nJets/F");
	miniTree->Branch("MuonP_isoR03_nJets", &MuonP_isoR03_nJets, "MuonP_isoR03_nJets/F");
	miniTree->Branch("MuonL_isoR03_nJets", &MuonL_isoR03_nJets, "MuonL_isoR03_nJets/F");
	miniTree->Branch("MuonS_isoR03_nJets", &MuonS_isoR03_nJets, "MuonS_isoR03_nJets/F");

	miniTree->Branch("MuonM_isoR03_nTracks", &MuonM_isoR03_nTracks, "MuonM_isoR03_nTracks/F");
	miniTree->Branch("MuonP_isoR03_nTracks", &MuonP_isoR03_nTracks, "MuonP_isoR03_nTracks/F");
	miniTree->Branch("MuonL_isoR03_nTracks", &MuonL_isoR03_nTracks, "MuonL_isoR03_nTracks/F");
	miniTree->Branch("MuonS_isoR03_nTracks", &MuonS_isoR03_nTracks, "MuonS_isoR03_nTracks/F");

	miniTree->Branch("MuonM_isoR03_sumPt", &MuonM_isoR03_sumPt, "MuonM_isoR03_sumPt/F");
	miniTree->Branch("MuonP_isoR03_sumPt", &MuonP_isoR03_sumPt, "MuonP_isoR03_sumPt/F");
	miniTree->Branch("MuonL_isoR03_sumPt", &MuonL_isoR03_sumPt, "MuonL_isoR03_sumPt/F");
	miniTree->Branch("MuonS_isoR03_sumPt", &MuonS_isoR03_sumPt, "MuonS_isoR03_sumPt/F");

	miniTree->Branch("MuonM_isoR05_emEt", &MuonM_isoR05_emEt, "MuonM_isoR05_emEt/F");
	miniTree->Branch("MuonP_isoR05_emEt", &MuonP_isoR05_emEt, "MuonP_isoR05_emEt/F");
	miniTree->Branch("MuonL_isoR05_emEt", &MuonL_isoR05_emEt, "MuonL_isoR05_emEt/F");
	miniTree->Branch("MuonS_isoR05_emEt", &MuonS_isoR05_emEt, "MuonS_isoR05_emEt/F");

	miniTree->Branch("MuonM_isoR05_hadEt", &MuonM_isoR05_hadEt, "MuonM_isoR05_hadEt/F");
	miniTree->Branch("MuonP_isoR05_hadEt", &MuonP_isoR05_hadEt, "MuonP_isoR05_hadEt/F");
	miniTree->Branch("MuonL_isoR05_hadEt", &MuonL_isoR05_hadEt, "MuonL_isoR05_hadEt/F");
	miniTree->Branch("MuonS_isoR05_hadEt", &MuonS_isoR05_hadEt, "MuonS_isoR05_hadEt/F");

	miniTree->Branch("MuonM_isoR05_hoEt", &MuonM_isoR05_hoEt, "MuonM_isoR05_hoEt/F");
	miniTree->Branch("MuonP_isoR05_hoEt", &MuonP_isoR05_hoEt, "MuonP_isoR05_hoEt/F");
	miniTree->Branch("MuonL_isoR05_hoEt", &MuonL_isoR05_hoEt, "MuonL_isoR05_hoEt/F");
	miniTree->Branch("MuonS_isoR05_hoEt", &MuonS_isoR05_hoEt, "MuonS_isoR05_hoEt/F");

	miniTree->Branch("MuonM_isoR05_nJets", &MuonM_isoR05_nJets, "MuonM_isoR05_nJets/F");
	miniTree->Branch("MuonP_isoR05_nJets", &MuonP_isoR05_nJets, "MuonP_isoR05_nJets/F");
	miniTree->Branch("MuonL_isoR05_nJets", &MuonL_isoR05_nJets, "MuonL_isoR05_nJets/F");
	miniTree->Branch("MuonS_isoR05_nJets", &MuonS_isoR05_nJets, "MuonS_isoR05_nJets/F");

	miniTree->Branch("MuonM_isoR05_nTracks", &MuonM_isoR05_nTracks, "MuonM_isoR05_nTracks/F");
	miniTree->Branch("MuonP_isoR05_nTracks", &MuonP_isoR05_nTracks, "MuonP_isoR05_nTracks/F");
	miniTree->Branch("MuonL_isoR05_nTracks", &MuonL_isoR05_nTracks, "MuonL_isoR05_nTracks/F");
	miniTree->Branch("MuonS_isoR05_nTracks", &MuonS_isoR05_nTracks, "MuonS_isoR05_nTracks/F");

	miniTree->Branch("MuonM_isoR05_sumPt", &MuonM_isoR05_sumPt, "MuonM_isoR05_sumPt/F");
	miniTree->Branch("MuonP_isoR05_sumPt", &MuonP_isoR05_sumPt, "MuonP_isoR05_sumPt/F");
	miniTree->Branch("MuonL_isoR05_sumPt", &MuonL_isoR05_sumPt, "MuonL_isoR05_sumPt/F");
	miniTree->Branch("MuonS_isoR05_sumPt", &MuonS_isoR05_sumPt, "MuonS_isoR05_sumPt/F");

	miniTree->Branch("MuonM_isoR05_nTracks", &MuonM_isoR05_nTracks, "MuonM_isoR05_nTracks/F");
	miniTree->Branch("MuonP_isoR05_nTracks", &MuonP_isoR05_nTracks, "MuonP_isoR05_nTracks/F");
	miniTree->Branch("MuonL_isoR05_nTracks", &MuonL_isoR05_nTracks, "MuonL_isoR05_nTracks/F");
	miniTree->Branch("MuonS_isoR05_nTracks", &MuonS_isoR05_nTracks, "MuonS_isoR05_nTracks/F");

	miniTree->Branch("MuonM_isoR05_sumPt", &MuonM_isoR05_sumPt, "MuonM_isoR05_sumPt/F");
	miniTree->Branch("MuonP_isoR05_sumPt", &MuonP_isoR05_sumPt, "MuonP_isoR05_sumPt/F");
	miniTree->Branch("MuonL_isoR05_sumPt", &MuonL_isoR05_sumPt, "MuonL_isoR05_sumPt/F");
	miniTree->Branch("MuonS_isoR05_sumPt", &MuonS_isoR05_sumPt, "MuonS_isoR05_sumPt/F");

	miniTree->Branch("MuonM_E", &MuonM_E, "MuonM_E/F");
	miniTree->Branch("MuonP_E", &MuonP_E, "MuonP_E/F");
	miniTree->Branch("MuonL_E", &MuonL_E, "MuonL_E/F");
	miniTree->Branch("MuonS_E", &MuonS_E, "MuonS_E/F");

	miniTree->Branch("MuonM_Px", &MuonM_Px, "MuonM_Px/F");
	miniTree->Branch("MuonP_Px", &MuonP_Px, "MuonP_Px/F");
	miniTree->Branch("MuonL_Px", &MuonL_Px, "MuonL_Px/F");
	miniTree->Branch("MuonS_Px", &MuonS_Px, "MuonS_Px/F");

	miniTree->Branch("MuonM_Py", &MuonM_Py, "MuonM_Py/F");
	miniTree->Branch("MuonP_Py", &MuonP_Py, "MuonP_Py/F");
	miniTree->Branch("MuonL_Py", &MuonL_Py, "MuonL_Py/F");
	miniTree->Branch("MuonS_Py", &MuonS_Py, "MuonS_Py/F");

	miniTree->Branch("MuonM_Pz", &MuonM_Pz, "MuonM_Pz/F");
	miniTree->Branch("MuonP_Pz", &MuonP_Pz, "MuonP_Pz/F");
	miniTree->Branch("MuonL_Pz", &MuonL_Pz, "MuonL_Pz/F");
	miniTree->Branch("MuonS_Pz", &MuonS_Pz, "MuonS_Pz/F");


	
	// ____________________________________________
	// mumu information
	// ____________________________________________

	miniTree->Branch("Mmumu", &Mmumu, "Mmumu/F");

	miniTree->Branch("Ptmumu", &Ptmumu, "Ptmumu/F");


	// ____________________________________________
	// MC Truth
	// ___________________________________________

	miniTree->Branch("MuonM_MC_E", &MuonM_MC_E, "MuonM_MC_E/F");
	miniTree->Branch("MuonM_MC_Px", &MuonM_MC_Px, "MuonM_MC_Px/F");
	miniTree->Branch("MuonM_MC_Py", &MuonM_MC_Py, "MuonM_MC_Py/F");
	miniTree->Branch("MuonM_MC_Pz", &MuonM_MC_Pz, "MuonM_MC_Pz/F");
	miniTree->Branch("MuonM_MC_Phi", &MuonM_MC_Phi, "MuonM_MC_Phi/F");
	miniTree->Branch("MuonM_MC_Eta", &MuonM_MC_Eta, "MuonM_MC_Eta/F");
	miniTree->Branch("MuonM_MC_Pt", &MuonM_MC_Pt, "MuonM_MC_Pt/F");
	miniTree->Branch("MuonP_MC_E", &MuonP_MC_E, "MuonP_MC_E/F");
	miniTree->Branch("MuonP_MC_Px", &MuonP_MC_Px, "MuonP_MC_Px/F");
	miniTree->Branch("MuonP_MC_Py", &MuonP_MC_Py, "MuonP_MC_Py/F");
	miniTree->Branch("MuonP_MC_Pz", &MuonP_MC_Pz, "MuonP_MC_Pz/F");
	miniTree->Branch("MuonP_MC_Phi", &MuonP_MC_Phi, "MuonP_MC_Phi/F");
	miniTree->Branch("MuonP_MC_Eta", &MuonP_MC_Eta, "MuonP_MC_Eta/F");
	miniTree->Branch("MuonP_MC_Pt", &MuonP_MC_Pt, "MuonP_MC_Pt/F");
	miniTree->Branch("MuonL_MC_E", &MuonL_MC_E, "MuonL_MC_E/F");
	miniTree->Branch("MuonL_MC_Px", &MuonL_MC_Px, "MuonL_MC_Px/F");
	miniTree->Branch("MuonL_MC_Py", &MuonL_MC_Py, "MuonL_MC_Py/F");
	miniTree->Branch("MuonL_MC_Pz", &MuonL_MC_Pz, "MuonL_MC_Pz/F");
	miniTree->Branch("MuonL_MC_Phi", &MuonL_MC_Phi, "MuonL_MC_Phi/F");
	miniTree->Branch("MuonL_MC_Eta", &MuonL_MC_Eta, "MuonL_MC_Eta/F");
	miniTree->Branch("MuonL_MC_Pt", &MuonL_MC_Pt, "MuonL_MC_Pt/F");
	miniTree->Branch("MuonS_MC_E", &MuonS_MC_E, "MuonS_MC_E/F");
	miniTree->Branch("MuonS_MC_Px", &MuonS_MC_Px, "MuonS_MC_Px/F");
	miniTree->Branch("MuonS_MC_Py", &MuonS_MC_Py, "MuonS_MC_Py/F");
	miniTree->Branch("MuonS_MC_Pz", &MuonS_MC_Pz, "MuonS_MC_Pz/F");
	miniTree->Branch("MuonS_MC_Phi", &MuonS_MC_Phi, "MuonS_MC_Phi/F");
	miniTree->Branch("MuonS_MC_Eta", &MuonS_MC_Eta, "MuonS_MC_Eta/F");
	miniTree->Branch("MuonS_MC_Pt", &MuonS_MC_Pt, "MuonS_MC_Pt/F");
	miniTree->Branch("Mmumu_Muons_MC", &Mmumu_Muons_MC, "Mmumu_Muons_MC/F");
	
	// SETUP PARAMETERS	
	unsigned int NbEvents = (int)inputEventTree->GetEntries();
//	unsigned int NbEvents = 40000;
//	bool powheg = false;
	bool powheg = true;
	bool signal = false;
	bool stew = false;
	bool zjet_veto = (bool)(isZgammaMC == 1 || isZgammaMC == 2);
	cout << "Nb of events : " << NbEvents << endl;
	cout << "Signal is: " << signal <<endl;
	cout << "Stew is: " << stew << endl;
	cout << "ZJet veto is: " << zjet_veto << endl;

	cout << "inputRunTree->GetEntries()= " << inputRunTree->GetEntries() << endl;
	cout << "inputEventTree->GetEntries()= " << inputEventTree->GetEntries() << endl;
	string lastFile = "";

//	double integratedLuminosity = 714.783728;
//	double integratedLuminosity = 1078.19387;
//	double integratedLuminosity = 1420.38;
//	double integratedLuminosity = 872.995;
//	double integratedLuminosity = 2.148*1000.0;
//	double integratedLuminosity = 1.131*1000.0 + 370.915 + 636.440;
//	double integratedLuminosity = 2.138 * 1000.0;
//	double integratedLuminosity = 4009.87400;
//	double XSectionDYToMuMu = 1300.0 * 1.2416;
	double XSectionDYToMuMu = 1626.0;
//	double XSectionTTJets = 94.0;
	double XSectionTTJets = 94.76;
	double XSectionWJetsToLNu = 7899.0;
	double XSectionQCDMu = 349988.0;

//	double InitialNumberDYToMuMu = 2148325.0;
	double InitialNumberDYToMuMu = 29743564.0;
//	double InitialNumberTTJets = 1089625.0;
	double InitialNumberTTJets = 3701947.0;
	double InitialNumberWJetsToLNu = 5413258.0;
	double InitialNumberQCDMu = 8797418.0;

	double minPtHat = -100;
	double maxPtHat = 1000000;
	int verbosity = 0;
	int Nb_events_outside_powheg_cuts = 0;
	int TOTALnbMuonsAfterID[12] = {0};
	int TOTALnbEventsAfterMuonID[12] = {0};
	int TOTALnbDimuonsAfterID[3] = {0};
	int TOTALnbEventsAfterDimuonID[3] = {0};
	int TOTALnbPhotonsAfterID[6] = {0};
	int TOTALnbEventsAfterPhotonID[6] = {0};
	int TOTALnbMuMuGammaAfterID[8] = {0};
	int TOTALnbEventsAfterMuMuGammaID[8] = {0};

	int NbEventsPerJob = NbEvents;
	int NbEventsBegin = 0;
	int NbEventsEnd = NbEvents;
	if( ntotjob == 9999 && ijob != -1 )
	{
		NbEventsPerJob = 200000;
		NbEventsBegin = ijob * NbEventsPerJob;
		NbEventsEnd = min( (ijob + 1)* NbEventsPerJob , (int)NbEvents);
		NbEvents = NbEventsEnd - NbEventsBegin;
	}
	cout << "NbEventsBegin= " << NbEventsBegin << "\tNbEventsEnd= " << NbEventsEnd << "\tNbEventsPerJob= " << NbEventsPerJob << endl;

	// LOOP over events
	for(unsigned int ievt=NbEventsBegin; ievt<NbEventsEnd; ievt++)
//	for(unsigned int ievt=0; ievt<NbEvents; ievt++)
//	for(unsigned int ievt=0; ievt<100; ievt++)
	{
		if(verbosity>4) cout << "analysing event ievt= " << ievt << endl;
		int nprint = (int)((double)NbEventsPerJob/(double)100.0);
		if( (NbEvents >= 100) && (ievt % nprint)==0 ){ cout<< ievt <<" events done over "<<NbEvents<<" ( "<<ceil((double)ievt/(double)NbEventsPerJob*100)<<" \% )"<<endl; }

		iEvent = ievt;
		inputEventTree->GetEvent(ievt);
//		if( lastFile == "" ){
//			lastFile = string(inputEventTree->GetCurrentFile()->GetName());
//			cout << ievt << "\t" << lastFile << endl;
//		}

		// ____________________________________________
		// Event information
		// ____________________________________________
		iEventID = event->eventId();
		iLumiID = event->luminosityBlock();
		iRunID = event->runId();
//		if(iRunID > 166967) continue;
		nVertices = vertices->GetEntries();
		nGenVertices = vertices->GetEntries();
		if( (isZgammaMC == 1) || (isZgammaMC == 2) )
		{
//			nGenVertices = event->nInTimePUVertices();
			nGenVertices = event->inTimePU_NumInteractions();
		}
		isMM = 0;

		weight_Xsection = weight_pileUp = 1.0;

		string sample_in(sample_char);

		if( sample_in.find("DYToMuMu") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionDYToMuMu) / (double)(InitialNumberDYToMuMu)) * (double)integratedLuminosity);
			weight_pileUp = weight_DYToMuMu(nGenVertices+1, lumi_set, pu_set);
		} else if( sample_in.find("QCD") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionQCDMu) / (double)(InitialNumberQCDMu)) * (double)integratedLuminosity);
			weight_pileUp = weight_QCDMu(nGenVertices+1);
		} else if( sample_in.find("TTJets") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionTTJets) / (double)(InitialNumberTTJets)) * (double)integratedLuminosity);
			weight_pileUp = weight_TTJets(nGenVertices+1);
		} else if( sample_in.find("WToMuNu") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionWJetsToLNu) / (double)(InitialNumberWJetsToLNu)) * (double)integratedLuminosity);
			weight_pileUp = weight_WJetsToLNu(nGenVertices+1);
		}
//cout << "nGenVertices= " << nGenVertices << "\tweight_pileUp= " << weight_pileUp << endl;

		// ____________________________________________
		// Muon variables
		// ____________________________________________
		NbMuons = muons->GetEntries();
		Pt_allMuons = Eta_allMuons = Phi_allMuons = Charge_allMuons = -99;
		MuonM_Pt = MuonP_Pt = MuonL_Pt = MuonS_Pt = -99;
		MuonM_Eta = MuonP_Eta = MuonL_Eta = MuonS_Eta = -99;
		MuonM_Phi = MuonP_Phi = MuonL_Phi = MuonS_Phi = -99;
		MuonL_Charge = MuonS_Charge = -99;
		MuonM_isoR03_emEt = MuonP_isoR03_emEt = MuonL_isoR03_emEt = MuonS_isoR03_emEt = -99;
		MuonM_isoR03_hadEt = MuonP_isoR03_hadEt = MuonL_isoR03_hadEt = MuonS_isoR03_hadEt = -99;
		MuonM_isoR03_hoEt = MuonP_isoR03_hoEt = MuonL_isoR03_hoEt = MuonS_isoR03_hoEt = -99;
		MuonM_isoR03_nJets = MuonP_isoR03_nJets = MuonL_isoR03_nJets = MuonS_isoR03_nJets = -99;
		MuonM_isoR03_nTracks = MuonP_isoR03_nTracks = MuonL_isoR03_nTracks = MuonS_isoR03_nTracks = -99;
		MuonM_isoR03_sumPt = MuonP_isoR03_sumPt = MuonL_isoR03_sumPt = MuonS_isoR03_sumPt = -99;
		MuonM_isoR05_emEt = MuonP_isoR05_emEt = MuonL_isoR05_emEt = MuonS_isoR05_emEt = -99;
		MuonM_isoR05_hadEt = MuonP_isoR05_hadEt = MuonL_isoR05_hadEt = MuonS_isoR05_hadEt = -99;
		MuonM_isoR05_hoEt = MuonP_isoR05_hoEt = MuonL_isoR05_hoEt = MuonS_isoR05_hoEt = -99;
		MuonM_isoR05_nJets = MuonP_isoR05_nJets = MuonL_isoR05_nJets = MuonS_isoR05_nJets = -99;
		MuonM_isoR05_nTracks = MuonP_isoR05_nTracks = MuonL_isoR05_nTracks = MuonS_isoR05_nTracks = -99;
		MuonM_isoR05_sumPt = MuonP_isoR05_sumPt = MuonL_isoR05_sumPt = MuonS_isoR05_sumPt = -99;

		MuonM_E = MuonP_E = MuonL_E = MuonS_E = -99;
		MuonM_Px = MuonP_Px = MuonL_Px = MuonS_Px = -99;
		MuonM_Py = MuonP_Py = MuonL_Py = MuonS_Py = -99;
		MuonM_Pz = MuonP_Pz = MuonL_Pz = MuonS_Pz = -99;

		// ____________________________________________
		// mumu information
		// ____________________________________________
		Mmumu = -99.0;
		Ptmumu = -99.0;
		// ____________________________________________
		// MC Truth
		// ___________________________________________
		MuonM_MC_E = MuonM_MC_Px = MuonM_MC_Py = MuonM_MC_Pz = MuonM_MC_Phi = MuonM_MC_Eta = MuonM_MC_Pt = -99.0;
		MuonP_MC_E = MuonP_MC_Px = MuonP_MC_Py = MuonP_MC_Pz = MuonP_MC_Phi = MuonP_MC_Eta = MuonP_MC_Pt = -99.0;
		MuonL_MC_E = MuonL_MC_Px = MuonL_MC_Py = MuonL_MC_Pz = MuonL_MC_Phi = MuonL_MC_Eta = MuonL_MC_Pt = -99.0;
		MuonS_MC_E = MuonS_MC_Px = MuonS_MC_Py = MuonS_MC_Pz = MuonS_MC_Phi = MuonS_MC_Eta = MuonS_MC_Pt = -99.0;
		Mmumu_Muons_MC = -99.0;
		// ____________________________________________
		// END OF INITIALIZATION
		// ____________________________________________

// ********************************************************************************************************************
// ***** MC-TRUTH SIGNAL MATCHING *****
// ********************************************************************************************************************


// ********************************************************************************************************************
// ***** HLT REQUIREMENT *****
// ********************************************************************************************************************

		// HLT information		

//		if(doHLT){
//			if( ievt==0 ){ inputRunTree->GetEvent(ievt); NumWantedHLTnames = InitializeHLTinfo(inputRunTree, runInfos, event->nHLTPaths(), ListWantedHLTnames, 1);	}
//			if ( string(inputEventTree->GetCurrentFile()->GetName()) != lastFile ){
//				inputRunTree->GetEntry(inputEventTree->GetTreeNumber());
//		 		lastFile = string(inputEventTree->GetCurrentFile()->GetName());
//		 		cout << ievt << "\t" << lastFile << endl;
//		 	 	NumWantedHLTnames = InitializeHLTinfo(inputRunTree, runInfos, event->nHLTPaths(), ListWantedHLTnames, 1);
//		 	}
//			doHLTInfo(event, runInfos, NumWantedHLTnames, 1, &Muon_eventPassHLT_Mu11);
//		}
//
//
//		if (!((event->ptHat()>=minPtHat)&&(event->ptHat()<maxPtHat)))
//		{
//			cerr << "CUT: event " << ievt << " ( " << iRunID << " , " << iLumiID << " , " << iEventID << " )" << " CUT for pthat filtering" << endl;
//			continue;
//		}
//		isAfterCutPthatFilter = 1;
//		nAfterCutPthatFilter++;

// ********************************************************************************************************************
// ***** MUON OBJECT SELECTION *****
// ********************************************************************************************************************

		if(verbosity>1) cout << "Start loop over muon objects" << endl;
		// Cleaning: muon quality
		vector<int> muonIsNotCommissioned;
		muonIsNotCommissioned.clear();
		vector<int> muonIdentified;
		muonIdentified.clear();
		vector<double> muonCorrectedPt;
		muonCorrectedPt.clear();
		int nbMuonsAfterID[12] = {0};
		
		if(verbosity>0) cerr << "\t\tThere is " << NbMuons << " muons in the muon collection" << endl;

		nbMuonsAfterID[0] = NbMuons;
		TOTALnbMuonsAfterID[0] += NbMuons;
		for(int imuon=0 ; imuon<NbMuons ; imuon++){
			TRootMuon *mymuon;
			mymuon = (TRootMuon*) muons->At(imuon);

			if(! (mymuon->isGlobalMuon()) ){
				muonIsNotCommissioned.push_back(1);
			if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because not global muon" << endl;
				continue;
			}
			nbMuonsAfterID[1]++;
			TOTALnbMuonsAfterID[1]++;

			if(! (mymuon->normalizedGlobalChi2()<10.0) ){// chi2/ndof of the global muon fit < 10
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because chi2/ndof of the global muon fit < 10 (" << mymuon->normalizedGlobalChi2() << ")" << endl;
				continue;
			}
			nbMuonsAfterID[2]++;
			TOTALnbMuonsAfterID[2]++;

			if(! (mymuon->numberOfValidGlobalHits()>0) ){// number of valid muon hits matched to the global fit
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because number of valid muon hits matched to the global fit too low (" << mymuon->numberOfValidGlobalHits() << ")" << endl;
				continue;
			}
			nbMuonsAfterID[3]++;
			TOTALnbMuonsAfterID[3]++;

			if(! (mymuon->isTrackerMuon()) ){// The muon must be identified as Tracker Muon.
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because not tracker muon" << endl;
				continue;
			}
			nbMuonsAfterID[4]++;
			TOTALnbMuonsAfterID[4]++;

			if(! (mymuon->numberOfMatches()>1) ){// number of muon stations with matched segments (global track: out-in fit)
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because number of muon stations with matched segments (global track: out-in fit) too low" << endl;
				continue;
			}
			nbMuonsAfterID[5]++;
			TOTALnbMuonsAfterID[5]++;

			if(! (mymuon->numberOfValidTrackerHits()>10) ){// number of tracker (pixels + strips) hits
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because number of tracker (pixels + strips) hits" << endl;
				continue;
			}
			nbMuonsAfterID[6]++;
			TOTALnbMuonsAfterID[6]++;

			if(! (mymuon->numberOfValidPixelHits()>0) ){// number of pixel hits
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because number of pixel hits" << endl;
				continue;
			}
			nbMuonsAfterID[7]++;
			TOTALnbMuonsAfterID[7]++;

			if(! (fabs(mymuon->GlobaldB())<0.2) )
			{// inner track transverse impact parameter w.r.t the beam spot |d_xy|
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because inner track transverse impact parameter w.r.t the beam spot |d_xy|" << endl;
				continue;
			}
			nbMuonsAfterID[8]++;
			TOTALnbMuonsAfterID[8]++;

			if(! (mymuon->isoR03_sumPt()<3.0) ){// sum of pT of tracks with pT >1.5 within a cone of DR < 0.3 around the muon direction, vetoing a cone of 0.015 around that direction
//			if(! ( (double)(mymuon->isoR03_sumPt() + mymuon->isoR03_emEt() + mymuon->isoR03_hadEt())/(double)(mymuon->Pt()) < 0.10 ) ){// combined isolation as described in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because sum of pT of tracks with pT >1.5 within a cone of DR < 0.3 around the muon direction, vetoing a cone of 0.015 around that direction" << endl;
				continue;
			}
			nbMuonsAfterID[9]++;
			TOTALnbMuonsAfterID[9]++;

		double corrected_Pt = mymuon->Pt();
		rochcor *rmcor = new rochcor();
		if( applyMuonScaleCorrection > 0 )
		{
			if( applyMuonScaleCorrection == 2 )
			{ // Apply SIDRA
				corrected_Pt = applySidra(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi(), generator);
			}
			if( applyMuonScaleCorrection == 1 )
			{ // Apply MuScleFit
				corrected_Pt = applyMuScleFit(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi());
			}
			if( applyMuonScaleCorrection == 21 )
			{ // Apply SIDRA then MuScleFit
				corrected_Pt = applySidra(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi(), generator);
				corrected_Pt = applyMuScleFit(corrected_Pt, mymuon->charge(), mymuon->Eta(), mymuon->Phi());
			}
			if( applyMuonScaleCorrection == 12 )
			{ // Apply SIDRA then MuScleFit
				corrected_Pt = applyMuScleFit(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi());
				corrected_Pt = applySidra(corrected_Pt, mymuon->charge(), mymuon->Eta(), mymuon->Phi(), generator);
			}
			if( applyMuonScaleCorrection == 3 )
			{	// Apply Rochester corrections
				TLorentzVector muonRochester(mymuon->Px(), mymuon->Py(), mymuon->Pz(), mymuon->E());
				TLorentzVector muonRochesterDummy(mymuon->Px(), mymuon->Py(), mymuon->Pz(), mymuon->E());
				// moption = 1	: recommended by the authors (better match in Z mass profile vs. eta/phi between the reconstructed and generated Z mass)
				// sysdev = 0 : no systematics yet
								if( isZgammaMC > 0 ) // If sample is MC
								{
										if( (lumi_set == "2011A") || (lumi_set == "2011A_rereco") ) runopt = 0;
										if( (lumi_set == "2011B") || (lumi_set == "2011B_rereco") ) runopt = 1;
										if( lumi_set == "2011A" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719;
										if( lumi_set == "2011A_rereco" ) integratedLuminosity = 2.221*1000.0;
										if( lumi_set == "2011B" ) integratedLuminosity = 2.714*1000.0;
										if( lumi_set == "2011" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0;
										if( lumi_set == "2011_rereco" ) integratedLuminosity = 2.221*1000.0 +	 2.714*1000.0;

										// event where to switch from Run2011A to Run2011B correction
										// We run over ntot=	DYToMuMu events
										// We run on average on avEventsPerJob= (ntot) / (ntotjob)= events per job
										// RunA represents Astat= (2.221) / (2.221 + 2.714) = 45.01 % of the total stat
										// RunB represents Bstat= (2.714) / (2.221 + 2.714) = 54.99 % of the total stat
										// The last event of 'RunA' is lastA= Astat*ntot
										// Which should happen in job ijobLastA= (int)(lastA)/(int)(avEventsPerJob)
										// and will be the event iEventLastA= (int)(lastA)%(int)(avEventsPerJob)
										// ***********************************
										// temp version
										// ***********************************
										if( (lumi_set == "2011" ) )
										{
												int ntot= 8950877;
												double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
												double Astat= (double)(215.552 + 951.716 + 389.876 + 706.719) / (double)(215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0);
												double Bstat= (double)(2.714*1000.0) / (double)(215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0);
												double lastA= Astat*ntot;
												int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
												int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
												if( ijob == (ijobLastA -1) )
												{
														if( ievt >= iEventLastA )
														{
																runopt = 1;
														} else runopt = 0;
												} else if( ijob >= ijobLastA ) runopt = 1;
										}
										if( (lumi_set == "2011_rereco" ) )
										{
												int ntot= 8950877;
												double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
												double Astat= (double)(2.221) / (double)(2.221 + 2.714);
												double Bstat= (double)(2.714) / (double)(2.221 + 2.714);
												double lastA= Astat*ntot;
												int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
												int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
												if( ijob == (ijobLastA -1) )
												{
														if( ievt >= iEventLastA )
														{
																runopt = 1;
														} else runopt = 0;
												} else if( ijob >= ijobLastA ) runopt = 1;
										}
										if( verbosity > 4) cerr << "### ievt= " << ievt << "\tijob= " << ijob << "\trunopt= " << runopt << endl;
										rmcor->momcor_mc(muonRochester, mymuon->charge(), sysdev, runopt);
//											if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
//											else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
								} else {
										rmcor->momcor_data(muonRochester, mymuon->charge(), sysdev, runopt);
//											if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev, runopt);
//											else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev, runopt);
								}

				/*if( isZgammaMC > 0 ) // If sample is MC
				{
					if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
					else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
				} else {
					if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev);
							else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev);
				}*/
				corrected_Pt = muonRochester.Pt();
			}
			if( applyMuonScaleCorrection == 31 )
			{	// Apply Rochester corrections on top of MuscleFit
				corrected_Pt = applyMuScleFit(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi());
				double corrected_Px_ = mymuon->Px() * (double)(corrected_Pt) / (double)(mymuon->Pt());
				double corrected_Py_ = mymuon->Py() * (double)(corrected_Pt) / (double)(mymuon->Pt());
				double corrected_Pz_ = mymuon->Pz();
				double m_mu = 105.658367e-3;
				double corrected_E_ = sqrt( m_mu * m_mu + (corrected_Pz_ * corrected_Pz_ + corrected_Pt * corrected_Pt) );
				TLorentzVector muonRochester(corrected_Px_, corrected_Py_, corrected_Pz_, corrected_E_);
				TLorentzVector muonRochesterDummy(corrected_Px_, corrected_Py_, corrected_Pz_, corrected_E_);

				// moption = 1	: recommended by the authors (better match in Z mass profile vs. eta/phi between the reconstructed and generated Z mass)
				// sysdev = 0 : no systematics yet
								if( isZgammaMC > 0 ) // If sample is MC
								{
										if( (lumi_set == "2011A") || (lumi_set == "2011A_rereco") ) runopt = 0;
										if( (lumi_set == "2011B") || (lumi_set == "2011B_rereco") ) runopt = 1;
										if( lumi_set == "2011A" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719;
										if( lumi_set == "2011A_rereco" ) integratedLuminosity = 2.221*1000.0;
										if( lumi_set == "2011B" ) integratedLuminosity = 2.714*1000.0;
										if( lumi_set == "2011" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0;
										if( lumi_set == "2011_rereco" ) integratedLuminosity = 2.221*1000.0 +	 2.714*1000.0;

										// event where to switch from Run2011A to Run2011B correction
										// We run over ntot=	DYToMuMu events
										// We run on average on avEventsPerJob= (ntot) / (ntotjob)= events per job
										// RunA represents Astat= (2.221) / (2.221 + 2.714) = 45.01 % of the total stat
										// RunB represents Bstat= (2.714) / (2.221 + 2.714) = 54.99 % of the total stat
										// The last event of 'RunA' is lastA= Astat*ntot
										// Which should happen in job ijobLastA= (int)(lastA)/(int)(avEventsPerJob)
										// and will be the event iEventLastA= (int)(lastA)%(int)(avEventsPerJob)
										// ***********************************
										// temp version
										// ***********************************
										if( (lumi_set == "2011" ) )
										{
												int ntot= 8950877;
												double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
												double Astat= (double)(215.552 + 951.716 + 389.876 + 706.719) / (double)(215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0);
												double Bstat= (double)(2.714*1000.0) / (double)(215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0);
												double lastA= Astat*ntot;
												int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
												int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
												if( ijob == (ijobLastA -1) )
												{
														if( ievt >= iEventLastA )
														{
																runopt = 1;
														} else runopt = 0;
												} else if( ijob >= ijobLastA ) runopt = 1;
										}
										if( (lumi_set == "2011_rereco" ) )
										{
												int ntot= 8950877;
												double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
												double Astat= (double)(2.221) / (double)(2.221 + 2.714);
												double Bstat= (double)(2.714) / (double)(2.221 + 2.714);
												double lastA= Astat*ntot;
												int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
												int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
												if( ijob == (ijobLastA -1) )
												{
														if( ievt >= iEventLastA )
														{
																runopt = 1;
														} else runopt = 0;
												} else if( ijob >= ijobLastA ) runopt = 1;
										}
										if( verbosity > 4) cerr << "### ievt= " << ievt << "\tijob= " << ijob << "\trunopt= " << runopt << endl;
										rmcor->momcor_mc(muonRochester, mymuon->charge(), sysdev, runopt);
//											if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
//											else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
								} else {
										rmcor->momcor_data(muonRochester, mymuon->charge(), sysdev, runopt);
//											if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev, runopt);
//											else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev, runopt);
								}
/*				if( isZgammaMC > 0 ) // If sample is MC
				{
					if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
					else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
				} else {
					if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev);
					else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev);
				}*/
				corrected_Pt = muonRochester.Pt();
			}
			if( applyMuonScaleCorrection == 32 )
			{	// Apply Rochester corrections on top of Sidra
				corrected_Pt = applySidra(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi(), generator);
				double corrected_Px_ = mymuon->Px() * (double)(corrected_Pt) / (double)(mymuon->Pt());
				double corrected_Py_ = mymuon->Py() * (double)(corrected_Pt) / (double)(mymuon->Pt());
				double corrected_Pz_ = mymuon->Pz();
				double m_mu = 105.658367e-3;
				double corrected_E_ = sqrt( m_mu * m_mu + (corrected_Pz_ * corrected_Pz_ + corrected_Pt * corrected_Pt) );
				TLorentzVector muonRochester(corrected_Px_, corrected_Py_, corrected_Pz_, corrected_E_);
				TLorentzVector muonRochesterDummy(corrected_Px_, corrected_Py_, corrected_Pz_, corrected_E_);

				// moption = 1	: recommended by the authors (better match in Z mass profile vs. eta/phi between the reconstructed and generated Z mass)
				// sysdev = 0 : no systematics yet
								if( isZgammaMC > 0 ) // If sample is MC
								{
										if( (lumi_set == "2011A") || (lumi_set == "2011A_rereco") ) runopt = 0;
										if( (lumi_set == "2011B") || (lumi_set == "2011B_rereco") ) runopt = 1;
										if( lumi_set == "2011A" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719;
										if( lumi_set == "2011A_rereco" ) integratedLuminosity = 2.221*1000.0;
										if( lumi_set == "2011B" ) integratedLuminosity = 2.714*1000.0;
										if( lumi_set == "2011" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0;
										if( lumi_set == "2011_rereco" ) integratedLuminosity = 2.221*1000.0 +	 2.714*1000.0;

										// event where to switch from Run2011A to Run2011B correction
										// We run over ntot=	DYToMuMu events
										// We run on average on avEventsPerJob= (ntot) / (ntotjob)= events per job
										// RunA represents Astat= (2.221) / (2.221 + 2.714) = 45.01 % of the total stat
										// RunB represents Bstat= (2.714) / (2.221 + 2.714) = 54.99 % of the total stat
										// The last event of 'RunA' is lastA= Astat*ntot
										// Which should happen in job ijobLastA= (int)(lastA)/(int)(avEventsPerJob)
										// and will be the event iEventLastA= (int)(lastA)%(int)(avEventsPerJob)
										// ***********************************
										// temp version
										// ***********************************
										if( (lumi_set == "2011" ) )
										{
												int ntot= 8950877;
												double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
												double Astat= (double)(215.552 + 951.716 + 389.876 + 706.719) / (double)(215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0);
												double Bstat= (double)(2.714*1000.0) / (double)(215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0);
												double lastA= Astat*ntot;
												int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
												int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
												if( ijob == (ijobLastA -1) )
												{
														if( ievt >= iEventLastA )
														{
																runopt = 1;
														} else runopt = 0;
												} else if( ijob >= ijobLastA ) runopt = 1;
										}
										if( (lumi_set == "2011_rereco" ) )
										{
												int ntot= 8950877;
												double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
												double Astat= (double)(2.221) / (double)(2.221 + 2.714);
												double Bstat= (double)(2.714) / (double)(2.221 + 2.714);
												double lastA= Astat*ntot;
												int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
												int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
												if( ijob == (ijobLastA -1) )
												{
														if( ievt >= iEventLastA )
														{
																runopt = 1;
														} else runopt = 0;
												} else if( ijob >= ijobLastA ) runopt = 1;
										}
										if( verbosity > 4) cerr << "### ievt= " << ievt << "\tijob= " << ijob << "\trunopt= " << runopt << endl;
										rmcor->momcor_mc(muonRochester, mymuon->charge(), sysdev, runopt);
//											if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
//											else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
								} else {
										rmcor->momcor_data(muonRochester, mymuon->charge(), sysdev, runopt);
//											if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev, runopt);
//											else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev, runopt);
								}
/*				if( isZgammaMC > 0 ) // If sample is MC
				{
					if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
					else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
				} else {
					if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev);
					else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev);
				}*/
				corrected_Pt = muonRochester.Pt();
			}
			if( applyMuonScaleCorrection == 99 )
			{ // Dummy muon correction to check code: it's NOT supposed to do anything
				corrected_Pt = mymuon->Pt();
			}
			// Sidra makes MC look like data
//			if( isZgammaMC > 0) corrected_Pt = applySidra(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi(), generator);
			// MuScleFit correct data absolute scale
//			corrected_Pt = applyMuScleFit(corrected_Pt, mymuon->charge(), mymuon->Eta(), mymuon->Phi());
		}
		double corrected_Px = applyMuonScaleCorrection > 0 ? mymuon->Px() * (double)(corrected_Pt) / (double)(mymuon->Pt()) : mymuon->Px();
		double corrected_Py = applyMuonScaleCorrection > 0 ? mymuon->Py() * (double)(corrected_Pt) / (double)(mymuon->Pt()) : mymuon->Py();
		double corrected_Pz = mymuon->Pz();
		double m_mu = 105.658367e-3;
		double corrected_E = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz * corrected_Pz + corrected_Pt * corrected_Pt) ) : mymuon->E();
//		cout << "px: " << mymuon->Px() << "\t" << corrected_Px << endl;
//		cout << "py: " << mymuon->Py() << "\t" << corrected_Py << endl;
//		cout << "pz: " << mymuon->Pz() << "\t" << corrected_Pz << endl;
//		cout << "E : " <<	mymuon->E() << "\t" << corrected_E << endl;
//		double corrected_E = mymuon->E();
		TLorentzVector correctedMuon(corrected_Px, corrected_Py, corrected_Pz, corrected_E);
//		cout << "###########" << endl;
//		cout << "mymuon->M()= " << mymuon->M() << endl;
//		cout << "correctedMuon.M()= " << correctedMuon.M() << endl;

			if(! (correctedMuon.Pt()>10.0) ){// transverse momentum
//			if(! (mymuon->Pt()>10.0) ){// transverse momentum
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because transverse momentum" << endl;
				continue;
			}
			nbMuonsAfterID[10]++;
			TOTALnbMuonsAfterID[10]++;

			if(! (fabs(mymuon->Eta())<2.4) ){// |eta_muon|< 2.1
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because high eta (" << mymuon->Eta() << ")" << endl;
				continue;
			}
			nbMuonsAfterID[11]++;
			TOTALnbMuonsAfterID[11]++;

//			if(! (mymuon->) ){// 
//				muonIsNotCommissioned.push_back(1);
//				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because" << endl;
//				continue;
//			}

			if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " accepted" << endl;
			muonIsNotCommissioned.push_back(0);
			muonIdentified.push_back(imuon);
			muonCorrectedPt.push_back(corrected_Pt);
////			mymuon->Clear();
		}
		unsigned int NbMuonsIdentified = muonIdentified.size();
		
		// Increasing counter
		for(int i = 0; i < 12 ; i++)
		{
			if(nbMuonsAfterID[i] >= 2){ TOTALnbEventsAfterMuonID[i]++;}
		}

		if(! (nbMuonsAfterID[11] >=2) )// Not enough dimuon candidates, skip the event
		{
				continue;
		}

// ********************************************************************************************************************
// ***** DIMUON OBJECT SELECTION *****
// ********************************************************************************************************************


		if(verbosity>1) cout << "Filling dimuon pair holder" << endl;


		// Making dimuon pairs holder
		int numberOfDimuons[3] = {0};
		numberOfDimuons[0] = factorial(nbMuonsAfterID[11] -1);

		if(! (numberOfDimuons[0] >= 1) )// Not enough dimuon candidates, skip the event. This cut is redundant with the previous one, should do nothing
		{
			continue;

		}

		TOTALnbDimuonsAfterID[0] += numberOfDimuons[0];
		TOTALnbEventsAfterDimuonID[0] += 1;

		pair <int, int> IDofMuons[3][numberOfDimuons[0]];
		pair <double, double> PTofMuons[3][numberOfDimuons[0]];

		if(verbosity>2) cout << "initializing dimuon pair object" << endl;
		// Initializing pair object
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < numberOfDimuons[0]; j++) IDofMuons[i][j] = make_pair(0.0, 0.0);
			for(int j = 0; j < numberOfDimuons[0]; j++) PTofMuons[i][j] = make_pair(0.0, 0.0);
		}

		if(verbosity>2) cout << "Filling pair object for dimuon pairs composed of ID'ed muons" << endl;
		// Filling pair object for dimuon pairs composed of ID'ed muons
//		for(int i_dimuons = 0; i_dimuons < numberOfDimuons[0]; i_dimuons++)
//		{
			int i_dimuons_ = 0;
			for(int muon_i = 0; muon_i < nbMuonsAfterID[11] ; muon_i++)
			{
				for(int muon_j = muon_i +1; muon_j < nbMuonsAfterID[11]; muon_j++)
				{
					IDofMuons[0][i_dimuons_] = make_pair(muonIdentified[muon_i], muonIdentified[muon_j]);
					PTofMuons[0][i_dimuons_] = make_pair(muonCorrectedPt[muon_i], muonCorrectedPt[muon_j]);
					i_dimuons_++;
				}
			}
//		}

// --------------------------------------------------------------------------------------------------------------------
// ----- dimuon candidate of opposite charge -----
// --------------------------------------------------------------------------------------------------------------------

		if(verbosity>2) cout << "loop over possible muon pairs, fill dimuons candidates with opposite charge" << endl;	
		// loop over possible muon pairs, fill dimuons candidates with opposite charge
		for(int i_dimuons = 0; i_dimuons < numberOfDimuons[0]; i_dimuons++)
		{
			TRootMuon *Muon1 = (TRootMuon*) muons->At(IDofMuons[0][i_dimuons].first);
			TRootMuon *Muon2 = (TRootMuon*) muons->At(IDofMuons[0][i_dimuons].second);
			double chargeproduct = ( Muon1->charge() ) * (Muon2->charge() );
			if( chargeproduct < 0.0 )
			{
				IDofMuons[1][numberOfDimuons[1]] = make_pair(IDofMuons[0][i_dimuons].first, IDofMuons[0][i_dimuons].second);
				PTofMuons[1][numberOfDimuons[1]] = make_pair(PTofMuons[0][i_dimuons].first, PTofMuons[0][i_dimuons].second);
				numberOfDimuons[1] += 1;
			}
////			Muon1->Clear();
////			Muon2->Clear();
		}
		if(! (numberOfDimuons[1] >= 1) ) continue; // Not enough dimuon candidates, skip the event
		TOTALnbDimuonsAfterID[1] += numberOfDimuons[1];
		TOTALnbEventsAfterDimuonID[1] += 1;

// --------------------------------------------------------------------------------------------------------------------
// ----- dimuon candidate of correct invariant mass -----
// --------------------------------------------------------------------------------------------------------------------

		if(verbosity>2) cout << "loop over possible muon pairs, fill dimuons candidates with valid invariant mass" << endl;
		// loop over possible muon pairs, fill dimuons candidates with valid invariant mass
		for(int i_dimuons = 0; i_dimuons < numberOfDimuons[1]; i_dimuons++)
		{
			TRootMuon *Muon1 = (TRootMuon*) muons->At(IDofMuons[1][i_dimuons].first);
			TRootMuon *Muon2 = (TRootMuon*) muons->At(IDofMuons[1][i_dimuons].second);

		double corrected_Pt1 = Muon1->Pt();
		double corrected_Pt2 = Muon2->Pt();
		if( applyMuonScaleCorrection > 0 )
		{
			corrected_Pt1 = PTofMuons[1][i_dimuons].first;
			corrected_Pt2 = PTofMuons[1][i_dimuons].second;
			// Sidra makes MC look like data
//			if( isZgammaMC > 0) corrected_Pt1 = applySidra(Muon1->Pt(), Muon1->charge(), Muon1->Eta(), Muon1->Phi(), generator);
//			if( isZgammaMC > 0) corrected_Pt2 = applySidra(Muon2->Pt(), Muon2->charge(), Muon2->Eta(), Muon2->Phi(), generator);
			// MuScleFit correct data absolute scale
//			corrected_Pt1 = applyMuScleFit(corrected_Pt1, Muon1->charge(), Muon1->Eta(), Muon1->Phi());
//			corrected_Pt2 = applyMuScleFit(corrected_Pt2, Muon2->charge(), Muon2->Eta(), Muon2->Phi());
		}
		double corrected_Pz1 = Muon1->Pz();
		double corrected_Pz2 = Muon2->Pz();
		double corrected_Px1 = applyMuonScaleCorrection > 0 ? Muon1->Px() * (double)(corrected_Pt1) / (double)(Muon1->Pt()) : Muon1->Px();
		double corrected_Px2 = applyMuonScaleCorrection > 0 ? Muon2->Px() * (double)(corrected_Pt2) / (double)(Muon2->Pt()) : Muon2->Px();
		double corrected_Py1 = applyMuonScaleCorrection > 0 ? Muon1->Py() * (double)(corrected_Pt1) / (double)(Muon1->Pt()) : Muon1->Py();
		double corrected_Py2 = applyMuonScaleCorrection > 0 ? Muon2->Py() * (double)(corrected_Pt2) / (double)(Muon2->Pt()) : Muon2->Py();
		double m_mu = 105.658367e-3;
//		double corrected_E1 = Muon1->E();
//		double corrected_E2 = Muon2->E();
		double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1)) : Muon1->E();
		double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2)) : Muon2->E();
/*
cout << "######" << endl;
		cout << "px1: " << Muon1->Px() << "\t" << corrected_Px1 << endl;
		cout << "py1: " << Muon1->Py() << "\t" << corrected_Py1 << endl;
		cout << "pz1: " << Muon1->Pz() << "\t" << corrected_Pz1 << endl;
		cout << "E1 : " <<	Muon1->E() << "\t" << corrected_E1 << endl;
		cout << "px2: " << Muon2->Px() << "\t" << corrected_Px2 << endl;
		cout << "py2: " << Muon2->Py() << "\t" << corrected_Py2 << endl;
		cout << "pz2: " << Muon2->Pz() << "\t" << corrected_Pz2 << endl;
		cout << "E2 : " <<	Muon2->E() << "\t" << corrected_E2 << endl;
*/
		TLorentzVector correctedMuon1(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
		TLorentzVector correctedMuon2(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);





			TLorentzVector mumu;
//			mumu = (*Muon1) + (*Muon2);
			mumu = (correctedMuon1) + (correctedMuon2);
			if( (low_m_mumu < mumu.M()) && (mumu.M() < high_m_mumu) )
			{
				IDofMuons[2][numberOfDimuons[2]] = make_pair(IDofMuons[1][i_dimuons].first, IDofMuons[1][i_dimuons].second);
				PTofMuons[2][numberOfDimuons[2]] = make_pair(PTofMuons[1][i_dimuons].first, PTofMuons[1][i_dimuons].second);
				numberOfDimuons[2] += 1;
			}
////			Muon1->Clear();
////			Muon2->Clear();
		}

		if(! (numberOfDimuons[2] >= 1) )// Not enough dimuon candidates, skip the event
		{
			 continue;
		}

		TOTALnbDimuonsAfterID[2] += numberOfDimuons[2];
		TOTALnbEventsAfterDimuonID[2] += 1;

		isMM = 1;
		for(int idimuon = 0 ; idimuon < numberOfDimuons[2]; idimuon++)
		{
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			mymuon1 = (TRootMuon*) muons->At( IDofMuons[2][idimuon].first );
			mymuon2 = (TRootMuon*) muons->At( IDofMuons[2][idimuon].second );
			
			TLorentzVector mumu;
			mumu = (*mymuon1) + (*mymuon2);
			Ptmumu = mumu.Pt();
			FillMM(mymuon1, mymuon2, doMC, applyMuonScaleCorrection, mcParticles, PTofMuons[2][idimuon].first, PTofMuons[2][idimuon].second);
			miniTree->Fill();
////			mumu.Clear();
////			mymuon1->Clear();
////			mymuon2->Clear();
		}


	} // fin boucle sur evts LOOP

	 	cout << "Nb_events_outside_powheg_cuts= " << Nb_events_outside_powheg_cuts << endl << endl;
		for(int i = 0; i < 12 ; i++)
		{
			cout << "TOTALnbMuonsAfterID["<<i<<"]=\t" << TOTALnbMuonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterMuonID["<<i<<"]=\t" << TOTALnbEventsAfterMuonID[i] << endl;
		}
		cout << endl;
		for(int i = 0; i < 3 ; i++)
		{
			cout << "TOTALnbDimuonsAfterID["<<i<<"]=\t" << TOTALnbDimuonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterDimuonID["<<i<<"]=\t" << TOTALnbEventsAfterDimuonID[i] << endl;
		}
		cout << endl;
// Writing stuff out
	OutputRootFile->Write();
	OutputRootFile->Close();

	// --- To check if everything is ok --- //
        system(Form("touch miniTree_muons_%d.done",ijob));


	delete OutputRootFile;
	delete inputEventTree;
	delete inputRunTree;



	return 0;

}

