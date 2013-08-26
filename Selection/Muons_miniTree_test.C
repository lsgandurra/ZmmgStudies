#include "Muons_miniTree_test.h"


	// ____________________________________________
	// Event information
	// ____________________________________________
	ULong64_t iEvent, iEventID, iLumiID, iRunID;
	Int_t isMM, isMM_nonFSR;
	
	Int_t isNotCommissionned;

	Int_t Muon_eventPassHLT_Mu11;
	Int_t nVertices;
	Int_t nGenVertices;
	Float_t weight_pileUp, weight_Xsection;

	Float_t rho;
	Float_t pu_TrueNumInteractions;
	Int_t pu_NumInteractions, inTimePU_NumInteractions, latePU_NumInteractions, earlyPU_NumInteractions, outOfTimePU_NumInteractions;
	Int_t pu_NumInteractions_inAcceptance, inTimePU_NumInteractions_inAcceptance, latePU_NumInteractions_inAcceptance, earlyPU_NumInteractions_inAcceptance, outOfTimePU_NumInteractions_inAcceptance;

	ULong64_t storeNumber, bunchCrossing, orbitNumber, collisionTimeStamp, microsecondCollisionTime;
	Float_t collisionTime;

	// ____________________________________________
	// Muon variables
	// ____________________________________________
	Int_t NbMuons;

	Float_t Pt_allMuons, Eta_allMuons, Phi_allMuons, Charge_allMuons;
	// (M minus charge, P plus charge), (F far, N near), (L leading, S subleading)
	Float_t MuonM_Pt, MuonP_Pt, MuonN_Pt, MuonF_Pt, MuonL_Pt, MuonS_Pt;
	Float_t MuonM_Eta, MuonP_Eta, MuonN_Eta, MuonF_Eta, MuonL_Eta, MuonS_Eta;
	Float_t MuonM_Phi, MuonP_Phi, MuonN_Phi, MuonF_Phi, MuonL_Phi, MuonS_Phi;
	Float_t MuonM_E, MuonP_E, MuonN_E, MuonF_E, MuonL_E, MuonS_E;
	Float_t MuonM_Px, MuonP_Px, MuonN_Px, MuonF_Px, MuonL_Px, MuonS_Px;
	Float_t MuonM_Py, MuonP_Py, MuonN_Py, MuonF_Py, MuonL_Py, MuonS_Py;
	Float_t MuonM_Pz, MuonP_Pz, MuonN_Pz, MuonF_Pz, MuonL_Pz, MuonS_Pz;
	Int_t MuonF_Charge, MuonN_Charge, MuonL_Charge, MuonS_Charge;

	Float_t MuonM_isoR03_emEt, MuonP_isoR03_emEt, MuonN_isoR03_emEt, MuonF_isoR03_emEt, MuonL_isoR03_emEt, MuonS_isoR03_emEt;
	Float_t MuonM_isoR03_hadEt, MuonP_isoR03_hadEt, MuonN_isoR03_hadEt, MuonF_isoR03_hadEt, MuonL_isoR03_hadEt, MuonS_isoR03_hadEt;
	Float_t MuonM_isoR03_hoEt, MuonP_isoR03_hoEt, MuonN_isoR03_hoEt, MuonF_isoR03_hoEt, MuonL_isoR03_hoEt, MuonS_isoR03_hoEt;
	Float_t MuonM_isoR03_nJets, MuonP_isoR03_nJets, MuonN_isoR03_nJets, MuonF_isoR03_nJets, MuonL_isoR03_nJets, MuonS_isoR03_nJets;
	Float_t MuonM_isoR03_nTracks, MuonP_isoR03_nTracks, MuonN_isoR03_nTracks, MuonF_isoR03_nTracks, MuonL_isoR03_nTracks, MuonS_isoR03_nTracks;
	Float_t MuonM_isoR03_sumPt, MuonP_isoR03_sumPt, MuonN_isoR03_sumPt, MuonF_isoR03_sumPt, MuonL_isoR03_sumPt, MuonS_isoR03_sumPt;

	Float_t MuonM_isoR05_emEt, MuonP_isoR05_emEt, MuonN_isoR05_emEt, MuonF_isoR05_emEt, MuonL_isoR05_emEt, MuonS_isoR05_emEt;
	Float_t MuonM_isoR05_hadEt, MuonP_isoR05_hadEt, MuonN_isoR05_hadEt, MuonF_isoR05_hadEt, MuonL_isoR05_hadEt, MuonS_isoR05_hadEt;
	Float_t MuonM_isoR05_hoEt, MuonP_isoR05_hoEt, MuonN_isoR05_hoEt, MuonF_isoR05_hoEt, MuonL_isoR05_hoEt, MuonS_isoR05_hoEt;
	Float_t MuonM_isoR05_nJets, MuonP_isoR05_nJets, MuonN_isoR05_nJets, MuonF_isoR05_nJets, MuonL_isoR05_nJets, MuonS_isoR05_nJets;
	Float_t MuonM_isoR05_nTracks, MuonP_isoR05_nTracks, MuonN_isoR05_nTracks, MuonF_isoR05_nTracks, MuonL_isoR05_nTracks, MuonS_isoR05_nTracks;
	Float_t MuonM_isoR05_sumPt, MuonP_isoR05_sumPt, MuonN_isoR05_sumPt, MuonF_isoR05_sumPt, MuonL_isoR05_sumPt, MuonS_isoR05_sumPt;

	// ____________________________________________
	// mugamma / mumu / mumugamma information
	// ____________________________________________

	Float_t Mmumu;
	Float_t Ptmumu;	

	// ____________________________________________
	// MC Truth
	// ___________________________________________

	Float_t MuonM_MC_E, MuonM_MC_Px, MuonM_MC_Py, MuonM_MC_Pz, MuonM_MC_Phi, MuonM_MC_Eta, MuonM_MC_Pt;
	Float_t MuonP_MC_E, MuonP_MC_Px, MuonP_MC_Py, MuonP_MC_Pz, MuonP_MC_Phi, MuonP_MC_Eta, MuonP_MC_Pt;
	Float_t MuonN_MC_E, MuonN_MC_Px, MuonN_MC_Py, MuonN_MC_Pz, MuonN_MC_Phi, MuonN_MC_Eta, MuonN_MC_Pt;
	Float_t MuonF_MC_E, MuonF_MC_Px, MuonF_MC_Py, MuonF_MC_Pz, MuonF_MC_Phi, MuonF_MC_Eta, MuonF_MC_Pt;
	Float_t MuonL_MC_E, MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_Phi, MuonL_MC_Eta, MuonL_MC_Pt;
	Float_t MuonS_MC_E, MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_Phi, MuonS_MC_Eta, MuonS_MC_Pt;
	Float_t Mmumu_Muons_MC;	

int main(int argc, char *argv[])
{

	cout << "argc= " << argc << endl;
	for(int iarg = 0 ; iarg < argc; iarg++)
	{
		cout << "argv[" << iarg << "]= " << argv[iarg] << endl;
	}
	
	if( argc == 1 )
	{
		cerr << "arguments should be passed !! sample (outputname) (ntotjob) (ijob) (isZgammaMC) (lumi_set) (pu_set) (low m_mumu cut) (high m_mumu cut) (photon energy correction scheme) (extra photon scale) (applyMuonScaleCorrection) (muon correction sys) (extra resolution) (itoy for mu sys)" << endl;
		return 1;
	}

	// ******************************************
	// First argument is sample
	// ******************************************
	char* sample_char = argv[1];

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
		if( lumi_set == "May10" ) integratedLuminosity = 215.552;
		if( lumi_set == "Promptv4" ) integratedLuminosity = 951.716;
		if( lumi_set == "July05" ) integratedLuminosity = 1.157*1000.0;
		if( lumi_set == "Aug05" ) integratedLuminosity = 389.876;
		if( lumi_set == "Oct03" ) integratedLuminosity = 706.719;
		if( lumi_set == "2011A" ){ integratedLuminosity = 706.370 + 385.819 + 1099; runopt = 0;}
		if( lumi_set == "2011A_rereco" ){ integratedLuminosity = 2.221*1000.0; runopt = 0;}
		if( lumi_set == "2011B" ){ integratedLuminosity = 2741; runopt = 1;}
		//if( lumi_set == "2011" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0;
		//if( lumi_set == "2011" ) integratedLuminosity = 775.208 + 606.941 + 3010 + 1220;
		if( lumi_set == "2011" ) integratedLuminosity = 706.370 + 385.819 + 2741 + 1099;
		if( lumi_set == "2011_rereco" ) integratedLuminosity = 2.221*1000.0 +	2.714*1000.0;
		if( lumi_set == "2012" ) { integratedLuminosity = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0; runopt = 0; }
		if( lumi_set == "2012ABC" ) { integratedLuminosity = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0; runopt = 0; } 
                if( lumi_set == "2012D" ) { integratedLuminosity = 7274.0; runopt = 1; }	
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
	double low_m_mumu = 40.0;
	if( argc > 8 )
	{
		std::stringstream ss ( argv[8] );
		ss >> low_m_mumu;
	}

	// ******************************************
	// Optional argument : high_m_mumu
	// ******************************************
	double high_m_mumu = 80.0;
	if( argc > 9 )
	{
		std::stringstream ss ( argv[9] );
		ss >> high_m_mumu;
	}
	
	// ******************************************
	// Optional argument photon correction scheme
	// ******************************************
	
	//ofstream outfile1(Form("mumugammaCutsTotal_%s_part%i_%i.txt",sample_char,ijob,isZgammaMC));	
	double EScale = 1.0;
	string correction = "";
	bool isManualCorrectionsApplied = false;
	if( argc > 10 )
	{
		correction = argv[10];
		if( (correction == "Louis") || (correction == "Anne-Fleur") || (correction == "START42_V11") || (correction == "ETHZ") || (correction == "Dynamic") || (correction == "MITregression"))
		{
			cout << correction << " correction set will be applied upstream" << endl;
			isManualCorrectionsApplied = true;
		} else {
			std::stringstream ss ( argv[10] );
			ss >> EScale;
		}
	}
	double EScale_inj = EScale;
	
	// ******************************************
	// Optional argument is extra injected scale
	// ******************************************
	double EScale_true_injected = 1.0;
	if( argc > 11 )
	{
		std::stringstream ss ( argv[11] );
		ss >> EScale_true_injected;
	}

	// ******************************************
	// Optional argument : applyMuonScaleCorrection : 0)nothing 1) MuScleFit 2)SIDRA 3)Rochester (21 & 31 also available)
	// ******************************************
	int applyMuonScaleCorrection = 0;
	if( argc > 12 )
	{
		std::stringstream ss ( argv[12] );
		ss >> applyMuonScaleCorrection;
		if( applyMuonScaleCorrection == 0 ) cout << "No muon correction will be applied" << endl;
		if( applyMuonScaleCorrection == 1 )	cout << "MuScleFit muon corrections will be applied upstream" << endl;
		if( applyMuonScaleCorrection == 2 )	cout << "SIDRA muon corrections will be applied upstream" << endl;
		if( applyMuonScaleCorrection == 3 )	cout << "Rochester muon corrections will be applied upstream" << endl;
	}

	// ******************************************
	// Optional argument is extra muon correction smearing (Rochester only)
	// ******************************************
	double sysdev = 0.0;
	if( argc > 13 )
	{
		std::stringstream ss ( argv[13] );
		ss >> sysdev;
	}

	// ******************************************
	// Optional argument is extra resolution
	// ******************************************
	double EResolution = 0.0;
	if( argc > 14 )
	{
		std::stringstream ss ( argv[14] );
		ss >> EResolution;
	}
	EResolution = (double)EResolution / (double)100.0;
	TTimeStamp *time = new TTimeStamp();
	TRandom3* generator = new TRandom3(time->GetNanoSec());
	delete time;
	time = 0;

	// ******************************************
	// Optional argument: random seed for muon corrections: toy number
	// ******************************************
	int itoy = 0;
	if( argc > 15 )
	{
		std::stringstream ss ( argv[15] );
		ss >> itoy;
	}
	UInt_t seeed;
	TTree *seedTree = new TTree();
	seedTree->ReadFile("seeds100toys_v2.txt", "seeed/i");
	seedTree->SetBranchAddress("seeed", &seeed);
	seedTree->GetEntry(itoy);
	cout << "Sanity check: before seed[1]= " << seeed << endl;
	delete seedTree;
	seedTree = 0;
	cout << "Sanity check: after seed[1]= " << seeed << endl;
	

		// TODO

//	char* inputfile = argv[6];

//	TProof * p = TProof::Open("ccaplmaster.in2p3.fr");
	gSystem->Load("libToto.so");
	gROOT->ProcessLine(".L rochcor2012v2v2.h+");
//	gSystem->Load("libFWCoreFWLite.so");
//	gSystem->Load("libDataFormatsFWLite.so");
//	AutoLibraryLoader::enable();
//	#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
//	#include "AutoLibraryLoader.h"
//	TruncatedPyramid a;

	bool doHLT = false;
	bool doMC = (bool)(isZgammaMC >= 1);
	bool doJetMC = false;
	bool doMETMC = false;
	bool doPDFInfo = false;
	bool doSignalMuMuGamma = (bool)(isZgammaMC == 1 || isZgammaMC == 2);
	bool doSignalTopTop = false;
	bool doPhotonConversionMC = false;
	bool doBeamSpot = false;
	bool doPrimaryVertex = true;
	bool doZeePrimaryVertex = false;
	bool doTrack = false;
	bool doJet = false;
	bool doMuon = true;
	bool doElectron = false;
	bool doPhoton = true;
	bool doCluster = true;
	bool doPhotonConversion = true;
	bool doMET = false;
	bool doBardak = false;
	bool doPhotonVertexCorrection = false;
	bool doPhotonIsolation = false;
	bool doR9Rescaling = true;

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
	if( ntotjob == 9999)
	{
		cout<<"Coucou"<<endl;
		inputEventTree->Add("test_DYtoMuMu2.root");
		inputRunTree->Add("test_DYtoMuMu2.root");
	
	} 
	else 
	{
		if (myfile.is_open())
		{
			while( myfile.good() )
			{
				iline++;
				getline(myfile, line);
				if( line == "") continue; // EOF !
				if( iline >= ilineBegin && iline < ilineEnd )
				{
					if( itoy > 0 )
					{ // if we are throwing toys to evaluate muon systematics, copy the files locally
						string localFile = exec(Form("echo %s | rev | cut -d / -f 1 | rev", line.c_str()));
						localFile = localFile.substr(0, localFile.size()-1); // removing end of line character
						string doesTheFileAlreadyExistLocally = exec(Form("if [[ -f %s ]]; then echo ok; else echo ko; fi", localFile.c_str()));
						doesTheFileAlreadyExistLocally = doesTheFileAlreadyExistLocally.substr(0, doesTheFileAlreadyExistLocally.size()-1); // removing end of line character
						if( doesTheFileAlreadyExistLocally == "ko" )
						{
							cout << "File " << localFile << " do not exists locally, copying from dcache" << endl;
							exec(Form("dccp %s%s .", protocol.c_str(), line.c_str()));
						} 
						else 
						{
							cout << "File " << localFile << " already exists locally, do not copy from dcache" << endl;
						}
						inputEventTree->Add(Form("%s", localFile.c_str()));
						inputRunTree->Add(Form("%s", localFile.c_str()));
					} 
					else 
					{ // if we are not throwing toys, don't bother copying locally, just run on the fly
						cout << "Adding file #" << iline << " ( / " << nlines << ") : " << line << endl;
						inputEventTree->Add(Form("%s%s", protocol.c_str(), line.c_str()));
						inputRunTree->Add(Form("%s%s", protocol.c_str(), line.c_str()));
					}
				}
			}
			myfile.close();
		} 
		else 
		{
			cout << "Unable to open file" << endl;
			return 987;
		}
	}
	
	cout<<endl<<"Coucou"<<endl;
	// INSERTFILES

	TFile* OutputRootFile = new TFile(Form("miniTree_%s_part%i.root", sample.c_str(), ijob), "RECREATE");
	TFile* OutputFriendFile = new TFile(Form("miniFriend_%s_part%i.root", sample.c_str(), ijob), "RECREATE");

	OutputRootFile->cd();
	
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

	string ListWantedHLTnames[1];
	ListWantedHLTnames[0] = "HLT_Mu11";
	
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

	vector<string> hltnames; //event->hltAcceptNames();
        miniTree->Branch ("hltnames", "vector<string>", &hltnames);

	miniTree->Branch("isMM", &isMM, "isMM/I");
	miniTree->Branch("isMM_nonFSR", &isMM_nonFSR, "isMM_nonFSR/I");

	miniTree->Branch("nVertices", &nVertices, "nVertices/I");
	miniTree->Branch("nGenVertices", &nGenVertices, "nGenVertices/I");
	miniTree->Branch("weight_pileUp", &weight_pileUp, "weight_pileUp/F");
	miniTree->Branch("weight_Xsection", &weight_Xsection, "weight_Xsection/F");

	miniTree->Branch("rho", &rho, "rho/F");
	miniTree->Branch("pu_TrueNumInteractions", &pu_TrueNumInteractions, "pu_TrueNumInteractions/F");
	miniTree->Branch("pu_NumInteractions", &pu_NumInteractions, "pu_NumInteractions/I");
	miniTree->Branch("inTimePU_NumInteractions", &inTimePU_NumInteractions, "inTimePU_NumInteractions/I");
	miniTree->Branch("latePU_NumInteractions", &latePU_NumInteractions, "latePU_NumInteractions/I");
	miniTree->Branch("earlyPU_NumInteractions", &earlyPU_NumInteractions, "earlyPU_NumInteractions/I");
	miniTree->Branch("outOfTimePU_NumInteractions", &outOfTimePU_NumInteractions, "outOfTimePU_NumInteractions/I");
	miniTree->Branch("pu_NumInteractions_inAcceptance", &pu_NumInteractions_inAcceptance, "pu_NumInteractions_inAcceptance/I");
	miniTree->Branch("inTimePU_NumInteractions_inAcceptance", &inTimePU_NumInteractions_inAcceptance, "inTimePU_NumInteractions_inAcceptance/I");
	miniTree->Branch("latePU_NumInteractions_inAcceptance", &latePU_NumInteractions_inAcceptance, "latePU_NumInteractions_inAcceptance/I");
	miniTree->Branch("earlyPU_NumInteractions_inAcceptance", &earlyPU_NumInteractions_inAcceptance, "earlyPU_NumInteractions_inAcceptance/I");
	miniTree->Branch("outOfTimePU_NumInteractions_inAcceptance", &outOfTimePU_NumInteractions_inAcceptance, "outOfTimePU_NumInteractions_inAcceptance/I");
	miniTree->Branch("storeNumber", &storeNumber, "storeNumber/l");
	miniTree->Branch("bunchCrossing", &bunchCrossing, "bunchCrossing/l");
	miniTree->Branch("orbitNumber", &orbitNumber, "orbitNumber/l");
	miniTree->Branch("collisionTimeStamp", &collisionTimeStamp, "collisionTimeStamp/l");
	miniTree->Branch("microsecondCollisionTime", &microsecondCollisionTime, "microsecondCollisionTime/l");
	miniTree->Branch("collisionTime", &collisionTime, "collisionTime/F");

	// ____________________________________________
	// Muon variables
	// ____________________________________________

	miniTree->Branch("NbMuons", &NbMuons, "NbMuons/I");

	miniTree->Branch("MuonM_Pt", &MuonM_Pt, "MuonM_Pt/F");
	miniTree->Branch("MuonP_Pt", &MuonP_Pt, "MuonP_Pt/F");
	miniTree->Branch("MuonF_Pt", &MuonF_Pt, "MuonF_Pt/F");
	miniTree->Branch("MuonN_Pt", &MuonN_Pt, "MuonN_Pt/F");
	miniTree->Branch("MuonL_Pt", &MuonL_Pt, "MuonL_Pt/F");
	miniTree->Branch("MuonS_Pt", &MuonS_Pt, "MuonS_Pt/F");

	miniTree->Branch("MuonM_Eta", &MuonM_Eta, "MuonM_Eta/F");
	miniTree->Branch("MuonP_Eta", &MuonP_Eta, "MuonP_Eta/F");
	miniTree->Branch("MuonF_Eta", &MuonF_Eta, "MuonF_Eta/F");
	miniTree->Branch("MuonN_Eta", &MuonN_Eta, "MuonN_Eta/F");
	miniTree->Branch("MuonL_Eta", &MuonL_Eta, "MuonL_Eta/F");
	miniTree->Branch("MuonS_Eta", &MuonS_Eta, "MuonS_Eta/F");

	miniTree->Branch("MuonM_Phi", &MuonM_Phi, "MuonM_Phi/F");
	miniTree->Branch("MuonP_Phi", &MuonP_Phi, "MuonP_Phi/F");
	miniTree->Branch("MuonF_Phi", &MuonF_Phi, "MuonF_Phi/F");
	miniTree->Branch("MuonN_Phi", &MuonN_Phi, "MuonN_Phi/F");
	miniTree->Branch("MuonL_Phi", &MuonL_Phi, "MuonL_Phi/F");
	miniTree->Branch("MuonS_Phi", &MuonS_Phi, "MuonS_Phi/F");

	miniTree->Branch("MuonM_E", &MuonM_E, "MuonM_E/F");
	miniTree->Branch("MuonP_E", &MuonP_E, "MuonP_E/F");
	miniTree->Branch("MuonF_E", &MuonF_E, "MuonF_E/F");
	miniTree->Branch("MuonN_E", &MuonN_E, "MuonN_E/F");
	miniTree->Branch("MuonL_E", &MuonL_E, "MuonL_E/F");
	miniTree->Branch("MuonS_E", &MuonS_E, "MuonS_E/F");

	miniTree->Branch("MuonM_Px", &MuonM_Px, "MuonM_Px/F");
	miniTree->Branch("MuonP_Px", &MuonP_Px, "MuonP_Px/F");
	miniTree->Branch("MuonF_Px", &MuonF_Px, "MuonF_Px/F");
	miniTree->Branch("MuonN_Px", &MuonN_Px, "MuonN_Px/F");
	miniTree->Branch("MuonL_Px", &MuonL_Px, "MuonL_Px/F");
	miniTree->Branch("MuonS_Px", &MuonS_Px, "MuonS_Px/F");

	miniTree->Branch("MuonM_Py", &MuonM_Py, "MuonM_Py/F");
	miniTree->Branch("MuonP_Py", &MuonP_Py, "MuonP_Py/F");
	miniTree->Branch("MuonF_Py", &MuonF_Py, "MuonF_Py/F");
	miniTree->Branch("MuonN_Py", &MuonN_Py, "MuonN_Py/F");
	miniTree->Branch("MuonL_Py", &MuonL_Py, "MuonL_Py/F");
	miniTree->Branch("MuonS_Py", &MuonS_Py, "MuonS_Py/F");

	miniTree->Branch("MuonM_Pz", &MuonM_Pz, "MuonM_Pz/F");
	miniTree->Branch("MuonP_Pz", &MuonP_Pz, "MuonP_Pz/F");
	miniTree->Branch("MuonF_Pz", &MuonF_Pz, "MuonF_Pz/F");
	miniTree->Branch("MuonN_Pz", &MuonN_Pz, "MuonN_Pz/F");
	miniTree->Branch("MuonL_Pz", &MuonL_Pz, "MuonL_Pz/F");
	miniTree->Branch("MuonS_Pz", &MuonS_Pz, "MuonS_Pz/F");

	miniTree->Branch("MuonF_Charge", &MuonF_Charge, "MuonF_Charge/I");
	miniTree->Branch("MuonN_Charge", &MuonN_Charge, "MuonN_Charge/I");
	miniTree->Branch("MuonL_Charge", &MuonL_Charge, "MuonL_Charge/I");
	miniTree->Branch("MuonS_Charge", &MuonS_Charge, "MuonS_Charge/I");
	
	miniTree->Branch("MuonM_isoR03_emEt", &MuonM_isoR03_emEt, "MuonM_isoR03_emEt/F");
	miniTree->Branch("MuonP_isoR03_emEt", &MuonP_isoR03_emEt, "MuonP_isoR03_emEt/F");
	miniTree->Branch("MuonF_isoR03_emEt", &MuonF_isoR03_emEt, "MuonF_isoR03_emEt/F");
	miniTree->Branch("MuonN_isoR03_emEt", &MuonN_isoR03_emEt, "MuonN_isoR03_emEt/F");
	miniTree->Branch("MuonL_isoR03_emEt", &MuonL_isoR03_emEt, "MuonL_isoR03_emEt/F");
	miniTree->Branch("MuonS_isoR03_emEt", &MuonS_isoR03_emEt, "MuonS_isoR03_emEt/F");

	miniTree->Branch("MuonM_isoR03_hadEt", &MuonM_isoR03_hadEt, "MuonM_isoR03_hadEt/F");
	miniTree->Branch("MuonP_isoR03_hadEt", &MuonP_isoR03_hadEt, "MuonP_isoR03_hadEt/F");
	miniTree->Branch("MuonF_isoR03_hadEt", &MuonF_isoR03_hadEt, "MuonF_isoR03_hadEt/F");
	miniTree->Branch("MuonN_isoR03_hadEt", &MuonN_isoR03_hadEt, "MuonN_isoR03_hadEt/F");
	miniTree->Branch("MuonL_isoR03_hadEt", &MuonL_isoR03_hadEt, "MuonL_isoR03_hadEt/F");
	miniTree->Branch("MuonS_isoR03_hadEt", &MuonS_isoR03_hadEt, "MuonS_isoR03_hadEt/F");

	miniTree->Branch("MuonM_isoR03_hoEt", &MuonM_isoR03_hoEt, "MuonM_isoR03_hoEt/F");
	miniTree->Branch("MuonP_isoR03_hoEt", &MuonP_isoR03_hoEt, "MuonP_isoR03_hoEt/F");
	miniTree->Branch("MuonF_isoR03_hoEt", &MuonF_isoR03_hoEt, "MuonF_isoR03_hoEt/F");
	miniTree->Branch("MuonN_isoR03_hoEt", &MuonN_isoR03_hoEt, "MuonN_isoR03_hoEt/F");
	miniTree->Branch("MuonL_isoR03_hoEt", &MuonL_isoR03_hoEt, "MuonL_isoR03_hoEt/F");
	miniTree->Branch("MuonS_isoR03_hoEt", &MuonS_isoR03_hoEt, "MuonS_isoR03_hoEt/F");

	miniTree->Branch("MuonM_isoR03_nJets", &MuonM_isoR03_nJets, "MuonM_isoR03_nJets/F");
	miniTree->Branch("MuonP_isoR03_nJets", &MuonP_isoR03_nJets, "MuonP_isoR03_nJets/F");
	miniTree->Branch("MuonF_isoR03_nJets", &MuonF_isoR03_nJets, "MuonF_isoR03_nJets/F");
	miniTree->Branch("MuonN_isoR03_nJets", &MuonN_isoR03_nJets, "MuonN_isoR03_nJets/F");
	miniTree->Branch("MuonL_isoR03_nJets", &MuonL_isoR03_nJets, "MuonL_isoR03_nJets/F");
	miniTree->Branch("MuonS_isoR03_nJets", &MuonS_isoR03_nJets, "MuonS_isoR03_nJets/F");

	miniTree->Branch("MuonM_isoR03_nTracks", &MuonM_isoR03_nTracks, "MuonM_isoR03_nTracks/F");
	miniTree->Branch("MuonP_isoR03_nTracks", &MuonP_isoR03_nTracks, "MuonP_isoR03_nTracks/F");
	miniTree->Branch("MuonF_isoR03_nTracks", &MuonF_isoR03_nTracks, "MuonF_isoR03_nTracks/F");
	miniTree->Branch("MuonN_isoR03_nTracks", &MuonN_isoR03_nTracks, "MuonN_isoR03_nTracks/F");
	miniTree->Branch("MuonL_isoR03_nTracks", &MuonL_isoR03_nTracks, "MuonL_isoR03_nTracks/F");
	miniTree->Branch("MuonS_isoR03_nTracks", &MuonS_isoR03_nTracks, "MuonS_isoR03_nTracks/F");

	miniTree->Branch("MuonM_isoR03_sumPt", &MuonM_isoR03_sumPt, "MuonM_isoR03_sumPt/F");
	miniTree->Branch("MuonP_isoR03_sumPt", &MuonP_isoR03_sumPt, "MuonP_isoR03_sumPt/F");
	miniTree->Branch("MuonF_isoR03_sumPt", &MuonF_isoR03_sumPt, "MuonF_isoR03_sumPt/F");
	miniTree->Branch("MuonN_isoR03_sumPt", &MuonN_isoR03_sumPt, "MuonN_isoR03_sumPt/F");
	miniTree->Branch("MuonL_isoR03_sumPt", &MuonL_isoR03_sumPt, "MuonL_isoR03_sumPt/F");
	miniTree->Branch("MuonS_isoR03_sumPt", &MuonS_isoR03_sumPt, "MuonS_isoR03_sumPt/F");

	miniTree->Branch("MuonM_isoR05_emEt", &MuonM_isoR05_emEt, "MuonM_isoR05_emEt/F");
	miniTree->Branch("MuonP_isoR05_emEt", &MuonP_isoR05_emEt, "MuonP_isoR05_emEt/F");
	miniTree->Branch("MuonF_isoR05_emEt", &MuonF_isoR05_emEt, "MuonF_isoR05_emEt/F");
	miniTree->Branch("MuonN_isoR05_emEt", &MuonN_isoR05_emEt, "MuonN_isoR05_emEt/F");
	miniTree->Branch("MuonL_isoR05_emEt", &MuonL_isoR05_emEt, "MuonL_isoR05_emEt/F");
	miniTree->Branch("MuonS_isoR05_emEt", &MuonS_isoR05_emEt, "MuonS_isoR05_emEt/F");

	miniTree->Branch("MuonM_isoR05_hadEt", &MuonM_isoR05_hadEt, "MuonM_isoR05_hadEt/F");
	miniTree->Branch("MuonP_isoR05_hadEt", &MuonP_isoR05_hadEt, "MuonP_isoR05_hadEt/F");
	miniTree->Branch("MuonF_isoR05_hadEt", &MuonF_isoR05_hadEt, "MuonF_isoR05_hadEt/F");
	miniTree->Branch("MuonN_isoR05_hadEt", &MuonN_isoR05_hadEt, "MuonN_isoR05_hadEt/F");
	miniTree->Branch("MuonL_isoR05_hadEt", &MuonL_isoR05_hadEt, "MuonL_isoR05_hadEt/F");
	miniTree->Branch("MuonS_isoR05_hadEt", &MuonS_isoR05_hadEt, "MuonS_isoR05_hadEt/F");

	miniTree->Branch("MuonM_isoR05_hoEt", &MuonM_isoR05_hoEt, "MuonM_isoR05_hoEt/F");
	miniTree->Branch("MuonP_isoR05_hoEt", &MuonP_isoR05_hoEt, "MuonP_isoR05_hoEt/F");
	miniTree->Branch("MuonF_isoR05_hoEt", &MuonF_isoR05_hoEt, "MuonF_isoR05_hoEt/F");
	miniTree->Branch("MuonN_isoR05_hoEt", &MuonN_isoR05_hoEt, "MuonN_isoR05_hoEt/F");
	miniTree->Branch("MuonL_isoR05_hoEt", &MuonL_isoR05_hoEt, "MuonL_isoR05_hoEt/F");
	miniTree->Branch("MuonS_isoR05_hoEt", &MuonS_isoR05_hoEt, "MuonS_isoR05_hoEt/F");

	miniTree->Branch("MuonM_isoR05_nJets", &MuonM_isoR05_nJets, "MuonM_isoR05_nJets/F");
	miniTree->Branch("MuonP_isoR05_nJets", &MuonP_isoR05_nJets, "MuonP_isoR05_nJets/F");
	miniTree->Branch("MuonF_isoR05_nJets", &MuonF_isoR05_nJets, "MuonF_isoR05_nJets/F");
	miniTree->Branch("MuonN_isoR05_nJets", &MuonN_isoR05_nJets, "MuonN_isoR05_nJets/F");
	miniTree->Branch("MuonL_isoR05_nJets", &MuonL_isoR05_nJets, "MuonL_isoR05_nJets/F");
	miniTree->Branch("MuonS_isoR05_nJets", &MuonS_isoR05_nJets, "MuonS_isoR05_nJets/F");

	miniTree->Branch("MuonM_isoR05_nTracks", &MuonM_isoR05_nTracks, "MuonM_isoR05_nTracks/F");
	miniTree->Branch("MuonP_isoR05_nTracks", &MuonP_isoR05_nTracks, "MuonP_isoR05_nTracks/F");
	miniTree->Branch("MuonF_isoR05_nTracks", &MuonF_isoR05_nTracks, "MuonF_isoR05_nTracks/F");
	miniTree->Branch("MuonN_isoR05_nTracks", &MuonN_isoR05_nTracks, "MuonN_isoR05_nTracks/F");
	miniTree->Branch("MuonL_isoR05_nTracks", &MuonL_isoR05_nTracks, "MuonL_isoR05_nTracks/F");
	miniTree->Branch("MuonS_isoR05_nTracks", &MuonS_isoR05_nTracks, "MuonS_isoR05_nTracks/F");

	miniTree->Branch("MuonM_isoR05_sumPt", &MuonM_isoR05_sumPt, "MuonM_isoR05_sumPt/F");
	miniTree->Branch("MuonP_isoR05_sumPt", &MuonP_isoR05_sumPt, "MuonP_isoR05_sumPt/F");
	miniTree->Branch("MuonF_isoR05_sumPt", &MuonF_isoR05_sumPt, "MuonF_isoR05_sumPt/F");
	miniTree->Branch("MuonN_isoR05_sumPt", &MuonN_isoR05_sumPt, "MuonN_isoR05_sumPt/F");
	miniTree->Branch("MuonL_isoR05_sumPt", &MuonL_isoR05_sumPt, "MuonL_isoR05_sumPt/F");
	miniTree->Branch("MuonS_isoR05_sumPt", &MuonS_isoR05_sumPt, "MuonS_isoR05_sumPt/F");
		
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
	miniTree->Branch("MuonN_MC_E", &MuonN_MC_E, "MuonN_MC_E/F");
	miniTree->Branch("MuonN_MC_Px", &MuonN_MC_Px, "MuonN_MC_Px/F");
	miniTree->Branch("MuonN_MC_Py", &MuonN_MC_Py, "MuonN_MC_Py/F");
	miniTree->Branch("MuonN_MC_Pz", &MuonN_MC_Pz, "MuonN_MC_Pz/F");
	miniTree->Branch("MuonN_MC_Phi", &MuonN_MC_Phi, "MuonN_MC_Phi/F");
	miniTree->Branch("MuonN_MC_Eta", &MuonN_MC_Eta, "MuonN_MC_Eta/F");
	miniTree->Branch("MuonN_MC_Pt", &MuonN_MC_Pt, "MuonN_MC_Pt/F");
	miniTree->Branch("MuonF_MC_E", &MuonF_MC_E, "MuonF_MC_E/F");
	miniTree->Branch("MuonF_MC_Px", &MuonF_MC_Px, "MuonF_MC_Px/F");
	miniTree->Branch("MuonF_MC_Py", &MuonF_MC_Py, "MuonF_MC_Py/F");
	miniTree->Branch("MuonF_MC_Pz", &MuonF_MC_Pz, "MuonF_MC_Pz/F");
	miniTree->Branch("MuonF_MC_Phi", &MuonF_MC_Phi, "MuonF_MC_Phi/F");
	miniTree->Branch("MuonF_MC_Eta", &MuonF_MC_Eta, "MuonF_MC_Eta/F");
	miniTree->Branch("MuonF_MC_Pt", &MuonF_MC_Pt, "MuonF_MC_Pt/F");
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
	bool powheg = true;
	bool signal = false;
	bool stew = false;
	bool zjet_veto = (bool)(isZgammaMC == 1 || isZgammaMC == 2);
	cout << "Nb of events : " << NbEvents << endl;
	cout << "Signal is: " << signal <<endl;
	cout << "Stew is: " << stew << endl;
	cout << "ZJet veto is: " << zjet_veto << endl;
	int nBeforeAllCuts, nAfterCutPthatFilter, nAfterCutCSA07ID, nAfterCutZJETVETO, nAfterVeryLooseMMG, nAfterJanLooseMMG, nAfterLooseMMG, nAfterCut1c, nAfterTightMMG, nAfterCut1e, nAfterCut2a, nAfterCut2b, nAfterCut2c, nAfterCut3, nAfterCut4, nAfterCut5, nAfterCut6, nAfterCut7, nAfterCut8, nAfterCut9, nAfterCut10, nSelected;
	nBeforeAllCuts = nAfterCutPthatFilter = nAfterCutCSA07ID = nAfterCutZJETVETO = nAfterVeryLooseMMG = nAfterJanLooseMMG = nAfterLooseMMG = nAfterCut1c = nAfterTightMMG = nAfterCut1e = nAfterCut2a = nAfterCut2b = nAfterCut2c = nAfterCut3 = nAfterCut4 = nAfterCut5 = nAfterCut6 = nAfterCut7 = nAfterCut8 = nAfterCut9 = nAfterCut10 = nSelected = 0;

	string lastFile = "";

	double XSectionDYToMuMu = 1914.894;
	double XSectionTTJets = 234.0;
	double XSectionWJetsToLNu = 37509.25;
	double XSectionQCDMu = 349988.0; //FIXME

	double InitialNumberDYToMuMu = 48819386.0;
	double InitialNumberTTJets = 6736135.0;
	double InitialNumberWJetsToLNu = 57709905.0;
	double InitialNumberQCDMu = 8797418.0; //FIXME

	double minPtHat = -100;
	double maxPtHat = 1000000;
	int verbosity = 0;
	int Nb_events_outside_powheg_cuts = 0;
	int TOTALnbMuonsAfterID[5] = {0};
	int TOTALnbEventsAfterMuonID[12] = {0};
	int TOTALnbDimuonsAfterID[3] = {0};
	int TOTALnbEventsAfterDimuonID[3] = {0};
	int TOTALnbPhotonsAfterID[6] = {0};
	int TOTALnbEventsAfterPhotonID[6] = {0};
	int TOTALnbMuMuGammaAfterID[10] = {0};
	int TOTALnbEventsAfterMuMuGammaID[10] = {0};

	int NbEventsPerJob = NbEvents;
	int NbEventsBegin = 0;
	int NbEventsEnd = NbEvents;
	if( ntotjob == 9999 && ijob != -1 )
	{
		NbEventsPerJob = 200000;
		NbEventsBegin = ijob * NbEventsPerJob;
		NbEventsEnd = min( (ijob + 1)* NbEventsPerJob , (int)NbEvents);
		NbEvents = NbEventsEnd - NbEventsBegin ;
	}
	cout << "NbEventsBegin= " << NbEventsBegin << "\tNbEventsEnd= " << NbEventsEnd << "\tNbEventsPerJob= " << NbEventsPerJob << endl;

/*	
	int runIDTab[55] = {0};
	int lumiIDTab[55] = {0};
	int eventIDTab[55] = {0};
	int tempNumber = 0;

	ifstream strangeFile("428_strange_events.txt");
        for(int i = 0; i < 55; i++)
	{     
		strangeFile>>runIDTab[i];
		strangeFile>>lumiIDTab[i];
		strangeFile>>eventIDTab[i];			
	}

	strangeFile.close();	
*/
	// ***************************************************************************************************
	// ***************************************************************************************************
	// ********************************** LOOP over events ***********************************************
	// ***************************************************************************************************
	// ***************************************************************************************************
	
	for(unsigned int ievt=NbEventsBegin; ievt<NbEventsEnd; ievt++)
	{
		if(verbosity>4) cout << "analysing event ievt= " << ievt << endl;
		nBeforeAllCuts++;
		int nprint = (int)((double)NbEventsPerJob/(double)100.0);
		if( (NbEvents >= 100) && (ievt % nprint)==0 ){ cout<< ievt <<" events done over "<<NbEvents<<" ( "<<ceil((double)ievt/(double)NbEventsPerJob*100)<<" \% )"<<endl; }

		iEvent = ievt;
		inputEventTree->GetEvent(ievt);
	
		hltnames = event->hltAcceptNames();

/*
		tempNumber = 0;
		for(int i = 0; i < 55; i++)
        	{
			if(event->runId() == runIDTab[i] && event->luminosityBlock() == lumiIDTab[i] && event->eventId() == eventIDTab[i])
			{
				tempNumber = 1;
			}	
		}

		if(tempNumber == 0) continue;	
*/
		// ____________________________________________
		// Event information
		// ____________________________________________
		iEventID = event->eventId();
		iLumiID = event->luminosityBlock();
		iRunID = event->runId();
		nVertices = vertices->GetEntries();
		nGenVertices = vertices->GetEntries();
		if( (isZgammaMC == 1) || (isZgammaMC == 2) )
		{
			nGenVertices = event->pu_TrueNumInteractions();
		}

		rho = event->rho();
/*
		// Can be uncommented for IpnTreeProducer versions > RECO_4_2_8_v3
		pu_TrueNumInteractions = event->pu_TrueNumInteractions();
		pu_NumInteractions = event->pu_NumInteractions();
		inTimePU_NumInteractions = event->inTimePU_NumInteractions();
		latePU_NumInteractions = event->latePU_NumInteractions();
		earlyPU_NumInteractions = event->earlyPU_NumInteractions();
		outOfTimePU_NumInteractions = event->outOfTimePU_NumInteractions();
		pu_NumInteractions_inAcceptance = event->pu_NumInteractions_inAcceptance();
		inTimePU_NumInteractions_inAcceptance = event->inTimePU_NumInteractions_inAcceptance();
		latePU_NumInteractions_inAcceptance = event->latePU_NumInteractions_inAcceptance();
		earlyPU_NumInteractions_inAcceptance = event->earlyPU_NumInteractions_inAcceptance();
		outOfTimePU_NumInteractions_inAcceptance = event->outOfTimePU_NumInteractions_inAcceptance();
*/
		storeNumber = event->storeNumber();
		bunchCrossing = event->bunchCrossing();
		orbitNumber = event->orbitNumber();
		collisionTimeStamp = event->collisionTimeStamp();
		microsecondCollisionTime = event->microsecondCollisionTime();
		collisionTime = event->collisionTime();


		weight_Xsection = weight_pileUp = 1.0;

		string sample_in(sample_char);
		
		if( sample_in.find("DYToMuMu") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionDYToMuMu) / (double)(InitialNumberDYToMuMu)) * (double)integratedLuminosity);
			weight_pileUp = weight_DYToMuMu(nGenVertices+1, lumi_set, pu_set);
		} 
		else if( sample_in.find("QCD") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionQCDMu) / (double)(InitialNumberQCDMu)) * (double)integratedLuminosity);
			weight_pileUp = weight_QCDMu(nGenVertices+1);
		} 
		else if( sample_in.find("TTJets") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionTTJets) / (double)(InitialNumberTTJets)) * (double)integratedLuminosity);
			weight_pileUp = weight_TTJets(nGenVertices+1, lumi_set, pu_set);
		} 
		else if( sample_in.find("WToMuNu") != string::npos	|| sample_in.find("WJetsToLNu") != string::npos)
		{
			weight_Xsection = (double)(	(double)((double)(XSectionWJetsToLNu) / (double)(InitialNumberWJetsToLNu)) * (double)integratedLuminosity);
			weight_pileUp = weight_WJetsToLNu(nGenVertices+1, lumi_set, pu_set);
		}


		isMM = isMM_nonFSR = 0;


		// ____________________________________________
		// Muon variables
		// ____________________________________________
		NbMuons = muons->GetEntries();
		Pt_allMuons = Eta_allMuons = Phi_allMuons = Charge_allMuons = -99;
		MuonM_Pt = MuonP_Pt = MuonN_Pt = MuonF_Pt = MuonL_Pt = MuonS_Pt = -99;
		MuonM_Eta = MuonP_Eta = MuonN_Eta = MuonF_Eta = MuonL_Eta = MuonS_Eta = -99;
		MuonM_Phi = MuonP_Phi = MuonN_Phi = MuonF_Phi = MuonL_Phi = MuonS_Phi = -99;
		MuonM_E = MuonP_E = MuonN_E = MuonF_E = MuonL_E = MuonS_E = -99;
		MuonM_Px = MuonP_Px = MuonN_Px = MuonF_Px = MuonL_Px = MuonS_Px = -99;
		MuonM_Py = MuonP_Py = MuonN_Py = MuonF_Py = MuonL_Py = MuonS_Py = -99;
		MuonM_Pz = MuonP_Pz = MuonN_Pz = MuonF_Pz = MuonL_Pz = MuonS_Pz = -99;
		MuonF_Charge = MuonN_Charge = MuonL_Charge = MuonS_Charge = -99;
		MuonM_isoR03_emEt = MuonP_isoR03_emEt = MuonN_isoR03_emEt = MuonF_isoR03_emEt = MuonL_isoR03_emEt = MuonS_isoR03_emEt = -99;
		MuonM_isoR03_hadEt = MuonP_isoR03_hadEt = MuonN_isoR03_hadEt = MuonF_isoR03_hadEt = MuonL_isoR03_hadEt = MuonS_isoR03_hadEt = -99;
		MuonM_isoR03_hoEt = MuonP_isoR03_hoEt = MuonN_isoR03_hoEt = MuonF_isoR03_hoEt = MuonL_isoR03_hoEt = MuonS_isoR03_hoEt = -99;
		MuonM_isoR03_nJets = MuonP_isoR03_nJets = MuonN_isoR03_nJets = MuonF_isoR03_nJets = MuonL_isoR03_nJets = MuonS_isoR03_nJets = -99;
		MuonM_isoR03_nTracks = MuonP_isoR03_nTracks = MuonN_isoR03_nTracks = MuonF_isoR03_nTracks = MuonL_isoR03_nTracks = MuonS_isoR03_nTracks = -99;
		MuonM_isoR03_sumPt = MuonP_isoR03_sumPt = MuonN_isoR03_sumPt = MuonF_isoR03_sumPt = MuonL_isoR03_sumPt = MuonS_isoR03_sumPt = -99;
		MuonM_isoR05_emEt = MuonP_isoR05_emEt = MuonN_isoR05_emEt = MuonF_isoR05_emEt = MuonL_isoR05_emEt = MuonS_isoR05_emEt = -99;
		MuonM_isoR05_hadEt = MuonP_isoR05_hadEt = MuonN_isoR05_hadEt = MuonF_isoR05_hadEt = MuonL_isoR05_hadEt = MuonS_isoR05_hadEt = -99;
		MuonM_isoR05_hoEt = MuonP_isoR05_hoEt = MuonN_isoR05_hoEt = MuonF_isoR05_hoEt = MuonL_isoR05_hoEt = MuonS_isoR05_hoEt = -99;
		MuonM_isoR05_nJets = MuonP_isoR05_nJets = MuonN_isoR05_nJets = MuonF_isoR05_nJets = MuonL_isoR05_nJets = MuonS_isoR05_nJets = -99;
		MuonM_isoR05_nTracks = MuonP_isoR05_nTracks = MuonN_isoR05_nTracks = MuonF_isoR05_nTracks = MuonL_isoR05_nTracks = MuonS_isoR05_nTracks = -99;
		MuonM_isoR05_sumPt = MuonP_isoR05_sumPt = MuonN_isoR05_sumPt = MuonF_isoR05_sumPt = MuonL_isoR05_sumPt = MuonS_isoR05_sumPt = -99;

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
		MuonN_MC_E = MuonN_MC_Px = MuonN_MC_Py = MuonN_MC_Pz = MuonN_MC_Phi = MuonN_MC_Eta = MuonN_MC_Pt = -99.0;
		MuonF_MC_E = MuonF_MC_Px = MuonF_MC_Py = MuonF_MC_Pz = MuonF_MC_Phi = MuonF_MC_Eta = MuonF_MC_Pt = -99.0;
		MuonL_MC_E = MuonL_MC_Px = MuonL_MC_Py = MuonL_MC_Pz = MuonL_MC_Phi = MuonL_MC_Eta = MuonL_MC_Pt = -99.0;
		MuonS_MC_E = MuonS_MC_Px = MuonS_MC_Py = MuonS_MC_Pz = MuonS_MC_Phi = MuonS_MC_Eta = MuonS_MC_Pt = -99.0;
		Mmumu_Muons_MC = -99.0;


		// ____________________________________________
		// END OF INITIALIZATION
		// ____________________________________________


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
		vector<double> muonIdentified_corrected_Pt;
		muonIdentified_corrected_Pt.clear();
		vector<double> muonIdentified_corrected_Pz;
                muonIdentified_corrected_Pz.clear();
		vector<double> muonIdentified_corrected_E;
                muonIdentified_corrected_E.clear();	
		int nbMuonsAfterID[5] = {0};
		
		if(verbosity>0) cerr << "\t\tThere is " << NbMuons << " muons in the muon collection" << endl;

		nbMuonsAfterID[0] = NbMuons;
		TOTALnbMuonsAfterID[0] += NbMuons;
		for(int imuon=0 ; imuon<NbMuons ; imuon++)
		{
			TRootMuon *mymuon;
			mymuon = (TRootMuon*) muons->At(imuon);
	
			/* // 2010 muon ID

			if(! (mymuon->isGlobalMuon()) )
			{
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because not global muon" << endl;
				continue;
			}
			nbMuonsAfterID[1]++;
			TOTALnbMuonsAfterID[1]++;

			if(! (mymuon->normalizedGlobalChi2()<10.0) )
			{// chi2/ndof of the global muon fit < 10
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because chi2/ndof of the global muon fit < 10 (" << mymuon->normalizedGlobalChi2() << ")" << endl;
				continue;
			}
			nbMuonsAfterID[2]++;
			TOTALnbMuonsAfterID[2]++;

			if(! (mymuon->numberOfValidGlobalHits()>0) )
			{// number of valid muon hits matched to the global fit
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because number of valid muon hits matched to the global fit too low (" << mymuon->numberOfValidGlobalHits() << ")" << endl;
				continue;
			}
			nbMuonsAfterID[3]++;
			TOTALnbMuonsAfterID[3]++;

			if(! (mymuon->isTrackerMuon()) )
			{// The muon must be identified as Tracker Muon.
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because not tracker muon" << endl;
				continue;
			}
			nbMuonsAfterID[4]++;
			TOTALnbMuonsAfterID[4]++;

			if(! (mymuon->numberOfMatches()>1) )
			{// number of muon stations with matched segments (global track: out-in fit)
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because number of muon stations with matched segments (global track: out-in fit) too low" << endl;
				continue;
			}
			nbMuonsAfterID[5]++;
			TOTALnbMuonsAfterID[5]++;

			if(! (mymuon->numberOfValidTrackerHits()>10) )
			{// number of tracker (pixels + strips) hits
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because number of tracker (pixels + strips) hits" << endl;
				continue;
			}
			nbMuonsAfterID[6]++;
			TOTALnbMuonsAfterID[6]++;

			if(! (mymuon->numberOfValidPixelHits()>0) )
			{// number of pixel hits
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


			if(! (mymuon->isoR03_sumPt()<3.0) )
			{// sum of pT of tracks with pT >1.5 within a cone of DR < 0.3 around the muon direction, vetoing a cone of 0.015 around that direction
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because sum of pT of tracks with pT >1.5 within a cone of DR < 0.3 around the muon direction, vetoing a cone of 0.015 around that direction" << endl;
				continue;
			}

			nbMuonsAfterID[9]++;
			TOTALnbMuonsAfterID[9]++;
			*/
			
			
			// 2012 Tight muon ID
			if(! (mymuon->isTightMuonPerso() == 1))
			{
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because isTightMuonPerso() != 1" << endl;
				continue;
			}	
			nbMuonsAfterID[1]++;
                        TOTALnbMuonsAfterID[1]++;	


			TLorentzVector correctedMuon(mymuon->Px(), mymuon->Py(), mymuon->Pz(), mymuon->E());
			float qter = 1.0;
			double corrected_Pt = mymuon->Pt();
			rochcor2012 *rmcor = new rochcor2012(seeed);
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
				{ // Apply Rochester corrections
					TLorentzVector muonRochester(mymuon->Px(), mymuon->Py(), mymuon->Pz(), mymuon->E());
					TLorentzVector muonRochesterDummy(mymuon->Px(), mymuon->Py(), mymuon->Pz(), mymuon->E());
					// moption = 1	: recommended by the authors (better match in Z mass profile vs. eta/phi between the reconstructed and generated Z mass)
					// sysdev = 0 : no systematics yet
		
					if( isZgammaMC > 0 ) // If sample is MC
					{
						if( (lumi_set == "2011A") || (lumi_set == "2011A_rereco") ) runopt = 0;
						if( (lumi_set == "2011B") || (lumi_set == "2011B_rereco") ) runopt = 1;
						if( lumi_set == "2011A" ) integratedLuminosity = 706.370 + 385.819 + 1099;
						if( lumi_set == "2011A_rereco" ) integratedLuminosity = 2.221*1000.0;
						if( lumi_set == "2011B" ) integratedLuminosity = 2741;
						//if( lumi_set == "2011" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0;
						if( lumi_set == "2011" ) integratedLuminosity = 706.370 + 385.819 + 2741 + 1099;
						if( lumi_set == "2011_rereco" ) integratedLuminosity = 2.221*1000.0 +	2.714*1000.0;
						if( lumi_set == "2012" ) integratedLuminosity = 808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0;
	
						// event where to switch from Run2011A to Run2011B correction
						// We run over ntot=  DYToMuMu events
						// We run on average on avEventsPerJob= (ntot) / (ntotjob)= events per job
						// RunA represents Astat= (2.221) / (2.221 + 2.714) = 45.01 % of the total stat
						// RunB represents Bstat= (2.714) / (2.221 + 2.714) = 54.99 % of the total stat
						// The last event of 'RunA' is lastA= Astat*ntot
						// Which should happen in job ijobLastA= (int)(lastA)/(int)(avEventsPerJob)
						// and will be the event iEventLastA= (int)(lastA)%(int)(avEventsPerJob)
						// ***********************************
						// temp version
						// ***********************************
						if( (lumi_set == "2011" ) && (sample_in.find("DYToMuMu") != string::npos))
						{
							int ntot= 9198125;
							double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
							double Astat= (double)(706.370 + 385.819 + 1099) / (double)(706.370 + 385.819 + 2741 + 1099);
							double Bstat= (double)(2741) / (double)(706.370 + 385.819 + 2741 + 1099);
							double lastA= Astat*ntot;
							int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
							int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
							if( ijob == (ijobLastA -1) )
							{
								if( ievt >= iEventLastA )
								{
									runopt = 1;
								} 
								else runopt = 0;
							} 
							else if( ijob >= ijobLastA ) runopt = 1;
						}
						if( (lumi_set == "2012" ) && (sample_in.find("DYToMuMu") != string::npos))
                                                {
                                                        int ntot= 3568581;
                                                        double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
                                                        double ABCstat= (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0) / (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0);
                                                        double Dstat= (double)(7274.0) / (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0);
                                                        double lastABC= ABCstat*ntot;
                                                        int ijobLastABC= (int)(lastABC)/(int)(avEventsPerJob);
                                                        int iEventLastABC= (int)(lastABC)%(int)(avEventsPerJob);
                                                        if( ijob == (ijobLastABC -1) )
                                                        {
                                                                if( ievt >= iEventLastABC )
                                                                {
                                                                        runopt = 1;
                                                                }
                                                                else runopt = 0;
                                                        }
                                                        else if( ijob >= ijobLastABC ) runopt = 1;
                                                }

						if( (lumi_set == "2011" ) && (sample_in.find("TTJets") != string::npos))
	                                        {
	                                                int ntot= 131388;
	                                                double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
	                                                double Astat= (double)(706.370 + 385.819 + 1099) / (double)(706.370 + 385.819 + 2741 + 1099);
	                                                double Bstat= (double)(2741) / (double)(706.370 + 385.819 + 2741 + 1099);
	                                                double lastA= Astat*ntot;
	                                                int ijobLastA= (int)(lastA)/(int)(avEventsPerJob);
	                                                int iEventLastA= (int)(lastA)%(int)(avEventsPerJob);
	                                                if( ijob == (ijobLastA -1) )
	                                                {
	                                                        if( ievt >= iEventLastA )
	                                                        {
	                                                                runopt = 1; 
	                                                        } 
								else runopt = 0; 
	                                                } 
							else if( ijob >= ijobLastA ) runopt = 1; 
	                                        }
						if( (lumi_set == "2012" ) && (sample_in.find("TTJets") != string::npos))
                                                {
                                                        int ntot= 110095;
                                                        double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
                                                        double ABCstat= (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0) / (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0);
                                                        double Dstat= (double)(7274.0) / (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0);
                                                        double lastABC= ABCstat*ntot;
                                                        int ijobLastABC= (int)(lastABC)/(int)(avEventsPerJob);
                                                        int iEventLastABC= (int)(lastABC)%(int)(avEventsPerJob);
                                                        if( ijob == (ijobLastABC -1) )
                                                        {
                                                                if( ievt >= iEventLastABC )
                                                                {
                                                                        runopt = 1;
                                                                }
                                                                else runopt = 0;
                                                        }
                                                        else if( ijob >= ijobLastABC ) runopt = 1;
                                                }
	
						if( (lumi_set == "2011" ) && (sample_in.find("WJetsToLNu") != string::npos))
	                                        {
	                                                int ntot= 14821;
	                                                double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
	                                                double Astat= (double)(706.370 + 385.819 + 1099) / (double)(706.370 + 385.819 + 2741 + 1099);
	                                                double Bstat= (double)(2741) / (double)(706.370 + 385.819 + 2741 + 1099);
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
						 if( (lumi_set == "2012" ) && (sample_in.find("WJetsToLNu") != string::npos))
                                                {
                                                        int ntot= 3353;
                                                        double avEventsPerJob= (double)(ntot)/(double)(ntotjob);
                                                        double ABCstat= (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0) / (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0);
                                                        double Dstat= (double)(7274.0) / (double)(808.472 + 82.136 + 4429.0 + 495.003 + 134.242 + 6397.0 + 7274.0);
                                                        double lastABC= ABCstat*ntot;
                                                        int ijobLastABC= (int)(lastABC)/(int)(avEventsPerJob);
                                                        int iEventLastABC= (int)(lastABC)%(int)(avEventsPerJob);
                                                        if( ijob == (ijobLastABC -1) )
                                                        {
                                                                if( ievt >= iEventLastABC )
                                                                {
                                                                        runopt = 1;
                                                                }
                                                                else runopt = 0;
                                                        }
                                                        else if( ijob >= ijobLastABC ) runopt = 1;
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
						rmcor->momcor_mc(muonRochester, mymuon->charge(), runopt, qter);
						//if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
						//else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
					} 
					else 
					{
						rmcor->momcor_data(muonRochester, mymuon->charge(), runopt, qter);
						//if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, runopt, qter);
						//else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, runopt, qter);
					}
					corrected_Pt = muonRochester.Pt();
					correctedMuon = muonRochester;
				}
				if( applyMuonScaleCorrection == 31 )
				{ // Apply Rochester corrections on top of MuscleFit
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
						rmcor->momcor_mc(muonRochester, mymuon->charge(), runopt, qter);
						//if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
						//else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
					} 
					else 
					{
						rmcor->momcor_data(muonRochester, mymuon->charge(), runopt, qter);
						//if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev);
						//else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev);
					}
					corrected_Pt = muonRochester.Pt();
				}
				if( applyMuonScaleCorrection == 32 )
				{ // Apply Rochester corrections on top of Sidra
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
						rmcor->momcor_mc(muonRochester, mymuon->charge(), runopt, qter);
						//if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
						//else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
					} 
					else 
					{
						rmcor->momcor_data(muonRochester,mymuon->charge(), runopt, qter);
						//if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev);
						//else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev);
					}
					corrected_Pt = muonRochester.Pt();
				}
				if( applyMuonScaleCorrection == 99 )
				{ // Dummy muon correction to check code: it's NOT supposed to do anything
					corrected_Pt = mymuon->Pt();
				}
				// Sidra makes MC look like data
				//if( isZgammaMC > 0) corrected_Pt = applySidra(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi(), generator);
				// MuScleFit correct data absolute scale
				//corrected_Pt = applyMuScleFit(corrected_Pt, mymuon->charge(), mymuon->Eta(), mymuon->Phi());
			}

			/*
			double corrected_Pz = mymuon->Pz();
			double corrected_Px = applyMuonScaleCorrection > 0 ? mymuon->Px() * (double)(corrected_Pt)/(double)(mymuon->Pt()) : mymuon->Px();
			double corrected_Py = applyMuonScaleCorrection > 0 ? mymuon->Py() * (double)(corrected_Pt)/(double)(mymuon->Pt()) : mymuon->Py();
			double m_mu = 105.658367e-3;
			double corrected_E = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz * corrected_Pz + corrected_Pt * corrected_Pt) ) : mymuon->E();
			//double corrected_E = mymuon->E();
			TLorentzVector correctedMuon(corrected_Px, corrected_Py, corrected_Pz, corrected_E);

			*/

			// 2012 muon ID cut
			if(! ( (mymuon->pfIsoChargedHadronPt04() / correctedMuon.Pt()) < 0.2 ) )
			{	
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because pfIsoChargedHadronPt04 / Pt > 0.2" << endl;
				continue;
			}
			nbMuonsAfterID[2]++;
                        TOTALnbMuonsAfterID[2]++;

		 	if(! (correctedMuon.Pt() > 10.0) )
			{// transverse momentum
				//if(! (mymuon->Pt() > 10.0) ){// transverse momentum
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because transverse momentum" << endl;
				continue;
			}
			nbMuonsAfterID[3]++;
			TOTALnbMuonsAfterID[3]++;

			if(! (fabs(mymuon->Eta())<2.4) )
			{// |eta_muon|< 2.1
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because high eta (" << mymuon->Eta() << ")" << endl;
				continue;
			}
			nbMuonsAfterID[4]++;
			TOTALnbMuonsAfterID[4]++;

			//if(! (mymuon->) ){// 
			//muonIsNotCommissioned.push_back(1);
			//if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because" << endl;
			//continue;
			//}

			if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " accepted" << endl;
			muonIsNotCommissioned.push_back(0);
			muonIdentified_corrected_Pt.push_back(correctedMuon.Pt());
			muonIdentified_corrected_Pz.push_back(correctedMuon.Pz());
			muonIdentified_corrected_E.push_back(correctedMuon.E());
			muonIdentified.push_back(imuon);
			//mymuon->Clear();
		
			delete rmcor;
			rmcor = 0;
		}
		unsigned int NbMuonsIdentified = muonIdentified.size();
		
		// Increasing counter
		for(int i = 0; i < 5 ; i++)
		{
			if(nbMuonsAfterID[i] >= 2)
			{ 
				TOTALnbEventsAfterMuonID[i]++;
			}
		}

		if(! (nbMuonsAfterID[4] >=2) )// Not enough dimuon candidates, skip the event
		{
				continue;
		}


		// ********************************************************************************************************************
		// ***** DIMUON OBJECT SELECTION *****
		// ********************************************************************************************************************

		if(verbosity>1) cout << "Filling dimuon pair holder" << endl;


		// Making dimuon pairs holder
		int numberOfDimuons[3] = {0};
		numberOfDimuons[0] = factorial(nbMuonsAfterID[4] -1);

		if(! (numberOfDimuons[0] >= 1) )// Not enough dimuon candidates, skip the event. This cut is redundant with the previous one, should do nothing
		{
			continue;

		}

		TOTALnbDimuonsAfterID[0] += numberOfDimuons[0];
		TOTALnbEventsAfterDimuonID[0] += 1;

		pair <int, int> IDofMuons[3][numberOfDimuons[0]];
		pair <double, double> PtofMuons[3][numberOfDimuons[0]];
		pair <double, double> PzofMuons[3][numberOfDimuons[0]];
		pair <double, double> EofMuons[3][numberOfDimuons[0]];

		if(verbosity>2) cout << "initializing dimuon pair object" << endl;
		// Initializing pair object
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < numberOfDimuons[0]; j++) IDofMuons[i][j] = make_pair(0, 0);
			for(int j = 0; j < numberOfDimuons[0]; j++) PtofMuons[i][j] = make_pair(0.0, 0.0);
		}

		if(verbosity>2) cout << "Filling pair object for dimuon pairs composed of ID'ed muons" << endl;
		// Filling pair object for dimuon pairs composed of ID'ed muons
		int i_dimuons_ = 0;
		for(int muon_i = 0; muon_i < nbMuonsAfterID[4] ; muon_i++)
		{
			for(int muon_j = muon_i +1; muon_j < nbMuonsAfterID[4]; muon_j++)
			{
				IDofMuons[0][i_dimuons_] = make_pair(muonIdentified[muon_i], muonIdentified[muon_j]);
				PtofMuons[0][i_dimuons_] = make_pair(muonIdentified_corrected_Pt[muon_i], muonIdentified_corrected_Pt[muon_j]);
				PzofMuons[0][i_dimuons_] = make_pair(muonIdentified_corrected_Pz[muon_i], muonIdentified_corrected_Pz[muon_j]);
				EofMuons[0][i_dimuons_] = make_pair(muonIdentified_corrected_E[muon_i], muonIdentified_corrected_E[muon_j]);		
				if(verbosity>5) cerr << "muonIdentified_corrected_Pt["<<muon_i<<"]= " << muonIdentified_corrected_Pt[muon_i] << endl;
				if(verbosity>5) cerr << "muonIdentified_corrected_Pt["<<muon_j<<"]= " << muonIdentified_corrected_Pt[muon_j] << endl;
				if(verbosity>5) cerr << "PtofMuons[0]["<<i_dimuons_<<"]=	(" << PtofMuons[0][i_dimuons_].first << " , " << PtofMuons[0][i_dimuons_].second << ")" << endl;
				i_dimuons_++;
			}
		}

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
				PtofMuons[1][numberOfDimuons[1]] = make_pair(PtofMuons[0][i_dimuons].first, PtofMuons[0][i_dimuons].second);
				PzofMuons[1][numberOfDimuons[1]] = make_pair(PzofMuons[0][i_dimuons].first, PzofMuons[0][i_dimuons].second);
				EofMuons[1][numberOfDimuons[1]] = make_pair(EofMuons[0][i_dimuons].first, EofMuons[0][i_dimuons].second);
				numberOfDimuons[1] += 1;
			}
			//Muon1->Clear();
			//Muon2->Clear();
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
			double corrected_Pt1 = PtofMuons[1][i_dimuons].first;
			double corrected_Pt2 = PtofMuons[1][i_dimuons].second;
			double corrected_Pz1 = PzofMuons[1][i_dimuons].first;
			double corrected_Pz2 = PzofMuons[1][i_dimuons].second;
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? Muon1->Px() * (double)(corrected_Pt1) / (double)(Muon1->Pt()) : Muon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? Muon2->Px() * (double)(corrected_Pt2) / (double)(Muon2->Pt()) : Muon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? Muon1->Py() * (double)(corrected_Pt1) / (double)(Muon1->Pt()) : Muon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? Muon2->Py() * (double)(corrected_Pt2) / (double)(Muon2->Pt()) : Muon2->Py();
			double m_mu = 105.658367e-3;
			//double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : Muon1->E();
			//double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : Muon2->E();
			double corrected_E1 = EofMuons[1][i_dimuons].first;;
			double corrected_E2 = EofMuons[1][i_dimuons].second;
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);

			TLorentzVector mumu;
			//mumu = (*Muon1) + (*Muon2);
			mumu = (*correctedMuon1) + (*correctedMuon2);
			if( (low_m_mumu < mumu.M()) && (mumu.M() < high_m_mumu) )
			{
				IDofMuons[2][numberOfDimuons[2]] = make_pair(IDofMuons[1][i_dimuons].first, IDofMuons[1][i_dimuons].second);
				PtofMuons[2][numberOfDimuons[2]] = make_pair(PtofMuons[1][i_dimuons].first, PtofMuons[1][i_dimuons].second);
				PzofMuons[2][numberOfDimuons[2]] = make_pair(PzofMuons[1][i_dimuons].first, PzofMuons[1][i_dimuons].second);
				EofMuons[2][numberOfDimuons[2]] = make_pair(EofMuons[1][i_dimuons].first, EofMuons[1][i_dimuons].second);
				numberOfDimuons[2] += 1;
			}
			//Muon1->Clear();
			//Muon2->Clear();
			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;	
			//correctedMuon1->Delete();
			//correctedMuon2->Delete();
		}

		if(! (numberOfDimuons[2] >= 1) )// Not enough dimuon candidates, skip the event
		{
			 continue;
		}

		TOTALnbDimuonsAfterID[2] += numberOfDimuons[2];
		TOTALnbEventsAfterDimuonID[2] += 1;

		isMM = 1;
		for(int i_dimuons = 0 ; i_dimuons < numberOfDimuons[2]; i_dimuons++)
		{
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			mymuon1 = (TRootMuon*) muons->At( IDofMuons[2][i_dimuons].first );
			mymuon2 = (TRootMuon*) muons->At( IDofMuons[2][i_dimuons].second );
		
			double corrected_Pt1 = PtofMuons[2][i_dimuons].first;
                        double corrected_Pt2 = PtofMuons[2][i_dimuons].second;
                        double corrected_Pz1 = PzofMuons[2][i_dimuons].first;
                        double corrected_Pz2 = PzofMuons[2][i_dimuons].second;
                        double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
                        double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
                        double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
                        double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
                        double m_mu = 105.658367e-3;
                        //double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
                        //double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : Muon2->E();
                        double corrected_E1 = EofMuons[2][i_dimuons].first;;
                        double corrected_E2 = EofMuons[2][i_dimuons].second;
                        TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
                        TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);

                        TLorentzVector mumu;
                        mumu = (*correctedMuon1) + (*correctedMuon2);
	
			Ptmumu = mumu.Pt();
			//FillMM(mymuon1, mymuon2, doMC, applyMuonScaleCorrection, mcParticles, PTofMuons[2][i_dimuons].first, PTofMuons[2][i_dimuons].second);
                        
			int fsr = 0; 
			int NbPhotons = photons->GetEntries();
			
			for(int iphoton=0 ; iphoton<NbPhotons ; iphoton++)
                        {
                                TRootPhoton *myphoton;
                                myphoton = (TRootPhoton*) photons->At(iphoton);
				
				double phiPhoton = myphoton->Phi();
                        	double etaPhoton = myphoton->Eta();
                        	double phiMuon1 = mymuon1->Phi();
                        	double etaMuon1 = mymuon1->Eta();
                        	double phiMuon2 = mymuon2->Phi();
                        	double etaMuon2 = mymuon2->Eta();
                        	double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
                        	double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
                        	double min_DeltaR = min(deltaRphomu1, deltaRphomu2);

				double far_muonPt = (deltaRphomu1 > deltaRphomu2) ? correctedMuon1->Pt() : correctedMuon2->Pt();

				if( min_DeltaR <= 0.8 || myphoton->Pt() > 10)
				{
					fsr = 1;
					break;
				}
		

                        }

			if(fsr == 0) 
			{
				isMM_nonFSR = 1;
				FillMM(mymuon1, mymuon2, correctedMuon1, correctedMuon2, doMC, doR9Rescaling, mcParticles);
				miniTree->Fill();
			}
			else
			{
				isMM_nonFSR = 0;
                                FillMM(mymuon1, mymuon2, correctedMuon1, correctedMuon2, doMC, doR9Rescaling, mcParticles);
                                miniTree->Fill();
			}
			//mumu.Clear();
			//mymuon1->Clear();
			//mymuon2->Clear();
		
		}


		// --------------------------------------------------------------------------------------------------------------------
                // ----- Cleanning of FSR events -----
                // --------------------------------------------------------------------------------------------------------------------	




	} // fin boucle sur evts LOOP


	cout << "Nb_events_outside_powheg_cuts= " << Nb_events_outside_powheg_cuts << endl << endl;
	for(int i = 0; i < 5 ; i++)
	{
		cout << "TOTALnbMuonsAfterID["<<i<<"]=\t" << TOTALnbMuonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterMuonID["<<i<<"]=\t" << TOTALnbEventsAfterMuonID[i] << endl;
	}
	cout << endl;
	for(int i = 0; i < 3 ; i++)
	{
		cout << "TOTALnbDimuonsAfterID["<<i<<"]=\t" << TOTALnbDimuonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterDimuonID["<<i<<"]=\t" << TOTALnbEventsAfterDimuonID[i] << endl;
	}
	cout << endl;

	cout << "Writing stuff out" << endl;
	// Writing stuff out
	OutputRootFile->Write();
	//OutputRootFile->Close();


	// --- To check if everything is ok --- //
        system(Form("touch miniTree_muons_%d.done",ijob));

	cout << "Cleaning" << endl;

	// Destruct the TRandom3	
	delete generator;
	generator = 0;

	// Destruct the branches / related TClonesArray of the input TChains
	delete event;
	event = 0;
	delete runInfos;
	runInfos = 0;
	delete mcParticles;
	mcParticles = 0;
	delete mcMuMuGammaEvent;
	mcMuMuGammaEvent = 0;
	delete mcPhotons;
	mcPhotons = 0;
	delete vertices;
	vertices = 0;
	delete muons;
	muons = 0;
	delete photons;
	photons = 0;
	delete superClusters;
	superClusters = 0;
	delete clusters;
	clusters = 0;

	// Deleting input trees	
	delete inputEventTree;
	inputEventTree = 0;
	delete inputRunTree;
	inputRunTree = 0;

	// Destructing the tree!
	// So its a bit of a mess: what needs to be done is: 
	// - SetDirectory(0) to tell the TTree to forget about the TFile that owns it
	miniTree->SetDirectory(0);
	// - delete the TTree
	delete miniTree;
	miniTree = 0;
	// - close the TFile
	OutputRootFile->Close();
	// - delete the TFile
	delete OutputRootFile;
	OutputRootFile = 0;

	// Exit!
	return 0;

}

