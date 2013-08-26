#include "SkimSelection_miniTree.h"


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
		cerr << "arguments should be passed !! sample (outputname) (ijob) (isZgammaMC) (applyMuonScaleCorrection)" << endl;
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
	// Optional argument : ijob
	// ******************************************
	int ijob = -1;
  if( argc > 3 )
  {
    std::stringstream ss ( argv[3] );
    ss >> ijob;
  }

  // ******************************************
  // Optional argument : isZgammaMC (1: FSR -- 2: nonFSR -- 3: MC info)
  // ******************************************
	int isZgammaMC = 0;
//	int n = -1;
	if( argc > 4 )
	{
		std::stringstream ss ( argv[4] );
		ss >> isZgammaMC;
//		ss >> n;
	}

	// ******************************************
  // Optional argument : applyMuonScaleCorrection
  // ******************************************
  int applyMuonScaleCorrection = 0;
  if( argc > 5 )
  {
    std::stringstream ss ( argv[5] );
    ss >> applyMuonScaleCorrection;
  }
	TTimeStamp *time = new TTimeStamp();
	TRandom3* generator = new TRandom3(time->GetNanoSec());




//	int n = 18;
//	TProof * p = TProof::Open("ccaplmaster.in2p3.fr");
	gSystem->Load("libToto.so");
	bool doHLT										= false;
	bool doMC										 = (bool)(isZgammaMC >= 1);
	bool doJetMC									= false;
	bool doMETMC									= false;
	bool doPDFInfo								= false;
	bool doSignalMuMuGamma				= (bool)(isZgammaMC >= 1);
	bool doSignalTopTop					 = false;
	bool doPhotonConversionMC		 = (bool)(isZgammaMC >= 1);
	bool doBeamSpot							 = false;
	bool doPrimaryVertex					= true;
	bool doZeePrimaryVertex			 = false;
	bool doTrack									= false;
	bool doJet										= false;
	bool doMuon									 = true;
	bool doElectron							 = false;
	bool doPhoton								 = true;
	bool doCluster								= false;
	bool doPhotonConversion			 = true;
	bool doMET										= false;
	bool doBardak								 = false;
	bool doPhotonVertexCorrection = false;
	bool doPhotonIsolation				= false;

	//cout<<endl<<"doMC du debut = "<<doMC<<endl;
	//cout<<endl<<"doSignalMuMuGamma = "<<doSignalMuMuGamma<<endl;
	// DATASET	
	TChain *inputEventTree = new TChain("eventTree");
	TChain *inputRunTree = new TChain("runTree");

	inputEventTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
	inputRunTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
//		inputRunTree->Add("DYToMuMu_*root");
//		inputRunTree->Add("DYToMuMu_*root");
//	inputEventTree->Add(inputfile);
//	inputRunTree->Add(inputfile);


	TFile* OutputRootFile = new TFile(Form("../SkimmedTotoSamples/%s/%s_part%i.root", sample.c_str(), sample.c_str(), ijob), "RECREATE");
//	TFile* OutputRootFile = new TFile(Form("../SkimmedTotoSamples/%s/%s_part%i.root", sample.c_str(), sample.c_str(), n), "RECREATE");
//	TFile* OutputRootFile = new TFile(Form("%s.root", sample.c_str()), "RECREATE");
//	TFile* OutputRootFile = new TFile(Form("%s_part%i.root", sample.c_str(), n), "RECREATE");
	TTree *outputEventTree = inputEventTree->CloneTree(0);	
	TTree *outputRunTree = inputRunTree->CloneTree(0);	
	
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

  int* NumWantedHLTnames;

  string ListWantedHLTnames[1];
	ListWantedHLTnames[0] = "HLT_Mu11";
	
	
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

	vector<int> SelectedEvent_RunNumber;
	vector<int> SelectedEvent_LumiNumber;
	vector<int> SelectedEvent_EventNumber;
	vector<double> SelectedEvent_mumugammaInvMass;
	vector<double> SelectedEvent_Eta_gamma;
	vector<double> SelectedEvent_Eta_muonNear;
	vector<double> SelectedEvent_Eta_muonFar;
	vector<double> SelectedEvent_Et_gamma;
	vector<double> SelectedEvent_Pt_muonNear;
	vector<double> SelectedEvent_Pt_muonFar;
	vector<double> SelectedEvent_DeltaRNear;
	vector<double> SelectedEvent_DeltaRFar;
	vector<double> SelectedEvent_mumuInvMass;
	SelectedEvent_RunNumber.clear();
	SelectedEvent_LumiNumber.clear();
	SelectedEvent_EventNumber.clear();
	SelectedEvent_mumugammaInvMass.clear();
	SelectedEvent_Eta_gamma.clear();
	SelectedEvent_Eta_muonNear.clear();
	SelectedEvent_Eta_muonFar.clear();
	SelectedEvent_Et_gamma.clear();
	SelectedEvent_Pt_muonNear.clear();
	SelectedEvent_Pt_muonFar.clear();
	SelectedEvent_DeltaRNear.clear();
	SelectedEvent_DeltaRFar.clear();
	SelectedEvent_mumuInvMass.clear();

//	cout << "inputRunTree->GetEntries()= " << inputRunTree->GetEntries() << endl;
//	cout << "inputRunTree->GetCurrentFile()->GetName()= " << inputRunTree->GetCurrentFile()->GetName() << endl;
//	cout << "inputRunTree->GetMaxTreeSize()= " << inputRunTree->GetMaxTreeSize() << endl;
//	inputRunTree->GetEvent(0);
//	inputRunTree->GetEntry(0);
//	cout << "inputRunTree->GetCurrentFile()->GetName()= " << inputRunTree->GetCurrentFile()->GetName() << endl;
//	inputRunTree->~TTree();
//	inputRunTree->Clear();
//	delete inputRunTree;
//	free(inputRunTree);
  string lastFile = "";

//	double integratedLuminosity = 714.783728;
//	double integratedLuminosity = 1078.19387;
	double integratedLuminosity = 1420.38;
//  double XSectionDYToMuMu = 1300.0 * 1.2416;
  double XSectionDYToMuMu = 1626.0;
//  double XSectionTTJets = 94.0;
  double XSectionTTJets = 94.76;
  double XSectionWJetsToLNu = 7899.0;
  double XSectionQCDMu = 349988.0;

//  double InitialNumberDYToMuMu = 2148325.0;
  double InitialNumberDYToMuMu = 29743564.0;
//  double InitialNumberTTJets = 1089625.0;
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
  if( ijob != -1 )
  {
    NbEventsPerJob = 250000;
    NbEventsBegin = ijob * NbEventsPerJob;
    NbEventsEnd = min( (ijob + 1)* NbEventsPerJob - 1 , (int)NbEvents);
    NbEvents = NbEventsEnd - NbEventsBegin + 1;
    cout << "NbEventsBegin= " << NbEventsBegin << "\tNbEventsEnd= " << NbEventsEnd << "\tNbEventsPerJob= " << NbEventsPerJob << endl;
  }

  // LOOP over events
  for(unsigned int ievt=NbEventsBegin; ievt<NbEventsEnd; ievt++)
//	for(unsigned int ievt=0; ievt<NbEvents; ievt++)
//	for(unsigned int ievt=0; ievt<100; ievt++)
	{
		if(verbosity>4) cout << "analysing event ievt= " << ievt << endl;
		nBeforeAllCuts++;
		int nprint = (int)((double)NbEventsPerJob/(double)100.0);
		if( (ievt % nprint)==0 ){ cout<< ievt <<" events done over "<<NbEvents<<" ( "<<ceil((double)ievt/(double)NbEventsPerJob*100)<<" \% )"<<endl; }

		Int_t iEvent = ievt;
		inputEventTree->GetEvent(ievt);
//    if( lastFile == "" ){
//      lastFile = string(inputEventTree->GetCurrentFile()->GetName());
//      cout << ievt << "\t" << lastFile << endl;
//    }

		// ____________________________________________
		// Event information
		// ____________________________________________
		Int_t iEventID = event->eventId();
		Int_t iLumiID = event->luminosityBlock();
		Int_t iRunID = event->runId();
//		if(iRunID != 149291) continue;
		Int_t nVertices = vertices->GetEntries();
		Int_t nGenVertices = vertices->GetEntries();
		if( (isZgammaMC >= 1) )
		{
			Int_t nGenVertices = event->nInTimePUVertices();
		}
		Int_t isSignalApplied = signal;
		Int_t isStewApplied = stew;
		Int_t isZJetsApplied = zjet_veto;
		string sample(sample_char);

// ********************************************************************************************************************
// ***** MC-TRUTH SIGNAL MATCHING *****
// ********************************************************************************************************************

    bool MCphotons_from_muons_from_Z = false;
    bool MC_first_muon_in_phase_space = false;
    bool MC_second_muon_in_phase_space = false;
    bool MCsignal_in_phase_space = false;

    if( zjet_veto && (!powheg) ){
      // ****
      // First loop: look for a photon
      for( int iMCparticle = 0 ; iMCparticle < mcParticles->GetEntries() ; iMCparticle++ ){
        TRootMCParticle *mcParticleCandidate = (TRootMCParticle *)mcParticles->At(iMCparticle);
        if( (mcParticleCandidate->status()==1) && (mcParticleCandidate->type() == 22) ){ // if the particle is a true MC photon
          if( abs(mcParticleCandidate->motherType()) == 13 ){// if the true MC photon origins from a muon
//							cout << "mother: " << abs(mcParticleCandidate->motherType()) << "\tgranny: " << abs(mcParticleCandidate->grannyType()) << "\toldgranny: " << abs(mcParticleCandidate->oldgrannyType()) << endl;
            if( abs(mcParticleCandidate->oldgrannyType()) == 23 ){// photon coming from a muon coming from a Z
//            if( abs(mcParticleCandidate->grannyType()) == 23 ){// photon coming from a muon coming from a Z
              if( (mcParticleCandidate->Pt()>8.0) && (abs(mcParticleCandidate->Eta())<3.0) ){
                MCphotons_from_muons_from_Z = true;
              }
            }
          }
        } // end of origin of the photon
      }// end of loop over MC particles
      // ****
      // Second loop: look for muons
      for( int iMCparticle = 0 ; iMCparticle < mcParticles->GetEntries() ; iMCparticle++ ){
        TRootMCParticle *mcParticleCandidate = (TRootMCParticle *)mcParticles->At(iMCparticle);
        if( MCphotons_from_muons_from_Z == true ){ // if there is a photon coming from a muon coming from a Z and photon in MC phase space, THEN, look for muons in phase space
          if( (MC_first_muon_in_phase_space == false) && (mcParticleCandidate->status()==1) && (abs(mcParticleCandidate->type()) == 13) ){// if the particle is a final state muon
            if( abs(mcParticleCandidate->motherType()) == 13 ){
              if( abs(mcParticleCandidate->oldgrannyType()) == 23 ){// muon is coming from a Z
//              if( abs(mcParticleCandidate->grannyType()) == 23 ){// muon is coming from a Z
                if( (mcParticleCandidate->Pt()>8.0) && (abs(mcParticleCandidate->Eta())<3.0) ){
                  MC_first_muon_in_phase_space = true;
                }
              }
            }
          } // end of selecting first muon
          if( (MC_first_muon_in_phase_space == true) && (mcParticleCandidate->status()==1) && (abs(mcParticleCandidate->type()) == 13) ){// if the particle is a final state muon
            if( abs(mcParticleCandidate->motherType()) == 13 ){
              if( abs(mcParticleCandidate->oldgrannyType()) == 23 ){// muon is coming from a Z
//              if( abs(mcParticleCandidate->grannyType()) == 23 ){// muon is coming from a Z
                if( (mcParticleCandidate->Pt()>8.0) && (abs(mcParticleCandidate->Eta())<3.0) ){
                  MC_second_muon_in_phase_space = true;
                  MCsignal_in_phase_space = true;
                }
              }
            }
          }// end of selecting second muon
        }
      }// end of loop over MC particles
 			if( ((isZgammaMC == 1) && (!MCsignal_in_phase_space)) || ((isZgammaMC == 2) && (MCsignal_in_phase_space)) )
      {
//        cerr<<"SAFE: photon(s) coming from muon, aborting event " << ievt << endl;
        continue;
      }
      nAfterCutZJETVETO++;
     }// end of if Z+Jets veto for anything but powheg

    if( zjet_veto && powheg )
    {
			// POWHEG GENERATOR CUTS
			TRootParticle *mymuplus;
      mymuplus = (mcMuMuGammaEvent->muplus());
      TRootParticle *mymuminus;
      mymuminus = (mcMuMuGammaEvent->muminus());
      TRootParticle dimuon;
      dimuon = *mymuplus + *mymuminus;
      //cout << "dimuon.M()= " << dimuon.M() << endl;
/*
      if( (dimuon.M() < 20.0) || (dimuon.M() > 500.0) )
      {
        cout << "This event " << ievt << "is outside powheg range" << endl;
        Nb_events_outside_powheg_cuts++;
//				continue;
      }
*/


			// FSR-tagging
      if( mcMuMuGammaEvent->nZ() == 1 )
      {// consistency check
//        cerr << "There is " << mcMuMuGammaEvent->nFSR() << " fsr photons in the event" << endl;
        if( mcMuMuGammaEvent->nFSR() < 1 )
        {// if there is no fsr photon, don't bother
          MCsignal_in_phase_space = false;
        }
        else
        {// if there is a fsr photon, check further
          for( int imcphoton = 0 ; imcphoton < mcMuMuGammaEvent->nFSR() ; imcphoton++ )
          {// loop over mc photons
            TRootParticle *myphoton;
            myphoton = (mcMuMuGammaEvent->photonFSR(imcphoton));
            // mc-acceptance cuts on photons
            if( (myphoton->Pt()>8.0) && (abs(myphoton->Eta())<3.0) ) MCphotons_from_muons_from_Z = true;
          }// end of loop over mc photons
          if( MCphotons_from_muons_from_Z )
          { // if there is a fsr photon passing acceptance cuts, then look at muons coming from the Z
            if( (mcMuMuGammaEvent->muplus()->Pt()>8.0) && (abs((mcMuMuGammaEvent->muplus()->Eta())<3.0)) )
            {
              MC_first_muon_in_phase_space = true;
              if( (mcMuMuGammaEvent->muminus()->Pt()>8.0) && (abs((mcMuMuGammaEvent->muminus()->Eta())<3.0)) )
              {
                MC_second_muon_in_phase_space = true;
                MCsignal_in_phase_space = true;
              }
            }
          } // end of cuts on muons coming from the Z
        }// end of check if the event is a mc fsr passing acceptance cuts

				if( ((isZgammaMC == 1) && (!MCsignal_in_phase_space)) || ((isZgammaMC == 2) && (MCsignal_in_phase_space)) )
        {
//          cerr<<"SAFE: photon(s) coming from muon, aborting event " << ievt << endl;

          continue;
        }
        nAfterCutZJETVETO++;
      } else {
        cout << "Failed POWHEG mumugamma consistency check" << endl;
      }// end of consistency check
    }// end of if Z+Jets veto for powheg

// ********************************************************************************************************************
// ***** HLT REQUIREMENT *****
// ********************************************************************************************************************

		// HLT information		

//		if(doHLT){
//			if( ievt==0 ){ inputRunTree->GetEvent(ievt); NumWantedHLTnames = InitializeHLTinfo(inputRunTree, runInfos, event->nHLTPaths(), ListWantedHLTnames, 1);  }
//			if ( string(inputEventTree->GetCurrentFile()->GetName()) != lastFile ){
//				inputRunTree->GetEntry(inputEventTree->GetTreeNumber());
//     	  lastFile = string(inputEventTree->GetCurrentFile()->GetName());
//     	  cout << ievt << "\t" << lastFile << endl;
//     	 	NumWantedHLTnames = InitializeHLTinfo(inputRunTree, runInfos, event->nHLTPaths(), ListWantedHLTnames, 1);
//     	}
//      doHLTInfo(event, runInfos, NumWantedHLTnames, 1, &Muon_eventPassHLT_Mu11);
//		}
//
//
//		if (!((event->ptHat()>=minPtHat)&&(event->ptHat()<maxPtHat)))
//		{
//      cerr << "CUT: event " << ievt << " ( " << iRunID << " , " << iLumiID << " , " << iEventID << " )" << " CUT for pthat filtering" << endl;
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
		int nbMuonsAfterID[12] = {0};
		Int_t NbMuons = muons->GetEntries();
		
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
        muonIsNotCommissioned.push_back(1);
        if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because sum of pT of tracks with pT >1.5 within a cone of DR < 0.3 around the muon direction, vetoing a cone of 0.015 around that direction" << endl;
        continue;
      }
			nbMuonsAfterID[9]++;
			TOTALnbMuonsAfterID[9]++;

		double corrected_Pt = mymuon->Pt();
 		if( applyMuonScaleCorrection > 0 )
		{
			// Sidra makes MC look like data
			if( isZgammaMC > 0) corrected_Pt = applySidra(mymuon->Pt(), mymuon->charge(), mymuon->Eta(), mymuon->Phi(), generator);
			// MuScleFit correct data absolute scale
			corrected_Pt = applyMuScleFit(corrected_Pt, mymuon->charge(), mymuon->Eta(), mymuon->Phi());
		}
		double corrected_Pz = mymuon->Pz();
		double corrected_Px = applyMuonScaleCorrection > 0 ? mymuon->Px() * corrected_Pt / mymuon->Pt() : mymuon->Px();
		double corrected_Py = applyMuonScaleCorrection > 0 ? mymuon->Py() * corrected_Pt / mymuon->Pt() : mymuon->Py();
		double m_mu = 105.658367e-3;
		double corrected_E = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz * corrected_Pz + corrected_Pt * corrected_Pt) ) : mymuon->E();
//		double corrected_E = mymuon->E();
		TLorentzVector correctedMuon(corrected_Px, corrected_Py, corrected_Pz, corrected_E);


     if(! (correctedMuon.Pt() > 10.0) ){// transverse momentum
//     if(! (mymuon->Pt() > 10.0) ){// transverse momentum
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

//      if(! (mymuon->) ){// 
//        muonIsNotCommissioned.push_back(1);
//        if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because" << endl;
//        continue;
//      }

      if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " accepted" << endl;
      muonIsNotCommissioned.push_back(0);
      muonIdentified.push_back(imuon);
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

		if(verbosity>2) cout << "initializing dimuon pair object" << endl;
		// Initializing pair object
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < numberOfDimuons[0]; j++) IDofMuons[i][j] = make_pair(0, 0);
		}

		if(verbosity>2) cout << "Filling pair object for dimuon pairs composed of ID'ed muons" << endl;
		// Filling pair object for dimuon pairs composed of ID'ed muons
		for(int i_dimuons = 0; i_dimuons < numberOfDimuons[0]; i_dimuons++)
    {
			for(int muon_i = 0; muon_i < nbMuonsAfterID[11] ; muon_i++)
			{
				for(int muon_j = muon_i +1; muon_j < nbMuonsAfterID[11]; muon_j++)
				{
					IDofMuons[0][i_dimuons] = make_pair(muonIdentified[muon_i], muonIdentified[muon_j]);
				}
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
				numberOfDimuons[1] += 1;
				IDofMuons[1][i_dimuons] = make_pair(IDofMuons[0][i_dimuons].first, IDofMuons[0][i_dimuons].second);
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
			// Sidra makes MC look like data
			if( isZgammaMC > 0) corrected_Pt1 = applySidra(Muon1->Pt(), Muon1->charge(), Muon1->Eta(), Muon1->Phi(), generator);
			if( isZgammaMC > 0) corrected_Pt2 = applySidra(Muon2->Pt(), Muon2->charge(), Muon2->Eta(), Muon2->Phi(), generator);
			// MuScleFit correct data absolute scale
			corrected_Pt1 = applyMuScleFit(corrected_Pt1, Muon1->charge(), Muon1->Eta(), Muon1->Phi());
			corrected_Pt2 = applyMuScleFit(corrected_Pt2, Muon2->charge(), Muon2->Eta(), Muon2->Phi());
		}
		double corrected_Pz1 = Muon1->Pz();
		double corrected_Pz2 = Muon2->Pz();
		double corrected_Px1 = applyMuonScaleCorrection > 0 ? Muon1->Px() * corrected_Pt1 / Muon1->Pt() : Muon1->Px();
		double corrected_Px2 = applyMuonScaleCorrection > 0 ? Muon2->Px() * corrected_Pt2 / Muon2->Pt() : Muon2->Px();
		double corrected_Py1 = applyMuonScaleCorrection > 0 ? Muon1->Py() * corrected_Pt1 / Muon1->Pt() : Muon1->Py();
		double corrected_Py2 = applyMuonScaleCorrection > 0 ? Muon2->Py() * corrected_Pt2 / Muon2->Pt() : Muon2->Py();
		double m_mu = 105.658367e-3;
//		double corrected_E1 = Muon1->E();
//		double corrected_E2 = Muon2->E();
		double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1)) : Muon1->E();
		double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2)) : Muon2->E();
		TLorentzVector correctedMuon1(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
		TLorentzVector correctedMuon2(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);



			TLorentzVector mumu;
//			mumu = (*Muon1) + (*Muon2);
			mumu = (correctedMuon1) + (correctedMuon2);
//			if( (40.0 < mumu.M()) && (mumu.M() < 80.0) )
			if( (30.0 < mumu.M()) && (mumu.M() < 90.0) )
			{
				numberOfDimuons[2] += 1;
				IDofMuons[2][i_dimuons] = make_pair(IDofMuons[1][i_dimuons].first, IDofMuons[1][i_dimuons].second);
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


		nAfterCut3++;

		nAfterCut4++;

		nAfterCut5++;

		nAfterCut6++;

		nAfterCut7++;

		nAfterCut8++;
		nAfterCut9++;
		nAfterCut10++;

		outputEventTree->Fill();
		outputRunTree->Fill();

	} // fin boucle sur evts LOOP

		cout << "Nb_events_outside_powheg_cuts= " << Nb_events_outside_powheg_cuts << endl << endl;
		for(int i = 0; i < 12 ; i++)
		{
			cout << "TOTALnbMuonsAfterID["<<i<<"]= " << TOTALnbMuonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterMuonID["<<i<<"]= " << TOTALnbEventsAfterMuonID[i] << endl;
		}
		cout << endl;
		for(int i = 0; i < 3 ; i++)
		{
			cout << "TOTALnbDimuonsAfterID["<<i<<"]= " << TOTALnbDimuonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterDimuonID["<<i<<"]= " << TOTALnbEventsAfterDimuonID[i] << endl;
		}
		cout << endl;
		for(int i = 0; i < 6 ; i++)
    {
      cout << "TOTALnbPhotonsAfterID["<<i<<"]= " << TOTALnbPhotonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterPhotonID["<<i<<"]= " << TOTALnbEventsAfterPhotonID[i] << endl;
		}
		cout << endl;
		for(int i = 0; i < 8 ; i++)
		{
			cout << "TOTALnbMuMuGammaAfterID["<<i<<"]= " << TOTALnbMuMuGammaAfterID[i] << "\t\t" << "TOTALnbEventsAfterMuMuGammaID["<<i<<"]= " << TOTALnbEventsAfterMuMuGammaID[i] << endl;
			if(i == 6) cout << endl;
		}

// Writing stuff out
	outputEventTree->AutoSave();
	outputRunTree->AutoSave();
	OutputRootFile->Write();
	OutputRootFile->Close();

	delete OutputRootFile;
	delete inputEventTree;
	delete inputRunTree;



	return 0;

}

