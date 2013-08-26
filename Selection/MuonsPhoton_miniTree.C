#include "Selection_miniTree.h"
#include "rochcor_v4_new.C"

	// ____________________________________________
	// Event information
	// ____________________________________________
	ULong64_t iEvent, iEventID, iLumiID, iRunID;
	Int_t isSignalApplied, isStewApplied, isZJetsApplied;

	Int_t isBeforeAllCuts, isAfterCutPthatFilter, isAfterCutZJETVETO;
	Int_t isVeryLooseMMG, isLooseMMG, isMM, isTightMMG, isMMGCandidate;
	Int_t isAfterFSRCut1, isAfterFSRCut2, isAfterFSRCut3;
	Int_t isAfterFSRCut4, isMultipleCandidate, isAfterCut5, isAfterCut6, isAfterCut7, isAfterCut8, isAfterCut9, isAfterCut10;
	Int_t isSelected;
	
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
	// Photon variables
	// ___________________________________________
	Int_t NbPhotons;
	Float_t Pt_allPhotons, Eta_allPhotons, Phi_allPhotons, Cross_allPhotons;
	Int_t isEBorEE_allPhotons, isEB_allPhotons, isEE_allPhotons, isEEP_allPhotons, isEEM_allPhotons;
	Float_t Photon_Eta, Photon_Phi;
	Float_t Photon_Px, Photon_Py, Photon_Pz;
	Int_t Photon_isEBorEE, Photon_isEB, Photon_isEE, Photon_isEEP, Photon_isEEM;

	Int_t Photon_hasPixelSeed, Photon_isAlsoElectron, Photon_Nclusters, Photon_nBasicClusters, Photon_nXtals;
	//Int_t Photon_isTightPhoton, Photon_isLoosePhoton;
	Int_t Photon_convNTracks, Photon_isConverted;
	Float_t Photon_convEoverP, Photon_convMass, Photon_convCotanTheta, Photon_convLikely, Photon_convVertexX, Photon_convVertexY, Photon_convVertexZ;
	Float_t Photon_E, Photon_Et, Photon_E2x2, Photon_E3x3, Photon_E5x5, Photon_Emax, Photon_E2nd;
	Float_t Photon_E_regression, Photon_E_regressionError, Photon_Et_regression;
	Float_t Photon_Ecorr_o_Ereco;
	Float_t Photon_r19, Photon_r9, Photon_cross;
	//Float_t Photon_caloConeSize, 
	Float_t Photon_PreshEnergy, Photon_HoE;
	Float_t Photon_sigmaEtaEta, Photon_sigmaIetaIeta;
	Float_t Photon_covEtaEta, Photon_covPhiPhi, Photon_covEtaPhi;
	Float_t Photon_etaWidth, Photon_phiWidth;
	Float_t Photon_dR03isoEcalRecHit, Photon_dR03isoHcalRecHit, Photon_dR03isoSolidTrkCone, Photon_dR03isoHollowTrkCone, Photon_dR03isoNTracksSolidCone, Photon_dR03isoNTracksHollowCone;
	Float_t Photon_dR04isoEcalRecHit, Photon_dR04isoHcalRecHit, Photon_dR04isoSolidTrkCone, Photon_dR04isoHollowTrkCone, Photon_dR04isoNTracksSolidCone, Photon_dR04isoNTracksHollowCone;
	Float_t Photon_seedTime, Photon_seedFlag;
	Int_t Photon_seedPosition1, Photon_seedPosition2;
	Float_t Photon_SC_Eta, Photon_SC_Phi, Photon_SC_brem;
	Float_t Photon_SC_E, Photon_SC_Et, Photon_SC_rawE, Photon_SC_rawEt;
	Float_t Photon_lambdaRatio, Photon_ratioSeed, Photon_ratioS4, Photon_lamdbaDivCov;
	Float_t Photon_ratioS4_corrected;
	Float_t Photon_SC_rawE_x_fEta;
	Float_t Photon_SC_rawE_x_fEta_x_fBrem, Photon_SC_rawE_x_fEta_x_fBrem_AF, Photon_SC_rawE_x_fEta_x_fBrem_L, Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta, Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta, Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta;
	Float_t Photon_secondMomentMaj, Photon_secondMomentMin, Photon_secondMomentAlpha;
	Float_t Photon_etaLAT, Photon_phiLAT, Photon_LAT, Photon_Zernike20, Photon_Zernike42, Photon_ESratio;

	// ____________________________________________
	// mugamma / mumu / mumugamma information
	// ____________________________________________
	Int_t iCandidate[8];
	Int_t iCandidate_temp[8][50];
	Int_t nCandidate[8];

	Float_t Mmumu, Mmumugamma, Mmumugamma_5x5, Mmumugamma_SC, Mmumugamma_SCraw, Mmumugamma_SCraw_fEta;
	Float_t Ptmumu;
	Float_t deltaRNear, deltaRFar, deltaRMinus, deltaRPlus, deltaRLeading, deltaRSubleading;
	Float_t mmg_k, mmg_ik, mmg_s, mmg_logk, mmg_logik, mmg_logs;
	Float_t mmg_k_5x5, mmg_ik_5x5, mmg_s_5x5, mmg_logk_5x5, mmg_logik_5x5, mmg_logs_5x5;
	Float_t mmg_k_SC, mmg_ik_SC, mmg_s_SC, mmg_logk_SC, mmg_logik_SC, mmg_logs_SC;
	Float_t mmg_k_SCraw, mmg_ik_SCraw, mmg_s_SCraw, mmg_logk_SCraw, mmg_logik_SCraw, mmg_logs_SCraw;
	Float_t mmg_k_SCraw_fEta, mmg_ik_SCraw_fEta, mmg_s_SCraw_fEta, mmg_logk_SCraw_fEta, mmg_logik_SCraw_fEta, mmg_logs_SCraw_fEta;
	Float_t Mmumugamma_SCraw_fEta_fBrem, Mmumugamma_SCraw_fEta_fBrem_AF, Mmumugamma_SCraw_fEta_fBrem_L, Mmumugamma_SCraw_fEta_fBrem_fEtEta, Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta, Mmumugamma_SCraw_fEta_fBrem_L_fEtEta;
	Float_t mmg_ik_SCraw_fEta_fBrem, mmg_ik_SCraw_fEta_fBrem_AF, mmg_ik_SCraw_fEta_fBrem_L, mmg_ik_SCraw_fEta_fBrem_fEtEta, mmg_ik_SCraw_fEta_fBrem_AF_fEtEta, mmg_ik_SCraw_fEta_fBrem_L_fEtEta; 

// MuonBeforeBrem = (MuonN+Photon) || MuonF
// (M minus charge, P plus charge), (F far, N near), (L leading, S subleading)
	Float_t MuonBeforeBremM_Pt, MuonBeforeBremP_Pt, MuonBeforeBremN_Pt, MuonBeforeBremF_Pt, MuonBeforeBremL_Pt, MuonBeforeBremS_Pt;
	Float_t MuonBeforeBremM_Eta, MuonBeforeBremP_Eta, MuonBeforeBremN_Eta, MuonBeforeBremF_Eta, MuonBeforeBremL_Eta, MuonBeforeBremS_Eta;
	Float_t MuonBeforeBremM_Phi, MuonBeforeBremP_Phi, MuonBeforeBremN_Phi, MuonBeforeBremF_Phi, MuonBeforeBremL_Phi, MuonBeforeBremS_Phi;
	Float_t MuonBeforeBremM_E, MuonBeforeBremP_E, MuonBeforeBremN_E, MuonBeforeBremF_E, MuonBeforeBremL_E, MuonBeforeBremS_E;
	Float_t MuonBeforeBremM_Px, MuonBeforeBremP_Px, MuonBeforeBremN_Px, MuonBeforeBremF_Px, MuonBeforeBremL_Px, MuonBeforeBremS_Px;
	Float_t MuonBeforeBremM_Py, MuonBeforeBremP_Py, MuonBeforeBremN_Py, MuonBeforeBremF_Py, MuonBeforeBremL_Py, MuonBeforeBremS_Py;
	Float_t MuonBeforeBremM_Pz, MuonBeforeBremP_Pz, MuonBeforeBremN_Pz, MuonBeforeBremF_Pz, MuonBeforeBremL_Pz, MuonBeforeBremS_Pz;
	Int_t MuonBeforeBremF_Charge, MuonBeforeBremN_Charge, MuonBeforeBremL_Charge, MuonBeforeBremS_Charge;



	// ____________________________________________
	// Neural Network variables
	// ____________________________________________

	Float_t Photon_NNshapeOutput;

	// ____________________________________________
	// Surface variables
	// ____________________________________________

				Float_t MZ_Surface;
	Float_t mmg_k_MZ_Surface, mmg_ik_MZ_Surface, mmg_s_MZ_Surface, mmg_logk_MZ_Surface, mmg_logik_MZ_Surface, mmg_logs_MZ_Surface;

	// ____________________________________________
	// MC Truth
	// ___________________________________________

	Float_t Photon_MC_E, Photon_MC_Px, Photon_MC_Py, Photon_MC_Pz, Photon_MC_Phi, Photon_MC_Eta, Photon_MC_Pt;
	Int_t Photon_MCisConverted;
	Float_t Photon_MCconvEoverP, Photon_MCconvMass, Photon_MCconvCotanTheta, Photon_MCconvVertexX, Photon_MCconvVertexY, Photon_MCconvVertexZ;
	Float_t MuonM_MC_E, MuonM_MC_Px, MuonM_MC_Py, MuonM_MC_Pz, MuonM_MC_Phi, MuonM_MC_Eta, MuonM_MC_Pt;
	Float_t MuonP_MC_E, MuonP_MC_Px, MuonP_MC_Py, MuonP_MC_Pz, MuonP_MC_Phi, MuonP_MC_Eta, MuonP_MC_Pt;
	Float_t MuonN_MC_E, MuonN_MC_Px, MuonN_MC_Py, MuonN_MC_Pz, MuonN_MC_Phi, MuonN_MC_Eta, MuonN_MC_Pt;
	Float_t MuonF_MC_E, MuonF_MC_Px, MuonF_MC_Py, MuonF_MC_Pz, MuonF_MC_Phi, MuonF_MC_Eta, MuonF_MC_Pt;
	Float_t MuonL_MC_E, MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_Phi, MuonL_MC_Eta, MuonL_MC_Pt;
	Float_t MuonS_MC_E, MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_Phi, MuonS_MC_Eta, MuonS_MC_Pt;
	Float_t Photon_SC_rawE_x_fEta_o_MC_E, Photon_E_o_MC_E;
	Float_t Photon_SC_rawE_x_fEta_x_fBrem_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_AF_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_L_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta_o_MC_E;	

	Float_t Mmumu_Photon_MC, Mmumugamma_Photon_MC, mmg_k_Photon_MC, mmg_ik_Photon_MC, mmg_s_Photon_MC, mmg_logk_Photon_MC, mmg_logik_Photon_MC, mmg_logs_Photon_MC;
	Float_t Mmumu_Muons_MC, Mmumugamma_Muons_MC, mmg_k_Muons_MC, mmg_ik_Muons_MC, mmg_s_Muons_MC, mmg_logk_Muons_MC, mmg_logik_Muons_MC, mmg_logs_Muons_MC;
	Float_t Mmumu_MMG_MC, Mmumugamma_MMG_MC, mmg_k_MMG_MC, mmg_ik_MMG_MC, mmg_s_MMG_MC, mmg_logk_MMG_MC, mmg_logik_MMG_MC, mmg_logs_MMG_MC;

	Float_t mmg_k_MZ, mmg_ik_MZ, mmg_s_MZ, mmg_logk_MZ, mmg_logik_MZ, mmg_logs_MZ;
	Float_t mmg_k_MZ_Photon_MC, mmg_ik_MZ_Photon_MC, mmg_s_MZ_Photon_MC, mmg_logk_MZ_Photon_MC, mmg_logik_MZ_Photon_MC, mmg_logs_MZ_Photon_MC;
	Float_t mmg_k_MZ_Muons_MC, mmg_ik_MZ_Muons_MC, mmg_s_MZ_Muons_MC, mmg_logk_MZ_Muons_MC, mmg_logik_MZ_Muons_MC, mmg_logs_MZ_Muons_MC;

	Float_t mmg_k_MZ_Muons_RECO_MC, mmg_ik_MZ_Muons_RECO_MC, mmg_s_MZ_Muons_RECO_MC, mmg_logk_MZ_Muons_RECO_MC, mmg_logik_MZ_Muons_RECO_MC, mmg_logs_MZ_Muons_RECO_MC;
	

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
		cerr << "arguments should be passed !! sample (outputname) (ntotjob) (ijob) (isZgammaMC) (lumi_set) (pu_set) (low m_mumu cut) (high m_mumu cut) (photon energy correction scheme) (extra photon scale) (applyMuonScaleCorrection) (muon correction sys) (extra resolution) (itoy for mu sys)" << endl;
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

	//cout<<endl<<"doMC du debut = "<<doMC<<endl;
	//cout<<endl<<"doSignalMuMuGamma = "<<doSignalMuMuGamma<<endl;
	// DATASET	
	TChain *inputEventTree = new TChain("eventTree");
	TChain *inputRunTree = new TChain("runTree");

	int MmumuLigneNumber = 0;	
				string sMmumu;
				ifstream FileMmumu("Mmumu.txt");
				while(std::getline(FileMmumu, sMmumu))
				{		
								MmumuLigneNumber++;
				}

	int binNumber = sqrt(MmumuLigneNumber);

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
	//inputEventTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
	//inputRunTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/TotoSamples/%s/%s*root", sample_char, sample_char));
	inputEventTree->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/UserCode/IpnTreeProducer/test/%s/*root", sample_char));
				inputRunTree->Add(Form("/sps/cms/sgandurr/CMSSW_4_2_8_patch7/src/UserCode/IpnTreeProducer/test/%s/*root", sample_char));
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

//	inputEventTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/SkimmedTotoSamples/%s/%s*root", sample_char, sample_char));
//	inputRunTree->Add(Form("/sps/cms/obondu/CMSSW_4_2_8__RECO_4_2_8_v2/src/Zmumugamma/SkimmedTotoSamples/%s/%s*root", sample_char, sample_char));


//		inputEventTree->Add("miniTree_TEST.root");
//		inputRunTree->Add("miniTree_TEST.root");
//	inputEventTree->Add(inputfile);
//	inputRunTree->Add(inputfile);
//	inputEventTree->Add("root://ccxroot.in2p3.fr:1999//hpss/in2p3.fr/group/cms/users/obondu/Zgamma/Fall11/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_part*.root");
//	inputRunTree->Add("root://ccxroot.in2p3.fr:1999//hpss/in2p3.fr/group/cms/users/obondu/Zgamma/Fall11/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_part*.root");
/*
	inputEventTree->Add("root://ccxroot.in2p3.fr:1999//hpss/in2p3.fr/group/cms/users/obondu/Zgamma/Fall11/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_part0.root");
	inputRunTree->Add("root://ccxroot.in2p3.fr:1999//hpss/in2p3.fr/group/cms/users/obondu/Zgamma/Fall11/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6/Skimv02_DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_PU_S6_part0.root");
*/

// INSERTFILES

	TFile* OutputRootFile = new TFile(Form("miniTreeMuonsPhoton_%s_part%i.root", sample.c_str(), ijob), "RECREATE");
	TFile* OutputFriendFile = new TFile(Form("miniFriendMuonsPhoton_%s_part%i.root", sample.c_str(), ijob), "RECREATE");

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

	TTree* miniTree = new TTree("miniTree","Mu Mu Gamma informations");
	TTree* miniTree_allmuons = new TTree("miniTree_allmuons","all muons informations");
	TTree* miniTree_allphotons = new TTree("miniTree_allphotons","all photons informations");
	OutputFriendFile->cd();
	TTree* miniFriend = new TTree("miniFriend","Mu Mu Gamma *bonus* informations");
	OutputRootFile->cd();
//	TTree *outputEventTree = inputEventTree->CloneTree(0);

	// ____________________________________________
	// Event information
	// ____________________________________________
	
	miniTree->Branch("iEvent", &iEvent, "iEvent/l");
	miniTree->Branch("iEventID", &iEventID, "iEventID/l");
	miniTree->Branch("iLumiID", &iLumiID, "iLumiID/l");
	miniTree->Branch("iRunID", &iRunID, "iRunID/l");

	miniTree->Branch("isSignalApplied", &isSignalApplied, "isSignalApplied/I");
	miniTree->Branch("isStewApplied", &isStewApplied, "isStewApplied/I");
	miniTree->Branch("isZJetsApplied", &isZJetsApplied, "isZJetsApplied/I");

	miniTree->Branch("isBeforeAllCuts", &isBeforeAllCuts, "isBeforeAllCuts/I");
	miniTree->Branch("isAfterCutPthatFilter", &isAfterCutPthatFilter, "isAfterCutPthatFilter/I");
	miniTree->Branch("isAfterCutZJETVETO", &isAfterCutZJETVETO, "isAfterCutZJETVETO/I");

	miniTree->Branch("isVeryLooseMMG", &isVeryLooseMMG, "isVeryLooseMMG/I");
	miniTree->Branch("isLooseMMG", &isLooseMMG, "isLooseMMG/I");
	miniTree->Branch("isMM", &isMM, "isMM/I");
	miniTree->Branch("isTightMMG", &isTightMMG, "isTightMMG/I");
	miniTree->Branch("isMMGCandidate", &isMMGCandidate, "isMMGCandidate/I");

	miniTree->Branch("isAfterFSRCut1", &isAfterFSRCut1, "isAfterFSRCut1/I");
	miniTree->Branch("isAfterFSRCut2", &isAfterFSRCut2, "isAfterFSRCut2/I");
	miniTree->Branch("isAfterFSRCut3", &isAfterFSRCut3, "isAfterFSRCut3/I");

	miniTree->Branch("isAfterFSRCut4", &isAfterFSRCut4, "isAfterFSRCut4/I");
	miniTree->Branch("isMultipleCandidate", &isMultipleCandidate, "isMultipleCandidate/I");
	miniTree->Branch("isAfterCut5", &isAfterCut5, "isAfterCut5/I");
	miniTree->Branch("isAfterCut6", &isAfterCut6, "isAfterCut6/I");
	miniTree->Branch("isAfterCut7", &isAfterCut7, "isAfterCut7/I");
	miniTree->Branch("isAfterCut8", &isAfterCut8, "isAfterCut8/I");
	miniTree->Branch("isAfterCut9", &isAfterCut9, "isAfterCut9/I");
	miniTree->Branch("isAfterCut10", &isAfterCut10, "isAfterCut10/I");

	miniTree->Branch("isSelected", &isSelected, "isSelected/I");
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

	miniTree_allmuons->Branch("iEvent", &iEvent, "iEvent/l");
	miniTree_allmuons->Branch("iEventID", &iEventID, "iEventID/l");
	miniTree_allmuons->Branch("iLumiID", &iLumiID, "iLumiID/l");
	miniTree_allmuons->Branch("iRunID", &iRunID, "iRunID/l");

	miniTree_allmuons->Branch("isSignalApplied", &isSignalApplied, "isSignalApplied/I");
	miniTree_allmuons->Branch("isStewApplied", &isStewApplied, "isStewApplied/I");
	miniTree_allmuons->Branch("isZJetsApplied", &isZJetsApplied, "isZJetsApplied/I");

	miniTree_allmuons->Branch("isBeforeAllCuts", &isBeforeAllCuts, "isBeforeAllCuts/I");
	miniTree_allmuons->Branch("isAfterCutPthatFilter", &isAfterCutPthatFilter, "isAfterCutPthatFilter/I");
	miniTree_allmuons->Branch("isAfterCutZJETVETO", &isAfterCutZJETVETO, "isAfterCutZJETVETO/I");

	miniTree_allmuons->Branch("isVeryLooseMMG", &isVeryLooseMMG, "isVeryLooseMMG/I");
	miniTree_allmuons->Branch("isLooseMMG", &isLooseMMG, "isLooseMMG/I");
	miniTree_allmuons->Branch("isMM", &isMM, "isMM/I");
	miniTree_allmuons->Branch("isTightMMG", &isTightMMG, "isTightMMG/I");
	miniTree_allmuons->Branch("isMMGCandidate", &isMMGCandidate, "isMMGCandidate/I");

	miniTree_allmuons->Branch("isAfterFSRCut1", &isAfterFSRCut1, "isAfterFSRCut1/I");
	miniTree_allmuons->Branch("isAfterFSRCut2", &isAfterFSRCut2, "isAfterFSRCut2/I");
	miniTree_allmuons->Branch("isAfterFSRCut3", &isAfterFSRCut3, "isAfterFSRCut3/I");

	miniTree_allmuons->Branch("isAfterFSRCut4", &isAfterFSRCut4, "isAfterFSRCut4/I");
	miniTree_allmuons->Branch("isMultipleCandidate", &isMultipleCandidate, "isMultipleCandidate/I");
	miniTree_allmuons->Branch("isAfterCut5", &isAfterCut5, "isAfterCut5/I");
	miniTree_allmuons->Branch("isAfterCut6", &isAfterCut6, "isAfterCut6/I");
	miniTree_allmuons->Branch("isAfterCut7", &isAfterCut7, "isAfterCut7/I");
	miniTree_allmuons->Branch("isAfterCut8", &isAfterCut8, "isAfterCut8/I");
	miniTree_allmuons->Branch("isAfterCut9", &isAfterCut9, "isAfterCut9/I");
	miniTree_allmuons->Branch("isAfterCut10", &isAfterCut10, "isAfterCut10/I");

	miniTree_allmuons->Branch("isSelected", &isSelected, "isSelected/I");

	miniTree_allphotons->Branch("iEvent", &iEvent, "iEvent/l");
	miniTree_allphotons->Branch("iEventID", &iEventID, "iEventID/l");
	miniTree_allphotons->Branch("iLumiID", &iLumiID, "iLumiID/l");
	miniTree_allphotons->Branch("iRunID", &iRunID, "iRunID/l");

	miniTree_allphotons->Branch("isSignalApplied", &isSignalApplied, "isSignalApplied/I");
	miniTree_allphotons->Branch("isStewApplied", &isStewApplied, "isStewApplied/I");
	miniTree_allphotons->Branch("isZJetsApplied", &isZJetsApplied, "isZJetsApplied/I");

	miniTree_allphotons->Branch("isBeforeAllCuts", &isBeforeAllCuts, "isBeforeAllCuts/I");
	miniTree_allphotons->Branch("isAfterCutPthatFilter", &isAfterCutPthatFilter, "isAfterCutPthatFilter/I");
	miniTree_allphotons->Branch("isAfterCutZJETVETO", &isAfterCutZJETVETO, "isAfterCutZJETVETO/I");

	miniTree_allphotons->Branch("isVeryLooseMMG", &isVeryLooseMMG, "isVeryLooseMMG/I");
	miniTree_allphotons->Branch("isLooseMMG", &isLooseMMG, "isLooseMMG/I");
	miniTree_allphotons->Branch("isMM", &isMM, "isMM/I");
	miniTree_allphotons->Branch("isTightMMG", &isTightMMG, "isTightMMG/I");
	miniTree_allphotons->Branch("isMMGCandidate", &isMMGCandidate, "isMMGCandidate/I");

	miniTree_allphotons->Branch("isAfterFSRCut1", &isAfterFSRCut1, "isAfterFSRCut1/I");
	miniTree_allphotons->Branch("isAfterFSRCut2", &isAfterFSRCut2, "isAfterFSRCut2/I");
	miniTree_allphotons->Branch("isAfterFSRCut3", &isAfterFSRCut3, "isAfterFSRCut3/I");

	miniTree_allphotons->Branch("isAfterFSRCut4", &isAfterFSRCut4, "isAfterFSRCut4/I");
	miniTree_allphotons->Branch("isMultipleCandidate", &isMultipleCandidate, "isMultipleCandidate/I");
	miniTree_allphotons->Branch("isAfterCut5", &isAfterCut5, "isAfterCut5/I");
	miniTree_allphotons->Branch("isAfterCut6", &isAfterCut6, "isAfterCut6/I");
	miniTree_allphotons->Branch("isAfterCut7", &isAfterCut7, "isAfterCut7/I");
	miniTree_allphotons->Branch("isAfterCut8", &isAfterCut8, "isAfterCut8/I");
	miniTree_allphotons->Branch("isAfterCut9", &isAfterCut9, "isAfterCut9/I");
	miniTree_allphotons->Branch("isAfterCut10", &isAfterCut10, "isAfterCut10/I");

	miniTree_allphotons->Branch("isSelected", &isSelected, "isSelected/I");
	miniTree_allphotons->Branch("isNotCommissionned", &isNotCommissionned, "isNotCommissionned/I");

	// ____________________________________________
	// Muon variables
	// ____________________________________________

	miniTree->Branch("NbMuons", &NbMuons, "NbMuons/I");

	miniTree_allmuons->Branch("Pt_allMuons", &Pt_allMuons, "Pt_allMuons/F");
	miniTree_allmuons->Branch("Eta_allMuons", &Eta_allMuons, "Eta_allMuons/F");
	miniTree_allmuons->Branch("Phi_allMuons", &Phi_allMuons, "Phi_allMuons/F");
	miniTree_allmuons->Branch("Charge_allMuons", &Charge_allMuons, "Charge_allMuons/F");

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
	// Photon variables
	// ___________________________________________

	miniTree->Branch("NbPhotons", &NbPhotons, "NbPhotons/I");

	miniTree_allphotons->Branch("Pt_allPhotons", &Pt_allPhotons, "Pt_allPhotons/F");
	miniTree_allphotons->Branch("Eta_allPhotons", &Eta_allPhotons, "Eta_allPhotons/F");
	miniTree_allphotons->Branch("Phi_allPhotons", &Phi_allPhotons, "Phi_allPhotons/F");
	miniTree_allphotons->Branch("Cross_allPhotons", &Cross_allPhotons, "Cross_allPhotons/F");
	
	miniTree_allphotons->Branch("isEBorEE_allPhotons", &isEBorEE_allPhotons, "isEBorEE_allPhotons/I");
	miniTree_allphotons->Branch("isEB_allPhotons", &isEB_allPhotons, "isEB_allPhotons/I");
	miniTree_allphotons->Branch("isEE_allPhotons", &isEE_allPhotons, "isEE_allPhotons/I");
	miniTree_allphotons->Branch("isEEM_allPhotons", &isEEM_allPhotons, "isEEM_allPhotons/I");
	miniTree_allphotons->Branch("isEEP_allPhotons", &isEEP_allPhotons, "isEEP_allPhotons/I");

	miniTree->Branch("Photon_Eta", &Photon_Eta, "Photon_Eta/F");
	miniTree->Branch("Photon_Phi", &Photon_Phi, "Photon_Phi/F");
	miniTree->Branch("Photon_Px", &Photon_Px, "Photon_Px/F");
	miniTree->Branch("Photon_Py", &Photon_Py, "Photon_Py/F");
	miniTree->Branch("Photon_Pz", &Photon_Pz, "Photon_Pz/F");

	miniTree->Branch("Photon_isEBorEE", &Photon_isEBorEE, "Photon_isEBorEE/I");
	miniTree->Branch("Photon_isEB", &Photon_isEB, "Photon_isEB/I");
	miniTree->Branch("Photon_isEE", &Photon_isEE, "Photon_isEE/I");
	miniTree->Branch("Photon_isEEP", &Photon_isEEP, "Photon_isEEP/I");
	miniTree->Branch("Photon_isEEM", &Photon_isEEM, "Photon_isEEM/I");

	miniTree->Branch("Photon_hasPixelSeed", &Photon_hasPixelSeed, "Photon_hasPixelSeed/I");
	miniTree->Branch("Photon_isAlsoElectron", &Photon_isAlsoElectron, "Photon_isAlsoElectron/I");
	miniTree->Branch("Photon_Nclusters", &Photon_Nclusters, "Photon_Nclusters/I");
	miniTree->Branch("Photon_nBasicClusters", &Photon_nBasicClusters, "Photon_nBasicClusters/I");
	miniTree->Branch("Photon_nXtals", &Photon_nXtals, "Photon_nXtals/I"); // Variable not filled in current version

	//miniTree->Branch("Photon_isTightPhoton", &Photon_isTightPhoton, "Photon_isTightPhoton/I");
	//miniTree->Branch("Photon_isLoosePhoton", &Photon_isLoosePhoton, "Photon_isLoosePhoton/I");

	miniTree->Branch("Photon_E", &Photon_E, "Photon_E/F");
	miniTree->Branch("Photon_Et", &Photon_Et, "Photon_Et/F");
	miniTree->Branch("Photon_E2x2", &Photon_E2x2, "Photon_E2x2/F");
	miniTree->Branch("Photon_E3x3", &Photon_E3x3, "Photon_E3x3/F");
	miniTree->Branch("Photon_E5x5", &Photon_E5x5, "Photon_E5x5/F");
	miniTree->Branch("Photon_Emax", &Photon_Emax, "Photon_Emax/F");
	miniTree->Branch("Photon_E2nd", &Photon_E2nd, "Photon_E2nd/F");
	miniTree->Branch("Photon_Ecorr_o_Ereco", &Photon_Ecorr_o_Ereco, "Photon_Ecorr_o_Ereco/F");

	miniTree->Branch("Photon_E_regression", &Photon_E_regression, "Photon_E_regression/F");
	miniTree->Branch("Photon_E_regressionError", &Photon_E_regressionError, "Photon_E_regressionError/F");
	miniTree->Branch("Photon_Et_regression", &Photon_Et_regression, "Photon_Et_regression/F");

	miniTree->Branch("Photon_r19", &Photon_r19, "Photon_r19/F");
	miniTree->Branch("Photon_r9", &Photon_r9, "Photon_r9/F");
	miniTree->Branch("Photon_cross", &Photon_cross, "Photon_cross/F");

	//miniTree->Branch("Photon_caloConeSize", &Photon_caloConeSize, "Photon_caloConeSize/F");
	miniTree->Branch("Photon_PreshEnergy", &Photon_PreshEnergy, "Photon_PreshEnergy/F");
	miniTree->Branch("Photon_HoE", &Photon_HoE, "Photon_HoE/F");
	miniTree->Branch("Photon_sigmaEtaEta", &Photon_sigmaEtaEta, "Photon_sigmaEtaEta/F");
	miniTree->Branch("Photon_sigmaIetaIeta", &Photon_sigmaIetaIeta, "Photon_sigmaIetaIeta/F");
	miniTree->Branch("Photon_covEtaEta", &Photon_covEtaEta, "Photon_covEtaEta/F");
	miniTree->Branch("Photon_covPhiPhi", &Photon_covPhiPhi, "Photon_covPhiPhi/F");
	miniTree->Branch("Photon_covEtaPhi", &Photon_covEtaPhi, "Photon_covEtaPhi/F");
	miniTree->Branch("Photon_convNTracks", &Photon_convNTracks, "Photon_convNTracks/I");
	miniTree->Branch("Photon_isConverted", &Photon_isConverted, "Photon_isConverted/I");
	miniTree->Branch("Photon_convEoverP", &Photon_convEoverP, "Photon_convEoverP/F");
	miniTree->Branch("Photon_convMass", &Photon_convMass, "Photon_convMass/F");
	miniTree->Branch("Photon_convCotanTheta", &Photon_convCotanTheta, "Photon_convCotanTheta/F");
	miniTree->Branch("Photon_convLikely", &Photon_convLikely, "Photon_convLikely/F");
	miniTree->Branch("Photon_convVertexX", &Photon_convVertexX, "Photon_convVertexX/F");
	miniTree->Branch("Photon_convVertexY", &Photon_convVertexY, "Photon_convVertexY/F");
	miniTree->Branch("Photon_convVertexZ", &Photon_convVertexZ, "Photon_convVertexZ/F");
		miniTree->Branch("Photon_etaWidth", &Photon_etaWidth, "Photon_etaWidth/F");
	miniTree->Branch("Photon_phiWidth", &Photon_phiWidth, "Photon_phiWidth/F");

	miniTree->Branch("Photon_dR03isoEcalRecHit", &Photon_dR03isoEcalRecHit, "Photon_dR03isoEcalRecHit/F");
	miniTree->Branch("Photon_dR03isoHcalRecHit", &Photon_dR03isoHcalRecHit, "Photon_dR03isoHcalRecHit/F");
	miniTree->Branch("Photon_dR03isoSolidTrkCone", &Photon_dR03isoSolidTrkCone, "Photon_dR03isoSolidTrkCone/F");
	miniTree->Branch("Photon_dR03isoHollowTrkCone", &Photon_dR03isoHollowTrkCone, "Photon_dR03isoHollowTrkCone/F");
	miniTree->Branch("Photon_dR03isoNTracksSolidCone", &Photon_dR03isoNTracksSolidCone, "Photon_dR03isoNTracksSolidCone/I");
	miniTree->Branch("Photon_dR03isoNTracksHollowCone", &Photon_dR03isoNTracksHollowCone, "Photon_dR03isoNTracksHollowCone/I");

	miniTree->Branch("Photon_dR04isoEcalRecHit", &Photon_dR04isoEcalRecHit, "Photon_dR04isoEcalRecHit/F");
	miniTree->Branch("Photon_dR04isoHcalRecHit", &Photon_dR04isoHcalRecHit, "Photon_dR04isoHcalRecHit/F");
	miniTree->Branch("Photon_dR04isoSolidTrkCone", &Photon_dR04isoSolidTrkCone, "Photon_dR04isoSolidTrkCone/F");
	miniTree->Branch("Photon_dR04isoHollowTrkCone", &Photon_dR04isoHollowTrkCone, "Photon_dR04isoHollowTrkCone/F");
	miniTree->Branch("Photon_dR04isoNTracksSolidCone", &Photon_dR04isoNTracksSolidCone, "Photon_dR04isoNTracksSolidCone/I");
	miniTree->Branch("Photon_dR04isoNTracksHollowCone", &Photon_dR04isoNTracksHollowCone, "Photon_dR04isoNTracksHollowCone/I");

	miniTree->Branch("Photon_seedTime", &Photon_seedTime, "Photon_seedTime/F");
	miniTree->Branch("Photon_seedFlag", &Photon_seedFlag, "Photon_seedFlag/F");
	miniTree->Branch("Photon_seedPosition1", &Photon_seedPosition1, "Photon_seedPosition1/I");
	miniTree->Branch("Photon_seedPosition2", &Photon_seedPosition2, "Photon_seedPosition2/I");
	miniTree->Branch("Photon_SC_Eta", &Photon_SC_Eta, "Photon_SC_Eta/F");
	miniTree->Branch("Photon_SC_Phi", &Photon_SC_Phi, "Photon_SC_Phi/F");
	miniTree->Branch("Photon_SC_brem", &Photon_SC_brem, "Photon_SC_brem/F");
	miniTree->Branch("Photon_SC_E", &Photon_SC_E, "Photon_SC_E/F");
	miniTree->Branch("Photon_SC_Et", &Photon_SC_Et, "Photon_SC_Et/F");
	miniTree->Branch("Photon_SC_rawEt", &Photon_SC_rawEt, "Photon_SC_rawEt/F");
	miniTree->Branch("Photon_SC_rawE", &Photon_SC_rawE, "Photon_SC_rawE/F");

	miniTree->Branch("Photon_lambdaRatio", &Photon_lambdaRatio, "Photon_lambdaRatio/F");
	miniTree->Branch("Photon_ratioSeed", &Photon_ratioSeed, "Photon_ratioSeed/F");
	miniTree->Branch("Photon_ratioS4", &Photon_ratioS4, "Photon_ratioS4/F");
	miniTree->Branch("Photon_lamdbaDivCov", &Photon_lamdbaDivCov, "Photon_lamdbaDivCov/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta", &Photon_SC_rawE_x_fEta, "Photon_SC_rawE_x_fEta/F");
	miniTree->Branch("Photon_ratioS4_corrected", &Photon_ratioS4_corrected, "Photon_ratioS4_corrected/F");

	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem", &Photon_SC_rawE_x_fEta_x_fBrem, "Photon_SC_rawE_x_fEta_x_fBrem/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_AF", &Photon_SC_rawE_x_fEta_x_fBrem_AF, "Photon_SC_rawE_x_fEta_x_fBrem_AF/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_L", &Photon_SC_rawE_x_fEta_x_fBrem_L, "Photon_SC_rawE_x_fEta_x_fBrem_L/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta", &Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta, "Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta", &Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta, "Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta", &Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta, "Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta/F");


	miniTree->Branch("Photon_secondMomentMaj", &Photon_secondMomentMaj, "Photon_secondMomentMaj/F");
	miniTree->Branch("Photon_secondMomentMin", &Photon_secondMomentMin, "Photon_secondMomentMin/F");
	miniTree->Branch("Photon_secondMomentAlpha", &Photon_secondMomentAlpha, "Photon_secondMomentAlpha/F");

	miniTree->Branch("Photon_etaLAT", &Photon_etaLAT, "Photon_etaLAT/F");
	miniTree->Branch("Photon_phiLAT", &Photon_phiLAT, "Photon_phiLAT/F");
	miniTree->Branch("Photon_LAT", &Photon_LAT, "Photon_LAT/F");
	miniTree->Branch("Photon_Zernike20", &Photon_Zernike20, "Photon_Zernike20/F");
	miniTree->Branch("Photon_Zernike42", &Photon_Zernike42, "Photon_Zernike42/F");
	miniTree->Branch("Photon_ESratio", &Photon_ESratio, "Photon_ESratio/F");
		
	// ____________________________________________
	// mugamma / mumu / mumugamma information
	// ____________________________________________

	miniTree->Branch("Mmumu", &Mmumu, "Mmumu/F");
	miniTree->Branch("Mmumugamma", &Mmumugamma, "Mmumugamma/F");
	miniTree->Branch("Mmumugamma_5x5", &Mmumugamma_5x5, "Mmumugamma_5x5/F");
	miniTree->Branch("Mmumugamma_SC", &Mmumugamma_SC, "Mmumugamma_SC/F");
	miniTree->Branch("Mmumugamma_SCraw", &Mmumugamma_SCraw, "Mmumugamma_SCraw/F");
	miniTree->Branch("Mmumugamma_SCraw_fEta", &Mmumugamma_SCraw_fEta, "Mmumugamma_SCraw_fEta/F");
	
	miniTree->Branch("Mmumugamma_SCraw_fEta_fBrem", &Mmumugamma_SCraw_fEta_fBrem, "Mmumugamma_SCraw_fEta_fBrem/F");
	miniTree->Branch("Mmumugamma_SCraw_fEta_fBrem_AF", &Mmumugamma_SCraw_fEta_fBrem_AF, "Mmumugamma_SCraw_fEta_fBrem_AF/F");	
	miniTree->Branch("Mmumugamma_SCraw_fEta_fBrem_L", &Mmumugamma_SCraw_fEta_fBrem_L, "Mmumugamma_SCraw_fEta_fBrem_L/F");
	miniTree->Branch("Mmumugamma_SCraw_fEta_fBrem_fEtEta", &Mmumugamma_SCraw_fEta_fBrem_fEtEta, "Mmumugamma_SCraw_fEta_fBrem_fEtEta/F");
	miniTree->Branch("Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta", &Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta, "Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta/F");
	miniTree->Branch("Mmumugamma_SCraw_fEta_fBrem_L_fEtEta", &Mmumugamma_SCraw_fEta_fBrem_L_fEtEta, "Mmumugamma_SCraw_fEta_fBrem_L_fEtEta/F");



	miniTree->Branch("Ptmumu", &Ptmumu, "Ptmumu/F");
	miniTree_allmuons->Branch("Ptmumu", &Ptmumu, "Ptmumu/F");

	miniTree->Branch("deltaRNear", &deltaRNear, "deltaRNear/F");
	miniTree->Branch("deltaRFar", &deltaRFar, "deltaRFar/F");
	miniTree->Branch("deltaRPlus", &deltaRPlus, "deltaRPlus/F");
	miniTree->Branch("deltaRMinus", &deltaRMinus, "deltaRMinus/F");
	miniTree->Branch("deltaRLeading", &deltaRLeading, "deltaRLeading/F");
	miniTree->Branch("deltaRSubleading", &deltaRSubleading, "deltaRSubleading/F");

	miniTree->Branch("mmg_k", &mmg_k, "mmg_k/F");
	miniTree->Branch("mmg_ik", &mmg_ik, "mmg_ik/F");
	miniTree->Branch("mmg_s", &mmg_s, "mmg_s/F");
	miniTree->Branch("mmg_logk", &mmg_logk, "mmg_logk/F");
	miniTree->Branch("mmg_logik", &mmg_logik, "mmg_logik/F");
	miniTree->Branch("mmg_logs", &mmg_logs, "mmg_logs/F");

	miniTree->Branch("mmg_k_5x5", &mmg_k_5x5, "mmg_k_5x5/F");
	miniTree->Branch("mmg_ik_5x5", &mmg_ik_5x5, "mmg_ik_5x5/F");
	miniTree->Branch("mmg_s_5x5", &mmg_s_5x5, "mmg_s_5x5/F");
	miniTree->Branch("mmg_logk_5x5", &mmg_logk_5x5, "mmg_logk_5x5/F");
	miniTree->Branch("mmg_logik_5x5", &mmg_logik_5x5, "mmg_logik_5x5/F");
	miniTree->Branch("mmg_logs_5x5", &mmg_logs_5x5, "mmg_logs_5x5/F");

	miniTree->Branch("mmg_k_SC", &mmg_k_SC, "mmg_k_SC/F");
	miniTree->Branch("mmg_ik_SC", &mmg_ik_SC, "mmg_ik_SC/F");
	miniTree->Branch("mmg_s_SC", &mmg_s_SC, "mmg_s_SC/F");
	miniTree->Branch("mmg_logk_SC", &mmg_logk_SC, "mmg_logk_SC/F");
	miniTree->Branch("mmg_logik_SC", &mmg_logik_SC, "mmg_logik_SC/F");
	miniTree->Branch("mmg_logs_SC", &mmg_logs_SC, "mmg_logs_SC/F");

	miniTree->Branch("mmg_k_SCraw", &mmg_k_SCraw, "mmg_k_SCraw/F");
	miniTree->Branch("mmg_ik_SCraw", &mmg_ik_SCraw, "mmg_ik_SCraw/F");
	miniTree->Branch("mmg_s_SCraw", &mmg_s_SCraw, "mmg_s_SCraw/F");
	miniTree->Branch("mmg_logk_SCraw", &mmg_logk_SCraw, "mmg_logk_SCraw/F");
	miniTree->Branch("mmg_logik_SCraw", &mmg_logik_SCraw, "mmg_logik_SCraw/F");
	miniTree->Branch("mmg_logs_SCraw", &mmg_logs_SCraw, "mmg_logs_SCraw/F");

	miniTree->Branch("mmg_k_SCraw_fEta", &mmg_k_SCraw_fEta, "mmg_k_SCraw_fEta/F");
	miniTree->Branch("mmg_ik_SCraw_fEta", &mmg_ik_SCraw_fEta, "mmg_ik_SCraw_fEta/F");
	miniTree->Branch("mmg_s_SCraw_fEta", &mmg_s_SCraw_fEta, "mmg_s_SCraw_fEta/F");
	miniTree->Branch("mmg_logk_SCraw_fEta", &mmg_logk_SCraw_fEta, "mmg_logk_SCraw_fEta/F");
	miniTree->Branch("mmg_logik_SCraw_fEta", &mmg_logik_SCraw_fEta, "mmg_logik_SCraw_fEta/F");
	miniTree->Branch("mmg_logs_SCraw_fEta", &mmg_logs_SCraw_fEta, "mmg_logs_SCraw_fEta/F");


	miniTree->Branch("mmg_ik_SCraw_fEta_fBrem", &mmg_ik_SCraw_fEta_fBrem, "mmg_ik_SCraw_fEta_fBrem/F");
	miniTree->Branch("mmg_ik_SCraw_fEta_fBrem_AF", &mmg_ik_SCraw_fEta_fBrem_AF, "mmg_ik_SCraw_fEta_fBrem_AF/F");
	miniTree->Branch("mmg_ik_SCraw_fEta_fBrem_L", &mmg_ik_SCraw_fEta_fBrem_L, "mmg_ik_SCraw_fEta_fBrem_L/F");
	miniTree->Branch("mmg_ik_SCraw_fEta_fBrem_fEtEta", &mmg_ik_SCraw_fEta_fBrem_fEtEta, "mmg_ik_SCraw_fEta_fBrem_fEtEta/F");
	miniTree->Branch("mmg_ik_SCraw_fEta_fBrem_AF_fEtEta", &mmg_ik_SCraw_fEta_fBrem_AF_fEtEta, "mmg_ik_SCraw_fEta_fBrem_AF_fEtEta/F");
	miniTree->Branch("mmg_ik_SCraw_fEta_fBrem_L_fEtEta", &mmg_ik_SCraw_fEta_fBrem_L_fEtEta, "mmg_ik_SCraw_fEta_fBrem_L_fEtEta/F");

	miniFriend->Branch("iCandidate", iCandidate, "iCandidate[8]/I");
	miniFriend->Branch("nCandidate", nCandidate, "nCandidate[8]/I");
	miniTree->Branch("MuonBeforeBremM_Pt", &MuonBeforeBremM_Pt, "MuonBeforeBremM_Pt/F");
	miniTree->Branch("MuonBeforeBremP_Pt", &MuonBeforeBremP_Pt, "MuonBeforeBremP_Pt/F");
	miniTree->Branch("MuonBeforeBremN_Pt", &MuonBeforeBremN_Pt, "MuonBeforeBremN_Pt/F");
	miniTree->Branch("MuonBeforeBremF_Pt", &MuonBeforeBremF_Pt, "MuonBeforeBremF_Pt/F");
	miniTree->Branch("MuonBeforeBremL_Pt", &MuonBeforeBremL_Pt, "MuonBeforeBremL_Pt/F");
	miniTree->Branch("MuonBeforeBremS_Pt", &MuonBeforeBremS_Pt, "MuonBeforeBremS_Pt/F");
	miniTree->Branch("MuonBeforeBremM_Eta", &MuonBeforeBremM_Eta, "MuonBeforeBremM_Eta/F");
	miniTree->Branch("MuonBeforeBremP_Eta", &MuonBeforeBremP_Eta, "MuonBeforeBremP_Eta/F");
	miniTree->Branch("MuonBeforeBremN_Eta", &MuonBeforeBremN_Eta, "MuonBeforeBremN_Eta/F");
	miniTree->Branch("MuonBeforeBremF_Eta", &MuonBeforeBremF_Eta, "MuonBeforeBremF_Eta/F");
	miniTree->Branch("MuonBeforeBremL_Eta", &MuonBeforeBremL_Eta, "MuonBeforeBremL_Eta/F");
	miniTree->Branch("MuonBeforeBremS_Eta", &MuonBeforeBremS_Eta, "MuonBeforeBremS_Eta/F");
	miniTree->Branch("MuonBeforeBremM_Phi", &MuonBeforeBremM_Phi, "MuonBeforeBremM_Phi/F");
	miniTree->Branch("MuonBeforeBremP_Phi", &MuonBeforeBremP_Phi, "MuonBeforeBremP_Phi/F");
	miniTree->Branch("MuonBeforeBremN_Phi", &MuonBeforeBremN_Phi, "MuonBeforeBremN_Phi/F");
	miniTree->Branch("MuonBeforeBremF_Phi", &MuonBeforeBremF_Phi, "MuonBeforeBremF_Phi/F");
	miniTree->Branch("MuonBeforeBremL_Phi", &MuonBeforeBremL_Phi, "MuonBeforeBremL_Phi/F");
	miniTree->Branch("MuonBeforeBremS_Phi", &MuonBeforeBremS_Phi, "MuonBeforeBremS_Phi/F");
	miniTree->Branch("MuonBeforeBremM_E", &MuonBeforeBremM_E, "MuonBeforeBremM_E/F");
	miniTree->Branch("MuonBeforeBremP_E", &MuonBeforeBremP_E, "MuonBeforeBremP_E/F");
	miniTree->Branch("MuonBeforeBremN_E", &MuonBeforeBremN_E, "MuonBeforeBremN_E/F");
	miniTree->Branch("MuonBeforeBremF_E", &MuonBeforeBremF_E, "MuonBeforeBremF_E/F");
	miniTree->Branch("MuonBeforeBremL_E", &MuonBeforeBremL_E, "MuonBeforeBremL_E/F");
	miniTree->Branch("MuonBeforeBremS_E", &MuonBeforeBremS_E, "MuonBeforeBremS_E/F");
	miniTree->Branch("MuonBeforeBremM_Px", &MuonBeforeBremM_Px, "MuonBeforeBremM_Px/F");
	miniTree->Branch("MuonBeforeBremP_Px", &MuonBeforeBremP_Px, "MuonBeforeBremP_Px/F");
	miniTree->Branch("MuonBeforeBremN_Px", &MuonBeforeBremN_Px, "MuonBeforeBremN_Px/F");
	miniTree->Branch("MuonBeforeBremF_Px", &MuonBeforeBremF_Px, "MuonBeforeBremF_Px/F");
	miniTree->Branch("MuonBeforeBremL_Px", &MuonBeforeBremL_Px, "MuonBeforeBremL_Px/F");
	miniTree->Branch("MuonBeforeBremS_Px", &MuonBeforeBremS_Px, "MuonBeforeBremS_Px/F");
	miniTree->Branch("MuonBeforeBremM_Py", &MuonBeforeBremM_Py, "MuonBeforeBremM_Py/F");
	miniTree->Branch("MuonBeforeBremP_Py", &MuonBeforeBremP_Py, "MuonBeforeBremP_Py/F");
	miniTree->Branch("MuonBeforeBremN_Py", &MuonBeforeBremN_Py, "MuonBeforeBremN_Py/F");
	miniTree->Branch("MuonBeforeBremF_Py", &MuonBeforeBremF_Py, "MuonBeforeBremF_Py/F");
	miniTree->Branch("MuonBeforeBremL_Py", &MuonBeforeBremL_Py, "MuonBeforeBremL_Py/F");
	miniTree->Branch("MuonBeforeBremS_Py", &MuonBeforeBremS_Py, "MuonBeforeBremS_Py/F");
	miniTree->Branch("MuonBeforeBremM_Pz", &MuonBeforeBremM_Pz, "MuonBeforeBremM_Pz/F");
	miniTree->Branch("MuonBeforeBremP_Pz", &MuonBeforeBremP_Pz, "MuonBeforeBremP_Pz/F");
	miniTree->Branch("MuonBeforeBremN_Pz", &MuonBeforeBremN_Pz, "MuonBeforeBremN_Pz/F");
	miniTree->Branch("MuonBeforeBremF_Pz", &MuonBeforeBremF_Pz, "MuonBeforeBremF_Pz/F");
	miniTree->Branch("MuonBeforeBremL_Pz", &MuonBeforeBremL_Pz, "MuonBeforeBremL_Pz/F");
	miniTree->Branch("MuonBeforeBremS_Pz", &MuonBeforeBremS_Pz, "MuonBeforeBremS_Pz/F");
	miniTree->Branch("MuonBeforeBremF_Charge", &MuonBeforeBremF_Charge, "MuonBeforeBremF_Charge/I");
	miniTree->Branch("MuonBeforeBremN_Charge", &MuonBeforeBremN_Charge, "MuonBeforeBremN_Charge/I");
	miniTree->Branch("MuonBeforeBremL_Charge", &MuonBeforeBremL_Charge, "MuonBeforeBremL_Charge/I");
	miniTree->Branch("MuonBeforeBremS_Charge", &MuonBeforeBremS_Charge, "MuonBeforeBremS_Charge/I");


 
	// ____________________________________________
	// Neural Network variables
	// ____________________________________________
	miniTree->Branch("Photon_NNshapeOutput", &Photon_NNshapeOutput, "Photon_NNshapeOutput/F");
 

	// ____________________________________________
	// Surface variables
	// ____________________________________________
				miniTree->Branch("MZ_Surface", &MZ_Surface, "MZ_Surface/F");
	miniTree->Branch("mmg_k_MZ_Surface", &mmg_k_MZ_Surface, "mmg_k_MZ_Surface/F");
				miniTree->Branch("mmg_ik_MZ_Surface", &mmg_ik_MZ_Surface, "mmg_ik_MZ_Surface/F");
				miniTree->Branch("mmg_s_MZ_Surface", &mmg_s_MZ_Surface, "mmg_s_MZ_Surface/F");
				miniTree->Branch("mmg_logk_MZ_Surface", &mmg_logk_MZ_Surface, "mmg_logk_MZ_Surface/F");
				miniTree->Branch("mmg_logik_MZ_Surface", &mmg_logik_MZ_Surface, "mmg_logik_MZ_Surface/F");
				miniTree->Branch("mmg_logs_MZ_Surface", &mmg_logs_MZ_Surface, "mmg_logs_MZ_Surface/F");


	// ____________________________________________
	// MC Truth
	// ___________________________________________

	miniTree->Branch("Photon_MC_E", &Photon_MC_E, "Photon_MC_E/F");
	miniTree->Branch("Photon_MC_Px", &Photon_MC_Px, "Photon_MC_Px/F");
	miniTree->Branch("Photon_MC_Py", &Photon_MC_Py, "Photon_MC_Py/F");
	miniTree->Branch("Photon_MC_Pz", &Photon_MC_Pz, "Photon_MC_Pz/F");
	miniTree->Branch("Photon_MC_Phi", &Photon_MC_Phi, "Photon_MC_Phi/F");
	miniTree->Branch("Photon_MC_Eta", &Photon_MC_Eta, "Photon_MC_Eta/F");
	miniTree->Branch("Photon_MC_Pt", &Photon_MC_Pt, "Photon_MC_Pt/F");
	miniTree->Branch("Photon_MCisConverted", &Photon_MCisConverted, "Photon_MCisConverted/F");
	miniTree->Branch("Photon_MCconvEoverP", &Photon_MCconvEoverP, "Photon_MCconvEoverP/F");
	miniTree->Branch("Photon_MCconvMass", &Photon_MCconvMass, "Photon_MCconvMass/F");
	miniTree->Branch("Photon_MCconvCotanTheta", &Photon_MCconvCotanTheta, "Photon_MCconvCotanTheta/F");
	miniTree->Branch("Photon_MCconvVertexX", &Photon_MCconvVertexX, "Photon_MCconvVertexX/F");
	miniTree->Branch("Photon_MCconvVertexY", &Photon_MCconvVertexY, "Photon_MCconvVertexY/F");
	miniTree->Branch("Photon_MCconvVertexZ", &Photon_MCconvVertexZ, "Photon_MCconvVertexZ/F");
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
	miniTree->Branch("Photon_SC_rawE_x_fEta_o_MC_E", &Photon_SC_rawE_x_fEta_o_MC_E, "Photon_SC_rawE_x_fEta_o_MC_E/F");
	miniTree->Branch("Photon_E_o_MC_E", &Photon_E_o_MC_E, "Photon_E_o_MC_E/F");

	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_o_MC_E", &Photon_SC_rawE_x_fEta_x_fBrem_o_MC_E, "Photon_SC_rawE_x_fEta_x_fBrem_o_MC_E/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_AF_o_MC_E", &Photon_SC_rawE_x_fEta_x_fBrem_AF_o_MC_E, "Photon_SC_rawE_x_fEta_x_fBrem_AF_o_MC_E/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_L_o_MC_E", &Photon_SC_rawE_x_fEta_x_fBrem_L_o_MC_E, "Photon_SC_rawE_x_fEta_x_fBrem_L_o_MC_E/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta_o_MC_E", &Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta_o_MC_E, "Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta_o_MC_E/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta_o_MC_E", &Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta_o_MC_E, "Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta_o_MC_E/F");
	miniTree->Branch("Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta_o_MC_E", &Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta_o_MC_E, "Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta_o_MC_E/F");


	miniTree->Branch("Mmumu_Photon_MC", &Mmumu_Photon_MC, "Mmumu_Photon_MC/F");
	miniTree->Branch("Mmumugamma_Photon_MC", &Mmumugamma_Photon_MC, "Mmumugamma_Photon_MC/F");
	miniTree->Branch("mmg_k_Photon_MC", &mmg_k_Photon_MC, "mmg_k_Photon_MC/F");
	miniTree->Branch("mmg_ik_Photon_MC", &mmg_ik_Photon_MC, "mmg_ik_Photon_MC/F");
	miniTree->Branch("mmg_s_Photon_MC", &mmg_s_Photon_MC, "mmg_s_Photon_MC/F");
	miniTree->Branch("mmg_logk_Photon_MC", &mmg_logk_Photon_MC, "mmg_logk_Photon_MC/F");
	miniTree->Branch("mmg_logik_Photon_MC", &mmg_logik_Photon_MC, "mmg_logik_Photon_MC/F");
	miniTree->Branch("mmg_logs_Photon_MC", &mmg_logs_Photon_MC, "mmg_logs_Photon_MC/F");
	miniTree->Branch("Mmumu_Muons_MC", &Mmumu_Muons_MC, "Mmumu_Muons_MC/F");
	miniTree->Branch("Mmumugamma_Muons_MC", &Mmumugamma_Muons_MC, "Mmumugamma_Muons_MC/F");
	miniTree->Branch("mmg_k_Muons_MC", &mmg_k_Muons_MC, "mmg_k_Muons_MC/F");
	miniTree->Branch("mmg_ik_Muons_MC", &mmg_ik_Muons_MC, "mmg_ik_Muons_MC/F");
	miniTree->Branch("mmg_s_Muons_MC", &mmg_s_Muons_MC, "mmg_s_Muons_MC/F");
	miniTree->Branch("mmg_logk_Muons_MC", &mmg_logk_Muons_MC, "mmg_logk_Muons_MC/F");
	miniTree->Branch("mmg_logik_Muons_MC", &mmg_logik_Muons_MC, "mmg_logik_Muons_MC/F");
	miniTree->Branch("mmg_logs_Muons_MC", &mmg_logs_Muons_MC, "mmg_logs_Muons_MC/F");
	miniTree->Branch("Mmumu_MMG_MC", &Mmumu_MMG_MC, "Mmumu_MMG_MC/F");
	miniTree->Branch("Mmumugamma_MMG_MC", &Mmumugamma_MMG_MC, "Mmumugamma_MMG_MC/F");
	miniTree->Branch("mmg_k_MMG_MC", &mmg_k_MMG_MC, "mmg_k_MMG_MC/F");
	miniTree->Branch("mmg_ik_MMG_MC", &mmg_ik_MMG_MC, "mmg_ik_MMG_MC/F");
	miniTree->Branch("mmg_s_MMG_MC", &mmg_s_MMG_MC, "mmg_s_MMG_MC/F");
	miniTree->Branch("mmg_logk_MMG_MC", &mmg_logk_MMG_MC, "mmg_logk_MMG_MC/F");
	miniTree->Branch("mmg_logik_MMG_MC", &mmg_logik_MMG_MC, "mmg_logik_MMG_MC/F");
	miniTree->Branch("mmg_logs_MMG_MC", &mmg_logs_MMG_MC, "mmg_logs_MMG_MC/F");

	miniTree->Branch("mmg_k_MZ",&mmg_k_MZ,"mmg_k_MZ/F");
	miniTree->Branch("mmg_ik_MZ",&mmg_ik_MZ,"mmg_ik_MZ/F");
	miniTree->Branch("mmg_s_MZ",&mmg_s_MZ,"mmg_s_MZ/F");
	miniTree->Branch("mmg_logk_MZ",&mmg_logk_MZ,"mmg_logk_MZ/F");
	miniTree->Branch("mmg_logik_MZ",&mmg_logik_MZ,"mmg_logik_MZ/F");
	miniTree->Branch("mmg_logs_MZ",&mmg_logs_MZ,"mmg_logs_MZ/F");
	miniTree->Branch("mmg_k_MZ_Photon_MC",&mmg_k_MZ_Photon_MC,"mmg_k_MZ_Photon_MC/F");
	miniTree->Branch("mmg_ik_MZ_Photon_MC",&mmg_ik_MZ_Photon_MC,"mmg_ik_MZ_Photon_MC/F");
	miniTree->Branch("mmg_s_MZ_Photon_MC",&mmg_s_MZ_Photon_MC,"mmg_s_MZ_Photon_MC/F");
	miniTree->Branch("mmg_logk_MZ_Photon_MC",&mmg_logk_MZ_Photon_MC,"mmg_logk_MZ_Photon_MC/F");
	miniTree->Branch("mmg_logik_MZ_Photon_MC",&mmg_logik_MZ_Photon_MC,"mmg_logik_MZ_Photon_MC/F");
	miniTree->Branch("mmg_logs_MZ_Photon_MC",&mmg_logs_MZ_Photon_MC,"mmg_logs_MZ_Photon_MC/F");
	miniTree->Branch("mmg_k_MZ_Muons_MC",&mmg_k_MZ_Muons_MC,"mmg_k_MZ_Muons_MC/F");
	miniTree->Branch("mmg_ik_MZ_Muons_MC",&mmg_ik_MZ_Muons_MC,"mmg_ik_MZ_Muons_MC/F");
	miniTree->Branch("mmg_s_MZ_Muons_MC",&mmg_s_MZ_Muons_MC,"mmg_s_MZ_Muons_MC/F");
	miniTree->Branch("mmg_logk_MZ_Muons_MC",&mmg_logk_MZ_Muons_MC,"mmg_logk_MZ_Muons_MC/F");
	miniTree->Branch("mmg_logik_MZ_Muons_MC",&mmg_logik_MZ_Muons_MC,"mmg_logik_MZ_Muons_MC/F");
	miniTree->Branch("mmg_logs_MZ_Muons_MC",&mmg_logs_MZ_Muons_MC,"mmg_logs_MZ_Muons_MC/F");
	miniTree->Branch("mmg_k_MZ_Muons_RECO_MC",&mmg_k_MZ_Muons_RECO_MC,"mmg_k_MZ_Muons_RECO_MC/F");
	miniTree->Branch("mmg_ik_MZ_Muons_RECO_MC",&mmg_ik_MZ_Muons_RECO_MC,"mmg_ik_MZ_Muons_RECO_MC/F");
	miniTree->Branch("mmg_s_MZ_Muons_RECO_MC",&mmg_s_MZ_Muons_RECO_MC,"mmg_s_MZ_Muons_RECO_MC/F");
	miniTree->Branch("mmg_logk_MZ_Muons_RECO_MC",&mmg_logk_MZ_Muons_RECO_MC,"mmg_logk_MZ_Muons_RECO_MC/F");
	miniTree->Branch("mmg_logik_MZ_Muons_RECO_MC",&mmg_logik_MZ_Muons_RECO_MC,"mmg_logik_MZ_Muons_RECO_MC/F");
	miniTree->Branch("mmg_logs_MZ_Muons_RECO_MC",&mmg_logs_MZ_Muons_RECO_MC,"mmg_logs_MZ_Muons_RECO_MC/F");

	TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );
	reader->AddVariable("pho_cEP",&Photon_covEtaPhi);
	reader->AddVariable("pho_SCbr",&Photon_SC_brem);
	reader->AddVariable("(pho_cEE+pho_cPP-sqrt((pho_cEE-pho_cPP)**2+4*pho_cEP**2))/(pho_cEE+pho_cPP+sqrt((pho_cEE-pho_cPP)**2+4*pho_cEP**2))",&Photon_lambdaRatio);
	reader->AddVariable("pho_eMax/pho_SCEraw",&Photon_ratioSeed);
	reader->AddVariable("pho_etawidth",&Photon_etaWidth);
	reader->AddVariable("pho_e2x2/pho_e5x5",&Photon_ratioS4_corrected);
	reader->AddVariable("(pho_cEE+pho_cPP-sqrt((pho_cEE-pho_cPP)**2+4*pho_cEP**2))/pho_cEE",&Photon_lamdbaDivCov);
	reader->AddVariable("pho_r19",&Photon_r19);
//	reader->BookMVA("MLP method","/sps/cms/hbrun/DiPhotons41X/diPhotonMC/weights/TMVAClassification_MLP.weights.xml");
	reader->BookMVA("MLP method","/sps/cms/hbrun/DiPhotons42X/diPhotonMC/weights/TMVAClassification_MLP.weights.xml");
//	reader->BookMVA("MLP method","/sps/cms/hbrun/DiPhotons42X/diPhotonMC/weights/TMVAClassification_MLP.weights.xml");


	
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
//	double XSectionWJetsToLNu = 7899.0;
	double XSectionWJetsToLNu = 27770.0;
	double XSectionQCDMu = 349988.0;

//	double InitialNumberDYToMuMu = 2148325.0;
	double InitialNumberDYToMuMu = 29743564.0;
//	double InitialNumberTTJets = 1089625.0;
	double InitialNumberTTJets = 3701947.0;
	double InitialNumberWJetsToLNu = 81345381.0;
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
//		NbEventsPerJob = 100;
		NbEventsBegin = ijob * NbEventsPerJob;
		NbEventsEnd = min( (ijob + 1)* NbEventsPerJob , (int)NbEvents);
		NbEvents = NbEventsEnd - NbEventsBegin ;
	}
	cout << "NbEventsBegin= " << NbEventsBegin << "\tNbEventsEnd= " << NbEventsEnd << "\tNbEventsPerJob= " << NbEventsPerJob << endl;



//return 1;


// ***************************************************************************************************
// ***************************************************************************************************
// ********************************** LOOP over events ***********************************************
// ***************************************************************************************************
// ***************************************************************************************************
	for(unsigned int ievt=NbEventsBegin; ievt<NbEventsEnd; ievt++)
//	for(unsigned int ievt=0; ievt<NbEvents; ievt++)
//	for(unsigned int ievt=0; ievt<100; ievt++)
	{
		if(verbosity>4) cout << "analysing event ievt= " << ievt << endl;
		nBeforeAllCuts++;
		isBeforeAllCuts = 1;
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
			//nGenVertices = event->nInTimePUVertices();
			//nGenVertices = event->pu_NumInteractions();
			//nGenVertices = event->inTimePU_NumInteractions();
			nGenVertices = event->pu_TrueNumInteractions();
		}
		isSignalApplied = signal;
		isStewApplied = stew;
		isZJetsApplied = zjet_veto;
		isAfterCutPthatFilter = isAfterCutZJETVETO = 0;
		isVeryLooseMMG = isLooseMMG = isMM = isTightMMG = isMMGCandidate = 0;
		isAfterFSRCut1 = isAfterFSRCut2 = isAfterFSRCut3 = 0;
		isAfterFSRCut4 = isMultipleCandidate = isAfterCut5 = isAfterCut6 = isAfterCut7 = isAfterCut8 = isAfterCut9 = isAfterCut10 = 0;
		isSelected = 0;

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
		} else if( sample_in.find("QCD") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionQCDMu) / (double)(InitialNumberQCDMu)) * (double)integratedLuminosity);
			weight_pileUp = weight_QCDMu(nGenVertices+1);
		} else if( sample_in.find("TTJets") != string::npos )
		{
			weight_Xsection = (double)(	(double)((double)(XSectionTTJets) / (double)(InitialNumberTTJets)) * (double)integratedLuminosity);
			weight_pileUp = weight_TTJets(nGenVertices+1);
		} else if( sample_in.find("WToMuNu") != string::npos	|| sample_in.find("WJetsToLNu") != string::npos)
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
		// Photon variables
		// ___________________________________________
		NbPhotons = photons->GetEntries();
		vector<double> Photon_scale;
		if( (argc > 5) && (EResolution != 0) )
		{ // If there is an extra resolution to smear the photon energy
			for(int iphoton = 0; iphoton < NbPhotons ; iphoton++)
			{
				Photon_scale.push_back(generator->Gaus(EScale_true_injected * EScale_inj, EResolution));
			}
//			cout << "Photon_scale = " << Photon_scale[Photon_scale.size() -1] << endl;
		} else {
			if( (argc > 4) && isManualCorrectionsApplied)
			{
				for(int iphoton = 0; iphoton < NbPhotons ; iphoton++)
				{
					TRootPhoton *myphotontocorrect;
					myphotontocorrect = (TRootPhoton*) photons->At(iphoton);
					Photon_scale.push_back(photonManualCorrectionFactor(myphotontocorrect, correction, clusters, superClusters, photons) * EScale_true_injected);
/*
					cout << "myphotontocorrect->isEB()= " << myphotontocorrect->isEB() << endl;
					cout << "myphotontocorrect->r9()= " << myphotontocorrect->r9() << endl;
					cout << "brem= " << (double)(myphotontocorrect->superCluster()->phiWidth()) / (double)(myphotontocorrect->superCluster()->etaWidth()) << endl;
					cout << "myphotontocorrect->superCluster()->rawEnergy()= " << myphotontocorrect->superCluster()->rawEnergy() << endl;
					cout << "myphotontocorrect->preshowerEnergy()= " << myphotontocorrect->preshowerEnergy() << endl;
*/
					if(verbosity>1) cout << "EScale_true_injected= " << EScale_true_injected	<< endl;
					if(verbosity>1) cout << "myphotontocorrect->Energy()= " << myphotontocorrect->Energy() << endl;
					if(verbosity>1) cout << "photonManualCorrectionFactor= " << Photon_scale[Photon_scale.size() - 1] << endl;
					if(verbosity>1) cout << "photonManualCorrectionFactor * myphotontocorrect->Energy()= " << Photon_scale[Photon_scale.size() - 1] * myphotontocorrect->Energy() << endl;
					if(verbosity>1) cout << endl << endl;

//					return 2;
				}
			} else {
				for(int iphoton = 0; iphoton < NbPhotons ; iphoton++)
				{
					Photon_scale.push_back(EScale_true_injected * EScale_inj);
				}
			}
		}
		Pt_allPhotons = Eta_allPhotons = Phi_allPhotons = Cross_allPhotons = -99;
		isNotCommissionned = -99;
		isEBorEE_allPhotons = 1;
		isEB_allPhotons, isEE_allPhotons, isEEM_allPhotons, isEEP_allPhotons = -99;
		Photon_Eta = Photon_Phi = -99;
		Photon_Px = Photon_Py = Photon_Pz = -99;
		Photon_isEBorEE = 1;
		Photon_isEB = Photon_isEE = Photon_isEEP = Photon_isEEM = -99;
		Photon_hasPixelSeed = Photon_isAlsoElectron = Photon_Nclusters = Photon_nBasicClusters = Photon_nXtals = -99;
		//Photon_isTightPhoton = Photon_isLoosePhoton = -99;
		Photon_convNTracks = Photon_isConverted = -99;
		Photon_convEoverP = Photon_convMass = Photon_convCotanTheta = Photon_convLikely = Photon_convVertexX = Photon_convVertexY = Photon_convVertexZ = -99.0;	
		Photon_E = Photon_Et = Photon_E2x2 = Photon_E3x3 = Photon_E5x5 = Photon_Emax = Photon_E2nd = -99.0;
		Photon_Ecorr_o_Ereco = -99.0;
		Photon_r19 = Photon_r9 = Photon_cross = -99.0;
		//Photon_caloConeSize = -99.0;
		Photon_PreshEnergy = Photon_HoE = Photon_sigmaEtaEta = Photon_sigmaIetaIeta = Photon_covEtaEta = Photon_covPhiPhi = Photon_covEtaPhi = Photon_etaWidth = Photon_phiWidth = -99.0;
		Photon_dR03isoEcalRecHit = Photon_dR03isoHcalRecHit = Photon_dR03isoSolidTrkCone = Photon_dR03isoHollowTrkCone = Photon_dR03isoNTracksSolidCone = Photon_dR03isoNTracksHollowCone = -99.0;
		Photon_dR04isoEcalRecHit = Photon_dR04isoHcalRecHit = Photon_dR04isoSolidTrkCone = Photon_dR04isoHollowTrkCone = Photon_dR04isoNTracksSolidCone = Photon_dR04isoNTracksHollowCone = -99.0;
		Photon_seedTime = Photon_seedFlag = -99.0;
		Photon_seedPosition1 = Photon_seedPosition2 = -99;
		Photon_SC_Eta = Photon_SC_Phi = Photon_SC_brem = -99.0;
		Photon_SC_E = Photon_SC_Et = Photon_SC_rawE = Photon_SC_rawEt = -99.0;
		Photon_lambdaRatio = Photon_ratioSeed = Photon_ratioS4 = Photon_lamdbaDivCov = -99.0;
		Photon_ratioS4_corrected = -99.0;
		Photon_SC_rawE_x_fEta = -99.0;
		Photon_SC_rawE_x_fEta_x_fBrem = Photon_SC_rawE_x_fEta_x_fBrem_AF = Photon_SC_rawE_x_fEta_x_fBrem_L = Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta = Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta = Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta = -99.0;
		Photon_secondMomentMaj = Photon_secondMomentMin = Photon_secondMomentAlpha = -99.0;
		Photon_etaLAT = Photon_phiLAT = Photon_LAT = Photon_Zernike20 = Photon_Zernike42 = Photon_ESratio = -99.0;
		Photon_E_regression = Photon_E_regressionError = Photon_Et_regression = -99.0;

		// ____________________________________________
		// mugamma / mumu / mumugamma information
		// ____________________________________________
		Mmumu = Mmumugamma = Mmumugamma_5x5 = Mmumugamma_SC = Mmumugamma_SCraw = Mmumugamma_SCraw_fEta = -99.0;
		Mmumugamma_SCraw_fEta_fBrem = Mmumugamma_SCraw_fEta_fBrem_AF = Mmumugamma_SCraw_fEta_fBrem_L = Mmumugamma_SCraw_fEta_fBrem_fEtEta = Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta = Mmumugamma_SCraw_fEta_fBrem_L_fEtEta = -99.0;
		Ptmumu = -99.0;
		deltaRNear = deltaRFar = deltaRPlus = deltaRMinus = deltaRLeading = deltaRSubleading = -99.0;
		mmg_k = mmg_ik = mmg_s = mmg_logk = mmg_logik = mmg_logs = -99.0;
		mmg_k_5x5 = mmg_ik_5x5 = mmg_s_5x5 = mmg_logk_5x5 = mmg_logik_5x5 = mmg_logs_5x5 = -99.0;
		mmg_k_SC = mmg_ik_SC = mmg_s_SC = mmg_logk_SC = mmg_logik_SC = mmg_logs_SC = -99.0;
		mmg_k_SCraw = mmg_ik_SCraw = mmg_s_SCraw = mmg_logk_SCraw = mmg_logik_SCraw = mmg_logs_SCraw = -99.0;
		mmg_k_SCraw_fEta = mmg_ik_SCraw_fEta = mmg_s_SCraw_fEta = mmg_logk_SCraw_fEta = mmg_logik_SCraw_fEta = mmg_logs_SCraw_fEta = -99.0;
		mmg_ik_SCraw_fEta_fBrem = mmg_ik_SCraw_fEta_fBrem_AF = mmg_ik_SCraw_fEta_fBrem_L = mmg_ik_SCraw_fEta_fBrem_fEtEta = mmg_ik_SCraw_fEta_fBrem_AF_fEtEta = mmg_ik_SCraw_fEta_fBrem_L_fEtEta = -99.0;

		MuonBeforeBremM_Pt = MuonBeforeBremP_Pt = MuonBeforeBremN_Pt = MuonBeforeBremF_Pt = MuonBeforeBremL_Pt = MuonBeforeBremS_Pt = -99.0;
		MuonBeforeBremM_Eta = MuonBeforeBremP_Eta = MuonBeforeBremN_Eta = MuonBeforeBremF_Eta = MuonBeforeBremL_Eta = MuonBeforeBremS_Eta = -99.0;
		MuonBeforeBremM_Phi = MuonBeforeBremP_Phi = MuonBeforeBremN_Phi = MuonBeforeBremF_Phi = MuonBeforeBremL_Phi = MuonBeforeBremS_Phi = -99.0;
		MuonBeforeBremM_E = MuonBeforeBremP_E = MuonBeforeBremN_E = MuonBeforeBremF_E = MuonBeforeBremL_E = MuonBeforeBremS_E = -99.0;
		MuonBeforeBremM_Px = MuonBeforeBremP_Px = MuonBeforeBremN_Px = MuonBeforeBremF_Px = MuonBeforeBremL_Px = MuonBeforeBremS_Px = -99.0;
		MuonBeforeBremM_Py = MuonBeforeBremP_Py = MuonBeforeBremN_Py = MuonBeforeBremF_Py = MuonBeforeBremL_Py = MuonBeforeBremS_Py = -99.0;
		MuonBeforeBremM_Pz = MuonBeforeBremP_Pz = MuonBeforeBremN_Pz = MuonBeforeBremF_Pz = MuonBeforeBremL_Pz = MuonBeforeBremS_Pz = -99.0;
		for(int ia= 0; ia < 8; ia ++) iCandidate[ia] = -99;
		for(int ia= 0; ia < 8; ia ++) nCandidate[ia] = -99;
		for(int ia= 0; ia < 8; ia ++) { for(int ib=0; ib < 50; ib++) iCandidate_temp[ia][ib] = -99; }
		MuonBeforeBremF_Charge = MuonBeforeBremN_Charge = MuonBeforeBremL_Charge = MuonBeforeBremS_Charge = -99;

	// ____________________________________________
	// Neural Network variables
		// ____________________________________________
		Photon_NNshapeOutput = -99.0;
	
	// ____________________________________________
				// Surface variables
				// ____________________________________________
								MZ_Surface = -99.0;
		mmg_k_MZ_Surface = mmg_ik_MZ_Surface = mmg_s_MZ_Surface = mmg_logk_MZ_Surface = mmg_logik_MZ_Surface = mmg_logs_MZ_Surface = -99.0;

		// ____________________________________________
		// MC Truth
		// ___________________________________________
		Photon_MC_E = Photon_MC_Px = Photon_MC_Py = Photon_MC_Pz = Photon_MC_Phi = Photon_MC_Eta = Photon_MC_Pt = -99.0;
		Photon_MCisConverted = -99;
		Photon_MCconvEoverP = Photon_MCconvMass = Photon_MCconvCotanTheta = Photon_MCconvVertexX = Photon_MCconvVertexY = Photon_MCconvVertexZ = -99.0;
		MuonM_MC_E = MuonM_MC_Px = MuonM_MC_Py = MuonM_MC_Pz = MuonM_MC_Phi = MuonM_MC_Eta = MuonM_MC_Pt = -99.0;
		MuonP_MC_E = MuonP_MC_Px = MuonP_MC_Py = MuonP_MC_Pz = MuonP_MC_Phi = MuonP_MC_Eta = MuonP_MC_Pt = -99.0;
		MuonN_MC_E = MuonN_MC_Px = MuonN_MC_Py = MuonN_MC_Pz = MuonN_MC_Phi = MuonN_MC_Eta = MuonN_MC_Pt = -99.0;
		MuonF_MC_E = MuonF_MC_Px = MuonF_MC_Py = MuonF_MC_Pz = MuonF_MC_Phi = MuonF_MC_Eta = MuonF_MC_Pt = -99.0;
		MuonL_MC_E = MuonL_MC_Px = MuonL_MC_Py = MuonL_MC_Pz = MuonL_MC_Phi = MuonL_MC_Eta = MuonL_MC_Pt = -99.0;
		MuonS_MC_E = MuonS_MC_Px = MuonS_MC_Py = MuonS_MC_Pz = MuonS_MC_Phi = MuonS_MC_Eta = MuonS_MC_Pt = -99.0;
		Photon_SC_rawE_x_fEta_o_MC_E = Photon_E_o_MC_E = -99.0;
		Photon_SC_rawE_x_fEta_x_fBrem_o_MC_E = Photon_SC_rawE_x_fEta_x_fBrem_AF_o_MC_E	= Photon_SC_rawE_x_fEta_x_fBrem_L_o_MC_E	= Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta_o_MC_E	= Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta_o_MC_E	= Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta_o_MC_E	= -99.0;

		Mmumu_Photon_MC = Mmumugamma_Photon_MC = mmg_k_Photon_MC = mmg_ik_Photon_MC = mmg_s_Photon_MC = mmg_logk_Photon_MC = mmg_logik_Photon_MC = mmg_logs_Photon_MC = -99.0;
		Mmumu_Muons_MC = Mmumugamma_Muons_MC = mmg_k_Muons_MC = mmg_ik_Muons_MC = mmg_s_Muons_MC = mmg_logk_Muons_MC = mmg_logik_Muons_MC = mmg_logs_Muons_MC = -99.0;
		Mmumu_MMG_MC = Mmumugamma_MMG_MC = mmg_k_MMG_MC = mmg_ik_MMG_MC = mmg_s_MMG_MC = mmg_logk_MMG_MC = mmg_logik_MMG_MC = mmg_logs_MMG_MC = -99.0;
	
		mmg_k_MZ = mmg_ik_MZ = mmg_s_MZ = mmg_logk_MZ = mmg_logik_MZ = mmg_logs_MZ = -99.0;					 
		mmg_k_MZ_Photon_MC = mmg_ik_MZ_Photon_MC = mmg_s_MZ_Photon_MC = mmg_logk_MZ_Photon_MC = mmg_logik_MZ_Photon_MC = mmg_logs_MZ_Photon_MC = -99.0;		
		mmg_k_MZ_Muons_MC = mmg_ik_MZ_Muons_MC = mmg_s_MZ_Muons_MC = mmg_logk_MZ_Muons_MC = mmg_logik_MZ_Muons_MC = mmg_logs_MZ_Muons_MC = -99.0;
		mmg_k_MZ_Muons_RECO_MC = mmg_ik_MZ_Muons_RECO_MC = mmg_s_MZ_Muons_RECO_MC = mmg_logk_MZ_Muons_RECO_MC = mmg_logik_MZ_Muons_RECO_MC = mmg_logs_MZ_Muons_RECO_MC = -99.0;	

		// ____________________________________________
		// END OF INITIALIZATION
		// ____________________________________________


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
//						if( abs(mcParticleCandidate->grannyType()) == 23 ){// photon coming from a muon coming from a Z
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
//							if( abs(mcParticleCandidate->grannyType()) == 23 ){// muon is coming from a Z
								if( (mcParticleCandidate->Pt()>8.0) && (abs(mcParticleCandidate->Eta())<3.0) ){
									MC_first_muon_in_phase_space = true;
								}
							}
						}
					} // end of selecting first muon
					if( (MC_first_muon_in_phase_space == true) && (mcParticleCandidate->status()==1) && (abs(mcParticleCandidate->type()) == 13) ){// if the particle is a final state muon
						if( abs(mcParticleCandidate->motherType()) == 13 ){
							if( abs(mcParticleCandidate->oldgrannyType()) == 23 ){// muon is coming from a Z
//							if( abs(mcParticleCandidate->grannyType()) == 23 ){// muon is coming from a Z
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
//				cerr<<"SAFE: photon(s) coming from muon, aborting event " << ievt << endl;
				continue;
			}
			isAfterCutZJETVETO = 1;
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
//				cerr << "There is " << mcMuMuGammaEvent->nFSR() << " fsr photons in the event" << endl;
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
//					cerr<<"SAFE: photon(s) coming from muon, aborting event " << ievt << endl;

					continue;
				}
				isAfterCutZJETVETO = 1;
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
				muonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tmuon " << imuon << " rejected because sum of pT of tracks with pT >1.5 within a cone of DR < 0.3 around the muon direction, vetoing a cone of 0.015 around that direction" << endl;
				continue;
			}

			nbMuonsAfterID[9]++;
			TOTALnbMuonsAfterID[9]++;

		double corrected_Pt = mymuon->Pt();
		rochcor *rmcor = new rochcor(seeed);
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
					if( lumi_set == "2011A" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719;
					if( lumi_set == "2011A_rereco" ) integratedLuminosity = 2.221*1000.0;
					if( lumi_set == "2011B" ) integratedLuminosity = 2.714*1000.0;
					if( lumi_set == "2011" ) integratedLuminosity = 215.552 + 951.716 + 389.876 + 706.719 + 2.714*1000.0;
					if( lumi_set == "2011_rereco" ) integratedLuminosity = 2.221*1000.0 +	2.714*1000.0;

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
//						if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
//						else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
				} else {
					rmcor->momcor_data(muonRochester, mymuon->charge(), sysdev, runopt);
//						if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev, runopt);
//						else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev, runopt);
				}
				corrected_Pt = muonRochester.Pt();
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
					rmcor->momcor_mc(muonRochester, mymuon->charge(), sysdev, runopt);
//					if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
//					else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
				} else {
					rmcor->momcor_data(muonRochester, mymuon->charge(), sysdev, runopt);
//					if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev);
//					else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev);
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
					rmcor->momcor_mc(muonRochester, mymuon->charge(), sysdev, runopt);
//					if( mymuon->charge() < 0 ) rmcor->momcor_mc(muonRochester, muonRochesterDummy, 1, sysdev);
//					else rmcor->momcor_mc(muonRochesterDummy, muonRochester, 1, sysdev);
				} else {
					rmcor->momcor_data(muonRochester,	mymuon->charge(), sysdev, runopt);
//					if( mymuon->charge() < 0 ) rmcor->momcor_data(muonRochester, muonRochesterDummy, 1, sysdev);
//					else rmcor->momcor_data(muonRochesterDummy, muonRochester, 1, sysdev);
				}
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
		double corrected_Pz = mymuon->Pz();
		double corrected_Px = applyMuonScaleCorrection > 0 ? mymuon->Px() * (double)(corrected_Pt)/(double)(mymuon->Pt()) : mymuon->Px();
		double corrected_Py = applyMuonScaleCorrection > 0 ? mymuon->Py() * (double)(corrected_Pt)/(double)(mymuon->Pt()) : mymuon->Py();
		double m_mu = 105.658367e-3;
		double corrected_E = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz * corrected_Pz + corrected_Pt * corrected_Pt) ) : mymuon->E();
//		double corrected_E = mymuon->E();
		TLorentzVector correctedMuon(corrected_Px, corrected_Py, corrected_Pz, corrected_E);


		 if(! (correctedMuon.Pt() > 10.0) ){// transverse momentum
//		 if(! (mymuon->Pt() > 10.0) ){// transverse momentum
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
			muonIdentified_corrected_Pt.push_back(correctedMuon.Pt());
			muonIdentified.push_back(imuon);
////			mymuon->Clear();
		
	delete rmcor;
	rmcor = 0;
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
		pair <double, double> PtofMuons[3][numberOfDimuons[0]];

		if(verbosity>2) cout << "initializing dimuon pair object" << endl;
		// Initializing pair object
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < numberOfDimuons[0]; j++) IDofMuons[i][j] = make_pair(0, 0);
			for(int j = 0; j < numberOfDimuons[0]; j++) PtofMuons[i][j] = make_pair(0.0, 0.0);
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
					PtofMuons[0][i_dimuons_] = make_pair(muonIdentified_corrected_Pt[muon_i], muonIdentified_corrected_Pt[muon_j]);
					if(verbosity>5) cerr << "muonIdentified_corrected_Pt["<<muon_i<<"]= " << muonIdentified_corrected_Pt[muon_i] << endl;
					if(verbosity>5) cerr << "muonIdentified_corrected_Pt["<<muon_j<<"]= " << muonIdentified_corrected_Pt[muon_j] << endl;
					if(verbosity>5) cerr << "PtofMuons[0]["<<i_dimuons_<<"]=	(" << PtofMuons[0][i_dimuons_].first << " , " << PtofMuons[0][i_dimuons_].second << ")" << endl;
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
				PtofMuons[1][numberOfDimuons[1]] = make_pair(PtofMuons[0][i_dimuons].first, PtofMuons[0][i_dimuons].second);
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
			double corrected_Pt1 = PtofMuons[1][i_dimuons].first;
			double corrected_Pt2 = PtofMuons[1][i_dimuons].second;
			double corrected_Pz1 = Muon1->Pz();
			double corrected_Pz2 = Muon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? Muon1->Px() * (double)(corrected_Pt1) / (double)(Muon1->Pt()) : Muon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? Muon2->Px() * (double)(corrected_Pt2) / (double)(Muon2->Pt()) : Muon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? Muon1->Py() * (double)(corrected_Pt1) / (double)(Muon1->Pt()) : Muon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? Muon2->Py() * (double)(corrected_Pt2) / (double)(Muon2->Pt()) : Muon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : Muon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : Muon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);

			TLorentzVector mumu;
//			mumu = (*Muon1) + (*Muon2);
			mumu = (*correctedMuon1) + (*correctedMuon2);
			if( (low_m_mumu < mumu.M()) && (mumu.M() < high_m_mumu) )
			{
				IDofMuons[2][numberOfDimuons[2]] = make_pair(IDofMuons[1][i_dimuons].first, IDofMuons[1][i_dimuons].second);
				PtofMuons[2][numberOfDimuons[2]] = make_pair(PtofMuons[1][i_dimuons].first, PtofMuons[1][i_dimuons].second);
				numberOfDimuons[2] += 1;
			}
////			Muon1->Clear();
////			Muon2->Clear();
			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;	
//			correctedMuon1->Delete();
//			correctedMuon2->Delete();
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
////			mumu.Clear();
////			mymuon1->Clear();
////			mymuon2->Clear();
			miniTree_allmuons->Fill();
		}

// ********************************************************************************************************************
// ***** PHOTON OBJECT SELECTION *****
// ********************************************************************************************************************

		if(verbosity>1) cout << "cleaning photon collection" << endl;
	 // Cleaning: removing not commissionned superclusters and spikes
		vector<int> photonsNoSpike;
		photonsNoSpike.clear();
		vector<int> photonIsNotCommissioned;
		if(verbosity>0) cerr << "\t\tThere is " << NbPhotons << " photons in the photon collection" << endl;
		photonIsNotCommissioned.clear();
		int nbPhotonsAfterID[6] = {0};
		nbPhotonsAfterID[0] = NbPhotons;
		TOTALnbPhotonsAfterID[0] += nbPhotonsAfterID[0];

		for(int iphoton=0 ; iphoton<NbPhotons ; iphoton++)
		{
			TRootPhoton *myphoton;
			myphoton = (TRootPhoton*) photons->At(iphoton);

			if( (fabs(myphoton->superCluster()->Eta()))>2.5 ){ // high eta clusters, not commisionned
				photonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tphoton " << iphoton << " rejected because high eta" << endl;
				continue;
			}
			nbPhotonsAfterID[1]++;
			if( nbPhotonsAfterID[1] > nbPhotonsAfterID[0] ) cout << "WHOOOPS!!" << endl;
			TOTALnbPhotonsAfterID[1] += 1;

			if( fabs(myphoton->superCluster()->Eta())>1.4442 && fabs(myphoton->superCluster()->Eta())<1.566 ){// eta gap clusters, not commissionned
				photonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tphoton " << iphoton << " rejected because eta gap" << endl;
				continue;
			}
			nbPhotonsAfterID[2]++;
			TOTALnbPhotonsAfterID[2] += 1;

/*
			if( myphoton->superCluster()->seedSeverity()==4 ){ // kWeird
				photonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tphoton " << iphoton << " rejected because kWeird" << endl;
				continue;
			}
*/
			nbPhotonsAfterID[3]++;
			TOTALnbPhotonsAfterID[3] += 1;
/*
			if( myphoton->superCluster()->seedRecoFlag()==2 ){ // kOutOfTime
				photonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tphoton " << iphoton << " rejected because kOutOfTime" << endl;
				continue;
			}
*/
			nbPhotonsAfterID[4]++;
			TOTALnbPhotonsAfterID[4] += 1;

			if( Photon_scale[iphoton]*(myphoton->Pt()) <= 10.0 ){ // Transverse momentum
				if(verbosity>2) cout << "\t\tmyphoton->Pt()= " << myphoton->Pt() << endl;
				if(verbosity>2) cout << "\t\tPhoton_scale["<<iphoton<<"]= " << Photon_scale[iphoton] << endl;
				photonIsNotCommissioned.push_back(1);
				if(verbosity>0) cerr << "\t\t\tphoton " << iphoton << " rejected because low pt" << endl;
				continue;
			}
			nbPhotonsAfterID[5]++;
			TOTALnbPhotonsAfterID[5] += 1;

// This is now implemented at clustering level
//			if( myphoton->superCluster()->seedSeverity()==5 ){ // kBad
//				photonIsNotCommissioned.push_back(1);
//				continue;
//				if(verbosity>0) cerr << "\t\t\tphoton " << iphoton << " rejected because kBad" << endl;
//			}


			if(verbosity>0) cerr << "\t\t\tphoton " << iphoton << " accepted" << endl;
			photonsNoSpike.push_back(iphoton);
			photonIsNotCommissioned.push_back(0);
		}
		unsigned int NbPhotonsNoSpike = photonsNoSpike.size();


		// Increasing counter
		for(int i = 0; i < 6 ; i++)
		{
			if(nbPhotonsAfterID[i] >= 1){ TOTALnbEventsAfterPhotonID[i]++;}
		}

		if(! (nbPhotonsAfterID[5] >=1) ) // Not enough photon candidates, skipping the event
		{
			if(verbosity>1) cerr << "Not enough photon candidates, skipping the event" << endl;
					continue;

		}
// ********************************************************************************************************************
// ***** MUMUGAMMA FSR TRIPLET OBJECT SELECTION *****
// ********************************************************************************************************************

// --------------------------------------------------------------------------------------------------------------------
// ----- initializing mumugamma triplet candidates holder -----
// --------------------------------------------------------------------------------------------------------------------

		if(verbosity>2) cout << "creating triplet pair object" << endl;
		int nbMuMuGammaAfterID[8] = {0};
		// Build mumugamma candidates
		int nbMuMuGammaCandidates = ( nbPhotonsAfterID[5] )*( numberOfDimuons[2] );
		nbMuMuGammaAfterID[0] = nbMuMuGammaCandidates;
		TOTALnbMuMuGammaAfterID[0] += nbMuMuGammaAfterID[0];
		TOTALnbEventsAfterMuMuGammaID[0]++ ;
		pair <int, pair<int, int> > MuMuGammaCandidates[8][nbMuMuGammaCandidates];
		pair <double, double > MuMuGammaCandidates_corrected[8][nbMuMuGammaCandidates];

		if(verbosity>2) cout << "initializing triplet objects" << endl;
		// Initializing triplet objects
		for(int i = 0; i < 8; i++)
		{
			for(int j = 0; j < nbMuMuGammaCandidates; j++)
			{
				MuMuGammaCandidates[i][j] = make_pair(0, make_pair(0,0));
				MuMuGammaCandidates_corrected[i][j] = make_pair(0.0 , 0.0);
			}
		}

		if(verbosity>2) cout << "Filling triplet objects with ID'ed dimuons and photons" << endl;
		// Filling triplet objects with ID'ed dimuons and photons
		int i_cand = 0;
		for(int i_dimuons = 0; i_dimuons < numberOfDimuons[2]; i_dimuons++)
		{
			for(int i_photon = 0; i_photon < nbPhotonsAfterID[5]; i_photon++)
			{
				if(verbosity>5) cerr << "photonsNoSpike[i_photon]= " << photonsNoSpike[i_photon] << endl;
				if(verbosity>5) cerr << "IDofMuons[2][i_dimuons].first= " << IDofMuons[2][i_dimuons].first << endl;
				if(verbosity>5) cerr << "muonIdentified[IDofMuons[2][i_dimuons].first]= " << muonIdentified[IDofMuons[2][i_dimuons].first] << endl;
				if(verbosity>5) cerr << "IDofMuons[2][i_dimuons].second= " << IDofMuons[2][i_dimuons].second << endl;
				if(verbosity>5) cerr << "muonIdentified[IDofMuons[2][i_dimuons].second]= " << muonIdentified[IDofMuons[2][i_dimuons].second] << endl;
				MuMuGammaCandidates[0][i_cand] = make_pair(photonsNoSpike[i_photon], make_pair( IDofMuons[2][i_dimuons].first, IDofMuons[2][i_dimuons].second));
				MuMuGammaCandidates_corrected[0][i_cand] = make_pair( PtofMuons[2][i_dimuons].first, PtofMuons[2][i_dimuons].second);
				if(verbosity>5) cerr << "PtofMuons[2][i_dimuons].first= " << PtofMuons[2][i_dimuons].first << endl;
				if(verbosity>5) cerr << "PtofMuons[2][i_dimuons].second= " << PtofMuons[2][i_dimuons].second << endl;
				
				i_cand++;
			}
		}	

		if(verbosity>3) cout << "start looping over triplet candidates" << endl;	
		if(verbosity>4) cout << "nbMuMuGammaAfterID[0]= " << nbMuMuGammaAfterID[0] << endl;

		nCandidate[0] = nbMuMuGammaAfterID[0];
		if(verbosity>4) cout << "nCandidate[0]= " << nCandidate[0] << endl;
//		for(int icand = 0 ; icand < nCandidate[0] ; icand++){iCandidate[0] = icand; miniFriend->Fill();}
// --------------------------------------------------------------------------------------------------------------------
// ----- cut on close_isoR03_hadEt -----
// --------------------------------------------------------------------------------------------------------------------
		isMMGCandidate = 1;
		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[0] ; i_mmg++)
		{
			iCandidate_temp[0][i_mmg] = i_mmg;
			if(verbosity>5) cerr << "examining i_mmg= " << i_mmg << " over nbMuMuGammaAfterID[0]= " << nbMuMuGammaAfterID[0] << endl;
			TRootPhoton* myphoton;
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			if(verbosity>5) cerr << "fetch photon photons->At(MuMuGammaCandidates[0][i_mmg].first)" << endl;
			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[0][i_mmg].first);
			if(verbosity>5) cerr << "fetch muon1 muons->At(MuMuGammaCandidates[0][i_mmg].second.first )" << endl;
			mymuon1 = (TRootMuon*) muons->At(MuMuGammaCandidates[0][i_mmg].second.first );
			if(verbosity>5) cerr << "fetch muon2 muons->At( MuMuGammaCandidates[0][i_mmg].second.second )" << endl;
			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[0][i_mmg].second.second );
			if(verbosity>5) cerr << "getting phi, eta" << endl;

			double corrected_Pt1 = MuMuGammaCandidates_corrected[0][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[0][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);

			if(verbosity>5) cerr << "check applyMuonScaleCorrection" << endl;
			if(verbosity>5) cerr << "correctedMuon1->Pt()= " << correctedMuon1->Pt() << endl;
			if(verbosity>5) cerr << "correctedMuon2->Pt()= " << correctedMuon2->Pt() << endl;



			double phiPhoton = myphoton->Phi();
			double etaPhoton = myphoton->Eta();
			double phiMuon1 = mymuon1->Phi();
			double etaMuon1 = mymuon1->Eta();
			double phiMuon2 = mymuon2->Phi();
			double etaMuon2 = mymuon2->Eta();
			if(verbosity>5) cerr << "compute deltaRs" << endl;
			double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
			double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
			if(verbosity>5) cerr << "check what is close_isoR03_hadEt" << endl;
			double close_isoR03_hadEt = (deltaRphomu1 < deltaRphomu2) ? (mymuon1->isoR03_hadEt()) : (mymuon2->isoR03_hadEt());

			if(verbosity>5) cerr << "close_isoR03_hadEt= " << close_isoR03_hadEt << endl;
			if( close_isoR03_hadEt >= 9999.0)
			{
				if(verbosity>3) cerr << "candidate thrown because close_isoR03_hadEt= " << close_isoR03_hadEt << endl;
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[0][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[0][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, Photon_scale[MuMuGammaCandidates[0][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
				iCandidate_temp[1][i_mmg] == -99;
				miniTree->Fill();
				continue;
			} else {
				if(verbosity>5) cerr << "candidate accepted: close_isoR03_hadEt= " << close_isoR03_hadEt << endl;
				iCandidate_temp[1][i_mmg] = nbMuMuGammaAfterID[1];
			}
		
			if(verbosity>5) cerr << " filling new pair " << endl;
			MuMuGammaCandidates[1][nbMuMuGammaAfterID[1]] = make_pair(MuMuGammaCandidates[0][i_mmg].first, make_pair(MuMuGammaCandidates[0][i_mmg].second.first, MuMuGammaCandidates[0][i_mmg].second.second) );
			MuMuGammaCandidates_corrected[1][nbMuMuGammaAfterID[1]] = make_pair(MuMuGammaCandidates_corrected[0][i_mmg].first, MuMuGammaCandidates_corrected[0][i_mmg].second);
			nbMuMuGammaAfterID[1]++;
			TOTALnbMuMuGammaAfterID[1]++;
		
			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;
//			correctedMuon1->Delete();
//			correctedMuon2->Delete();

		}
		if(verbosity>4) cout << "nbMuMuGammaAfterID[1]= " << nbMuMuGammaAfterID[1] << endl;
		nCandidate[1] = nbMuMuGammaAfterID[1];
		if(verbosity>4) cout << "nCandidate[1]= " << nCandidate[1] << endl;
		for(int icand = 0 ; icand < nCandidate[0] ; icand++ )
		{
			if(	iCandidate_temp[1][icand] == -99 )
			{
				iCandidate[0] = icand;
				if(verbosity>4) cout << "\t## event cut: " << "\tiCandidate[0]= " << iCandidate[0] << endl;
				miniFriend->Fill();
			} else {
				iCandidate_temp[1][iCandidate_temp[1][icand]] = iCandidate_temp[1][icand];
				iCandidate_temp[0][iCandidate_temp[1][icand]] = iCandidate_temp[0][icand];
			}
		}
		if(verbosity>4) for( int ia = 0 ; ia < 8 ; ia++ ){for( int ib = 0 ; ib < nCandidate[0] ; ib ++){ cout << "\t" << iCandidate_temp[ia][ib];} cout << endl;}
//		for(int icand = 0 ; icand < nCandidate[0] - nCandidate[1] ; icand++){iCandidate[1] = icand; miniFriend->Fill();}
		if(! (nbMuMuGammaAfterID[1] > 0) )
		{
					if(verbosity>2) cerr << "not enough triplet candidate passing close_isoR03_hadEt cut " << endl;
				continue;

		}

		TOTALnbEventsAfterMuMuGammaID[1]++ ;
		isAfterFSRCut1 = 1;
// --------------------------------------------------------------------------------------------------------------------
// ----- cut on far_isoR03_emEt -----
// --------------------------------------------------------------------------------------------------------------------

//		cout << "loop over nbMuMuGammaAfterID[1]= " << nbMuMuGammaAfterID[1] << " candidates" << endl;
		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[1] ; i_mmg++)
		{
			TRootPhoton* myphoton;
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[1][i_mmg].first);
			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[1][i_mmg].second.first );
			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[1][i_mmg].second.second );


			double corrected_Pt1 = MuMuGammaCandidates_corrected[1][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[1][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);

			if(verbosity > 5) cerr << "corrected_Pt1= " << corrected_Pt1 << endl;
			if(verbosity > 5) cerr << "corrected_Pt2= " << corrected_Pt2 << endl;
			if(verbosity > 5) cerr << "correctedMuon1->Pt()= " << correctedMuon1->Pt() << endl;
			if(verbosity > 5) cerr << "correctedMuon2->Pt()= " << correctedMuon2->Pt() << endl;

			double phiPhoton = myphoton->Phi();
			double etaPhoton = myphoton->Eta();
			double phiMuon1 = mymuon1->Phi();
			double etaMuon1 = mymuon1->Eta();
			double phiMuon2 = mymuon2->Phi();
			double etaMuon2 = mymuon2->Eta();
			double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
			double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
			double far_isoR03_emEt = (deltaRphomu1 > deltaRphomu2) ? mymuon1->isoR03_emEt() : mymuon2->isoR03_emEt();
			if( far_isoR03_emEt >= 9999.0) 
			{
				if(verbosity>4) cerr << "candidate thrown because far_isoR03_emEt= " << far_isoR03_emEt << endl;
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[1][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[1][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, Photon_scale[MuMuGammaCandidates[1][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
				iCandidate_temp[2][i_mmg] == -99;
				miniTree->Fill();
				continue;
			}
			iCandidate_temp[2][i_mmg] = nbMuMuGammaAfterID[2];
			MuMuGammaCandidates[2][nbMuMuGammaAfterID[2]] = make_pair(MuMuGammaCandidates[1][i_mmg].first, make_pair(MuMuGammaCandidates[1][i_mmg].second.first, MuMuGammaCandidates[1][i_mmg].second.second) );
			MuMuGammaCandidates_corrected[2][nbMuMuGammaAfterID[2]] = make_pair(MuMuGammaCandidates_corrected[1][i_mmg].first, MuMuGammaCandidates_corrected[1][i_mmg].second);
			nbMuMuGammaAfterID[2]++;
			

			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;
//			correctedMuon1->Delete();
//			correctedMuon2->Delete();


		 }
		if(verbosity>4) cout << "nbMuMuGammaAfterID[2]= " << nbMuMuGammaAfterID[2] << endl;
		nCandidate[2] = nbMuMuGammaAfterID[2];
		if(verbosity>4) cout << "nCandidate[2]= " << nCandidate[2] << endl;
		for(int icand = 0 ; icand < nCandidate[1] ; icand++ )
		{
			if(	iCandidate_temp[2][icand] == -99 )
			{
				iCandidate[1] = icand;
				iCandidate[0] = iCandidate_temp[0][icand];
				if(verbosity>4) cout << "\t## event cut: ";
				if(verbosity>4) cout << "\tiCandidate[1]= " << iCandidate[1];
				if(verbosity>4) cout << "\tiCandidate[0]= " << iCandidate[0];
				if(verbosity>4) cout << endl;
				miniFriend->Fill();
			} else {
				iCandidate_temp[2][iCandidate_temp[2][icand]] = iCandidate_temp[2][icand];
				iCandidate_temp[1][iCandidate_temp[2][icand]] = iCandidate_temp[1][icand];
				iCandidate_temp[0][iCandidate_temp[2][icand]] = iCandidate_temp[0][icand];
			}
		}
		if(verbosity>4) for( int ia = 0 ; ia < 8 ; ia++ ){for( int ib = 0 ; ib < nCandidate[0] ; ib ++){ cout << "\t" << iCandidate_temp[ia][ib];} cout << endl;}
//		for(int icand = 0 ; icand < nCandidate[1] - nCandidate[2] ; icand++){iCandidate[2] = icand; miniFriend->Fill();}
		if(! (nbMuMuGammaAfterID[2] > 0) )
		{
			continue;
		}

		TOTALnbMuMuGammaAfterID[2] += nbMuMuGammaAfterID[2];
		TOTALnbEventsAfterMuMuGammaID[2]++ ;
		isAfterFSRCut2 = 1;
// --------------------------------------------------------------------------------------------------------------------
// ----- cut on min_DeltaR -----
// --------------------------------------------------------------------------------------------------------------------

//		cout << "loop over nbMuMuGammaAfterID[2]= " << nbMuMuGammaAfterID[2] << " candidates" << endl;
		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[2] ; i_mmg++)
		{
			TRootPhoton* myphoton;
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[2][i_mmg].first);
			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[2][i_mmg].second.first );
			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[2][i_mmg].second.second );

			double corrected_Pt1 = MuMuGammaCandidates_corrected[2][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[2][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);




			double phiPhoton = myphoton->Phi();
			double etaPhoton = myphoton->Eta();
			double phiMuon1 = mymuon1->Phi();
			double etaMuon1 = mymuon1->Eta();
			double phiMuon2 = mymuon2->Phi();
			double etaMuon2 = mymuon2->Eta();
			double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
			double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
			double min_DeltaR = min(deltaRphomu1, deltaRphomu2);
			if( min_DeltaR >= 9999.0 )
			{
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[2][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[2][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, Photon_scale[MuMuGammaCandidates[2][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
				iCandidate_temp[3][i_mmg] == -99;
				miniTree->Fill();
				continue;
			}
			iCandidate_temp[3][i_mmg] = nbMuMuGammaAfterID[3];
			MuMuGammaCandidates[3][nbMuMuGammaAfterID[3]] = make_pair(MuMuGammaCandidates[2][i_mmg].first, make_pair(MuMuGammaCandidates[2][i_mmg].second.first, MuMuGammaCandidates[2][i_mmg].second.second) );
			MuMuGammaCandidates_corrected[3][nbMuMuGammaAfterID[3]] = make_pair(MuMuGammaCandidates_corrected[2][i_mmg].first, MuMuGammaCandidates_corrected[2][i_mmg].second);
			nbMuMuGammaAfterID[3]++;

			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;
//			correctedMuon1->Delete();
//			correctedMuon2->Delete();


		}
		nCandidate[3] = nbMuMuGammaAfterID[3];
		if(verbosity>4) cout << "nCandidate[3]= " << nCandidate[3] << endl;
		for(int icand = 0 ; icand < nCandidate[2] ; icand++ )
		{
			if(	iCandidate_temp[3][icand] == -99 )
			{
				iCandidate[2] = icand;
				iCandidate[1] = iCandidate_temp[1][icand];
				iCandidate[0] = iCandidate_temp[0][icand];
				if(verbosity>4) cout << "\t## event cut: ";
				if(verbosity>4) cout << "\tiCandidate[2]= " << iCandidate[2];
				if(verbosity>4) cout << "\tiCandidate[1]= " << iCandidate[1];
				if(verbosity>4) cout << "\tiCandidate[0]= " << iCandidate[0];
				if(verbosity>4) cout << endl;
				miniFriend->Fill();
			} else {
				iCandidate_temp[3][iCandidate_temp[3][icand]] = iCandidate_temp[3][icand];
				iCandidate_temp[2][iCandidate_temp[3][icand]] = iCandidate_temp[2][icand];
				iCandidate_temp[1][iCandidate_temp[3][icand]] = iCandidate_temp[1][icand];
				iCandidate_temp[0][iCandidate_temp[3][icand]] = iCandidate_temp[0][icand];
			}
		}
		if(verbosity>4) for( int ia = 0 ; ia < 8 ; ia++ ){for( int ib = 0 ; ib < nCandidate[0] ; ib ++){ cout << "\t" << iCandidate_temp[ia][ib];} cout << endl;}
//		for(int icand = 0 ; icand < nCandidate[2] - nCandidate[3] ; icand++){iCandidate[3] = icand; miniFriend->Fill();}
		if(! (nbMuMuGammaAfterID[3] > 0) )
		{
				 continue;
		}

		TOTALnbMuMuGammaAfterID[3] += nbMuMuGammaAfterID[3];
		TOTALnbEventsAfterMuMuGammaID[3]++ ;
		isAfterFSRCut3 = 1;
// --------------------------------------------------------------------------------------------------------------------
// ----- cut on far_muonPt -----
// --------------------------------------------------------------------------------------------------------------------

//		cout << "loop over nbMuMuGammaAfterID[3]= " << nbMuMuGammaAfterID[3] << " candidates" << endl;
		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[3] ; i_mmg++)
		{
			TRootPhoton* myphoton;
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[3][i_mmg].first);
			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[3][i_mmg].second.first );
			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[3][i_mmg].second.second );

			double corrected_Pt1 = MuMuGammaCandidates_corrected[3][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[3][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);



			double phiPhoton = myphoton->Phi();
			double etaPhoton = myphoton->Eta();
			double phiMuon1 = mymuon1->Phi();
			double etaMuon1 = mymuon1->Eta();
			double phiMuon2 = mymuon2->Phi();
			double etaMuon2 = mymuon2->Eta();
			double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
			double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
			double far_muonPt = (deltaRphomu1 > deltaRphomu2) ? mymuon1->Pt() : mymuon2->Pt();
			if( far_muonPt <= -5.0 )
			{
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[3][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[3][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, Photon_scale[MuMuGammaCandidates[3][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
				iCandidate_temp[4][i_mmg] == -99;
				miniTree->Fill();
				continue;
			}
			iCandidate_temp[4][i_mmg] = nbMuMuGammaAfterID[4];
			MuMuGammaCandidates[4][nbMuMuGammaAfterID[4]] = make_pair(MuMuGammaCandidates[3][i_mmg].first, make_pair(MuMuGammaCandidates[3][i_mmg].second.first, MuMuGammaCandidates[3][i_mmg].second.second) );
			MuMuGammaCandidates_corrected[4][nbMuMuGammaAfterID[4]] = make_pair(MuMuGammaCandidates_corrected[3][i_mmg].first, MuMuGammaCandidates_corrected[3][i_mmg].second);
			nbMuMuGammaAfterID[4]++;

			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;
//			correctedMuon1->Delete();
//			correctedMuon2->Delete();

		}
		nCandidate[4] = nbMuMuGammaAfterID[4];
		if(verbosity>4) cout << "nCandidate[4]= " << nCandidate[4] << endl;
		for(int icand = 0 ; icand < nCandidate[3] ; icand++ )
		{
			if(	iCandidate_temp[4][icand] == -99 )
			{
				iCandidate[3] = icand;
				iCandidate[2] = iCandidate_temp[2][icand];
				iCandidate[1] = iCandidate_temp[1][icand];
				iCandidate[0] = iCandidate_temp[0][icand];
				if(verbosity>4) cout << "\t## event cut: ";
				if(verbosity>4) cout << "\tiCandidate[3]= " << iCandidate[3];
				if(verbosity>4) cout << "\tiCandidate[2]= " << iCandidate[2];
				if(verbosity>4) cout << "\tiCandidate[1]= " << iCandidate[1];
				if(verbosity>4) cout << "\tiCandidate[0]= " << iCandidate[0];
				if(verbosity>4) cout << endl;
				miniFriend->Fill();
			} else {
				iCandidate_temp[4][iCandidate_temp[4][icand]] = iCandidate_temp[4][icand];
				iCandidate_temp[3][iCandidate_temp[4][icand]] = iCandidate_temp[3][icand];
				iCandidate_temp[2][iCandidate_temp[4][icand]] = iCandidate_temp[2][icand];
				iCandidate_temp[1][iCandidate_temp[4][icand]] = iCandidate_temp[1][icand];
				iCandidate_temp[0][iCandidate_temp[4][icand]] = iCandidate_temp[0][icand];
			}
		}
		if(verbosity>4) for( int ia = 0 ; ia < 8 ; ia++ ){for( int ib = 0 ; ib < nCandidate[0] ; ib ++){ cout << "\t" << iCandidate_temp[ia][ib];} cout << endl;}
//		for(int icand = 0 ; icand < nCandidate[3] - nCandidate[4] ; icand++){iCandidate[4] = icand; miniFriend->Fill();}
		if(! (nbMuMuGammaAfterID[4] > 0) )
		{
			continue;

		}

		TOTALnbMuMuGammaAfterID[4] += nbMuMuGammaAfterID[4];
		TOTALnbEventsAfterMuMuGammaID[4]++ ;
		isAfterFSRCut4 = 1;
// --------------------------------------------------------------------------------------------------------------------
// ----- cut on mumugamma.M() (very loose window) -----
// --------------------------------------------------------------------------------------------------------------------

		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[4] ; i_mmg++)
		{
			TRootPhoton* myphoton;
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[4][i_mmg].first);
			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[4][i_mmg].second.first );
			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[4][i_mmg].second.second );

			double corrected_Pt1 = MuMuGammaCandidates_corrected[4][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[4][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);


			double phiPhoton = myphoton->Phi();
			double etaPhoton = myphoton->Eta();
			double phiMuon1 = mymuon1->Phi();
			double etaMuon1 = mymuon1->Eta();
			double phiMuon2 = mymuon2->Phi();
			double etaMuon2 = mymuon2->Eta();
			double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
			double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
			TLorentzVector mumugamma;
			EScale = Photon_scale[MuMuGammaCandidates[4][i_mmg].first];
			TLorentzVector *PhotonEScale = new TLorentzVector( EScale*(myphoton->Px()), EScale*(myphoton->Py()), EScale*(myphoton->Pz()), EScale*(myphoton->Energy()));
//			mumugamma = (*PhotonEScale) + (*mymuon1) + (*mymuon2);
			mumugamma = (*PhotonEScale) + (*correctedMuon1) + (*correctedMuon2);
			if( (mumugamma.M() < 0.0) || (9999.0 < mumugamma.M())	)
			{
//				FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, EScale, doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
				iCandidate_temp[5][i_mmg] == -99;
				miniTree->Fill();
				continue;
			}
			iCandidate_temp[5][i_mmg] = nbMuMuGammaAfterID[5];
			MuMuGammaCandidates[5][nbMuMuGammaAfterID[5]] = make_pair(MuMuGammaCandidates[4][i_mmg].first, make_pair(MuMuGammaCandidates[4][i_mmg].second.first, MuMuGammaCandidates[4][i_mmg].second.second) );
			MuMuGammaCandidates_corrected[5][nbMuMuGammaAfterID[5]] = make_pair(MuMuGammaCandidates_corrected[4][i_mmg].first, MuMuGammaCandidates_corrected[4][i_mmg].second);
			nbMuMuGammaAfterID[5]++;
		
			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;
			delete PhotonEScale;
			PhotonEScale = 0;
//			correctedMuon1->Delete();
//			correctedMuon2->Delete();
//			PhotonEScale->Delete();
			
		}
		nCandidate[5] = nbMuMuGammaAfterID[5];
		if(verbosity>4) cout << "nCandidate[5]= " << nCandidate[5] << endl;
		for(int icand = 0 ; icand < nCandidate[4] ; icand++ )
		{
			if(	iCandidate_temp[5][icand] == -99 )
			{
				iCandidate[4] = icand;
				iCandidate[3] = iCandidate_temp[3][icand];
				iCandidate[2] = iCandidate_temp[2][icand];
				iCandidate[1] = iCandidate_temp[1][icand];
				iCandidate[0] = iCandidate_temp[0][icand];
				if(verbosity>4) cout << "\t## event cut: ";
				if(verbosity>4) cout << "\tiCandidate[4]= " << iCandidate[4];
				if(verbosity>4) cout << "\tiCandidate[3]= " << iCandidate[3];
				if(verbosity>4) cout << "\tiCandidate[2]= " << iCandidate[2];
				if(verbosity>4) cout << "\tiCandidate[1]= " << iCandidate[1];
				if(verbosity>4) cout << "\tiCandidate[0]= " << iCandidate[0];
				if(verbosity>4) cout << endl;
				miniFriend->Fill();
			} else {
				iCandidate_temp[5][iCandidate_temp[5][icand]] = iCandidate_temp[5][icand];
				iCandidate_temp[4][iCandidate_temp[5][icand]] = iCandidate_temp[4][icand];
				iCandidate_temp[3][iCandidate_temp[5][icand]] = iCandidate_temp[3][icand];
				iCandidate_temp[2][iCandidate_temp[5][icand]] = iCandidate_temp[2][icand];
				iCandidate_temp[1][iCandidate_temp[5][icand]] = iCandidate_temp[1][icand];
				iCandidate_temp[0][iCandidate_temp[5][icand]] = iCandidate_temp[0][icand];
			}
		}
		if(verbosity>4) for( int ia = 0 ; ia < 8 ; ia++ ){for( int ib = 0 ; ib < nCandidate[0] ; ib ++){ cout << "\t" << iCandidate_temp[ia][ib];} cout << endl;}
//		for(int icand = 0 ; icand < nCandidate[4] - nCandidate[5] ; icand++){iCandidate[5] = icand; miniFriend->Fill();}
		if(! (nbMuMuGammaAfterID[5] > 0) )
		{
//			miniTree->Fill();
			continue;
		}
//		for(int i_mmg = 0; i_mmg < 1 ; i_mmg++)
//		{
//			TRootPhoton* myphoton;
//			TRootMuon* mymuon1;
//			TRootMuon* mymuon2;
//			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[5][i_mmg].first);
//			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[5][i_mmg].second.first );
//			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[5][i_mmg].second.second );
//			FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//		}
		TOTALnbMuMuGammaAfterID[5] += nbMuMuGammaAfterID[5];
		TOTALnbEventsAfterMuMuGammaID[5]++ ;
			isVeryLooseMMG = 1;
			nAfterVeryLooseMMG++;
// --------------------------------------------------------------------------------------------------------------------
// ----- cut on mumugamma.M() (loose window) -----
// --------------------------------------------------------------------------------------------------------------------


		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[5] ; i_mmg++)
		{
			TRootPhoton* myphoton;
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[5][i_mmg].first);
			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[5][i_mmg].second.first );
			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[5][i_mmg].second.second );

			double corrected_Pt1 = MuMuGammaCandidates_corrected[5][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[5][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);


			double phiPhoton = myphoton->Phi();
			double etaPhoton = myphoton->Eta();
			double phiMuon1 = mymuon1->Phi();
			double etaMuon1 = mymuon1->Eta();
			double phiMuon2 = mymuon2->Phi();
			double etaMuon2 = mymuon2->Eta();
			double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
			double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
			TLorentzVector mumugamma;
			EScale = Photon_scale[MuMuGammaCandidates[5][i_mmg].first];
			TLorentzVector *PhotonEScale = new TLorentzVector( EScale*(myphoton->Px()), EScale*(myphoton->Py()), EScale*(myphoton->Pz()), EScale*(myphoton->Energy()));
			mumugamma = (*PhotonEScale) + (*mymuon1) + (*mymuon2);
			if( (mumugamma.M() < 0.0) || (9999.0 < mumugamma.M())	)
//			if( (mumugamma.M() < 70.0) || (110.0 < mumugamma.M())	)
			{
//				cout << "non-loose event: rejected: mumugamma.M()= " << mumugamma.M() << endl;
//				FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, EScale, doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
				iCandidate_temp[6][i_mmg] == -99;
				miniTree->Fill();
				continue;
			}
			iCandidate_temp[6][i_mmg] = nbMuMuGammaAfterID[6];
			MuMuGammaCandidates[6][nbMuMuGammaAfterID[6]] = make_pair(MuMuGammaCandidates[5][i_mmg].first, make_pair(MuMuGammaCandidates[5][i_mmg].second.first, MuMuGammaCandidates[5][i_mmg].second.second) );
			MuMuGammaCandidates_corrected[6][nbMuMuGammaAfterID[6]] = make_pair(MuMuGammaCandidates_corrected[5][i_mmg].first, MuMuGammaCandidates_corrected[5][i_mmg].second);
			nbMuMuGammaAfterID[6]++;
 
	 
			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;
			delete PhotonEScale;
			PhotonEScale = 0;
//			correctedMuon1->Delete();
//			correctedMuon2->Delete();
//			PhotonEScale->Delete();



		}
		nCandidate[6] = nbMuMuGammaAfterID[6];
		if(verbosity>4) cout << "nCandidate[6]= " << nCandidate[6] << endl;
		for(int icand = 0 ; icand < nCandidate[5] ; icand++ )
		{
			if(	iCandidate_temp[6][icand] == -99 )
			{
				iCandidate[5] = icand;
				iCandidate[4] = iCandidate_temp[4][icand];
				iCandidate[3] = iCandidate_temp[3][icand];
				iCandidate[2] = iCandidate_temp[2][icand];
				iCandidate[1] = iCandidate_temp[1][icand];
				iCandidate[0] = iCandidate_temp[0][icand];
				if(verbosity>4) cout << "\t## event cut: ";
				if(verbosity>4) cout << "\tiCandidate[5]= " << iCandidate[5];
				if(verbosity>4) cout << "\tiCandidate[4]= " << iCandidate[4];
				if(verbosity>4) cout << "\tiCandidate[3]= " << iCandidate[3];
				if(verbosity>4) cout << "\tiCandidate[2]= " << iCandidate[2];
				if(verbosity>4) cout << "\tiCandidate[1]= " << iCandidate[1];
				if(verbosity>4) cout << "\tiCandidate[0]= " << iCandidate[0];
				if(verbosity>4) cout << endl;
				miniFriend->Fill();
			} else {
				iCandidate_temp[6][iCandidate_temp[6][icand]] = iCandidate_temp[6][icand];
				iCandidate_temp[5][iCandidate_temp[6][icand]] = iCandidate_temp[5][icand];
				iCandidate_temp[4][iCandidate_temp[6][icand]] = iCandidate_temp[4][icand];
				iCandidate_temp[3][iCandidate_temp[6][icand]] = iCandidate_temp[3][icand];
				iCandidate_temp[2][iCandidate_temp[6][icand]] = iCandidate_temp[2][icand];
				iCandidate_temp[1][iCandidate_temp[6][icand]] = iCandidate_temp[1][icand];
				iCandidate_temp[0][iCandidate_temp[6][icand]] = iCandidate_temp[0][icand];
			}
		}
		if(verbosity>4) for( int ia = 0 ; ia < 8 ; ia++ ){for( int ib = 0 ; ib < nCandidate[0] ; ib ++){ cout << "\t" << iCandidate_temp[ia][ib];} cout << endl;}
//		for(int icand = 0 ; icand <	nCandidate[5] - nCandidate[6] ; icand++){iCandidate[6] = icand; miniFriend->Fill();}
		if(! (nbMuMuGammaAfterID[6] > 0) )
		{
//			miniTree->Fill();
			continue;
		}
			isLooseMMG = 1;
			nAfterLooseMMG++;

		TOTALnbMuMuGammaAfterID[6] += nbMuMuGammaAfterID[6];
		TOTALnbEventsAfterMuMuGammaID[6]++ ;

//		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[6] ; i_mmg++)
//		for(int i_mmg = 0; i_mmg < 1 ; i_mmg++)
//		{
//			TRootPhoton* myphoton;
//			TRootMuon* mymuon1;
//			TRootMuon* mymuon2;
//			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[6][i_mmg].first);
//			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[6][i_mmg].second.first );
//			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[6][i_mmg].second.second );
//			FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//		}
// --------------------------------------------------------------------------------------------------------------------
// ----- cut on mumugamma.M() (tight window) -----
// --------------------------------------------------------------------------------------------------------------------


		for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[6] ; i_mmg++)
		{
			TRootPhoton* myphoton;
			TRootMuon* mymuon1;
			TRootMuon* mymuon2;
			myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[6][i_mmg].first);
			mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[6][i_mmg].second.first );
			mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[6][i_mmg].second.second );

			double corrected_Pt1 = MuMuGammaCandidates_corrected[6][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[6][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);



			double phiPhoton = myphoton->Phi();
			double etaPhoton = myphoton->Eta();
			double phiMuon1 = mymuon1->Phi();
			double etaMuon1 = mymuon1->Eta();
			double phiMuon2 = mymuon2->Phi();
			double etaMuon2 = mymuon2->Eta();
			double deltaRphomu1 = DeltaR(etaPhoton, phiPhoton, etaMuon1, phiMuon1);
			double deltaRphomu2 = DeltaR(etaPhoton, phiPhoton, etaMuon2, phiMuon2);
			TLorentzVector mumugamma;
			EScale = Photon_scale[MuMuGammaCandidates[6][i_mmg].first];
			TLorentzVector *PhotonEScale = new TLorentzVector( EScale*(myphoton->Px()), EScale*(myphoton->Py()), EScale*(myphoton->Pz()), EScale*(myphoton->Energy()));
			mumugamma = (*PhotonEScale) + (*mymuon1) + (*mymuon2);
			if( (mumugamma.M() < 0.0) || (9999.0 < mumugamma.M())	)
//			if( (mumugamma.M() < 87.2) || (95.2 < mumugamma.M())	)
			{
//				cout << "non-tight event: rejected: mumugamma.M()= " << mumugamma.M() << endl;
//				cout << "*** isVeryLooseMMG:isLooseMMG:isTightMMG= " << isVeryLooseMMG << isLooseMMG << isTightMMG << endl;
//				FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, EScale, doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
//				FillMMG(myphoton, mymuon1, mymuon2, EScale, doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				cout << "*** isVeryLooseMMG:isLooseMMG:isTightMMG= " << isVeryLooseMMG << isLooseMMG << isTightMMG << endl;
				iCandidate_temp[7][i_mmg] == -99;
				miniTree->Fill();
				continue;
			}
			iCandidate_temp[7][i_mmg] = nbMuMuGammaAfterID[7];
			MuMuGammaCandidates[7][nbMuMuGammaAfterID[7]] = make_pair(MuMuGammaCandidates[6][i_mmg].first, make_pair(MuMuGammaCandidates[6][i_mmg].second.first, MuMuGammaCandidates[6][i_mmg].second.second) );
			MuMuGammaCandidates_corrected[7][nbMuMuGammaAfterID[7]] = make_pair(MuMuGammaCandidates_corrected[6][i_mmg].first, MuMuGammaCandidates_corrected[6][i_mmg].second);
			nbMuMuGammaAfterID[7]++;
 
	
			delete correctedMuon1;
			correctedMuon1 = 0;
			delete correctedMuon2;
			correctedMuon2 = 0;
			delete PhotonEScale;
			PhotonEScale = 0;
/*
			correctedMuon1->Delete();
			correctedMuon2->Delete();
			PhotonEScale->Delete();
*/


		}
		nCandidate[7] = nbMuMuGammaAfterID[7];
		if(verbosity>4) cout << "nCandidate[7]= " << nCandidate[7] << endl;
		for(int icand = 0 ; icand < nCandidate[6] ; icand++ )
		{
			if(	iCandidate_temp[7][icand] == -99 )
			{
				iCandidate[6] = icand;
				iCandidate[5] = iCandidate_temp[5][icand];
				iCandidate[4] = iCandidate_temp[4][icand];
				iCandidate[3] = iCandidate_temp[3][icand];
				iCandidate[2] = iCandidate_temp[2][icand];
				iCandidate[1] = iCandidate_temp[1][icand];
				iCandidate[0] = iCandidate_temp[0][icand];
				if(verbosity>4) cout << "\t## event cut: ";
				if(verbosity>4) cout << "\tiCandidate[6]= " << iCandidate[6];
				if(verbosity>4) cout << "\tiCandidate[5]= " << iCandidate[5];
				if(verbosity>4) cout << "\tiCandidate[4]= " << iCandidate[4];
				if(verbosity>4) cout << "\tiCandidate[3]= " << iCandidate[3];
				if(verbosity>4) cout << "\tiCandidate[2]= " << iCandidate[2];
				if(verbosity>4) cout << "\tiCandidate[1]= " << iCandidate[1];
				if(verbosity>4) cout << "\tiCandidate[0]= " << iCandidate[0];
				if(verbosity>4) cout << endl;
				miniFriend->Fill();
			} else {
/*
				iCandidate_temp[7][iCandidate_temp[7][icand]] = iCandidate_temp[7][icand];
				iCandidate_temp[6][iCandidate_temp[7][icand]] = iCandidate_temp[6][icand];
				iCandidate_temp[5][iCandidate_temp[7][icand]] = iCandidate_temp[5][icand];
				iCandidate_temp[4][iCandidate_temp[7][icand]] = iCandidate_temp[4][icand];
				iCandidate_temp[3][iCandidate_temp[7][icand]] = iCandidate_temp[3][icand];
				iCandidate_temp[2][iCandidate_temp[7][icand]] = iCandidate_temp[2][icand];
				iCandidate_temp[1][iCandidate_temp[7][icand]] = iCandidate_temp[1][icand];
				iCandidate_temp[0][iCandidate_temp[7][icand]] = iCandidate_temp[0][icand];
*/
				iCandidate[7] = iCandidate_temp[7][icand];
				iCandidate[6] = iCandidate_temp[6][icand];
				iCandidate[5] = iCandidate_temp[5][icand];
				iCandidate[4] = iCandidate_temp[4][icand];
				iCandidate[3] = iCandidate_temp[3][icand];
				iCandidate[2] = iCandidate_temp[2][icand];
				iCandidate[1] = iCandidate_temp[1][icand];
				iCandidate[0] = iCandidate_temp[0][icand];
				if(verbosity>4) cout << "\t## event PASS: ";
				if(verbosity>4) cout << "\tiCandidate[7]= " << iCandidate[7];
				if(verbosity>4) cout << "\tiCandidate[6]= " << iCandidate[6];
				if(verbosity>4) cout << "\tiCandidate[5]= " << iCandidate[5];
				if(verbosity>4) cout << "\tiCandidate[4]= " << iCandidate[4];
				if(verbosity>4) cout << "\tiCandidate[3]= " << iCandidate[3];
				if(verbosity>4) cout << "\tiCandidate[2]= " << iCandidate[2];
				if(verbosity>4) cout << "\tiCandidate[1]= " << iCandidate[1];
				if(verbosity>4) cout << "\tiCandidate[0]= " << iCandidate[0];
				if(verbosity>4) cout << endl;
				miniFriend->Fill();
			}
		}
		if(verbosity>4) for( int ia = 0 ; ia < 8 ; ia++ ){for( int ib = 0 ; ib < nCandidate[0] ; ib ++){ cout << "\t" << iCandidate_temp[ia][ib];} cout << endl;}
//		for(int icand = 0 ; icand <	nCandidate[6] - nCandidate[7] ; icand++){iCandidate[7] = icand; miniFriend->Fill();}
		 if(! (nbMuMuGammaAfterID[7] > 0) )
		{
			continue;
		}

		TOTALnbMuMuGammaAfterID[7] += nbMuMuGammaAfterID[7];
		TOTALnbEventsAfterMuMuGammaID[7]++ ;
//		nCandidate[7] = nbMuMuGammaAfterID[7];
//		cout << "nCandidate[7]= " << nCandidate[7] << endl;
//		for(int icand = 0 ; icand < nCandidate[7] ; icand++){iCandidate[7] = icand; miniFriend->Fill();}

		isTightMMG = 1;
		nAfterTightMMG++;

		nAfterCut1e++;

		nAfterCut2a++;

		nAfterCut2b++;

		nAfterCut2c++;
		if( nbMuMuGammaAfterID[7] != 1 ){
			isMultipleCandidate = 1;
			for(int i_mmg = 0; i_mmg < nbMuMuGammaAfterID[7] ; i_mmg++)
			{
				TRootPhoton* myphoton;
				TRootMuon* mymuon1;
				TRootMuon* mymuon2;
				myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[7][i_mmg].first);
				mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[7][i_mmg].second.first );
				mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[7][i_mmg].second.second );

			double corrected_Pt1 = MuMuGammaCandidates_corrected[7][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[7][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);



//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[7][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//				FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[7][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader);
				FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, Photon_scale[MuMuGammaCandidates[7][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
				miniTree->Fill();
	
			delete correctedMuon1;
		 	correctedMuon1 = 0;
		 	delete correctedMuon2;
		 	correctedMuon2 = 0;
/*
			correctedMuon1->Delete();
			correctedMuon2->Delete();
*/

			}
			continue;
		}

		int i_mmg = 0;
		isMultipleCandidate = 0;
		TRootPhoton* myphoton;
		TRootMuon* mymuon1;
		TRootMuon* mymuon2;
		myphoton = (TRootPhoton*) photons->At(MuMuGammaCandidates[7][i_mmg].first);
		mymuon1 = (TRootMuon*) muons->At( MuMuGammaCandidates[7][i_mmg].second.first );
		mymuon2 = (TRootMuon*) muons->At( MuMuGammaCandidates[7][i_mmg].second.second );

			double corrected_Pt1 = MuMuGammaCandidates_corrected[7][i_mmg].first;
			double corrected_Pt2 = MuMuGammaCandidates_corrected[7][i_mmg].second;
			double corrected_Pz1 = mymuon1->Pz();
			double corrected_Pz2 = mymuon2->Pz();
			double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
			double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
			double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
			double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
			double m_mu = 105.658367e-3;
			double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1) ) : mymuon1->E();
			double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2) ) : mymuon2->E();
			TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
			TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);


//		FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[7][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, mcPhotons, reader);
//		FillMMG(myphoton, mymuon1, mymuon2, Photon_scale[MuMuGammaCandidates[7][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader);
		FillMMG(myphoton, mymuon1, mymuon2, correctedMuon1, correctedMuon2, Photon_scale[MuMuGammaCandidates[7][i_mmg].first], doMC, doPhotonConversionMC, mcParticles, reader, binNumber);
		isSelected = 1;
		nSelected++;
		cerr << "OK: Surviving veto event: "<< ievt << " ( " << iRunID << " , " << iLumiID << " , " << iEventID << " )"	<< endl;
		miniTree->Fill();

		nAfterCut3++;

		nAfterCut4++;

		isAfterCut5 = 1;
		nAfterCut5++;

		isAfterCut6 = 1;
		nAfterCut6++;

		isAfterCut7 = 1;
		nAfterCut7++;

		isAfterCut8 = 1;
		nAfterCut8++;
		isAfterCut9 = 1;
		nAfterCut9++;
		isAfterCut10 = 1;
		nAfterCut10++;
	
	
		delete correctedMuon1;
	 	correctedMuon1 = 0;
	 	delete correctedMuon2;
	 	correctedMuon2 = 0;
/*
			correctedMuon1->Delete();
			correctedMuon2->Delete();
*/

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
		for(int i = 0; i < 6 ; i++)
		{
			cout << "TOTALnbPhotonsAfterID["<<i<<"]=\t" << TOTALnbPhotonsAfterID[i] << "\t\t" << "TOTALnbEventsAfterPhotonID["<<i<<"]=\t" << TOTALnbEventsAfterPhotonID[i] << endl;
		}
		cout << endl;
		for(int i = 0; i < 8 ; i++)
		{
			cout << "TOTALnbMuMuGammaAfterID["<<i<<"]=\t" << TOTALnbMuMuGammaAfterID[i] << "\t\t" << "TOTALnbEventsAfterMuMuGammaID["<<i<<"]=\t" << TOTALnbEventsAfterMuMuGammaID[i] << endl;
			if(i == 6) cout << endl;
		}


	cout << "Writing stuff out" << endl;
	// Writing stuff out
	OutputRootFile->Write();
	OutputFriendFile->Write();

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

	// Deleting TMVA Reader
	delete reader;
	reader = 0;

	// Deleting input trees	
	delete inputEventTree;
	inputEventTree = 0;
	delete inputRunTree;
	inputRunTree = 0;

	// Destructing the tree!
	// So its a bit of a mess: what needs to be done is: 
	// - SetDirectory(0) to tell the TTree to forget about the TFile that owns it
	miniTree->SetDirectory(0);
	miniTree_allmuons->SetDirectory(0);
	miniTree_allphotons->SetDirectory(0);
	// - delete the TTree
	delete miniFriend;
	miniFriend = 0;
	delete miniTree;
	miniTree = 0;
	delete miniTree_allmuons;
	miniTree_allmuons = 0;
	delete miniTree_allphotons;
	miniTree_allphotons = 0;
	// - close the TFile
	OutputRootFile->Close();
	OutputFriendFile->Close();
	// - delete the TFile
	delete OutputFriendFile;
	OutputFriendFile = 0;
	delete OutputRootFile;
	OutputRootFile = 0;


	// Exit!
	return 0;
}

