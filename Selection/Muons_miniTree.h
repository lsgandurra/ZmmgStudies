// ROOT headers
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
#include "TRandom3.h"
#include "TTimeStamp.h"
#include "TSystem.h"
//#include "TProof.h"

// C++ headers
#include <sstream>
#include <iostream>
#include <fstream>
#include <utility>

// IpnTreeProducer headers
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

// Muon corrections headers


using namespace std;

// *****************************************************************************************************
// ******************* Compute pile-up weights
// *****************************************************************************************************

double weight_DYToMuMu(int nGenVertices, string lumi_set, string pu_set)
{
	if( lumi_set == "May10" )
	{
		if( pu_set == "PU_S6" )
		{
			double weight[36] = {0, 0.721319, 1.51119, 2.04427, 2.48331, 2.63479, 2.56301, 2.28339, 1.88574, 1.4547, 1.05749, 0.722332, 0.469077, 0.287616, 0.167336, 0.0926727, 0.0496952, 0.0256116, 0.012822, 0.00632371, 0.00302901, 0.00143204, 0.000671079, 0.000309772, 0.000141109, 6.43521e-05, 2.90263e-05, 1.29881e-05, 5.80744e-06, 2.57824e-06, 1.13896e-06, 5.00778e-07, 2.16283e-07, 9.38425e-08, 3.9934e-08, 2.51262e-08
//			double weight[25] = {0, 0.00557876, 0.0720262, 0.259919, 1.00233, 4.0137, 4.95845, 3.81446, 1.65209, 0.43272, 0.106763, 0.028188, 0.00586923, 0.00130453, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0575699, 0.456986, 0.938683, 1.55675, 2.01064, 2.1879, 2.08267, 1.7888, 1.40147, 1.03889, 0.72754, 0.494383, 0.325802, 0.20656, 0.129838, 0.0802393, 0.0489262, 0.0291108, 0.0178913, 0.0106616, 0.00642673, 0.00377871, 0.00232346, 0.00140085, 0.000857247, 0.000509539, 0.000302549, 0.000194001, 0.000110133, 6.37239e-05, 4.23769e-05, 2.58207e-05, 1.73478e-05, 9.19466e-06, 7.75149e-06
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "Promptv4" )
	{
		if( pu_set == "PU_S6" )
		{
			double weight[36] = {0, 0.793115, 1.67919, 2.24498, 2.67134, 2.76649, 2.62244, 2.27402, 1.82562, 1.36703, 0.963039, 0.63633, 0.398949, 0.235656, 0.13176, 0.0699228, 0.0358045, 0.0175461, 0.00831046, 0.00385472, 0.00172482, 0.000756113, 0.000325951, 0.000137284, 5.65972e-05, 2.31767e-05, 9.31827e-06, 3.69178e-06, 1.45294e-06, 5.64848e-07, 2.17563e-07, 8.31074e-08, 3.1095e-08, 1.166e-08, 4.28295e-09, 2.17473e-09
//			double weight[25] = {0, 0.000675729, 0.0326698, 0.187647, 1.62785, 4.26306, 5.40223, 3.489, 1.17513, 0.284859, 0.0481176, 0.00569589, 0.0002741, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0633001, 0.507791, 1.03085, 1.67462, 2.11114, 2.23863, 2.07412, 1.73176, 1.31701, 0.946101, 0.640917, 0.420471, 0.266944, 0.162645, 0.0979641, 0.057811, 0.0335185, 0.0188679, 0.0109059, 0.00607106, 0.0033933, 0.00183537, 0.0010297, 0.000561865, 0.000308742, 0.000163576, 8.59978e-05, 4.85363e-05, 2.41283e-05, 1.21725e-05, 7.03273e-06, 3.71224e-06, 2.15547e-06, 9.86134e-07, 6.7091e-07
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "July05" )
	{
		if( pu_set == "PU_S6" )
		{
			double weight[36] = {0, 0.779616, 1.6476, 2.20724, 2.63598, 2.74173, 2.61127, 2.27578, 1.83692, 1.38352, 0.980796, 0.6525, 0.412134, 0.245425, 0.138449, 0.0742001, 0.0384162, 0.0190625, 0.00915869, 0.00431893, 0.00197003, 0.000883198, 0.000390841, 0.000169714, 7.24868e-05, 3.09184e-05, 1.30237e-05, 5.43963e-06, 2.27166e-06, 9.43399e-07, 3.908e-07, 1.61636e-07, 6.59134e-08, 2.71116e-08, 1.09859e-08, 6.48997e-09
//			double weight[25] = {0, 0.00159758, 0.0400694, 0.201236, 1.51024, 4.21618, 5.31879, 3.5502, 1.26481, 0.312659, 0.0591438, 0.00992477, 0.00132608, 0.000245274, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0622227, 0.498239, 1.01352, 1.65246, 2.09225, 2.22909, 2.07573, 1.74249, 1.33289, 0.963547, 0.657204, 0.434368, 0.27801, 0.170902, 0.103957, 0.0620279, 0.0364154, 0.0207938, 0.0122193, 0.00693416, 0.00396364, 0.00220075, 0.00127295, 0.000719609, 0.00041187, 0.000228623, 0.000126713, 7.58861e-05, 4.02985e-05, 2.1865e-05, 1.3678e-05, 7.86899e-06, 5.01186e-06, 2.52947e-06, 2.00217e-06
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "Aug05" )
	{
		if( pu_set == "PU_S6" )
		{
			double weight[36] = {0, 0.621884, 1.23155, 1.60508, 1.93665, 2.10963, 2.17317, 2.10514, 1.93023, 1.6785, 1.38919, 1.08615, 0.808613, 0.567516, 0.376386, 0.236125, 0.142266, 0.0815734, 0.0449254, 0.0240704, 0.0123565, 0.00617209, 0.00301131, 0.00142604, 0.000656869, 0.000298747, 0.000132648, 5.77309e-05, 2.4836e-05, 1.05061e-05, 4.38461e-06, 1.80779e-06, 7.27497e-07, 2.92517e-07, 1.1483e-07, 6.38265e-08
//			double weight[25] = {0, 0.00272668, 0.0470127, 0.428861, 1.198, 2.57645, 2.63078, 3.09, 3.21441, 2.12765, 0.765809, 0.164888, 0.017926, 0.000221654, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0496338, 0.372423, 0.737016, 1.21406, 1.60989, 1.85511, 1.92009, 1.831, 1.61708, 1.36475, 1.09398, 0.852235, 0.642864, 0.464613, 0.330819, 0.229706, 0.155831, 0.101998, 0.0681009, 0.043493, 0.0276993, 0.0169561, 0.0106961, 0.00652103, 0.00397967, 0.00232855, 0.00134481, 0.000829663, 0.000448781, 0.000245316, 0.000152979, 8.68513e-05, 5.40748e-05, 2.64391e-05, 1.96906e-05
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "Oct03" )
	{
		if( pu_set == "PU_S6" )
		{
			double weight[36] = {0, 0.455005, 0.961839, 1.32142, 1.65395, 1.85799, 1.97655, 1.98843, 1.90714, 1.747, 1.53216, 1.27521, 1.0139, 0.761588, 0.541268, 0.364095, 0.235236, 0.144596, 0.0853228, 0.0489447, 0.0268788, 0.0143505, 0.00747748, 0.00377898, 0.00185653, 0.000900163, 0.000426006, 0.000197621, 9.06457e-05, 4.09073e-05, 1.82296e-05, 8.0353e-06, 3.46224e-06, 1.49329e-06, 6.30138e-07, 3.90769e-07
//			double weight[25] = {0, 0.000183741, 0.00193635, 0.080368, 1.15505, 2.10757, 2.29994, 2.47695, 2.77971, 2.75815, 1.7542, 0.605981, 0.106182, 0.0107263, 0.000959857, 9.10193e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0363149, 0.290862, 0.606765, 1.03684, 1.41786, 1.68726, 1.81364, 1.8091, 1.68307, 1.50521, 1.28441, 1.06859, 0.862703, 0.668144, 0.51011, 0.379819, 0.276225, 0.193716, 0.138476, 0.094609, 0.0644025, 0.0421042, 0.0283444, 0.0184306, 0.0119912, 0.00747826, 0.00460347, 0.00302808, 0.00174741, 0.00101993, 0.000679965, 0.000413336, 0.00027605, 0.000145087, 0.000120553
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "2011A" )
	{
		if( pu_set == "PU_S6" )
		{
			double weight[36] = {0, 0.653879, 1.36786, 1.83486, 2.21774, 2.36488, 2.34331, 2.15919, 1.87421, 1.54448, 1.2183, 0.91604, 0.66296, 0.457507, 0.301672, 0.190079, 0.116053, 0.0679454, 0.038449, 0.0212753, 0.0113258, 0.00588586, 0.00299559, 0.00148298, 0.000715391, 0.000341282, 0.00015918, 7.28784e-05, 3.30306e-05, 1.47437e-05, 6.50407e-06, 2.84003e-06, 1.21299e-06, 5.18855e-07, 2.17241e-07, 1.33291e-07
//			double weight[25] = {0, 0.00135944, 0.0296401, 0.203218, 1.34886, 3.29462, 3.94141, 3.14494, 2.05837, 1.36667, 0.695758, 0.217848, 0.0360863, 0.00343302, 0.000292303, 2.77179e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0521874, 0.413644, 0.842525, 1.39027, 1.80467, 2.00035, 1.96939, 1.77786, 1.48797, 1.19688, 0.922644, 0.698725, 0.518249, 0.372385, 0.266307, 0.187383, 0.129797, 0.0872939, 0.0601927, 0.0398649, 0.0264147, 0.0168676, 0.0111232, 0.007102, 0.00454628, 0.0027943, 0.00169766, 0.00110341, 0.000629799, 0.000363898, 0.00024033, 0.000144811, 9.59158e-05, 5.0019e-05, 4.11205e-05
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "2011B" )
	{
		if( pu_set == "PU_S6" )
		{
			double weight[36] = {0, 0.0210887, 0.0684082, 0.142701, 0.265012, 0.428769, 0.637971, 0.875684, 1.12559, 1.3677, 1.58606, 1.74977, 1.8559, 1.87615, 1.81248, 1.67434, 1.50053, 1.29156, 1.07654, 0.879335, 0.692634, 0.533926, 0.404124, 0.298322, 0.215172, 0.153903, 0.107925, 0.0745003, 0.0510545, 0.0345561, 0.0231817, 0.0154372, 0.0100837, 0.00661527, 0.00425963, 0.00538563
//			double weight[25] = {0, 0.000574776, 0.00134631, 0.000388427, 0.00292684, 0.0195094, 0.133842, 0.669048, 1.21, 1.54708, 1.93797, 2.24333, 2.45257, 2.52476, 2.35355, 1.8487, 1.16996, 0.590205, 0.243815, 0.0920991, 0.0329105, 0.0108703, 0.00286849, 0.000388309, 1.21386e-07
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.00168313, 0.0206868, 0.065525, 0.166132, 0.327199, 0.5446, 0.798706, 1.06772, 1.31765, 1.55816, 1.76238, 1.95602, 2.12525, 2.23733, 2.3458, 2.4228, 2.46728, 2.44417, 2.48785, 2.43795, 2.39617, 2.27554, 2.23758, 2.13611, 2.05017, 1.89456, 1.73544, 1.70551, 1.47611, 1.297, 1.30633, 1.20383, 1.2229, 0.980765, 1.66148
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "2011" || lumi_set == "2011ff" || lumi_set == "2011f"	|| lumi_set == "2011_rereco")
	{
		if( pu_set == "PU_S6" )
		{
//			double weight[25] = {0, 0.000886659, 0.0125924, 0.0810079, 0.537903, 1.32128, 1.64725, 1.65315, 1.54721, 1.47538, 1.44422, 1.43825, 1.49208, 1.5226, 1.41819, 1.11389, 0.70493, 0.355613, 0.146905, 0.055492, 0.0198294, 0.00654963, 0.00172834, 0.000233966, 7.31382e-08
//			double weight[36] = {0, 0.315323, 0.672626, 0.929519, 1.17299, 1.32902, 1.43092, 1.47249, 1.47368, 1.4499, 1.41506, 1.3621, 1.30121, 1.21651, 1.10998, 0.984188, 0.856777, 0.722602, 0.59385, 0.480355, 0.375839, 0.288398, 0.217607, 0.160298, 0.115454, 0.0824999, 0.0578162, 0.039893, 0.0273306, 0.018495, 0.0124057, 0.00826054, 0.00539555, 0.00353955, 0.00227909, 0.00288148};
//			double weight[36] = {0, 0, 0.00199742, 0.000274429, 0.000639867, 0.00723113, 0.0362718, 0.260005, 0.789205, 1.06541, 1.39698, 1.87352, 2.60395, 2.8074, 2.95221, 2.81277, 2.29894, 1.76985, 0.674663, 0.309738, 0.111751, 0.0383114, 0.0114226, 0.00312103, 0.000120836, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
			double weight[50] = {0, 0.526829, 0.0199659, 0.041029, 0.297575, 0.945121, 1.25478, 1.37131, 1.33925, 1.21619, 1.31131, 1.3953, 1.58576, 1.44483, 1.39461, 1.35017, 1.31175, 1.39853, 0.846608, 0.656049, 0.373876, 0.188962, 0.0829103, 0.0353935, 0.0118849, 0.00415899, 0.00152407, 0.000559973, 0.000338563, 8.18194e-05, 2.39062e-05, 8.35835e-06, 6.12925e-06, 6.12198e-06, 3.45724e-05, 1.58819e-05, 6.8505e-05, 4.91191e-05, 0.00286974, 0, 0.000113381, 0, 0, 0, 0, 0, 0, 0, 0, 0};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0251666, 0.203404, 0.426814, 0.735332, 1.01419, 1.22149, 1.34305, 1.39792, 1.39684, 1.39017, 1.37192, 1.37141, 1.37803, 1.37017, 1.37888, 1.38338, 1.3804, 1.34827, 1.35904, 1.32289, 1.29428, 1.2253, 1.20232, 1.14617, 1.099, 1.01493, 0.929284, 0.912994, 0.790041, 0.69409, 0.699025, 0.644141, 0.654322, 0.524752, 0.888945};
			return weight[nGenVertices];
		}
	}
//			double weight[36] = {0.0};
			double weight[36] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
//			double weight[25] = {0.0};
	return weight[nGenVertices];
}




double weight_TTJets(int nGenVertices)
{
	double weight[51] = {
0, 0, 0.130953, 0.478315, 0.814199, 1.11549, 1.36463, 1.4534, 1.47222, 1.38149, 1.30462, 1.19369, 1.08147, 1.00757, 0.93383, 0.83356, 0.748923, 0.691019, 0.699203, 0.586132, 0.499935, 0.442856, 0.377828, 0.310338, 0.176515, 0.201486, 0.292975, 0.0715264, 0.0572582, 0.0292919, 0.061007, 0.0248481, 0, 0, 0.00147735, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0// May10_Promptv4_Aug05_Promptv6 2.15 fb-1
};
	return weight[nGenVertices];
}
double weight_WJetsToLNu(int nGenVertices)
{
	double weight[51] = {
0, 0, 0.207242, 0.445507, 0.859512, 1.07573, 1.26513, 1.46919, 1.21893, 1.2368, 0.996704, 0.898286, 1.02384, 1.10968, 1.15656, 1.41378, 1.19411, 1.00502, 0.973797, 1.64394, 0, 0.947574, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0// May10_Promptv4_Aug05_Promptv6 2.15 fb-1
};
	return weight[nGenVertices];
}
double	weight_QCDMu(int nGenVertices)
{
	double weight[51] = {
0, 0, 0.216306, 0.484129, 0.782602, 1.0295, 1.30368, 1.31239, 1.26043, 1.10026, 1.07655, 1.03941, 1.0667, 1.00138, 1.15552, 0.968367, 1.46147, 1.95245, 3.40521, 1.9162, 2.08934, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 // May10_Promptv4_Aug05_Promptv6 2.15 fb-1
};
	return weight[nGenVertices];
}

// *****************************************************************************************************
// ******************* factorial for combinatorics
// *****************************************************************************************************
int factorial(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * factorial(number - 1);
	return temp;
}

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

// *****************************************************************************************************
// ******************* Compute DeltaR between two four-momentum vectors
// *****************************************************************************************************
double DeltaR( double eta1, double phi1, double eta2, double phi2)
{
	double DeltaEta = fabs( eta1-eta2 );
	double DeltaPhi = fabs( phi1-phi2 );
	// Returning DeltaPhi in the correct range (0, 2pi)
	while (DeltaPhi >	 TMath::Pi()) DeltaPhi -= 2*TMath::Pi();
	while (DeltaPhi <= -TMath::Pi()) DeltaPhi += 2*TMath::Pi();
	return sqrt(DeltaEta*DeltaEta + DeltaPhi*DeltaPhi);
}

// *****************************************************************************************************
// ******************* Process gen info : get out four momentum of gen particle
// *****************************************************************************************************
void doGenInfo(TRootParticle* myparticle, TClonesArray* mcParticles, float* particule_trueE, float* particule_truePx, float* particule_truePy, float* particule_truePz, float* particule_truePhi, float* particule_trueEta, float* particule_truePt, int particle_pdgId = 0)
{
	TRootMCParticle* mygenparticle;
	int NbMCpartInCone=0;
	double bestPtdiff=500.0;
	int igpsl=-1;
	for (int igp=0; igp<mcParticles->GetEntriesFast(); igp++) {
		mygenparticle = (TRootMCParticle*) mcParticles->At(igp);
		//cout<<endl<<"deltaR = "<<DeltaR(mygenparticle->Eta(), mygenparticle->Phi(), myparticle->Eta(), myparticle->Phi())<<endl;
		//cout<<"mygenparticle->Mag() = "<<mygenparticle->Mag()<<endl;
		if (DeltaR(mygenparticle->Eta(), mygenparticle->Phi(), myparticle->Eta(), myparticle->Phi())<0.3){
			if ( (mygenparticle->status()==1) && ( (particle_pdgId==0)?true:((mygenparticle->type())==particle_pdgId) ) ){
				NbMCpartInCone++;
				if (fabs( (mygenparticle->Pt()) - (myparticle->Pt()) )<bestPtdiff){
					bestPtdiff=fabs(mygenparticle->Pt()-myparticle->Pt());
					igpsl=igp;
				}
			}
		}
	}
	if (igpsl!=-1){
		mygenparticle = (TRootMCParticle*) mcParticles->At(igpsl);
		*particule_trueE = mygenparticle->Energy();
		*particule_truePx = mygenparticle->Px();
		*particule_truePy = mygenparticle->Py();
		*particule_truePz = mygenparticle->Pz();
		*particule_truePhi = mygenparticle->Phi();
		*particule_trueEta = mygenparticle->Eta();
		*particule_truePt = mygenparticle->Pt();

	}
	return;
}

// *****************************************************************************************************
// ******************* Declaration of main function
// *****************************************************************************************************
int main(int argc, char *argv[]);


	// ____________________________________________
	// Event information
	// ____________________________________________
	extern ULong64_t iEvent, iEventID, iLumiID, iRunID;
	extern Int_t isMM;
	extern Int_t nVertices;
	extern Int_t nGenVertices;
	extern Float_t weight_pileUp, weight_Xsection;

	// ____________________________________________
	// Muon variables 
	// ____________________________________________
	extern Int_t NbMuons;
	
	extern Float_t Pt_allMuons, Eta_allMuons, Phi_allMuons, Charge_allMuons;
// (M minus charge, P plus charge), (L leading, S subleading)
	extern Float_t MuonM_Pt, MuonP_Pt, MuonL_Pt, MuonS_Pt;
	extern Float_t MuonM_Eta, MuonP_Eta, MuonL_Eta, MuonS_Eta;
	extern Float_t MuonM_Phi, MuonP_Phi, MuonL_Phi, MuonS_Phi;
	extern Int_t MuonL_Charge, MuonS_Charge;

	extern Float_t MuonM_isoR03_emEt, MuonP_isoR03_emEt, MuonL_isoR03_emEt, MuonS_isoR03_emEt;
	extern Float_t MuonM_isoR03_hadEt, MuonP_isoR03_hadEt, MuonL_isoR03_hadEt, MuonS_isoR03_hadEt;
	extern Float_t MuonM_isoR03_hoEt, MuonP_isoR03_hoEt, MuonL_isoR03_hoEt, MuonS_isoR03_hoEt;
	extern Float_t MuonM_isoR03_nJets, MuonP_isoR03_nJets, MuonL_isoR03_nJets, MuonS_isoR03_nJets;
	extern Float_t MuonM_isoR03_nTracks, MuonP_isoR03_nTracks, MuonL_isoR03_nTracks, MuonS_isoR03_nTracks;
	extern Float_t MuonM_isoR03_sumPt, MuonP_isoR03_sumPt, MuonL_isoR03_sumPt, MuonS_isoR03_sumPt;

	extern Float_t MuonM_isoR05_emEt, MuonP_isoR05_emEt, MuonL_isoR05_emEt, MuonS_isoR05_emEt;
	extern Float_t MuonM_isoR05_hadEt, MuonP_isoR05_hadEt, MuonL_isoR05_hadEt, MuonS_isoR05_hadEt;
	extern Float_t MuonM_isoR05_hoEt, MuonP_isoR05_hoEt, MuonL_isoR05_hoEt, MuonS_isoR05_hoEt;
	extern Float_t MuonM_isoR05_nJets, MuonP_isoR05_nJets, MuonL_isoR05_nJets, MuonS_isoR05_nJets;
	extern Float_t MuonM_isoR05_nTracks, MuonP_isoR05_nTracks, MuonL_isoR05_nTracks, MuonS_isoR05_nTracks;
	extern Float_t MuonM_isoR05_sumPt, MuonP_isoR05_sumPt, MuonL_isoR05_sumPt, MuonS_isoR05_sumPt;

	extern Float_t MuonM_E, MuonP_E, MuonL_E, MuonS_E;
	extern Float_t MuonM_Px, MuonP_Px, MuonL_Px, MuonS_Px;
	extern Float_t MuonM_Py, MuonP_Py, MuonL_Py, MuonS_Py;
	extern Float_t MuonM_Pz, MuonP_Pz, MuonL_Pz, MuonS_Pz;

	// ____________________________________________
	// mumu information
	// ____________________________________________

	extern Float_t Mmumu;
	extern Float_t Ptmumu;
	// ____________________________________________
	// MC Truth
	// ___________________________________________

	extern Float_t MuonM_MC_E, MuonM_MC_Px, MuonM_MC_Py, MuonM_MC_Pz, MuonM_MC_Phi, MuonM_MC_Eta, MuonM_MC_Pt;
	extern Float_t MuonP_MC_E, MuonP_MC_Px, MuonP_MC_Py, MuonP_MC_Pz, MuonP_MC_Phi, MuonP_MC_Eta, MuonP_MC_Pt;
	extern Float_t MuonL_MC_E, MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_Phi, MuonL_MC_Eta, MuonL_MC_Pt;
	extern Float_t MuonS_MC_E, MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_Phi, MuonS_MC_Eta, MuonS_MC_Pt;
	extern Float_t Mmumu_Muons_MC;

// *****************************************************************************************************
// ******************* SIDRA muon corrections
// *****************************************************************************************************
double applySidra( double _pt, double charge, double eta, double phi, TRandom3* generator)
{ 
	double pt = _pt;
	double sgn_charge = charge >= 0.0 ? 1.0 : -1.0;
// Correct MC
	double a = 0.0650687e-3;
	double b = 0.212987e-3;
	double c = 1.53414;
	pt = (double)(1.0)/(double)(pt);
	pt -= a - b * sgn_charge * sin( phi + c );
// Apply Corrections
	double A = 0.143812;
	double B = 0.0404834;
	double Ap = 0.0995898;
	double Bp = -0.0692569;
	double Cp = 0.0952174;
	double phi0 = -1.08881;
	pt += (double)(( A + B * eta * eta ) * (generator->Gaus(0,1))) / (double)(1000.) + (double)(Ap + Cp * sgn_charge * sin( phi + phi0 )+ Bp * sgn_charge * eta )/(double)(1000.);
	pt = (double)(1.0)/(double)(pt);
	return pt;
} 

// *****************************************************************************************************
// ******************* MuScleFit muon correction
// *****************************************************************************************************
double applyMuScleFit(double _pt, double charge, double eta, double phi)
{
	double b = -5.03313e-6;
	double c = -4.41463e-5;
	double d0 = -0.000148871;
	double e0 = 1.59501;
	double d1 = 7.95495e-5;
	double e1 = -0.364823;
	double d2 = 0.000152032;
	double e2 = 0.410195;
	double d = eta > .9 ? d1 : (eta < -.9 ? d2 : d0);
	double e = eta > .9 ? e1 : (eta < -.9 ? e2 : e0);
	double pt = _pt;
	double sgn_eta = eta >= 0.0 ? 1.0 : -1.0;
	double sgn_charge = charge >= 0.0 ? 1.0 : -1.0;
	pt = _pt * (1.0 + b * _pt	+ c * sgn_charge * _pt * sgn_eta * eta * eta + sgn_charge * d * _pt * sin( phi + e ));
	return pt;

}




// *****************************************************************************************************
// ******************* Fill minitree
// *****************************************************************************************************
//int FillMM(TRootMuon* mymuon1, TRootMuon* mymuon2, bool doMC, bool applyMuonScaleCorrection, TClonesArray* mcParticles, TRandom3* generator){
int FillMM(TRootMuon* mymuon1, TRootMuon* mymuon2, bool doMC, bool applyMuonScaleCorrection, TClonesArray* mcParticles, double corrected_Pt1_, double corrected_Pt2_){

		double corrected_Pt1 = mymuon1->Pt();
		double corrected_Pt2 = mymuon2->Pt();
		if( applyMuonScaleCorrection > 0 )
		{
/*
			// Sidra makes MC look like data
			if( doMC > 0) corrected_Pt1 = applySidra(mymuon1->Pt(), mymuon1->charge(), mymuon1->Eta(), mymuon1->Phi(), generator);
			if( doMC > 0) corrected_Pt2 = applySidra(mymuon2->Pt(), mymuon2->charge(), mymuon2->Eta(), mymuon2->Phi(), generator);
			// MuScleFit correct data absolute scale
			corrected_Pt1 = applyMuScleFit(corrected_Pt1, mymuon1->charge(), mymuon1->Eta(), mymuon1->Phi());
			corrected_Pt2 = applyMuScleFit(corrected_Pt2, mymuon2->charge(), mymuon2->Eta(), mymuon2->Phi());
*/
			corrected_Pt1 = corrected_Pt1_;
			corrected_Pt2 = corrected_Pt2_;
		}
		double corrected_Pz1 = mymuon1->Pz();
		double corrected_Pz2 = mymuon2->Pz();
		double corrected_Px1 = applyMuonScaleCorrection > 0 ? mymuon1->Px() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Px();
		double corrected_Px2 = applyMuonScaleCorrection > 0 ? mymuon2->Px() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Px();
		double corrected_Py1 = applyMuonScaleCorrection > 0 ? mymuon1->Py() * (double)(corrected_Pt1) / (double)(mymuon1->Pt()) : mymuon1->Py();
		double corrected_Py2 = applyMuonScaleCorrection > 0 ? mymuon2->Py() * (double)(corrected_Pt2) / (double)(mymuon2->Pt()) : mymuon2->Py();
		double m_mu = 105.658367e-3;
//		double corrected_E1 = mymuon1->E();
//		double corrected_E2 = mymuon2->E();
		double corrected_E1 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz1 * corrected_Pz1 + corrected_Pt1 * corrected_Pt1)) : mymuon1->E();
		double corrected_E2 = applyMuonScaleCorrection > 0 ? sqrt( m_mu * m_mu + (corrected_Pz2 * corrected_Pz2 + corrected_Pt2 * corrected_Pt2)) : mymuon2->E();
		TLorentzVector *correctedMuon1 = new TLorentzVector(corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1);
		TLorentzVector *correctedMuon2 = new TLorentzVector(corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2);

		TRootMuon *mycorrectedMuon1 = new TRootMuon( corrected_Px1, corrected_Py1, corrected_Pz1, corrected_E1, mymuon1->vx(), mymuon1->vy(), mymuon1->vz(), 1, mymuon1->charge() );
		TRootMuon *mycorrectedMuon2 = new TRootMuon( corrected_Px2, corrected_Py2, corrected_Pz2, corrected_E2, mymuon2->vx(), mymuon2->vy(), mymuon2->vz(), 1, mymuon2->charge() );



			// Fill muons stuff
			TRootMuon *leadingMuon;
			TRootMuon *subleadingMuon;
			if( (mycorrectedMuon1->Pt()) > (mycorrectedMuon2->Pt()) )			{
				leadingMuon = mycorrectedMuon1;
				subleadingMuon = mycorrectedMuon2;
			} else {
				leadingMuon = mycorrectedMuon2;
				subleadingMuon = mycorrectedMuon1;
			}
			MuonL_Pt = leadingMuon->Pt();
			MuonL_Eta = leadingMuon->Eta();
			MuonL_Phi = leadingMuon->Phi();
			MuonL_Charge = leadingMuon->charge();
			MuonL_isoR03_emEt = leadingMuon->isoR03_emEt();
			MuonL_isoR03_hadEt = leadingMuon->isoR03_hadEt();
			MuonL_isoR03_hoEt = leadingMuon->isoR03_hoEt();
			MuonL_isoR03_nJets = leadingMuon->isoR03_nJets();
			MuonL_isoR03_nTracks = leadingMuon->isoR03_nTracks();
			MuonL_isoR03_sumPt = leadingMuon->isoR03_sumPt();
			MuonL_isoR05_emEt = leadingMuon->isoR05_emEt();
			MuonL_isoR05_hadEt = leadingMuon->isoR05_hadEt();
			MuonL_isoR05_hoEt = leadingMuon->isoR05_hoEt();
			MuonL_isoR05_nJets = leadingMuon->isoR05_nJets();
			MuonL_isoR05_nTracks = leadingMuon->isoR05_nTracks();
			MuonL_isoR05_sumPt = leadingMuon->isoR05_sumPt();
			MuonL_E = leadingMuon->Energy();
			MuonL_Px = leadingMuon->Px();
			MuonL_Py = leadingMuon->Py();
			MuonL_Pz = leadingMuon->Pz();

			MuonS_Pt = subleadingMuon->Pt();
			MuonS_Eta = subleadingMuon->Eta();
			MuonS_Phi = subleadingMuon->Phi();
			MuonS_Charge = subleadingMuon->charge();
			MuonS_isoR03_emEt = subleadingMuon->isoR03_emEt();
			MuonS_isoR03_hadEt = subleadingMuon->isoR03_hadEt();
			MuonS_isoR03_hoEt = subleadingMuon->isoR03_hoEt();
			MuonS_isoR03_nJets = subleadingMuon->isoR03_nJets();
			MuonS_isoR03_nTracks = subleadingMuon->isoR03_nTracks();
			MuonS_isoR03_sumPt = subleadingMuon->isoR03_sumPt();
			MuonS_isoR05_emEt = subleadingMuon->isoR05_emEt();
			MuonS_isoR05_hadEt = subleadingMuon->isoR05_hadEt();
			MuonS_isoR05_hoEt = subleadingMuon->isoR05_hoEt();
			MuonS_isoR05_nJets = subleadingMuon->isoR05_nJets();
			MuonS_isoR05_nTracks = subleadingMuon->isoR05_nTracks();
			MuonS_isoR05_sumPt = subleadingMuon->isoR05_sumPt();
			MuonS_E = subleadingMuon->Energy();
			MuonS_Px = subleadingMuon->Px();
			MuonS_Py = subleadingMuon->Py();
			MuonS_Pz = subleadingMuon->Pz();


			TLorentzVector mumu;
			mumu = (*leadingMuon) + (*subleadingMuon);
//			cout << "mumu.M()= " << mumu.M() << endl;
//			TLorentzVector mumu2;
//			mumu2 = (*correctedMuon1) + (*correctedMuon2);
//			cout << "mumu2.M()= " << mumu2.M() << endl;

			Ptmumu = mumu.Pt();
			double mumuInvMass = mumu.M();
//			cerr << "\t\tINFO: Dimuon invariant mass : Mmumu = " << mumuInvMass << endl;
			mumu.Clear();
			Mmumu = mumuInvMass;

		double phiMuon = mycorrectedMuon1->Phi();
		double etaMuon = mycorrectedMuon1->Eta();
		double phiMuon_oppositeCharge = mycorrectedMuon2->Phi();
		double etaMuon_oppositeCharge = mycorrectedMuon2->Eta();

		TRootMuon *minusMuon;
		TRootMuon *plusMuon;
		if( mycorrectedMuon1->charge()>0 ){
			plusMuon	= (TRootMuon*) mycorrectedMuon1;
			minusMuon = (TRootMuon*) mycorrectedMuon2;
		} else {
			minusMuon = (TRootMuon*) mycorrectedMuon1;
			plusMuon	= (TRootMuon*) mycorrectedMuon2;
		}


		// FILLING MINITREE INFORMATION
		MuonM_Pt = minusMuon->Pt();
		MuonM_Eta = minusMuon->Eta();
		MuonM_Phi = minusMuon->Phi();
		MuonM_E = minusMuon->Energy();
		MuonM_Px = minusMuon->Px();
		MuonM_Py = minusMuon->Py();
		MuonM_Pz = minusMuon->Pz();
		MuonP_Pt = plusMuon->Pt();
		MuonP_Eta = plusMuon->Eta();
		MuonP_Phi = plusMuon->Phi();
		MuonP_E = plusMuon->Energy();
		MuonP_Px = plusMuon->Px();
		MuonP_Py = plusMuon->Py();
		MuonP_Pz = plusMuon->Pz();

		if( doMC )
		{ 
	// Compute Stuff, with MC truth information
			doGenInfo( (TRootParticle*) minusMuon, mcParticles, &MuonM_MC_E, &MuonM_MC_Px, &MuonM_MC_Py, &MuonM_MC_Pz, &MuonM_MC_Phi, &MuonM_MC_Eta, &MuonM_MC_Pt, 13 );
			doGenInfo( (TRootParticle*) plusMuon, mcParticles, &MuonP_MC_E, &MuonP_MC_Px, &MuonP_MC_Py, &MuonP_MC_Pz, &MuonP_MC_Phi, &MuonP_MC_Eta, &MuonP_MC_Pt, -13 );
			doGenInfo( (TRootParticle*) leadingMuon, mcParticles, &MuonL_MC_E, &MuonL_MC_Px, &MuonL_MC_Py, &MuonL_MC_Pz, &MuonL_MC_Phi, &MuonL_MC_Eta, &MuonL_MC_Pt, (-1)*(leadingMuon->charge())*13 );
			doGenInfo( (TRootParticle*) subleadingMuon, mcParticles, &MuonS_MC_E, &MuonS_MC_Px, &MuonS_MC_Py, &MuonS_MC_Pz, &MuonS_MC_Phi, &MuonS_MC_Eta, &MuonS_MC_Pt, (-1)*(subleadingMuon->charge())*13 );

			TLorentzVector mumu_Muons_MC;

			TLorentzVector *MuonLMC = new TLorentzVector( MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_E);
			TLorentzVector *MuonSMC = new TLorentzVector( MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_E);

			mumu_Muons_MC = (*MuonLMC) + (*MuonSMC);
			Mmumu_Muons_MC = mumu_Muons_MC.M();
		MuonLMC->Clear();
		MuonSMC->Clear();
		}// end doMC

return 0;
}

