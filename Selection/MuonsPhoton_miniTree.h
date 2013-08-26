// TMVA headers
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

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
#include <string>
#include <stdio.h>
#include <stdlib.h>

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

// personal headers
#include "corrections.h"


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
	} else if( lumi_set == "2011A" || lumi_set == "2011A_rereco")
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
			//double weight[36] = {0, 0.0210887, 0.0684082, 0.142701, 0.265012, 0.428769, 0.637971, 0.875684, 1.12559, 1.3677, 1.58606, 1.74977, 1.8559, 1.87615, 1.81248, 1.67434, 1.50053, 1.29156, 1.07654, 0.879335, 0.692634, 0.533926, 0.404124, 0.298322, 0.215172, 0.153903, 0.107925, 0.0745003, 0.0510545, 0.0345561, 0.0231817, 0.0154372, 0.0100837, 0.00661527, 0.00425963, 0.00538563
//			double weight[25] = {0, 0.000574776, 0.00134631, 0.000388427, 0.00292684, 0.0195094, 0.133842, 0.669048, 1.21, 1.54708, 1.93797, 2.24333, 2.45257, 2.52476, 2.35355, 1.8487, 1.16996, 0.590205, 0.243815, 0.0920991, 0.0329105, 0.0108703, 0.00286849, 0.000388309, 1.21386e-07
			double weight[36] = {0, 0, 0.00199742, 0.000274429, 0.000639867, 0.00723113, 0.0362718, 0.260005, 0.789205, 1.06541, 1.39698, 1.87352, 2.60395, 2.8074, 2.95221, 2.81277, 2.29894, 1.76985, 0.674663, 0.309738, 0.111751, 0.0383114, 0.0114226, 0.00312103, 0.000120836, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
    {
      double weight[36] = {0, 0.00168313, 0.0206868, 0.065525, 0.166132, 0.327199, 0.5446, 0.798706, 1.06772, 1.31765, 1.55816, 1.76238, 1.95602, 2.12525, 2.23733, 2.3458, 2.4228, 2.46728, 2.44417, 2.48785, 2.43795, 2.39617, 2.27554, 2.23758, 2.13611, 2.05017, 1.89456, 1.73544, 1.70551, 1.47611, 1.297, 1.30633, 1.20383, 1.2229, 0.980765, 1.66148
			};
			return weight[nGenVertices];
		}
	} else if( lumi_set == "2011" || lumi_set == "2011ff" || lumi_set == "2011f"  || lumi_set == "2011_rereco" )
	{
		if( pu_set == "PU_S6" )
		{
//			double weight[36] = {0, 0.315323, 0.672626, 0.929519, 1.17299, 1.32902, 1.43092, 1.47249, 1.47368, 1.4499, 1.41506, 1.3621, 1.30121, 1.21651, 1.10998, 0.984188, 0.856777, 0.722602, 0.59385, 0.480355, 0.375839, 0.288398, 0.217607, 0.160298, 0.115454, 0.0824999, 0.0578162, 0.039893, 0.0273306, 0.018495, 0.0124057, 0.00826054, 0.00539555, 0.00353955, 0.00227909, 0.00288148
			//double weight[36] = {0, 0.315325, 0.672629, 0.929519, 1.17299, 1.32902, 1.43091, 1.47248, 1.47368, 1.4499, 1.41506, 1.3621, 1.30121, 1.21651, 1.10999, 0.984189, 0.856772, 0.722602, 0.593852, 0.480355, 0.375841, 0.288399, 0.217607, 0.160298, 0.115453, 0.0825002, 0.0578159, 0.0398932, 0.0273307, 0.0184951, 0.0124057, 0.00826039, 0.00539557, 0.00353944, 0.0022791, 0.00288149
//			double weight[25] = {0, 0.000886659, 0.0125924, 0.0810079, 0.537903, 1.32128, 1.64725, 1.65315, 1.54721, 1.47538, 1.44422, 1.43825, 1.49208, 1.5226, 1.41819, 1.11389, 0.70493, 0.355613, 0.146905, 0.055492, 0.0198294, 0.00654963, 0.00172834, 0.000233966, 7.31382e-08
//			double weight[36] = {0, 0, 0.00199742, 0.000274429, 0.000639867, 0.00723113, 0.0362718, 0.260005, 0.789205, 1.06541, 1.39698, 1.87352, 2.60395, 2.8074, 2.95221, 2.81277, 2.29894, 1.76985, 0.674663, 0.309738, 0.111751, 0.0383114, 0.0114226, 0.00312103, 0.000120836, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//			double weight[36] = {0, 0.278681, 0.0195797, 0.0557011, 0.518554, 1.33704, 1.47464, 1.45241, 1.36472, 1.2217, 1.33893, 1.40766, 1.60477, 1.47503, 1.38344, 1.23752, 1.04787, 0.923383, 0.442392, 0.261088, 0.109806, 0.040314, 0.012964, 0.00420289, 0.00113162, 0.000336742, 0.000109707, 3.63179e-05, 1.96727e-05, 4.52912e-06, 1.83561e-06, 1.87659e-06, 5.32118e-06, 1.26775e-05, 9.1875e-05, 4.19224e-05
			double weight[50] = {0, 0.526829, 0.0199659, 0.041029, 0.297575, 0.945121, 1.25478, 1.37131, 1.33925, 1.21619, 1.31131, 1.3953, 1.58576, 1.44483, 1.39461, 1.35017, 1.31175, 1.39853, 0.846608, 0.656049, 0.373876, 0.188962, 0.0829103, 0.0353935, 0.0118849, 0.00415899, 0.00152407, 0.000559973, 0.000338563, 8.18194e-05, 2.39062e-05, 8.35835e-06, 6.12925e-06, 6.12198e-06, 3.45724e-05, 1.58819e-05, 6.8505e-05, 4.91191e-05, 0.00286974, 0, 0.000113381, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
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
0, 0.309682, 0.665237, 0.943737, 1.17551, 1.36737, 1.44975, 1.49434, 1.45382, 1.491, 1.42228, 1.38831, 1.30475, 1.21827, 1.12734, 0.992364, 0.866169, 0.714571, 0.58875, 0.467808, 0.370044, 0.275162, 0.205839, 0.154965, 0.108118, 0.0817004, 0.0547621, 0.039171, 0.0267462, 0.0189005, 0.012636, 0.00799086, 0.00525691, 0.00373897, 0.00217098, 0.00276357// 2011ff
//0, 0, 0.130953, 0.478315, 0.814199, 1.11549, 1.36463, 1.4534, 1.47222, 1.38149, 1.30462, 1.19369, 1.08147, 1.00757, 0.93383, 0.83356, 0.748923, 0.691019, 0.699203, 0.586132, 0.499935, 0.442856, 0.377828, 0.310338, 0.176515, 0.201486, 0.292975, 0.0715264, 0.0572582, 0.0292919, 0.061007, 0.0248481, 0, 0, 0.00147735, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0// May10_Promptv4_Aug05_Promptv6 2.15 fb-1
};
	return weight[nGenVertices];
}
double weight_WJetsToLNu(int nGenVertices)
{
	double weight[51] = {
0, 0.327797, 0.711539, 1.02162, 1.29027, 1.4698, 1.61995, 1.50805, 1.57083, 1.42989, 1.46513, 1.27826, 1.32338, 1.19892, 1.0591, 0.971407, 0.867411, 0.748671, 0.536506, 0.466712, 0.360607, 0.271065, 0.200123, 0.139836, 0.110706, 0.0784664, 0.049712, 0.0348633, 0.0259981, 0.0157856, 0.0102124, 0.00841392, 0.00394399, 0.00302661, 0.00198323, 0.00247728// 2011ff
//0, 0, 0.207242, 0.445507, 0.859512, 1.07573, 1.26513, 1.46919, 1.21893, 1.2368, 0.996704, 0.898286, 1.02384, 1.10968, 1.15656, 1.41378, 1.19411, 1.00502, 0.973797, 1.64394, 0, 0.947574, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0// May10_Promptv4_Aug05_Promptv6 2.15 fb-1
};
  return weight[nGenVertices];
}
double  weight_QCDMu(int nGenVertices)
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
// ******************* electromagnetic calorimeter cracks
// *****************************************************************************************************
int photonIsInCrack(double sc_abs_eta, double sc_abs_phi)
{
	double phiCrackSize = (double)(21.5) / (double) (1290.0);
	double phiCrackPosition = (double)(TMath::Pi()) / (double)(9.0);
	double phiOffset = -(double)(10.0 * TMath::Pi()) / (double)(180.0);
	bool isInEtaCrack = (sc_abs_eta < 0.018 ||
		(sc_abs_eta >0.423 && sc_abs_eta <0.461) ||
		(sc_abs_eta >0.770 && sc_abs_eta <0.806) ||
		(sc_abs_eta >1.127 && sc_abs_eta <1.163) ||
		(sc_abs_eta >1.460 && sc_abs_eta <1.558));
	bool isInPhiCrack = (	(sc_abs_phi + phiOffset) < phiCrackSize ||
		( (1.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (1.0*phiCrackPosition + phiCrackSize) ) ||
		( (2.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (2.0*phiCrackPosition + phiCrackSize) ) ||
		( (3.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (3.0*phiCrackPosition + phiCrackSize) ) ||
		( (4.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (4.0*phiCrackPosition + phiCrackSize) ) ||
		( (5.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (5.0*phiCrackPosition + phiCrackSize) ) ||
		( (6.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (6.0*phiCrackPosition + phiCrackSize) ) ||
		( (7.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (7.0*phiCrackPosition + phiCrackSize) ) ||
		( (8.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) && (sc_abs_phi + phiOffset) < (8.0*phiCrackPosition + phiCrackSize) ) ||
		( (9.0*phiCrackPosition - phiCrackSize) < (sc_abs_phi + phiOffset) ) );
//	isInPhiCrack = false;
	if (isInEtaCrack) return 1;
	if (isInPhiCrack) return 2;
	return 0;
}

// *****************************************************************************************************
// ******************* ETHZ photon corrections
// *****************************************************************************************************
//float ETHZ_getValue(bool isEB, float SC_rawE, float SC_Eta, float ES_E, float phiWidth, float etaWidth, float f_eta)
float ETHZ_getValue(TRootPhoton *myphoton, float f_eta)
{
	int mode = 1; // Photons !
  float corr = 1.;
  float corr2 = 1.;
  float energy = 0;

  if (myphoton->isEBPho()){
//    float cetacorr = fEta(myphoton->superCluster()->rawEnergy(), myphoton->superCluster()->Eta(), 0)/myphoton->superCluster()->rawEnergy();
		float cetacorr = f_eta;

    energy = myphoton->superCluster()->rawEnergy()*cetacorr; //previously in CMSSW
  }
  else if (myphoton->isEEPho()){
    energy = myphoton->superCluster()->rawEnergy()+myphoton->preshowerEnergy();
  }

  float newEnergy = energy;

  if (mode==0){ //Electron

    corr = ETHZ_fBremEta(myphoton->superCluster()->phiWidth()/myphoton->superCluster()->etaWidth(), myphoton->superCluster()->Eta(), 0);

    float et = energy*TMath::Sin(2*TMath::ATan(TMath::Exp(-myphoton->superCluster()->Eta())))/corr;

    if (myphoton->isEBPho()) corr2 = corr * ETHZ_fEt(et, 0);
    if (myphoton->isEEPho()) corr2 = corr * ETHZ_fEnergy(energy/corr, 1);

    newEnergy = energy/corr2; 

  }

  if (mode==1){ //low R9 Photons

    corr = ETHZ_fBremEta(myphoton->superCluster()->phiWidth()/myphoton->superCluster()->etaWidth(), myphoton->superCluster()->Eta(), 1);

    float et = energy*TMath::Sin(2*TMath::ATan(TMath::Exp(-myphoton->superCluster()->Eta())))/corr;

    if (myphoton->isEBPho()) corr2 = corr * ETHZ_fEt(et, 2);
    if (myphoton->isEEPho()) corr2 = corr * ETHZ_fEnergy(energy/corr, 3);

    newEnergy = energy/corr2; 

  }

// WITHOUT CRACK CORRECTIONS
//  return newEnergy;
// WITH CRACK CORRECTIONS
  return newEnergy * myphoton->superCluster()->crackCorrectionEtaPhi();

}

// *****************************************************************************************************
// ******************* by-hand photon corrections
// *****************************************************************************************************
double photonManualCorrectionFactor(TRootPhoton *myphoton, string correctionSet, TClonesArray* clusters, TClonesArray* superClusters, TClonesArray* photons)
{
	int verbositybis = 0;
	int photonIsInCrack_ = photonIsInCrack(    fabs(myphoton->superCluster()->Eta()), fabs(myphoton->superCluster()->Phi())   );
	if( photonIsInCrack_ > 0 )
	{
		if( verbositybis > 0) cout << "Photon is in crack : " << ((photonIsInCrack_ == 1) ? "Eta" : "Phi") << endl;
//		return 1.0; // IF COMMENTED, BE SURE CRACK CORRECTIONS ARE TURNED ON BELOW
	}

	vector<double> param_Ceta;
	parameters_Ceta(param_Ceta, correctionSet);
	double f_eta = fEta(param_Ceta, myphoton->superCluster()->Eta());

	if( correctionSet == "ETHZ" )
	{
		if( verbositybis > 1) cout << "ETHZ_getValue(myphoton, f_eta)= " << ETHZ_getValue(myphoton, f_eta) << endl;
		if( ((myphoton->isEBPho()) && (myphoton->r9()>0.94)) || ((myphoton->isEEPho()) && (myphoton->r9()>0.95)) ) return 1.0;
		return ETHZ_getValue(myphoton, f_eta) / (double)(myphoton->Energy());
	} else if( correctionSet == "Dynamic" )
	{
		cout << "We have: " << endl;
		cout << "photons->GetEntries()= " << photons->GetEntries() << endl;
		cout << "superClusters->GetEntries()= " << superClusters->GetEntries() << endl;
		cout << "clusters->GetEntries()= " << clusters->GetEntries() << endl;
	}
	else if( correctionSet == "MITregression" )
	{
		//if( verbositybis > 1) cout << "myphoton->energyRegression() = " << myphoton->energyRegression() << endl;
		//return myphoton->energyRegression() / (double)(myphoton->Energy());
	
		///// NEW VARIABLE IN TOTO FOR PHOTON REGRESSION /////

		if( verbositybis > 1) cout << "myphoton->energyRegression() = " << myphoton->energyRegression() << endl;
		return myphoton->energyRegression() / (double)(myphoton->Energy());

	}
	 else {
		vector<double> param_fbrem;
		vector<double> param_feteta;
		parameters_fbrem(param_fbrem, correctionSet, myphoton->isEBPho());
		parameters_feteta(param_feteta, correctionSet, myphoton->isEBPho());
		double brem = (double)(myphoton->superCluster()->phiWidth()) / (double)(myphoton->superCluster()->etaWidth());
		double f_brem = BremCor(param_fbrem, brem);
		double sc_e = (myphoton->isEBPho()==1) ? (f_eta * myphoton->superCluster()->rawEnergy()) : (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy());
		double sc_e_noCrack = (myphoton->isEBPho()==1) ? (f_eta * myphoton->superCluster()->rawEnergy()) : (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy());
		double f_crack = myphoton->superCluster()->crackCorrectionEtaPhi();
		double sc_et = sc_e * (sin(myphoton->superCluster()->Theta()));
	  double sc_et_noCrack = sc_e_noCrack * (sin(myphoton->superCluster()->Theta()));
		double f_et_eta = EtEtaCor(param_feteta, f_brem * sc_et, myphoton->superCluster()->Eta(), myphoton->isEBPho());
	  double f_et_eta_noCrack = EtEtaCor(param_feteta, f_brem * sc_et_noCrack, myphoton->superCluster()->Eta(), myphoton->isEBPho());
	  if( verbositybis > 1)
		{
			cout << "###\tmyphoton->superCluster()->crackCorrectionEta()= " << myphoton->superCluster()->crackCorrectionEta() << endl;
    	cout << "###\tmyphoton->superCluster()->crackCorrectionPhi()= " << myphoton->superCluster()->crackCorrectionPhi() << endl;
    	cout << "###\tmyphoton->superCluster()->crackCorrectionEtaPhi()= " << myphoton->superCluster()->crackCorrectionEtaPhi() << endl;
		}
		if( (myphoton->isEBPho()) && (myphoton->r9()<0.94) )
		{
	    if( verbositybis > 1)
			{
				cout << "f_et_eta * f_brem * f_eta * myphoton->superCluster()->rawEnergy()= " << f_et_eta * f_brem * f_eta * myphoton->superCluster()->rawEnergy() << endl;
    		cout << "###\tf_et_eta_noCrack * f_brem * f_eta * sc_e_noCrack * f_crack= " << f_et_eta_noCrack * f_brem * f_eta * sc_e_noCrack * f_crack << endl;
			}
// WITH CRACK CORRECTIONS
    return (double)(f_et_eta_noCrack * f_brem * sc_e_noCrack * f_crack) / (double)(myphoton->Energy());
// WITHOUT CRACK CORRECTIONS
//		return (double)(f_et_eta * f_brem * f_eta * myphoton->superCluster()->rawEnergy()) / (double)(myphoton->Energy());
		}
		if( (myphoton->isEBPho()) && (myphoton->r9()>0.94) )
		{
			if( verbositybis > 1) cout << "f_eta * myphoton->e5x5()= " << f_eta * myphoton->e5x5() << endl;
			return (double)(f_eta * myphoton->e5x5()) / (double)(myphoton->Energy());
		}
		if( (myphoton->isEEPho()) && (myphoton->r9()<0.95) )
		{
			if( verbositybis > 1) cout << "f_et_eta * f_brem * (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy())= " << f_et_eta * f_brem * (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy()) << endl;
			return (double)(f_et_eta * f_brem * (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy())) / (double)(myphoton->Energy());
		}
		if( (myphoton->isEEPho()) && (myphoton->r9()>0.95) )
		{
			if( verbositybis > 1) cout << "(myphoton->e5x5() + myphoton->preshowerEnergy())= " << (myphoton->e5x5() + myphoton->preshowerEnergy()) <<  endl;
			return (double)(myphoton->e5x5() + myphoton->preshowerEnergy()) / (double)(myphoton->Energy());
		}
	} // end if classical corrections
}



// *****************************************************************************************************
// ******************* Find conversion MC Truth
// *****************************************************************************************************
void findConversionMCtruth(TRootPhoton *myPhoton, TClonesArray *theMCphotons, int &Photon_MCisConverted, float &Photon_MCconvEoverP, float &Photon_MCconvMass, float &Photon_MCconvCotanTheta, float &Photon_MCconvVertexX, float &Photon_MCconvVertexY, float &Photon_MCconvVertexZ)
{
float dr = 0;
int theIteMin = -1000;
float theDiff;
float theMinDiff = 100000;
for (unsigned int i =0 ; i < theMCphotons->GetEntriesFast() ; i++){
	TRootMCPhoton *theMCphoton = (TRootMCPhoton*) theMCphotons->At(i);
	//dr = DeltaR(myPhoton->Phi(),theMCphoton->Phi(),myPhoton->Eta(),theMCphoton->Eta());
	dr = DeltaR(myPhoton->Eta(),myPhoton->Phi(),theMCphoton->Eta(),theMCphoton->Phi());
	cout<<endl<<"dr = "<<dr<<endl;
	if (dr < 0.3){
		theDiff = fabs(theMCphoton->Pt()-myPhoton->Pt());
		if (theDiff < theMinDiff){
			theMinDiff = theDiff;
			theIteMin = i;	
		}
	}
}
if (theIteMin == -1000){
	Photon_MCisConverted = 0;
	Photon_MCconvEoverP = -10000;
	Photon_MCconvMass = -10000;
	Photon_MCconvCotanTheta = -10000;
	Photon_MCconvVertexX = -10000;
	Photon_MCconvVertexY = -10000;
	Photon_MCconvVertexZ = -10000;
}
else {
	TRootMCPhoton *theMCphoton = (TRootMCPhoton*) theMCphotons->At(theIteMin);
	Photon_MCisConverted = 1;
	Photon_MCconvEoverP = theMCphoton->convEoverP();
	Photon_MCconvMass = theMCphoton->convMass();
	Photon_MCconvCotanTheta = theMCphoton->convDeltaCotanTheta();
	Photon_MCconvVertexX = theMCphoton->conv_vx();
	Photon_MCconvVertexY = theMCphoton->conv_vy();
	Photon_MCconvVertexZ = theMCphoton->conv_vz();

}
}



// *****************************************************************************************************
// ******************* Tentative to remove the muon from photon tracker isolation cone
// *****************************************************************************************************
double PtMuonsConefunction(TRootMuon *mymuon, TRootPhoton *myphoton, double deltaR)
{
        double PtMuonsCone = 0;
    
        double dzCut = 0.0;
        TRootTrack trackbis = (TRootTrack)mymuon->innerTrack();
        dzCut = trackbis.vz() - myphoton->vz();
        double etLow_ = 0.0;
        double lip_ = 0.2;
        double drb_ = 0.1;
        if(dzCut < lip_ && mymuon->innerTrack().Et() > 0 && fabs(mymuon->dB()) < drb_)
        {
                double dr = DeltaR(mymuon->innerTrack().Eta(), mymuon->innerTrack().Phi(), myphoton->Eta(), myphoton->Phi());
                double deta = mymuon->innerTrack().Eta() - myphoton->Eta();
                if(fabs(myphoton->Eta()) < 1.479)
                {
                        double extRadius_ = deltaR;
                        double intRadiusBarrel_ = 0.04;
                        double stripBarrel_ = 0.015;
                        if(fabs(dr) < extRadius_ && fabs(dr) >= intRadiusBarrel_ && fabs(deta) >= stripBarrel_)    
                        {
                            PtMuonsCone = mymuon->innerTrack().Pt();
                            //cout<<"mymuon->innerTrack().Pt() = "<<mymuon->innerTrack().Pt()<<endl;
                        }
                }
                else
                {
                        double extRadius_ = deltaR;
                        double intRadiusBarrel_ = 0.04;
                        double stripBarrel_ = 0.015;
                        if(fabs(dr) < extRadius_ && fabs(dr) >= intRadiusBarrel_ && fabs(deta) >= stripBarrel_)
                        {
                            PtMuonsCone = mymuon->innerTrack().Pt();
                            //cout<<"mymuon->innerTrack().Pt() = "<<mymuon->innerTrack().Pt()<<endl;
                        }
                }

        }

        return PtMuonsCone;
}

string DoubleToString(double x)
{

        std::string s;
        {
                std::ostringstream oss;
                oss << x;
                s = oss.str();
        }
        std::cout << "x = " << x << " s = " << s << std::endl;

        return s;
}

double StringToDouble(string ligne)
{
        std::istringstream stm;
        stm.str(ligne.c_str());
        double d;
        stm >>d;

        return d;
}


// *****************************************************************************************************
// ******************* Declaration of main function
// *****************************************************************************************************
int main(int argc, char *argv[]);


  // ____________________________________________
  // Event information
  // ____________________________________________
  extern ULong64_t iEvent, iEventID, iLumiID, iRunID;
  extern Int_t isSignalApplied, isStewApplied, isZJetsApplied;

  extern Int_t isBeforeAllCuts, isAfterCutPthatFilter, isAfterCutZJETVETO;
  extern Int_t isVeryLooseMMG, isLooseMMG, isMM, isTightMMG, isMMGCandidate;
  extern Int_t isAfterFSRCut1, isAfterFSRCut2, isAfterFSRCut3;
  extern Int_t isAfterFSRCut4, isMultipleCandidate, isAfterCut5, isAfterCut6, isAfterCut7, isAfterCut8, isAfterCut9, isAfterCut10;
  extern Int_t isSelected;

  extern Int_t isNotCommissionned;
  
  extern Int_t Muon_eventPassHLT_Mu11;
  extern Int_t nVertices;
  extern Int_t nGenVertices;
  extern Float_t weight_pileUp, weight_Xsection;

  extern Float_t rho;
  extern Float_t pu_TrueNumInteractions;
  extern Int_t pu_NumInteractions, inTimePU_NumInteractions, latePU_NumInteractions, earlyPU_NumInteractions, outOfTimePU_NumInteractions;
  extern Int_t pu_NumInteractions_inAcceptance, inTimePU_NumInteractions_inAcceptance, latePU_NumInteractions_inAcceptance, earlyPU_NumInteractions_inAcceptance, outOfTimePU_NumInteractions_inAcceptance;
  extern ULong64_t storeNumber, bunchCrossing, orbitNumber, collisionTimeStamp, microsecondCollisionTime;
  extern Float_t collisionTime;


  // ____________________________________________
  // Muon variables 
  // ____________________________________________
  extern Int_t NbMuons;
  
  extern Float_t Pt_allMuons, Eta_allMuons, Phi_allMuons, Charge_allMuons;
// (M minus charge, P plus charge), (F far, N near), (L leading, S subleading)
  extern Float_t MuonM_Pt, MuonP_Pt, MuonN_Pt, MuonF_Pt, MuonL_Pt, MuonS_Pt;
  extern Float_t MuonM_Eta, MuonP_Eta, MuonN_Eta, MuonF_Eta, MuonL_Eta, MuonS_Eta;
  extern Float_t MuonM_Phi, MuonP_Phi, MuonN_Phi, MuonF_Phi, MuonL_Phi, MuonS_Phi;
  extern Float_t MuonM_E, MuonP_E, MuonN_E, MuonF_E, MuonL_E, MuonS_E;
  extern Float_t MuonM_Px, MuonP_Px, MuonN_Px, MuonF_Px, MuonL_Px, MuonS_Px;
  extern Float_t MuonM_Py, MuonP_Py, MuonN_Py, MuonF_Py, MuonL_Py, MuonS_Py;
  extern Float_t MuonM_Pz, MuonP_Pz, MuonN_Pz, MuonF_Pz, MuonL_Pz, MuonS_Pz;
  extern Int_t MuonF_Charge, MuonN_Charge, MuonL_Charge, MuonS_Charge;

  extern Float_t MuonM_isoR03_emEt, MuonP_isoR03_emEt, MuonN_isoR03_emEt, MuonF_isoR03_emEt, MuonL_isoR03_emEt, MuonS_isoR03_emEt; 
  extern Float_t MuonM_isoR03_hadEt, MuonP_isoR03_hadEt, MuonN_isoR03_hadEt, MuonF_isoR03_hadEt, MuonL_isoR03_hadEt, MuonS_isoR03_hadEt;
  extern Float_t MuonM_isoR03_hoEt, MuonP_isoR03_hoEt, MuonN_isoR03_hoEt, MuonF_isoR03_hoEt, MuonL_isoR03_hoEt, MuonS_isoR03_hoEt;   
  extern Float_t MuonM_isoR03_nJets, MuonP_isoR03_nJets, MuonN_isoR03_nJets, MuonF_isoR03_nJets, MuonL_isoR03_nJets, MuonS_isoR03_nJets;
  extern Float_t MuonM_isoR03_nTracks, MuonP_isoR03_nTracks, MuonN_isoR03_nTracks, MuonF_isoR03_nTracks, MuonL_isoR03_nTracks, MuonS_isoR03_nTracks;
  extern Float_t MuonM_isoR03_sumPt, MuonP_isoR03_sumPt, MuonN_isoR03_sumPt, MuonF_isoR03_sumPt, MuonL_isoR03_sumPt, MuonS_isoR03_sumPt;
    
  extern Float_t MuonM_isoR05_emEt, MuonP_isoR05_emEt, MuonN_isoR05_emEt, MuonF_isoR05_emEt, MuonL_isoR05_emEt, MuonS_isoR05_emEt;
  extern Float_t MuonM_isoR05_hadEt, MuonP_isoR05_hadEt, MuonN_isoR05_hadEt, MuonF_isoR05_hadEt, MuonL_isoR05_hadEt, MuonS_isoR05_hadEt;
  extern Float_t MuonM_isoR05_hoEt, MuonP_isoR05_hoEt, MuonN_isoR05_hoEt, MuonF_isoR05_hoEt, MuonL_isoR05_hoEt, MuonS_isoR05_hoEt;
  extern Float_t MuonM_isoR05_nJets, MuonP_isoR05_nJets, MuonN_isoR05_nJets, MuonF_isoR05_nJets, MuonL_isoR05_nJets, MuonS_isoR05_nJets;
  extern Float_t MuonM_isoR05_nTracks, MuonP_isoR05_nTracks, MuonN_isoR05_nTracks, MuonF_isoR05_nTracks, MuonL_isoR05_nTracks, MuonS_isoR05_nTracks;
  extern Float_t MuonM_isoR05_sumPt, MuonP_isoR05_sumPt, MuonN_isoR05_sumPt, MuonF_isoR05_sumPt, MuonL_isoR05_sumPt, MuonS_isoR05_sumPt;

  // ____________________________________________
  // Photon variables
  // ___________________________________________
  extern Int_t NbPhotons;
  extern Float_t Pt_allPhotons, Eta_allPhotons, Phi_allPhotons, Cross_allPhotons;
  extern Int_t isEBorEE_allPhotons, isEB_allPhotons, isEE_allPhotons, isEEP_allPhotons, isEEM_allPhotons;
  extern Float_t Photon_Eta, Photon_Phi;
	extern Float_t Photon_Px, Photon_Py, Photon_Pz;
  extern Int_t Photon_isEBorEE, Photon_isEB, Photon_isEE, Photon_isEEP, Photon_isEEM;

  extern Int_t Photon_hasPixelSeed, Photon_isAlsoElectron, Photon_Nclusters, Photon_nBasicClusters, Photon_nXtals;
  //extern Int_t Photon_isTightPhoton, Photon_isLoosePhoton;
  extern Int_t Photon_convNTracks, Photon_isConverted;
  extern Float_t Photon_convEoverP, Photon_convMass, Photon_convCotanTheta, Photon_convLikely, Photon_convVertexX, Photon_convVertexY, Photon_convVertexZ;
  extern Float_t Photon_E, Photon_Et, Photon_E2x2, Photon_E3x3, Photon_E5x5, Photon_Emax, Photon_E2nd;
	extern Float_t Photon_Ecorr_o_Ereco;
  extern Float_t Photon_r19, Photon_r9, Photon_cross;
  //extern Float_t Photon_caloConeSize; 
  extern Float_t Photon_PreshEnergy, Photon_HoE;
  extern Float_t Photon_sigmaEtaEta, Photon_sigmaIetaIeta;
  extern Float_t Photon_covEtaEta, Photon_covPhiPhi, Photon_covEtaPhi;
  extern Float_t Photon_etaWidth, Photon_phiWidth;
  extern Float_t Photon_dR03isoEcalRecHit, Photon_dR03isoHcalRecHit, Photon_dR03isoSolidTrkCone, Photon_dR03isoHollowTrkCone, Photon_dR03isoNTracksSolidCone, Photon_dR03isoNTracksHollowCone;
  extern Float_t Photon_dR04isoEcalRecHit, Photon_dR04isoHcalRecHit, Photon_dR04isoSolidTrkCone, Photon_dR04isoHollowTrkCone, Photon_dR04isoNTracksSolidCone, Photon_dR04isoNTracksHollowCone;
  extern Float_t Photon_seedTime, Photon_seedFlag;
  extern Int_t Photon_seedPosition1, Photon_seedPosition2;
  extern Float_t Photon_SC_Eta, Photon_SC_Phi, Photon_SC_brem;
  extern Float_t Photon_SC_E, Photon_SC_Et, Photon_SC_rawE, Photon_SC_rawEt;
	extern Float_t Photon_lambdaRatio, Photon_ratioSeed, Photon_ratioS4, Photon_lamdbaDivCov;
	extern Float_t Photon_ratioS4_corrected;
	extern Float_t Photon_SC_rawE_x_fEta, Photon_SC_rawE_x_fEta_x_fBrem, Photon_SC_rawE_x_fEta_x_fBrem_AF, Photon_SC_rawE_x_fEta_x_fBrem_L, Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta, Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta, Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta;
	extern Float_t Photon_secondMomentMaj, Photon_secondMomentMin, Photon_secondMomentAlpha;
	extern Float_t Photon_etaLAT, Photon_phiLAT, Photon_LAT, Photon_Zernike20, Photon_Zernike42, Photon_ESratio;

	extern Float_t Photon_E_regression, Photon_E_regressionError, Photon_Et_regression;

  // ____________________________________________
  // mugamma / mumu / mumugamma information
  // ____________________________________________

  extern Float_t Mmumu, Mmumugamma, Mmumugamma_5x5, Mmumugamma_SC, Mmumugamma_SCraw, Mmumugamma_SCraw_fEta, Mmumugamma_SCraw_fEta_fBrem, Mmumugamma_SCraw_fEta_fBrem_AF, Mmumugamma_SCraw_fEta_fBrem_L, Mmumugamma_SCraw_fEta_fBrem_fEtEta, Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta, Mmumugamma_SCraw_fEta_fBrem_L_fEtEta;
  extern Float_t Ptmumu;
  extern Float_t deltaRNear, deltaRFar, deltaRMinus, deltaRPlus, deltaRLeading, deltaRSubleading;
  extern Float_t mmg_k, mmg_ik, mmg_s, mmg_logk, mmg_logik, mmg_logs;
  extern Float_t mmg_k_5x5, mmg_ik_5x5, mmg_s_5x5, mmg_logk_5x5, mmg_logik_5x5, mmg_logs_5x5;
  extern Float_t mmg_k_SC, mmg_ik_SC, mmg_s_SC, mmg_logk_SC, mmg_logik_SC, mmg_logs_SC;
  extern Float_t mmg_k_SCraw, mmg_ik_SCraw, mmg_s_SCraw, mmg_logk_SCraw, mmg_logik_SCraw, mmg_logs_SCraw;
  extern Float_t mmg_k_SCraw_fEta, mmg_ik_SCraw_fEta, mmg_s_SCraw_fEta, mmg_logk_SCraw_fEta, mmg_logik_SCraw_fEta, mmg_logs_SCraw_fEta;
  extern Float_t mmg_ik_SCraw_fEta_fBrem, mmg_ik_SCraw_fEta_fBrem_AF, mmg_ik_SCraw_fEta_fBrem_L, mmg_ik_SCraw_fEta_fBrem_fEtEta, mmg_ik_SCraw_fEta_fBrem_AF_fEtEta, mmg_ik_SCraw_fEta_fBrem_L_fEtEta;

	extern Float_t MuonBeforeBremM_Pt, MuonBeforeBremP_Pt, MuonBeforeBremN_Pt, MuonBeforeBremF_Pt, MuonBeforeBremL_Pt, MuonBeforeBremS_Pt;
	extern Float_t MuonBeforeBremM_Eta, MuonBeforeBremP_Eta, MuonBeforeBremN_Eta, MuonBeforeBremF_Eta, MuonBeforeBremL_Eta, MuonBeforeBremS_Eta;
	extern Float_t MuonBeforeBremM_Phi, MuonBeforeBremP_Phi, MuonBeforeBremN_Phi, MuonBeforeBremF_Phi, MuonBeforeBremL_Phi, MuonBeforeBremS_Phi;
	extern Float_t MuonBeforeBremM_E, MuonBeforeBremP_E, MuonBeforeBremN_E, MuonBeforeBremF_E, MuonBeforeBremL_E, MuonBeforeBremS_E;
	extern Float_t MuonBeforeBremM_Px, MuonBeforeBremP_Px, MuonBeforeBremN_Px, MuonBeforeBremF_Px, MuonBeforeBremL_Px, MuonBeforeBremS_Px;
	extern Float_t MuonBeforeBremM_Py, MuonBeforeBremP_Py, MuonBeforeBremN_Py, MuonBeforeBremF_Py, MuonBeforeBremL_Py, MuonBeforeBremS_Py;
	extern Float_t MuonBeforeBremM_Pz, MuonBeforeBremP_Pz, MuonBeforeBremN_Pz, MuonBeforeBremF_Pz, MuonBeforeBremL_Pz, MuonBeforeBremS_Pz;
	extern Int_t MuonBeforeBremF_Charge, MuonBeforeBremN_Charge, MuonBeforeBremL_Charge, MuonBeforeBremS_Charge;

  // ____________________________________________
  // Neural Network variables
  // ____________________________________________
	extern Float_t Photon_NNshapeOutput;

  // ____________________________________________
  // Surface variables
  // ____________________________________________
        extern Float_t MZ_Surface;
	extern Float_t mmg_k_MZ_Surface, mmg_ik_MZ_Surface, mmg_s_MZ_Surface, mmg_logk_MZ_Surface, mmg_logik_MZ_Surface, mmg_logs_MZ_Surface;	

  // ____________________________________________
  // MC Truth
  // ___________________________________________

  extern Float_t Photon_MC_E, Photon_MC_Px, Photon_MC_Py, Photon_MC_Pz, Photon_MC_Phi, Photon_MC_Eta, Photon_MC_Pt;
  extern Int_t Photon_MCisConverted;
  extern Float_t Photon_MCconvEoverP, Photon_MCconvMass, Photon_MCconvCotanTheta, Photon_MCconvVertexX, Photon_MCconvVertexY, Photon_MCconvVertexZ;
  extern Float_t MuonM_MC_E, MuonM_MC_Px, MuonM_MC_Py, MuonM_MC_Pz, MuonM_MC_Phi, MuonM_MC_Eta, MuonM_MC_Pt;
  extern Float_t MuonP_MC_E, MuonP_MC_Px, MuonP_MC_Py, MuonP_MC_Pz, MuonP_MC_Phi, MuonP_MC_Eta, MuonP_MC_Pt;
  extern Float_t MuonN_MC_E, MuonN_MC_Px, MuonN_MC_Py, MuonN_MC_Pz, MuonN_MC_Phi, MuonN_MC_Eta, MuonN_MC_Pt;
  extern Float_t MuonF_MC_E, MuonF_MC_Px, MuonF_MC_Py, MuonF_MC_Pz, MuonF_MC_Phi, MuonF_MC_Eta, MuonF_MC_Pt;
  extern Float_t MuonL_MC_E, MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_Phi, MuonL_MC_Eta, MuonL_MC_Pt;
  extern Float_t MuonS_MC_E, MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_Phi, MuonS_MC_Eta, MuonS_MC_Pt;
  extern Float_t Photon_SC_rawE_x_fEta_o_MC_E, Photon_E_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_AF_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_L_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta_o_MC_E, Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta_o_MC_E;

  extern Float_t Mmumu_Photon_MC, Mmumugamma_Photon_MC, mmg_k_Photon_MC, mmg_ik_Photon_MC, mmg_s_Photon_MC, mmg_logk_Photon_MC, mmg_logik_Photon_MC, mmg_logs_Photon_MC;
  extern Float_t Mmumu_Muons_MC, Mmumugamma_Muons_MC, mmg_k_Muons_MC, mmg_ik_Muons_MC, mmg_s_Muons_MC, mmg_logk_Muons_MC, mmg_logik_Muons_MC, mmg_logs_Muons_MC;
  extern Float_t Mmumu_MMG_MC, Mmumugamma_MMG_MC, mmg_k_MMG_MC, mmg_ik_MMG_MC, mmg_s_MMG_MC, mmg_logk_MMG_MC, mmg_logik_MMG_MC, mmg_logs_MMG_MC;

  extern Float_t mmg_k_MZ, mmg_ik_MZ, mmg_s_MZ, mmg_logk_MZ, mmg_logik_MZ, mmg_logs_MZ;
  extern Float_t mmg_k_MZ_Photon_MC, mmg_ik_MZ_Photon_MC, mmg_s_MZ_Photon_MC, mmg_logk_MZ_Photon_MC, mmg_logik_MZ_Photon_MC, mmg_logs_MZ_Photon_MC;
  extern Float_t mmg_k_MZ_Muons_MC, mmg_ik_MZ_Muons_MC, mmg_s_MZ_Muons_MC, mmg_logk_MZ_Muons_MC, mmg_logik_MZ_Muons_MC, mmg_logs_MZ_Muons_MC;
  extern Float_t mmg_k_MZ_Muons_RECO_MC, mmg_ik_MZ_Muons_RECO_MC, mmg_s_MZ_Muons_RECO_MC, mmg_logk_MZ_Muons_RECO_MC, mmg_logik_MZ_Muons_RECO_MC, mmg_logs_MZ_Muons_RECO_MC;

//int FillMMG(TRootPhoton* myphoton, TRootMuon* mymuon1, TRootMuon* mymuon2, double EScale, bool doMC, TClonesArray* mcParticles);

// *****************************************************************************************************
// ******************* SIDRA muon corrections
// *****************************************************************************************************
double applySidra( double _pt, double charge, double eta, double phi, TRandom3* generator)
{
	double pt = _pt;
// Correct MC
	double a = 0.0650687e-3;
	double b = 0.212987e-3;
	double c = 1.53414;
	pt = (double)1.0/(double)pt;
	pt -= a - b * charge * sin( phi + c );
// Apply Corrections
	double A = 0.143812;
	double B = 0.0404834;
	double Ap = 0.0995898;
	double Bp = -0.0692569;
	double Cp = 0.0952174;
	double phi0 = -1.08881;
	pt += (double)(( A + B * eta * eta ) * (generator->Gaus(0,1))) /(double)(1000.) + (double)(Ap + Cp * charge * sin( phi + phi0 )+ Bp * charge * eta )/(double)(1000.);
	pt = (double)1.0/(double)pt;
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
	pt = _pt * (1.0 + b * _pt  + c * sgn_charge * _pt * sgn_eta * eta * eta + sgn_charge * d * _pt * sin( phi + e ));
	return pt;

}



// *****************************************************************************************************
// ******************* Fill minitree
// *****************************************************************************************************
//int FillMMG(TRootPhoton* myphoton, TRootMuon* mymuon1, TRootMuon* mymuon2, double EScale, bool doMC, bool doPhotonConversionMC, TClonesArray* mcParticles, TClonesArray* mcPhotons, TMVA::Reader* reader){
//int FillMMG(TRootPhoton* myphoton, TRootMuon* mymuon1, TRootMuon* mymuon2, double EScale, bool doMC, bool doPhotonConversionMC, TClonesArray* mcParticles, TMVA::Reader* reader){
int FillMMG(TRootPhoton* myphoton, TRootMuon* mymuon1, TRootMuon* mymuon2, TLorentzVector* correctedmymuon1, TLorentzVector* correctedmymuon2, double EScale, bool doMC, bool doPhotonConversionMC, TClonesArray* mcParticles, TMVA::Reader* reader, int binNumber){

//			cout << "correctedmymuon1->Pt()= " << correctedmymuon1->Pt() << endl;
//			cout << "correctedmymuon2->Pt()= " << correctedmymuon2->Pt() << endl;

      // Fill photon stuff
      Photon_Eta = myphoton->Eta();
      Photon_Phi = myphoton->Phi();
      Photon_Px = myphoton->Px();
      Photon_Py = myphoton->Py();
      Photon_Pz = myphoton->Pz();
      if( myphoton->isEBPho() ){ Photon_isEB=1; } else { Photon_isEB=0; }
      if( myphoton->isEEPho() ){ Photon_isEE=1; } else { Photon_isEE=0; }
      if( myphoton->isEEPho() && myphoton->Eta()<0 ){ Photon_isEEM=1; } else { Photon_isEEM=0; }
      if( myphoton->isEEPho() && myphoton->Eta()>0 ){ Photon_isEEP=1; } else { Photon_isEEP=0; }
      Photon_hasPixelSeed = myphoton->hasPixelSeed();
      Photon_isAlsoElectron = myphoton->isAlsoElectron();
      Photon_Nclusters = myphoton->nbClusters();
  //    cout << "myphoton->superCluster()=" << myphoton->superCluster() << endl; // FIXME
      Photon_nBasicClusters = myphoton->superCluster()->nBasicClusters();
      //Photon_isTightPhoton = myphoton->isTightPhoton();
      //Photon_isLoosePhoton = myphoton->isLoosePhoton();

      Photon_isConverted = 0;
      if (myphoton->convNTracks() > 0 ) Photon_isConverted = 1;	
      Photon_convNTracks = myphoton->convNTracks();
      Photon_convEoverP = myphoton->convEoverP();
      Photon_convMass = myphoton->convMass();
      Photon_convCotanTheta = myphoton->convCotanTheta();
      Photon_convLikely = myphoton->convLikely();
      Photon_convVertexX = myphoton->convVertex().x();
      Photon_convVertexY = myphoton->convVertex().y();
      Photon_convVertexZ = myphoton->convVertex().z();

/*
      if (doPhotonConversionMC)
      {
				findConversionMCtruth(myphoton, mcPhotons, Photon_MCisConverted, Photon_MCconvEoverP, Photon_MCconvMass, Photon_MCconvCotanTheta, Photon_MCconvVertexX, Photon_MCconvVertexY, Photon_MCconvVertexZ);
      }
*/


 
      Photon_E = EScale*(myphoton->Energy());
			Photon_Ecorr_o_Ereco = EScale;
      Photon_Et = EScale*(myphoton->Et());
      Photon_E2x2 = EScale*(myphoton->e2x2());
      Photon_E3x3 = EScale*(myphoton->e3x3());
      Photon_E5x5 = EScale*(myphoton->e5x5());
      Photon_Emax = EScale*(myphoton->eMax());
      Photon_E2nd = EScale*(myphoton->e2nd());
      Photon_r19 = myphoton->r19();
      Photon_r9 = myphoton->r9();
  //FIXME   Photon_cross = 1-((myphoton->superCluster()->s4())/(myphoton->superCluster()->eMax()));   
      //Photon_caloConeSize = myphoton->caloConeSize();
      Photon_PreshEnergy = myphoton->preshowerEnergy();
      Photon_HoE = myphoton->hoe();
      Photon_sigmaEtaEta = myphoton->sigmaEtaEta();
      Photon_sigmaIetaIeta = myphoton->sigmaIetaIeta();
      Photon_covEtaEta = myphoton->covEtaEta();
      Photon_covPhiPhi = myphoton->covPhiPhi();
      Photon_covEtaPhi = myphoton->covEtaPhi();
      Photon_etaWidth = myphoton->superCluster()->etaWidth();
      Photon_phiWidth = myphoton->superCluster()->phiWidth();
      Photon_dR03isoEcalRecHit = myphoton->dR03IsolationEcalRecHit();
      Photon_dR03isoHcalRecHit = myphoton->dR03IsolationHcalRecHit();
      Photon_dR03isoSolidTrkCone = myphoton->dR03IsolationSolidTrkCone();
      Photon_dR03isoHollowTrkCone = myphoton->dR03IsolationHollowTrkCone();
      Photon_dR03isoNTracksSolidCone = myphoton->dR03IsolationNTracksSolidCone();
      Photon_dR03isoNTracksHollowCone = myphoton->dR03IsolationNTracksHollowCone();
      Photon_dR04isoEcalRecHit = myphoton->dR04IsolationEcalRecHit();
      Photon_dR04isoHcalRecHit = myphoton->dR04IsolationHcalRecHit();
      Photon_dR04isoSolidTrkCone = myphoton->dR04IsolationSolidTrkCone();
      Photon_dR04isoHollowTrkCone = myphoton->dR04IsolationHollowTrkCone();
      Photon_dR04isoNTracksSolidCone = myphoton->dR04IsolationNTracksSolidCone();
      Photon_dR04isoNTracksHollowCone = myphoton->dR04IsolationNTracksHollowCone();
      Photon_seedTime = myphoton->superCluster()->seedTime();
      Photon_seedFlag = myphoton->superCluster()->seedRecoFlag();
      Photon_seedPosition1 = myphoton->superCluster()->seedPosition1();
      Photon_seedPosition2 = myphoton->superCluster()->seedPosition2();
      Photon_SC_Eta = myphoton->superCluster()->Eta();
      Photon_SC_Phi = myphoton->superCluster()->Phi();
      Photon_SC_brem = (double)(Photon_phiWidth) / (double)(Photon_etaWidth);
      Photon_SC_E = myphoton->superCluster()->Mag();
      Photon_SC_Et = Photon_SC_E * (sin(myphoton->superCluster()->Theta()));
      Photon_SC_rawE = myphoton->superCluster()->rawEnergy();
      Photon_SC_rawEt = Photon_SC_rawE * (sin(myphoton->superCluster()->Theta()));
			Photon_lambdaRatio = 0.0;
	    if ((Photon_covEtaEta+Photon_covPhiPhi+sqrt((Photon_covEtaEta-Photon_covPhiPhi)*(Photon_covEtaEta-Photon_covPhiPhi)+4*Photon_covEtaPhi*Photon_covEtaPhi))!=0) Photon_lambdaRatio = (Photon_covEtaEta+Photon_covPhiPhi-sqrt((Photon_covEtaEta-Photon_covPhiPhi)*(Photon_covEtaEta-Photon_covPhiPhi)+4*Photon_covEtaPhi*Photon_covEtaPhi))/(Photon_covEtaEta+Photon_covPhiPhi+sqrt((Photon_covEtaEta-Photon_covPhiPhi)*(Photon_covEtaEta-Photon_covPhiPhi)+4*Photon_covEtaPhi*Photon_covEtaPhi));
			Photon_ratioSeed = 0.0;
			if (Photon_SC_rawE != 0) Photon_ratioSeed = Photon_Emax/Photon_SC_rawE; 
			Photon_ratioS4 = 0.0;
			if (Photon_E5x5 != 0)  Photon_ratioS4 = Photon_E2x2/Photon_E5x5;
			if(Photon_isEB) Photon_ratioS4_corrected = 1.008 * Photon_ratioS4;
			if(Photon_isEE) Photon_ratioS4_corrected = 1.008074 * Photon_ratioS4;
			Photon_lamdbaDivCov = 0.0;
			if (Photon_covEtaEta != 0) Photon_lamdbaDivCov = (Photon_covEtaEta+Photon_covPhiPhi-sqrt((Photon_covEtaEta-Photon_covPhiPhi)*(Photon_covEtaEta-Photon_covPhiPhi)+4*Photon_covEtaPhi*Photon_covEtaPhi))/Photon_covEtaEta;
//			Photon_SC_rawE_x_fEta = Photon_SC_rawE * fEta(Photon_SC_Eta);
			
//			Photon_SC_rawE_x_fEta_x_fBrem = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), 0);
//			Photon_SC_rawE_x_fEta_x_fBrem_AF = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), 1);
//			Photon_SC_rawE_x_fEta_x_fBrem_L = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), 2);
/*
			Photon_SC_rawE_x_fEta_x_fBrem = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11");
			Photon_SC_rawE_x_fEta_x_fBrem_AF = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur");
			Photon_SC_rawE_x_fEta_x_fBrem_L = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis");

			float etraw_fEta_fBrem = Photon_SC_rawE_x_fEta_x_fBrem/cosh(Photon_SC_Eta);
			float etraw_fEta_fBrem_AF = Photon_SC_rawE_x_fEta_x_fBrem_AF/cosh(Photon_SC_Eta);
			float etraw_fEta_fBrem_L = Photon_SC_rawE_x_fEta_x_fBrem_L/cosh(Photon_SC_Eta);
*/
//			Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), 0) * EtEtaCor(etraw_fEta_fBrem, abs(Photon_SC_Eta), myphoton->isEEPho(), 0);
//			Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), 0) * EtEtaCor(etraw_fEta_fBrem_AF, abs(Photon_SC_Eta), myphoton->isEEPho(), 1);
//			Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), 0) * EtEtaCor(etraw_fEta_fBrem_L, abs(Photon_SC_Eta), myphoton->isEEPho(), 1);
/*
			Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11") * EtEtaCor(etraw_fEta_fBrem, abs(Photon_SC_Eta), myphoton->isEEPho(), "START42_V11");
			Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11") * EtEtaCor(etraw_fEta_fBrem_AF, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur");
			Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta = Photon_SC_rawE * fEta(Photon_SC_Eta) * BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11") * EtEtaCor(etraw_fEta_fBrem_L, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur");
*/

			Photon_secondMomentMaj = myphoton->secondMomentMaj();
			Photon_secondMomentMin = myphoton->secondMomentMin();
			Photon_secondMomentAlpha = myphoton->secondMomentAlpha();

			Photon_etaLAT = myphoton->etaLAT();
			Photon_phiLAT = myphoton->phiLAT();
			Photon_LAT = myphoton->lat();
			Photon_Zernike20 = myphoton->zernike20();
			Photon_Zernike42 = myphoton->zernike42();
			Photon_ESratio = myphoton->superCluster()->esRatio();
	

// Read NN output from weight file
			Photon_NNshapeOutput = reader->EvaluateMVA("MLP method");

      // Fill muons stuff
      TRootMuon *leadingMuon;
      TRootMuon *subleadingMuon;
      TLorentzVector *correctedleadingMuon;
      TLorentzVector *correctedsubleadingMuon;
      if( (correctedmymuon1->Pt()) > (correctedmymuon2->Pt()) )      {
        leadingMuon = mymuon1;
        correctedleadingMuon = correctedmymuon1;
        subleadingMuon = mymuon2;
        correctedsubleadingMuon = correctedmymuon2;
      } else {
        leadingMuon = mymuon2;
        correctedleadingMuon = correctedmymuon2;
        subleadingMuon = mymuon1;
        correctedsubleadingMuon = correctedmymuon1;
      }
      MuonL_Pt = correctedleadingMuon->Pt();
      MuonL_Eta = leadingMuon->Eta();
      MuonL_Phi = leadingMuon->Phi();
      MuonL_E = leadingMuon->Energy();
      MuonL_Px = leadingMuon->Px();
      MuonL_Py = leadingMuon->Py();
      MuonL_Pz = leadingMuon->Pz();
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
      MuonS_Pt = correctedsubleadingMuon->Pt();
      MuonS_Eta = subleadingMuon->Eta();
      MuonS_Phi = subleadingMuon->Phi();
      MuonS_E = subleadingMuon->Energy();
      MuonS_Px = subleadingMuon->Px();
      MuonS_Py = subleadingMuon->Py();
      MuonS_Pz = subleadingMuon->Pz();
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

      TLorentzVector mumu;
      mumu = (*correctedleadingMuon) + (*correctedsubleadingMuon);
      Ptmumu = mumu.Pt();
      double mumuInvMass = mumu.M();
//      cerr << "\t\tINFO: Dimuon invariant mass : Mmumu = " << mumuInvMass << endl;
      mumu.Clear();
      Mmumu = mumuInvMass;

    // ********************************************************************
    // *** Compute mumugamma invariant mass ***
    // ********************************************************************
    // Re-compute modified photon momentum so as to keep m_photon = 0
    Double_t Px_SC, Px_5x5, Px_SCraw, Px_SCraw_fEta;
    Double_t Py_SC, Py_5x5, Py_SCraw, Py_SCraw_fEta;
    Double_t Pz_SC, Pz_5x5, Pz_SCraw, Pz_SCraw_fEta;
	
    Double_t Px_SCraw_fEta_fBrem, Px_SCraw_fEta_fBrem_AF, Px_SCraw_fEta_fBrem_L, Px_SCraw_fEta_fBrem_fEtEta, Px_SCraw_fEta_fBrem_AF_fEtEta, Px_SCraw_fEta_fBrem_L_fEtEta;
    Double_t Py_SCraw_fEta_fBrem, Py_SCraw_fEta_fBrem_AF, Py_SCraw_fEta_fBrem_L, Py_SCraw_fEta_fBrem_fEtEta, Py_SCraw_fEta_fBrem_AF_fEtEta, Py_SCraw_fEta_fBrem_L_fEtEta;
    Double_t Pz_SCraw_fEta_fBrem, Pz_SCraw_fEta_fBrem_AF, Pz_SCraw_fEta_fBrem_L, Pz_SCraw_fEta_fBrem_fEtEta, Pz_SCraw_fEta_fBrem_AF_fEtEta, Pz_SCraw_fEta_fBrem_L_fEtEta;



    Px_SC = Px_5x5 = Px_SCraw = Px_SCraw_fEta = 0.0;
    Py_SC = Py_5x5 = Py_SCraw = Py_SCraw_fEta = 0.0;
    Pz_SC = Pz_5x5 = Pz_SCraw = Pz_SCraw_fEta = 0.0;

    Px_SC = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->Mag());
    Py_SC = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->Mag());
    Pz_SC = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->Mag());

    Px_5x5 = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->e5x5());
    Py_5x5 = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->e5x5());
    Pz_5x5 = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->e5x5());

    Px_SCraw = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy());
    Py_SCraw = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy());
    Pz_SCraw = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy());

/*
    Px_SCraw_fEta = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta));
    Py_SCraw_fEta = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta));
    Pz_SCraw_fEta = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta));


    Px_SCraw_fEta_fBrem = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11"));
    Px_SCraw_fEta_fBrem_AF = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur"));
    Px_SCraw_fEta_fBrem_L = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis"));

    Py_SCraw_fEta_fBrem = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11"));
    Py_SCraw_fEta_fBrem_AF = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur"));
    Py_SCraw_fEta_fBrem_L = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis"));

    Pz_SCraw_fEta_fBrem = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11"));
    Pz_SCraw_fEta_fBrem_AF = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur"));
    Pz_SCraw_fEta_fBrem_L = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis"));


    Px_SCraw_fEta_fBrem_fEtEta = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11")) * (EtEtaCor(etraw_fEta_fBrem, abs(Photon_SC_Eta), myphoton->isEEPho(), "START42_V11"));
    Px_SCraw_fEta_fBrem_AF_fEtEta  = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur")) * (EtEtaCor(etraw_fEta_fBrem_AF, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur"));
    Px_SCraw_fEta_fBrem_L_fEtEta  = (myphoton->Px()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis")) * (EtEtaCor(etraw_fEta_fBrem_L, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur"));

    Py_SCraw_fEta_fBrem_fEtEta = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11")) * (EtEtaCor(etraw_fEta_fBrem, abs(Photon_SC_Eta), myphoton->isEEPho(), "START42_V11"));
    Py_SCraw_fEta_fBrem_AF_fEtEta  = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur")) * (EtEtaCor(etraw_fEta_fBrem_AF, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur"));
    Py_SCraw_fEta_fBrem_L_fEtEta  = (myphoton->Py()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis")) * (EtEtaCor(etraw_fEta_fBrem_L, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur"));

    Pz_SCraw_fEta_fBrem_fEtEta = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11")) * (EtEtaCor(etraw_fEta_fBrem, abs(Photon_SC_Eta), myphoton->isEEPho(), "START42_V11"));
    Pz_SCraw_fEta_fBrem_AF_fEtEta  = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur")) * (EtEtaCor(etraw_fEta_fBrem_AF, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur"));
    Pz_SCraw_fEta_fBrem_L_fEtEta  = (myphoton->Pz()) / (myphoton->Energy()) * (myphoton->superCluster()->rawEnergy()) * (fEta(Photon_SC_Eta)) * (BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis")) * (EtEtaCor(etraw_fEta_fBrem_L, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur"));
*/


    TLorentzVector mumugamma;
    TLorentzVector mumuSC;
    TLorentzVector mumu5x5;
    TLorentzVector mumuSC_raw;
/*    TLorentzVector mumuSC_raw_fEta;
	
    TLorentzVector mumuSC_raw_fEta_fBrem, mumuSC_raw_fEta_fBrem_AF, mumuSC_raw_fEta_fBrem_L, mumuSC_raw_fEta_fBrem_fEtEta, mumuSC_raw_fEta_fBrem_AF_fEtEta, mumuSC_raw_fEta_fBrem_L_fEtEta;
*/

    TLorentzVector *PhotonEScale = new TLorentzVector( EScale*(myphoton->Px()), EScale*(myphoton->Py()), EScale*(myphoton->Pz()), EScale*(myphoton->Energy()));
    TLorentzVector *PhotonSC = new TLorentzVector( EScale*Px_SC, EScale*Py_SC, EScale*Pz_SC, EScale*(myphoton->superCluster()->Mag()) );
    TLorentzVector *Photon5x5 = new TLorentzVector(EScale*Px_5x5 , EScale*Py_5x5, EScale*Pz_5x5, EScale*(myphoton->e5x5()));
    TLorentzVector *PhotonSC_raw = new TLorentzVector(EScale*Px_SCraw , EScale*Py_SCraw, EScale*Pz_SCraw, EScale*(myphoton->superCluster()->rawEnergy()));

/*
    TLorentzVector *PhotonSC_raw_fEta = new TLorentzVector(EScale*Px_SCraw_fEta , EScale*Py_SCraw_fEta, EScale*Pz_SCraw_fEta, EScale*(myphoton->superCluster()->rawEnergy())*(fEta(Photon_SC_Eta)));

    TLorentzVector *PhotonSC_raw_fEta_fBrem = new TLorentzVector(EScale*Px_SCraw_fEta_fBrem, EScale*Py_SCraw_fEta_fBrem, EScale*Pz_SCraw_fEta_fBrem, EScale*(myphoton->superCluster()->rawEnergy())*(fEta(Photon_SC_Eta))*(BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11")));

    TLorentzVector *PhotonSC_raw_fEta_fBrem_AF = new TLorentzVector(EScale*Px_SCraw_fEta_fBrem_AF, EScale*Py_SCraw_fEta_fBrem_AF, EScale*Pz_SCraw_fEta_fBrem_AF, EScale*(myphoton->superCluster()->rawEnergy())*(fEta(Photon_SC_Eta))*(BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur")));

    TLorentzVector *PhotonSC_raw_fEta_fBrem_L = new TLorentzVector(EScale*Px_SCraw_fEta_fBrem_L, EScale*Py_SCraw_fEta_fBrem_L, EScale*Pz_SCraw_fEta_fBrem_L, EScale*(myphoton->superCluster()->rawEnergy())*(fEta(Photon_SC_Eta))*(BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis")));

    
    TLorentzVector *PhotonSC_raw_fEta_fBrem_fEtEta = new TLorentzVector(EScale*Px_SCraw_fEta_fBrem_fEtEta, EScale*Py_SCraw_fEta_fBrem_fEtEta, EScale*Pz_SCraw_fEta_fBrem_fEtEta, EScale*(myphoton->superCluster()->rawEnergy())*(fEta(Photon_SC_Eta))*(BremCor(Photon_SC_brem, myphoton->isEEPho(), "START42_V11"))*(EtEtaCor(etraw_fEta_fBrem, abs(Photon_SC_Eta), myphoton->isEEPho(), "START42_V11")));

    TLorentzVector *PhotonSC_raw_fEta_fBrem_AF_fEtEta = new TLorentzVector(EScale*Px_SCraw_fEta_fBrem_AF_fEtEta, EScale*Py_SCraw_fEta_fBrem_AF_fEtEta, EScale*Pz_SCraw_fEta_fBrem_AF_fEtEta, EScale*(myphoton->superCluster()->rawEnergy())*(fEta(Photon_SC_Eta))*(BremCor(Photon_SC_brem, myphoton->isEEPho(), "Anne-Fleur"))*(EtEtaCor(etraw_fEta_fBrem_AF, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur")));

    TLorentzVector *PhotonSC_raw_fEta_fBrem_L_fEtEta = new TLorentzVector(EScale*Px_SCraw_fEta_fBrem_L_fEtEta, EScale*Py_SCraw_fEta_fBrem_L_fEtEta, EScale*Pz_SCraw_fEta_fBrem_L_fEtEta, EScale*(myphoton->superCluster()->rawEnergy())*(fEta(Photon_SC_Eta))*(BremCor(Photon_SC_brem, myphoton->isEEPho(), "Louis"))*(EtEtaCor(etraw_fEta_fBrem_L, abs(Photon_SC_Eta), myphoton->isEEPho(), "Anne-Fleur")));
*/



//    mumugamma = (*leadingMuon) + (*subleadingMuon) + (*myphoton);
    mumugamma = (*correctedleadingMuon) + (*correctedsubleadingMuon) + (*PhotonEScale);
    mumuSC = (*correctedleadingMuon) + (*correctedsubleadingMuon) + (*PhotonSC);
    mumu5x5 = (*correctedleadingMuon) + (*correctedsubleadingMuon) + (*Photon5x5);
    mumuSC_raw = (*correctedleadingMuon) + (*correctedsubleadingMuon) + (*PhotonSC_raw);
/*
    mumuSC_raw_fEta = (*leadingMuon) + (*subleadingMuon) + (*PhotonSC_raw_fEta);

    mumuSC_raw_fEta_fBrem = (*leadingMuon) + (*subleadingMuon) + (*PhotonSC_raw_fEta_fBrem);
    mumuSC_raw_fEta_fBrem_AF = (*leadingMuon) + (*subleadingMuon) + (*PhotonSC_raw_fEta_fBrem_AF);
    mumuSC_raw_fEta_fBrem_L = (*leadingMuon) + (*subleadingMuon) + (*PhotonSC_raw_fEta_fBrem_L);
    mumuSC_raw_fEta_fBrem_fEtEta = (*leadingMuon) + (*subleadingMuon) + (*PhotonSC_raw_fEta_fBrem_fEtEta);
    mumuSC_raw_fEta_fBrem_AF_fEtEta = (*leadingMuon) + (*subleadingMuon) + (*PhotonSC_raw_fEta_fBrem_AF_fEtEta);
    mumuSC_raw_fEta_fBrem_L_fEtEta = (*leadingMuon) + (*subleadingMuon) + (*PhotonSC_raw_fEta_fBrem_L_fEtEta);
*/
    double mumugammaInvMass = mumugamma.M();
    double mumuSCInvMass = mumuSC.M();
    double mumu5x5InvMass = mumu5x5.M();
    double mumuSC_rawInvMass = mumuSC_raw.M();
/*
    double mumuSC_rawInvMass_fEta = mumuSC_raw_fEta.M();
    
    double mumuSC_rawInvMass_fEta_fBrem = mumuSC_raw_fEta_fBrem.M();
    double mumuSC_rawInvMass_fEta_fBrem_AF = mumuSC_raw_fEta_fBrem_AF.M();
    double mumuSC_rawInvMass_fEta_fBrem_L = mumuSC_raw_fEta_fBrem_L.M();
    double mumuSC_rawInvMass_fEta_fBrem_fEtEta = mumuSC_raw_fEta_fBrem_fEtEta.M();
    double mumuSC_rawInvMass_fEta_fBrem_AF_fEtEta = mumuSC_raw_fEta_fBrem_AF_fEtEta.M();
    double mumuSC_rawInvMass_fEta_fBrem_L_fEtEta = mumuSC_raw_fEta_fBrem_L_fEtEta.M();
*/
    Mmumugamma = mumugammaInvMass;
//		if( ((Mmumugamma>97.2) && (Mmumugamma<110.0)) || ((Mmumugamma>70.0) && (Mmumugamma<87.2)) ) cout << "isVeryLooseMMG:isLooseMMG:isTightMMG " << isVeryLooseMMG << isLooseMMG << isTightMMG << "\t\t\tMmumugamma= " << Mmumugamma << endl;
//		cout << "ievt= " << iEvent << "\t\tVL:L:T:M " << isVeryLooseMMG << isLooseMMG << isTightMMG << isMultipleCandidate << "\t\t\tMmumugamma= " << Mmumugamma << endl;
//		cout << "ievt= " << iEvent << "\t\t" << isMMGCandidate << isAfterFSRCut1 << isAfterFSRCut2 << isAfterFSRCut3 << isAfterFSRCut4 << isVeryLooseMMG << isLooseMMG << isTightMMG << isMultipleCandidate << "\t\t\tMmumugamma= " << Mmumugamma << endl;
    Mmumugamma_SC = mumuSCInvMass;
    Mmumugamma_5x5 = mumu5x5InvMass;
    Mmumugamma_SCraw = mumuSC_rawInvMass;
/*
    Mmumugamma_SCraw_fEta = mumuSC_rawInvMass_fEta;
   
    Mmumugamma_SCraw_fEta_fBrem = mumuSC_rawInvMass_fEta_fBrem;
    Mmumugamma_SCraw_fEta_fBrem_AF = mumuSC_rawInvMass_fEta_fBrem_AF;
    Mmumugamma_SCraw_fEta_fBrem_L = mumuSC_rawInvMass_fEta_fBrem_L;
    Mmumugamma_SCraw_fEta_fBrem_fEtEta = mumuSC_rawInvMass_fEta_fBrem_fEtEta;
    Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta = mumuSC_rawInvMass_fEta_fBrem_AF_fEtEta;
    Mmumugamma_SCraw_fEta_fBrem_L_fEtEta = mumuSC_rawInvMass_fEta_fBrem_L_fEtEta;
*/
    mumugamma.Clear();
    mumuSC.Clear();
    mumu5x5.Clear();
    mumuSC_raw.Clear();
    PhotonEScale->Clear();
    PhotonSC->Clear();
    Photon5x5->Clear();
    PhotonSC_raw->Clear();
//    PhotonSC_raw_fEta->Clear();

    mmg_k = (double)(pow(91.1876,2) - pow(Mmumu,2) ) / (double)(pow(Mmumugamma,2) - pow(Mmumu,2));
    mmg_ik = (double)(pow(Mmumugamma,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_s = mmg_ik - 1.0;
    mmg_logk = log(mmg_k);
    mmg_logik = log(mmg_ik);
    mmg_logs = log(mmg_s);

    mmg_k_5x5 = (double)(pow(91.1876,2) - pow(Mmumu,2) ) / (double)(pow(Mmumugamma_5x5,2) - pow(Mmumu,2));
    mmg_ik_5x5 = (double)(pow(Mmumugamma_5x5,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_s_5x5 = mmg_ik_5x5 - 1.0;
    mmg_logk_5x5 = log(mmg_k_5x5);
    mmg_logik_5x5 = log(mmg_ik_5x5);
    mmg_logs_5x5 = log(mmg_s_5x5);

    mmg_k_SC = (double)(pow(91.1876,2) - pow(Mmumu,2) ) / (double)(pow(Mmumugamma_SC,2) - pow(Mmumu,2));
    mmg_ik_SC = (double)(pow(Mmumugamma_SC,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_s_SC = mmg_ik_SC - 1.0;
    mmg_logk_SC = log(mmg_k_SC);
    mmg_logik_SC = log(mmg_ik_SC);
    mmg_logs_SC = log(mmg_s_SC);

    mmg_k_SCraw = (double)(pow(91.1876,2) - pow(Mmumu,2) ) / (double)(pow(Mmumugamma_SCraw,2) - pow(Mmumu,2));
    mmg_ik_SCraw = (double)(pow(Mmumugamma_SCraw,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_s_SCraw = mmg_ik_SCraw - 1.0;
    mmg_logk_SCraw = log(mmg_k_SCraw);
    mmg_logik_SCraw = log(mmg_ik_SCraw);
    mmg_logs_SCraw = log(mmg_s_SCraw);

    mmg_k_SCraw_fEta = (double)(pow(91.1876,2) - pow(Mmumu,2) ) / (double)(pow(Mmumugamma_SCraw_fEta,2) - pow(Mmumu,2));
    mmg_ik_SCraw_fEta = (double)(pow(Mmumugamma_SCraw_fEta,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_s_SCraw_fEta = mmg_ik_SCraw_fEta - 1.0;
    mmg_logk_SCraw_fEta = log(mmg_k_SCraw_fEta);
    mmg_logik_SCraw_fEta = log(mmg_ik_SCraw_fEta);
    mmg_logs_SCraw_fEta = log(mmg_s_SCraw_fEta);


    mmg_ik_SCraw_fEta_fBrem = (double)(pow(Mmumugamma_SCraw_fEta_fBrem,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_ik_SCraw_fEta_fBrem_AF = (double)(pow(Mmumugamma_SCraw_fEta_fBrem_AF,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_ik_SCraw_fEta_fBrem_L = (double)(pow(Mmumugamma_SCraw_fEta_fBrem_L,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_ik_SCraw_fEta_fBrem_fEtEta = (double)(pow(Mmumugamma_SCraw_fEta_fBrem_fEtEta,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_ik_SCraw_fEta_fBrem_AF_fEtEta = (double)(pow(Mmumugamma_SCraw_fEta_fBrem_AF_fEtEta,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );
    mmg_ik_SCraw_fEta_fBrem_L_fEtEta = (double)(pow(Mmumugamma_SCraw_fEta_fBrem_L_fEtEta,2) - pow(Mmumu,2)) / (double)(pow(91.1876,2) - pow(Mmumu,2) );


//    cerr << "\t\tINFO: mumugamma invariant mass : Mmumugamma = " << mumugammaInvMass << endl;

    double phiPhoton = myphoton->Phi();
    double etaPhoton = myphoton->Eta();
    double phiMuon = mymuon1->Phi();
    double etaMuon = mymuon1->Eta();
    double phiMuon_oppositeCharge = mymuon2->Phi();
    double etaMuon_oppositeCharge = mymuon2->Eta();
    double deltaRPM = DeltaR(etaPhoton, phiPhoton, etaMuon, phiMuon);
    double deltaRPAM = DeltaR(etaPhoton, phiPhoton, etaMuon_oppositeCharge, phiMuon_oppositeCharge);

    double deltaRmin;
    TRootMuon *farMuon;
    TRootMuon *nearMuon;
    TRootMuon *minusMuon;
    TRootMuon *plusMuon;
    TLorentzVector *correctedfarMuon;
    TLorentzVector *correctednearMuon;
    TLorentzVector *correctedminusMuon;
    TLorentzVector *correctedplusMuon;

    if(deltaRPM < deltaRPAM){
      deltaRmin = deltaRPM;
      farMuon = (TRootMuon*) mymuon2;
      nearMuon = (TRootMuon*) mymuon1;
      correctedfarMuon = (TLorentzVector*) correctedmymuon2;
      correctednearMuon = (TLorentzVector*) correctedmymuon1;
    } else {
      deltaRmin = deltaRPAM;
      farMuon = (TRootMuon*) mymuon1;
      nearMuon = (TRootMuon*) mymuon2;
      correctedfarMuon = (TLorentzVector*) correctedmymuon1;
      correctednearMuon = (TLorentzVector*) correctedmymuon2;
    }
    if( mymuon1->charge()>0 ){
      plusMuon  = (TRootMuon*) mymuon1;
      minusMuon = (TRootMuon*) mymuon2;
      correctedplusMuon  = (TLorentzVector*) correctedmymuon1;
      correctedminusMuon = (TLorentzVector*) correctedmymuon2;
    } else {
      minusMuon = (TRootMuon*) mymuon1;
      plusMuon  = (TRootMuon*) mymuon2;
      correctedminusMuon = (TLorentzVector*) correctedmymuon1;
      correctedplusMuon  = (TLorentzVector*) correctedmymuon2;
    }

    //-----Implementation of s with MZ_Surface -----//


    //TLorentzVector * nearMuonPlusPhoton;
    //nearMuonPlusPhoton = correctednearMuon + myphoton;

		TLorentzVector MuonBeforeBremN;
		TLorentzVector MuonBeforeBremF;
		TLorentzVector MuonBeforeBremP;
		TLorentzVector MuonBeforeBremM;
		TLorentzVector MuonBeforeBremL;
		TLorentzVector MuonBeforeBremS;
		MuonBeforeBremN = (*correctednearMuon) + (*PhotonEScale);
		MuonBeforeBremF = (*correctedfarMuon);
		if( MuonBeforeBremF.Pt() > MuonBeforeBremN.Pt() )
		{
			MuonBeforeBremL = MuonBeforeBremF;
			MuonBeforeBremS = MuonBeforeBremN;
		} else {
			MuonBeforeBremL = MuonBeforeBremN;
			MuonBeforeBremS = MuonBeforeBremF;
		}
		if( nearMuon->charge() > 0 )
		{
			MuonBeforeBremM = MuonBeforeBremF;
			MuonBeforeBremP = MuonBeforeBremN;
		} else {
			MuonBeforeBremM = MuonBeforeBremN;
			MuonBeforeBremP = MuonBeforeBremF;
		}

		MuonBeforeBremM_Pt = MuonBeforeBremM.Pt();
		MuonBeforeBremP_Pt = MuonBeforeBremP.Pt();
		MuonBeforeBremN_Pt = MuonBeforeBremN.Pt();
		MuonBeforeBremF_Pt = MuonBeforeBremF.Pt();
		MuonBeforeBremL_Pt = MuonBeforeBremL.Pt();
		MuonBeforeBremS_Pt = MuonBeforeBremS.Pt();
		MuonBeforeBremM_Eta = MuonBeforeBremM.Eta();
		MuonBeforeBremP_Eta = MuonBeforeBremP.Eta();
		MuonBeforeBremN_Eta = MuonBeforeBremN.Eta();
		MuonBeforeBremF_Eta = MuonBeforeBremF.Eta();
		MuonBeforeBremL_Eta = MuonBeforeBremL.Eta();
		MuonBeforeBremS_Eta = MuonBeforeBremS.Eta();

		MuonBeforeBremM_Phi = MuonBeforeBremM.Phi();
		MuonBeforeBremP_Phi = MuonBeforeBremP.Phi();
		MuonBeforeBremN_Phi = MuonBeforeBremN.Phi();
		MuonBeforeBremF_Phi = MuonBeforeBremF.Phi();
		MuonBeforeBremL_Phi = MuonBeforeBremL.Phi();
		MuonBeforeBremS_Phi = MuonBeforeBremS.Phi();

		MuonBeforeBremM_E = MuonBeforeBremM.Energy();
		MuonBeforeBremP_E = MuonBeforeBremP.Energy();
		MuonBeforeBremN_E = MuonBeforeBremN.Energy();
		MuonBeforeBremF_E = MuonBeforeBremF.Energy();
		MuonBeforeBremL_E = MuonBeforeBremL.Energy();
		MuonBeforeBremS_E = MuonBeforeBremS.Energy();

		MuonBeforeBremM_Px = MuonBeforeBremM.Px();
		MuonBeforeBremP_Px = MuonBeforeBremP.Px();
		MuonBeforeBremN_Px = MuonBeforeBremN.Px();
		MuonBeforeBremF_Px = MuonBeforeBremF.Px();
		MuonBeforeBremL_Px = MuonBeforeBremL.Px();
		MuonBeforeBremS_Px = MuonBeforeBremS.Px();

		MuonBeforeBremM_Py = MuonBeforeBremM.Py();
		MuonBeforeBremP_Py = MuonBeforeBremP.Py();
		MuonBeforeBremN_Py = MuonBeforeBremN.Py();
		MuonBeforeBremF_Py = MuonBeforeBremF.Py();
		MuonBeforeBremL_Py = MuonBeforeBremL.Py();
		MuonBeforeBremS_Py = MuonBeforeBremS.Py();

		MuonBeforeBremM_Pz = MuonBeforeBremM.Pz();
		MuonBeforeBremP_Pz = MuonBeforeBremP.Pz();
		MuonBeforeBremN_Pz = MuonBeforeBremN.Pz();
		MuonBeforeBremF_Pz = MuonBeforeBremF.Pz();
		MuonBeforeBremL_Pz = MuonBeforeBremL.Pz();
		MuonBeforeBremS_Pz = MuonBeforeBremS.Pz();

		MuonBeforeBremF_Charge = farMuon->charge();
		MuonBeforeBremN_Charge = nearMuon->charge();
		MuonBeforeBremL_Charge = leadingMuon->charge();
		MuonBeforeBremS_Charge = subleadingMuon->charge();
		
	
//    double nearMuonPlusPhoton_Pt = sqrt(pow(correctednearMuon->Px() + myphoton->Px(),2) + pow(correctednearMuon->Py() + myphoton->Py(),2));
    double nearMuonPlusPhoton_Pt = MuonBeforeBremN_Pt;
 
    //double nearMuonPlusPhoton_Pt = nearMuonPlusPhoton->Pt();
//    double farMuon_Pt = correctedfarMuon->Pt();
    double farMuon_Pt = MuonBeforeBremF_Pt;

    double lead_Pt = MuonBeforeBremL_Pt;
    double trail_Pt = MuonBeforeBremS_Pt;

/*
    if(nearMuonPlusPhoton_Pt > farMuon_Pt) 
    {
	lead_Pt = nearMuonPlusPhoton_Pt;
	trail_Pt = farMuon_Pt;
    }
    if(nearMuonPlusPhoton_Pt < farMuon_Pt)
    {
        lead_Pt = farMuon_Pt;
        trail_Pt = nearMuonPlusPhoton_Pt;
    }
*/
    //double bin1 = lead_Pt / 2.0 + 2.0; //modify if we change the binning of the surface 
    //double bin2 = trail_Pt / 2.0 + 2.0;

    double bin1 = lead_Pt / (200.0 / binNumber) + 1.0; //ok if we use a surface with Pt between 0 and 200 GeV >> bin de 200.0 / binNumber GeV
    double bin2 = trail_Pt / (200.0 / binNumber) + 1.0;


    int bin1Int = (int) bin1;
    int bin2Int = (int) bin2;
  
    int multi = (bin1Int-1) * binNumber + bin2Int;
    string ligne;
    ifstream monFlux("Mmumu.txt");
    for(int k = 0; k <multi; k++)
    {
	getline(monFlux, ligne);
    }

    MZ_Surface = StringToDouble(ligne);

    monFlux.close();


    mmg_k_MZ_Surface = (double)(pow(MZ_Surface,2) - pow(Mmumu,2) ) / (double)(pow(Mmumugamma,2) - pow(Mmumu,2));
    mmg_ik_MZ_Surface = (double)(pow(Mmumugamma,2) - pow(Mmumu,2)) / (double)(pow(MZ_Surface,2) - pow(Mmumu,2) );
    mmg_s_MZ_Surface = mmg_ik_MZ_Surface - 1.0;
    mmg_logk_MZ_Surface = log(mmg_k_MZ_Surface);
    mmg_logik_MZ_Surface = log(mmg_ik_MZ_Surface);
    mmg_logs_MZ_Surface = log(mmg_s_MZ_Surface);



    //-----END of : Implementation of s with MZ_Surface -----//




    // FILLING MINITREE INFORMATION
    MuonM_Pt = correctedminusMuon->Pt();
    MuonM_Eta = minusMuon->Eta();
    MuonM_Phi = minusMuon->Phi();
    MuonM_E = minusMuon->Energy();
    MuonM_Px = minusMuon->Px();
    MuonM_Py = minusMuon->Py();
    MuonM_Pz = minusMuon->Pz();
    MuonM_isoR03_emEt = minusMuon->isoR03_emEt();
    MuonM_isoR03_hadEt = minusMuon->isoR03_hadEt();
    MuonM_isoR03_hoEt = minusMuon->isoR03_hoEt();
    MuonM_isoR03_nJets = minusMuon->isoR03_nJets();
    MuonM_isoR03_nTracks = minusMuon->isoR03_nTracks();
    MuonM_isoR03_sumPt = minusMuon->isoR03_sumPt();
    MuonM_isoR05_emEt = minusMuon->isoR05_emEt();
    MuonM_isoR05_hadEt = minusMuon->isoR05_hadEt();
    MuonM_isoR05_hoEt = minusMuon->isoR05_hoEt();
    MuonM_isoR05_nJets = minusMuon->isoR05_nJets();
    MuonM_isoR05_nTracks = minusMuon->isoR05_nTracks();
    MuonM_isoR05_sumPt = minusMuon->isoR05_sumPt();
    MuonP_Pt = correctedplusMuon->Pt();
    MuonP_Eta = plusMuon->Eta();
    MuonP_Phi = plusMuon->Phi();
    MuonP_E = plusMuon->Energy();
    MuonP_Px = plusMuon->Px();
    MuonP_Py = plusMuon->Py();
    MuonP_Pz = plusMuon->Pz();
    MuonP_isoR03_emEt = plusMuon->isoR03_emEt();
    MuonP_isoR03_hadEt = plusMuon->isoR03_hadEt();
    MuonP_isoR03_hoEt = plusMuon->isoR03_hoEt();
    MuonP_isoR03_nJets = plusMuon->isoR03_nJets();
    MuonP_isoR03_nTracks = plusMuon->isoR03_nTracks();
    MuonP_isoR03_sumPt = plusMuon->isoR03_sumPt();
    MuonP_isoR05_emEt = plusMuon->isoR05_emEt();
    MuonP_isoR05_hadEt = plusMuon->isoR05_hadEt();
    MuonP_isoR05_hoEt = plusMuon->isoR05_hoEt();
    MuonP_isoR05_nJets = plusMuon->isoR05_nJets();
    MuonP_isoR05_nTracks = plusMuon->isoR05_nTracks();
    MuonP_isoR05_sumPt = plusMuon->isoR05_sumPt();
    MuonF_Pt = correctedfarMuon->Pt();
    MuonF_Eta = farMuon->Eta();
    MuonF_Phi = farMuon->Phi();
    MuonF_E = farMuon->Energy();
    MuonF_Px = farMuon->Px();
    MuonF_Py = farMuon->Py();
    MuonF_Pz = farMuon->Pz();
    MuonF_Charge = farMuon->charge();
    MuonF_isoR03_emEt = farMuon->isoR03_emEt();
    MuonF_isoR03_hadEt = farMuon->isoR03_hadEt();
    MuonF_isoR03_hoEt = farMuon->isoR03_hoEt();
    MuonF_isoR03_nJets = farMuon->isoR03_nJets();
    MuonF_isoR03_nTracks = farMuon->isoR03_nTracks();
    MuonF_isoR03_sumPt = farMuon->isoR03_sumPt();
    MuonF_isoR05_emEt = farMuon->isoR05_emEt();
    MuonF_isoR05_hadEt = farMuon->isoR05_hadEt();
    MuonF_isoR05_hoEt = farMuon->isoR05_hoEt();
    MuonF_isoR05_nJets = farMuon->isoR05_nJets();
    MuonF_isoR05_nTracks = farMuon->isoR05_nTracks();
    MuonF_isoR05_sumPt = farMuon->isoR05_sumPt();
    MuonN_Pt = correctednearMuon->Pt();
    MuonN_Eta = nearMuon->Eta();
    MuonN_Phi = nearMuon->Phi();
    MuonN_E = nearMuon->Energy();
    MuonN_Px = nearMuon->Px();
    MuonN_Py = nearMuon->Py();
    MuonN_Pz = nearMuon->Pz();
    MuonN_Charge = nearMuon->charge();
    MuonN_isoR03_emEt = nearMuon->isoR03_emEt();
    MuonN_isoR03_hadEt = nearMuon->isoR03_hadEt();
    MuonN_isoR03_hoEt = nearMuon->isoR03_hoEt();
    MuonN_isoR03_nJets = nearMuon->isoR03_nJets();
    MuonN_isoR03_nTracks = nearMuon->isoR03_nTracks();
    MuonN_isoR03_sumPt = nearMuon->isoR03_sumPt();
    MuonN_isoR05_emEt = nearMuon->isoR05_emEt();
    MuonN_isoR05_hadEt = nearMuon->isoR05_hadEt();
    MuonN_isoR05_hoEt = nearMuon->isoR05_hoEt();
    MuonN_isoR05_nJets = nearMuon->isoR05_nJets();
    MuonN_isoR05_nTracks = nearMuon->isoR05_nTracks();
    MuonN_isoR05_sumPt = nearMuon->isoR05_sumPt();

    deltaRNear = DeltaR(etaPhoton, phiPhoton, nearMuon->Eta(), nearMuon->Phi());
    deltaRFar = DeltaR(etaPhoton, phiPhoton, farMuon->Eta(), farMuon->Phi());
    deltaRMinus = DeltaR(etaPhoton, phiPhoton, minusMuon->Eta(), minusMuon->Phi());
    deltaRPlus = DeltaR(etaPhoton, phiPhoton, plusMuon->Eta(), plusMuon->Phi());
    deltaRLeading = DeltaR(etaPhoton, phiPhoton, leadingMuon->Eta(), leadingMuon->Phi());
    deltaRSubleading = DeltaR(etaPhoton, phiPhoton, subleadingMuon->Eta(), subleadingMuon->Phi());

    if( doMC )
    { 
	//cout<<endl<<"doMC dans le if = "<<doMC<<endl; 
	// Compute Stuff, with MC truth information
      doGenInfo( (TRootParticle*) myphoton, mcParticles, &Photon_MC_E, &Photon_MC_Px, &Photon_MC_Py, &Photon_MC_Pz, &Photon_MC_Phi, &Photon_MC_Eta, &Photon_MC_Pt, 22 );
      doGenInfo( (TRootParticle*) correctedminusMuon, mcParticles, &MuonM_MC_E, &MuonM_MC_Px, &MuonM_MC_Py, &MuonM_MC_Pz, &MuonM_MC_Phi, &MuonM_MC_Eta, &MuonM_MC_Pt, 13 );
      doGenInfo( (TRootParticle*) correctedplusMuon, mcParticles, &MuonP_MC_E, &MuonP_MC_Px, &MuonP_MC_Py, &MuonP_MC_Pz, &MuonP_MC_Phi, &MuonP_MC_Eta, &MuonP_MC_Pt, -13 );
      doGenInfo( (TRootParticle*) correctednearMuon, mcParticles, &MuonN_MC_E, &MuonN_MC_Px, &MuonN_MC_Py, &MuonN_MC_Pz, &MuonN_MC_Phi, &MuonN_MC_Eta, &MuonN_MC_Pt, (-1)*(nearMuon->charge())*13 );
      doGenInfo( (TRootParticle*) correctedfarMuon, mcParticles, &MuonF_MC_E, &MuonF_MC_Px, &MuonF_MC_Py, &MuonF_MC_Pz, &MuonF_MC_Phi, &MuonF_MC_Eta, &MuonF_MC_Pt, (-1)*(farMuon->charge())*13 );
      doGenInfo( (TRootParticle*) correctedleadingMuon, mcParticles, &MuonL_MC_E, &MuonL_MC_Px, &MuonL_MC_Py, &MuonL_MC_Pz, &MuonL_MC_Phi, &MuonL_MC_Eta, &MuonL_MC_Pt, (-1)*(leadingMuon->charge())*13 );
      doGenInfo( (TRootParticle*) correctedsubleadingMuon, mcParticles, &MuonS_MC_E, &MuonS_MC_Px, &MuonS_MC_Py, &MuonS_MC_Pz, &MuonS_MC_Phi, &MuonS_MC_Eta, &MuonS_MC_Pt, (-1)*(subleadingMuon->charge())*13 );

			if(Photon_MC_E != 0.0) Photon_SC_rawE_x_fEta_o_MC_E = (double)(Photon_SC_rawE_x_fEta) / (double)(Photon_MC_E);
			if(Photon_MC_E != 0.0) Photon_E_o_MC_E = (double)(Photon_E) / (double)(Photon_MC_E);
			
			if(Photon_MC_E != 0.0) Photon_SC_rawE_x_fEta_x_fBrem_o_MC_E  = (double)(Photon_SC_rawE_x_fEta_x_fBrem) / (double)(Photon_MC_E);
			if(Photon_MC_E != 0.0) Photon_SC_rawE_x_fEta_x_fBrem_AF_o_MC_E  = (double)(Photon_SC_rawE_x_fEta_x_fBrem_AF) / (double)(Photon_MC_E);
			if(Photon_MC_E != 0.0) Photon_SC_rawE_x_fEta_x_fBrem_L_o_MC_E  = (double)(Photon_SC_rawE_x_fEta_x_fBrem_L) / (double)(Photon_MC_E);

			if(Photon_MC_E != 0.0) Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta_o_MC_E  = (double)(Photon_SC_rawE_x_fEta_x_fBrem_x_fEtEta) / (double)(Photon_MC_E);
			if(Photon_MC_E != 0.0) Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta_o_MC_E  = (double)(Photon_SC_rawE_x_fEta_x_fBrem_AF_x_fEtEta) / (double)(Photon_MC_E);
			if(Photon_MC_E != 0.0) Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta_o_MC_E  = (double)(Photon_SC_rawE_x_fEta_x_fBrem_L_x_fEtEta) / (double)(Photon_MC_E);

	//cout<<endl<<"Photon_E_o_MC_E = "<<Photon_E_o_MC_E<<endl;

      TLorentzVector mumu_Photon_MC;
      TLorentzVector mumugamma_Photon_MC;
      TLorentzVector mumu_Muons_MC;
      TLorentzVector mumugamma_Muons_MC;
      TLorentzVector mumu_MMG_MC;
      TLorentzVector mumugamma_MMG_MC;

      TLorentzVector *PhotonMC = new TLorentzVector( Photon_MC_Px, Photon_MC_Py, Photon_MC_Pz, Photon_MC_E);
      TLorentzVector *MuonLMC = new TLorentzVector( MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_E);
      TLorentzVector *MuonSMC = new TLorentzVector( MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_E);

      TLorentzVector *PhotonEScale = new TLorentzVector( EScale*(myphoton->Px()), EScale*(myphoton->Py()), EScale*(myphoton->Pz()), EScale*(myphoton->Energy()));
      TLorentzVector *PhotonSC = new TLorentzVector( myphoton->Px(), myphoton->Py(), myphoton->Pz(), EScale*(myphoton->superCluster()->Mag()));
      TLorentzVector *Photon5x5 = new TLorentzVector( myphoton->Px(), myphoton->Py(), myphoton->Pz(), EScale*(myphoton->e5x5()));
      TLorentzVector *PhotonSC_raw = new TLorentzVector( myphoton->Px(), myphoton->Py(), myphoton->Pz(), EScale*(myphoton->superCluster()->rawEnergy()));

      mumu_Photon_MC = (*correctedleadingMuon) + (*correctedsubleadingMuon);
      mumugamma_Photon_MC = (*correctedleadingMuon) + (*correctedsubleadingMuon) + (*PhotonMC);
      mumu_Muons_MC = (*MuonLMC) + (*MuonSMC);
      mumugamma_Muons_MC = (*MuonLMC) + (*MuonSMC) + (*PhotonEScale);
      mumu_MMG_MC = (*MuonLMC) + (*MuonSMC);
      mumugamma_MMG_MC = (*MuonLMC) + (*MuonSMC) + (*PhotonMC);
      Mmumu_Photon_MC = mumu_Photon_MC.M();
      Mmumugamma_Photon_MC = mumugamma_Photon_MC.M();
      Mmumu_Muons_MC = mumu_Muons_MC.M();
      Mmumugamma_Muons_MC = mumugamma_Muons_MC.M();
      Mmumu_MMG_MC = mumu_MMG_MC.M();
      Mmumugamma_MMG_MC = mumugamma_MMG_MC.M();

      mmg_k_Photon_MC = (double)(pow(91.1876,2) - pow(Mmumu_Photon_MC,2) ) / (double)(pow(Mmumugamma_Photon_MC,2) - pow(Mmumu_Photon_MC,2));
      mmg_ik_Photon_MC = (double)(pow(Mmumugamma_Photon_MC,2) - pow(Mmumu_Photon_MC,2)) / (double)(pow(91.1876,2) - pow(Mmumu_Photon_MC,2) );
      mmg_s_Photon_MC = mmg_ik_Photon_MC -1.0;
      mmg_logk_Photon_MC = log(mmg_k_Photon_MC);
      mmg_logik_Photon_MC = log(mmg_ik_Photon_MC);
      mmg_logs_Photon_MC = log(mmg_s_Photon_MC);


      mmg_k_Muons_MC = (double)(pow(91.1876,2) - pow(Mmumu_Muons_MC,2) ) / (double)(pow(Mmumugamma_Muons_MC,2) - pow(Mmumu_Muons_MC,2));
      mmg_ik_Muons_MC = (double)(pow(Mmumugamma_Muons_MC,2) - pow(Mmumu_Muons_MC,2)) / (double)(pow(91.1876,2) - pow(Mmumu_Muons_MC,2) );
      mmg_s_Muons_MC = mmg_ik_Muons_MC -1.0;
      mmg_logk_Muons_MC = log(mmg_k_Muons_MC);
      mmg_logik_Muons_MC = log(mmg_ik_Muons_MC);
      mmg_logs_Muons_MC = log(mmg_s_Muons_MC);

      mmg_k_MMG_MC = (double)(pow(91.1876,2) - pow(Mmumu_MMG_MC,2) ) / (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu_MMG_MC,2));
      mmg_ik_MMG_MC = (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu_MMG_MC,2)) / (double)(pow(91.1876,2) - pow(Mmumu_MMG_MC,2) );
      mmg_s_MMG_MC = mmg_ik_MMG_MC -1.0;
      mmg_logk_MMG_MC = log(mmg_k_MMG_MC);
      mmg_logik_MMG_MC = log(mmg_ik_MMG_MC);
      mmg_logs_MMG_MC = log(mmg_s_MMG_MC);

	
      mmg_k_MZ = (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu,2) ) / (double)(pow(Mmumugamma,2) - pow(Mmumu,2));
      mmg_ik_MZ = (double)(pow(Mmumugamma,2) - pow(Mmumu,2)) / (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu,2) );
      mmg_s_MZ = mmg_ik_MZ -1.0;
      mmg_logk_MZ = log(mmg_k_MZ);
      mmg_logik_MZ = log(mmg_ik_MZ);
      mmg_logs_MZ = log(mmg_s_MZ);     

 
      mmg_k_MZ_Photon_MC = (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu_Photon_MC,2) ) / (double)(pow(Mmumugamma_Photon_MC,2) - pow(Mmumu_Photon_MC,2));
      mmg_ik_MZ_Photon_MC = (double)(pow(Mmumugamma_Photon_MC,2) - pow(Mmumu_Photon_MC,2)) / (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu_Photon_MC,2) );
      mmg_s_MZ_Photon_MC = mmg_ik_MZ_Photon_MC -1.0;
      mmg_logk_MZ_Photon_MC = log(mmg_k_MZ_Photon_MC);
      mmg_logik_MZ_Photon_MC = log(mmg_ik_MZ_Photon_MC);
      mmg_logs_MZ_Photon_MC = log(mmg_s_MZ_Photon_MC);
      
      mmg_k_MZ_Muons_MC = (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu_Muons_MC,2) ) / (double)(pow(Mmumugamma_Muons_MC,2) - pow(Mmumu_Muons_MC,2));
      mmg_ik_MZ_Muons_MC = (double)(pow(Mmumugamma_Muons_MC,2) - pow(Mmumu_Muons_MC,2)) / (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu_Muons_MC,2) );
      mmg_s_MZ_Muons_MC = mmg_ik_MZ_Muons_MC -1.0;
      mmg_logk_MZ_Muons_MC = log(mmg_k_MZ_Muons_MC);
      mmg_logik_MZ_Muons_MC = log(mmg_ik_MZ_Muons_MC);
      mmg_logs_MZ_Muons_MC = log(mmg_s_MZ_Muons_MC);

      mmg_k_MZ_Muons_RECO_MC = (double)(pow(Mmumugamma_MMG_MC,2) - pow(Mmumu_Muons_MC,2) ) / (double)(pow(Mmumugamma,2) - pow(Mmumu,2));
      mmg_ik_MZ_Muons_RECO_MC = (double)(pow(Mmumugamma_Muons_MC,2) - pow(Mmumu_Muons_MC,2)) / (double)(pow(Mmumugamma,2) - pow(Mmumu,2) );
      mmg_s_MZ_Muons_RECO_MC = mmg_ik_MZ_Muons_RECO_MC -1.0;
      mmg_logk_MZ_Muons_RECO_MC = log(mmg_k_MZ_Muons_RECO_MC);
      mmg_logik_MZ_Muons_RECO_MC = log(mmg_ik_MZ_Muons_RECO_MC);
      mmg_logs_MZ_Muons_RECO_MC = log(mmg_s_MZ_Muons_RECO_MC); 	


//      Mmumu_Photon_MC = Mmumugamma_Photon_MC = mmg_k_Photon_MC = mmg_ik_Photon_MC = mmg_s_Photon_MC = mmg_logk_Photon_MC = mmg_logik_Photon_MC = mmg_logs_Photon_MC = -99.0;
//      Mmumu_Muons_MC = Mmumugamma_Muons_MC = mmg_k_Muons_MC = mmg_ik_Muons_MC = mmg_s_Muons_MC = mmg_logk_Muons_MC = mmg_logik_Muons_MC = mmg_logs_Muons_MC = -99.0;
//      Mmumu_MMG_MC = Mmumugamma_MMG_MC = mmg_k_MMG_MC = mmg_ik_MMG_MC = mmg_s_MMG_MC = mmg_logk_MMG_MC = mmg_logik_MMG_MC = mmg_logs_MMG_MC = -99.0;
    PhotonMC->Clear();
    MuonLMC->Clear();
    MuonSMC->Clear();
    PhotonEScale->Clear();
    PhotonSC->Clear();
    Photon5x5->Clear();
    PhotonSC_raw->Clear();
    
			delete PhotonMC;
			PhotonMC = 0;
			delete MuonLMC;
			MuonLMC = 0;
			delete MuonSMC;
			MuonSMC = 0;

			delete PhotonEScale;
			PhotonEScale = 0;
			delete PhotonSC;
			PhotonSC = 0;  
			delete Photon5x5;
			Photon5x5 = 0;
			delete PhotonSC_raw;
			PhotonSC_raw = 0;

    }// end doMC

		delete PhotonEScale;
		PhotonEScale = 0;
		delete PhotonSC;
		PhotonSC = 0;
		delete Photon5x5;
		Photon5x5 = 0;
		delete PhotonSC_raw;
		PhotonSC_raw = 0;



return 0;
}


