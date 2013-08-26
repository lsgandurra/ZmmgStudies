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
#include "rochcor2012v2.h"


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
			//double weight[50] = {0, 0.526829, 0.0199659, 0.041029, 0.297575, 0.945121, 1.25478, 1.37131, 1.33925, 1.21619, 1.31131, 1.3953, 1.58576, 1.44483, 1.39461, 1.35017, 1.31175, 1.39853, 0.846608, 0.656049, 0.373876, 0.188962, 0.0829103, 0.0353935, 0.0118849, 0.00415899, 0.00152407, 0.000559973, 0.000338563, 8.18194e-05, 2.39062e-05, 8.35835e-06, 6.12925e-06, 6.12198e-06, 3.45724e-05, 1.58819e-05, 6.8505e-05, 4.91191e-05, 0.00286974, 0, 0.000113381, 0, 0, 0, 0, 0, 0, 0, 0, 0
			double weight[60] = {0, 1.29035, 0.969051, 1.10542, 1.01565, 1.07791, 1.0295, 1.1186, 1.17847, 1.15068, 1.29173, 1.36484, 1.52146, 1.32135, 1.18026, 1.0568, 0.984744, 1.07815, 0.735267, 0.713699, 0.572953, 0.459397, 0.359919, 0.305021, 0.220799, 0.17518, 0.148288, 0.124672, 0.168702, 0.0893422, 0.0562933, 0.0405267, 0.0493557, 0.0452175, 0.127515, 0.0228427, 0.0374335, 0.0102364, 0.229358, 0, 0.00158044, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			};
			return weight[nGenVertices];
		} else if ( pu_set == "PU_S4" )
		{
			double weight[36] = {0, 0.0251666, 0.203404, 0.426814, 0.735332, 1.01419, 1.22149, 1.34305, 1.39792, 1.39684, 1.39017, 1.37192, 1.37141, 1.37803, 1.37017, 1.37888, 1.38338, 1.3804, 1.34827, 1.35904, 1.32289, 1.29428, 1.2253, 1.20232, 1.14617, 1.099, 1.01493, 0.929284, 0.912994, 0.790041, 0.69409, 0.699025, 0.644141, 0.654322, 0.524752, 0.888945};
			return weight[nGenVertices];
		}
	} else if(lumi_set == "2012")
        {
                if( pu_set == "PU_S10" )
                {
                        double weight[60] = {0, 10.4923, 24.4296, 37.9208, 31.8, 40.0539, 28.9985, 6.42726, 3.29256, 2.72439, 2.74116, 2.76669, 2.64701, 2.37805, 2.0256, 1.67002, 1.4141, 1.2635, 1.19499, 1.175, 1.16815, 1.15999, 1.12558, 1.06777, 1.00448, 0.933594, 0.861446, 0.801237, 0.739057, 0.683081, 0.626923, 0.578676, 0.533231, 0.489387, 0.449333, 0.413669, 0.379853, 0.347238, 0.31829, 0.290759, 0.266329, 0.243079, 0.221235, 0.203534, 0.185281, 0.16678, 0.156065, 0.137128, 0.128382, 0.118529, 0.106961, 0.0966754, 0.0900311, 0.0779896, 0.0726763, 0.0621801, 0.0612443, 0.0643741, 0.0467334, 0.047643};//minBias = 69400 mb
                        //double weight[60] = {0, 6.97794, 15.9626, 25.1337, 21.5631, 27.8706, 20.7291, 4.72141, 2.48527, 2.11211, 2.18107, 2.25708, 2.21129, 2.03136, 1.76656, 1.48469, 1.27975, 1.16261, 1.11701, 1.11517, 1.12553, 1.13497, 1.11901, 1.07958, 1.03399, 0.979687, 0.922786, 0.877374, 0.828428, 0.784851, 0.739304, 0.701227, 0.664716, 0.628224, 0.594535, 0.56465, 0.535298, 0.505554, 0.479074, 0.452693, 0.429154, 0.405581, 0.382396, 0.364585, 0.344078, 0.321201, 0.311801, 0.284283, 0.276236, 0.26475, 0.248048, 0.232794, 0.225124, 0.202509, 0.195954, 0.174067, 0.177969, 0.194123, 0.146186, 0.154512}; //minBias = 73500 mb
                        return weight[nGenVertices];

                }
        }

        double weight[60] = {0, 10.4923, 24.4296, 37.9208, 31.8, 40.0539, 28.9985, 6.42726, 3.29256, 2.72439, 2.74116, 2.76669, 2.64701, 2.37805, 2.0256, 1.67002, 1.4141, 1.2635, 1.19499, 1.175, 1.16815, 1.15999, 1.12558, 1.06777, 1.00448, 0.933594, 0.861446, 0.801237, 0.739057, 0.683081, 0.626923, 0.578676, 0.533231, 0.489387, 0.449333, 0.413669, 0.379853, 0.347238, 0.31829, 0.290759, 0.266329, 0.243079, 0.221235, 0.203534, 0.185281, 0.16678, 0.156065, 0.137128, 0.128382, 0.118529, 0.106961, 0.0966754, 0.0900311, 0.0779896, 0.0726763, 0.0621801, 0.0612443, 0.0643741, 0.0467334, 0.047643};//minBias = 69400 mb        

        return weight[nGenVertices];

}


double weight_TTJets(int nGenVertices, string lumi_set, string pu_set)
{
	//double weight[60] = {0, 1.04204, 0.953355, 1.18724, 1.08439, 1.07778, 1.03545, 1.05117, 1.17321, 1.27515, 1.25051, 1.3289, 1.62975, 1.43617, 1.16916, 1.07712, 1.00101, 1.16663, 0.659773, 0.62577, 0.556534, 0.480728, 0.30211, 0.288529, 0.221898, 0.184643, 0.161076, 0.145077, 0.122756, 0.0903329, 0.0414645, 0.0415935, 0.0427404, 0.0941495, 0.0595874, 0.0297779, 0, 0.00829203, 0, 0, 0.00155194, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	//return weight[nGenVertices];

	if(lumi_set == "2012" && pu_set == "PU_S7")
        {
                double weight[60] = {0, 0.323701, 0, 9.35923, 16.8744, 4.13585, 7.21989, 9.76902, 11.9687, 11.3925, 13.9582, 11.3853, 11.3239, 7.78842, 6.62434, 5.79726, 4.69908, 3.55859, 2.6404, 2.11875, 1.69048, 1.43388, 1.22031, 1.08963, 0.992514, 0.922732, 0.794686, 0.711696, 0.612313, 0.560809, 0.471973, 0.40871, 0.358024, 0.31271, 0.267614, 0.232174, 0.209486, 0.1874, 0.172988, 0.157301, 0.14804, 0.134552, 0.129954, 0.125054, 0.128653, 0.130938, 0.149274, 0.165733, 0.184317, 0.310683, 0.263024, 0.558884, 0.605857, 0.457497, 0, 0, 0, 0, 0, 0}; //minBias = 69400 mb
                //double weight[60] = {0, 0.215278, 0, 6.20324, 11.4423, 2.87784, 5.16102, 7.17624, 9.03417, 8.83212, 11.1062, 9.28817, 9.45992, 6.65298, 5.77719, 5.15393, 4.25264, 3.27444, 2.46809, 2.01086, 1.62881, 1.40294, 1.21319, 1.10168, 1.02168, 0.968288, 0.851272, 0.779324, 0.686358, 0.644362, 0.556578, 0.495267, 0.446306, 0.401425, 0.354094, 0.316913, 0.295214, 0.272841, 0.260372, 0.244908, 0.238548, 0.224502, 0.22462, 0.224006, 0.238916, 0.252173, 0.298233, 0.343583, 0.396589, 0.693947, 0.609966, 1.34579, 1.51496, 1.18794, 0, 0, 0, 0, 0, 0}; //minBias = 73500 mb
                return weight[nGenVertices];

        }

        double weight[60] = {0, 0.323701, 0, 9.35923, 16.8744, 4.13585, 7.21989, 9.76902, 11.9687, 11.3925, 13.9582, 11.3853, 11.3239, 7.78842, 6.62434, 5.79726, 4.69908, 3.55859, 2.6404, 2.11875, 1.69048, 1.43388, 1.22031, 1.08963, 0.992514, 0.922732, 0.794686, 0.711696, 0.612313, 0.560809, 0.471973, 0.40871, 0.358024, 0.31271, 0.267614, 0.232174, 0.209486, 0.1874, 0.172988, 0.157301, 0.14804, 0.134552, 0.129954, 0.125054, 0.128653, 0.130938, 0.149274, 0.165733, 0.184317, 0.310683, 0.263024, 0.558884, 0.605857, 0.457497, 0, 0, 0, 0, 0, 0}; //minBias = 69400 mb
        return weight[nGenVertices];

}

double weight_WJetsToLNu(int nGenVertices, string lumi_set, string pu_set)
{
  	//double weight[60] = {0, 1.81661, 1.06742, 1.24536, 1.07591, 1.14997, 1.11121, 1.14298, 1.29881, 1.21755, 1.27235, 1.2782, 1.56192, 1.29055, 1.10867, 1.13326, 0.976632, 1.08336, 0.693502, 0.682135, 0.530099, 0.425931, 0.354297, 0.272527, 0.192103, 0.147848, 0.123977, 0.115278, 0.132297, 0.0819375, 0.0464697, 0.03329, 0.0344375, 0.0323669, 0.0941032, 0.0227669, 0.0548884, 0.00644365, 0, 0, 0.000587713, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	//return weight[nGenVertices];

	if(lumi_set == "2012" && pu_set == "PU_S10")
        {
                double weight[60] = {0, 0, 0, 0, 0, 0, 0, 0, 1.6478, 1.81253, 1.79692, 1.40752, 2.5692, 1.89001, 1.39015, 1.09703, 0.969595, 1.04985, 1.14386, 0.930022, 0.91867, 1.15763, 0.993279, 1.31673, 1.13406, 0.881241, 0.914166, 0.827539, 0.802942, 0.669971, 0.980831, 0.475241, 0.604143, 0.772618, 0.903048, 0.802019, 0.545218, 0.624473, 0.549896, 0.339962, 0.518853, 0.325879, 0.291156, 1.07095, 0.389309, 0.279784, 0.397604, 0, 0.194203, 0, 0.0908993, 0.0612415, 0.0204273, 0, 0, 0, 0, 0, 0, 0.0018708}; //minBias = 69400 mb
                //double weight[60] = {0, 0, 0, 0, 0, 0, 0, 0, 1.24378, 1.40518, 1.42976, 1.14826, 2.14629, 1.61447, 1.21237, 0.975294, 0.877477, 0.96602, 1.06922, 0.882665, 0.885155, 1.13266, 0.987485, 1.3313, 1.16738, 0.924749, 0.979259, 0.906174, 0.900038, 0.769788, 1.15665, 0.575887, 0.753114, 0.991807, 1.19487, 1.09474, 0.768335, 0.909189, 0.827675, 0.5293, 0.836065, 0.543734, 0.503251, 1.91837, 0.722971, 0.538836, 0.79437, 0, 0.41786, 0, 0.2108, 0.147469, 0.0510788, 0, 0, 0, 0, 0, 0, 0.00606726}; //minBias = 73500 mb
                return weight[nGenVertices];

        }

        double weight[60] = {0, 0, 0, 0, 0, 0, 0, 0, 1.6478, 1.81253, 1.79692, 1.40752, 2.5692, 1.89001, 1.39015, 1.09703, 0.969595, 1.04985, 1.14386, 0.930022, 0.91867, 1.15763, 0.993279, 1.31673, 1.13406, 0.881241, 0.914166, 0.827539, 0.802942, 0.669971, 0.980831, 0.475241, 0.604143, 0.772618, 0.903048, 0.802019, 0.545218, 0.624473, 0.549896, 0.339962, 0.518853, 0.325879, 0.291156, 1.07095, 0.389309, 0.279784, 0.397604, 0, 0.194203, 0, 0.0908993, 0.0612415, 0.0204273, 0, 0, 0, 0, 0, 0, 0.0018708}; //minBias = 69400 mb
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

  if (myphoton->isEB()){
//    float cetacorr = fEta(myphoton->superCluster()->rawEnergy(), myphoton->superCluster()->Eta(), 0)/myphoton->superCluster()->rawEnergy();
		float cetacorr = f_eta;

    energy = myphoton->superCluster()->rawEnergy()*cetacorr; //previously in CMSSW
  }
  else if (myphoton->isEE()){
    energy = myphoton->superCluster()->rawEnergy()+myphoton->preshowerEnergy();
  }

  float newEnergy = energy;

  if (mode==0){ //Electron

    corr = ETHZ_fBremEta(myphoton->superCluster()->phiWidth()/myphoton->superCluster()->etaWidth(), myphoton->superCluster()->Eta(), 0);

    float et = energy*TMath::Sin(2*TMath::ATan(TMath::Exp(-myphoton->superCluster()->Eta())))/corr;

    if (myphoton->isEB()) corr2 = corr * ETHZ_fEt(et, 0);
    if (myphoton->isEE()) corr2 = corr * ETHZ_fEnergy(energy/corr, 1);

    newEnergy = energy/corr2; 

  }

  if (mode==1){ //low R9 Photons

    corr = ETHZ_fBremEta(myphoton->superCluster()->phiWidth()/myphoton->superCluster()->etaWidth(), myphoton->superCluster()->Eta(), 1);

    float et = energy*TMath::Sin(2*TMath::ATan(TMath::Exp(-myphoton->superCluster()->Eta())))/corr;

    if (myphoton->isEB()) corr2 = corr * ETHZ_fEt(et, 2);
    if (myphoton->isEE()) corr2 = corr * ETHZ_fEnergy(energy/corr, 3);

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
		if( ((myphoton->isEB()) && (myphoton->r9()>0.94)) || ((myphoton->isEE()) && (myphoton->r9()>0.95)) ) return 1.0;
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
		//if( verbositybis > 1) cout << "myphoton->joshEnergyRegression() = " << myphoton->joshEnergyRegression() << endl;
		//return myphoton->joshEnergyRegression() / (double)(myphoton->Energy());
	
		///// NEW VARIABLE IN TOTO FOR PHOTON REGRESSION /////

		if( verbositybis > 1) cout << "myphoton->joshEnergyRegression() = " << myphoton->joshEnergyRegression() << endl;
		return myphoton->joshEnergyRegression() / (double)(myphoton->Energy());

	}
	 else {
		vector<double> param_fbrem;
		vector<double> param_feteta;
		parameters_fbrem(param_fbrem, correctionSet, myphoton->isEB());
		parameters_feteta(param_feteta, correctionSet, myphoton->isEB());
		double brem = (double)(myphoton->superCluster()->phiWidth()) / (double)(myphoton->superCluster()->etaWidth());
		double f_brem = BremCor(param_fbrem, brem);
		double sc_e = (myphoton->isEB()==1) ? (f_eta * myphoton->superCluster()->rawEnergy()) : (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy());
		double sc_e_noCrack = (myphoton->isEB()==1) ? (f_eta * myphoton->superCluster()->rawEnergy()) : (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy());
		double f_crack = myphoton->superCluster()->crackCorrectionEtaPhi();
		double sc_et = sc_e * (sin(myphoton->superCluster()->Theta()));
	  double sc_et_noCrack = sc_e_noCrack * (sin(myphoton->superCluster()->Theta()));
		double f_et_eta = EtEtaCor(param_feteta, f_brem * sc_et, myphoton->superCluster()->Eta(), myphoton->isEB());
	  double f_et_eta_noCrack = EtEtaCor(param_feteta, f_brem * sc_et_noCrack, myphoton->superCluster()->Eta(), myphoton->isEB());
	  if( verbositybis > 1)
		{
			cout << "###\tmyphoton->superCluster()->crackCorrectionEta()= " << myphoton->superCluster()->crackCorrectionEta() << endl;
    	cout << "###\tmyphoton->superCluster()->crackCorrectionPhi()= " << myphoton->superCluster()->crackCorrectionPhi() << endl;
    	cout << "###\tmyphoton->superCluster()->crackCorrectionEtaPhi()= " << myphoton->superCluster()->crackCorrectionEtaPhi() << endl;
		}
		if( (myphoton->isEB()) && (myphoton->r9()<0.94) )
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
		if( (myphoton->isEB()) && (myphoton->r9()>0.94) )
		{
			if( verbositybis > 1) cout << "f_eta * myphoton->e5x5()= " << f_eta * myphoton->e5x5() << endl;
			return (double)(f_eta * myphoton->e5x5()) / (double)(myphoton->Energy());
		}
		if( (myphoton->isEE()) && (myphoton->r9()<0.95) )
		{
			if( verbositybis > 1) cout << "f_et_eta * f_brem * (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy())= " << f_et_eta * f_brem * (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy()) << endl;
			return (double)(f_et_eta * f_brem * (myphoton->superCluster()->rawEnergy() + myphoton->preshowerEnergy())) / (double)(myphoton->Energy());
		}
		if( (myphoton->isEE()) && (myphoton->r9()>0.95) )
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
/*
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
*/
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
  extern vector<string> hltnames;
  extern Int_t nVertices;
  extern Int_t nGenVertices;
  extern Float_t weight_pileUp, weight_Xsection;
  extern Int_t isMM, isMM_nonFSR;

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
  // mumu information
  // ____________________________________________

  extern Float_t Mmumu;
  extern Float_t Ptmumu;


  // ____________________________________________
  // MC Truth
  // ___________________________________________

  extern Float_t MuonM_MC_E, MuonM_MC_Px, MuonM_MC_Py, MuonM_MC_Pz, MuonM_MC_Phi, MuonM_MC_Eta, MuonM_MC_Pt;
  extern Float_t MuonP_MC_E, MuonP_MC_Px, MuonP_MC_Py, MuonP_MC_Pz, MuonP_MC_Phi, MuonP_MC_Eta, MuonP_MC_Pt;
  extern Float_t MuonN_MC_E, MuonN_MC_Px, MuonN_MC_Py, MuonN_MC_Pz, MuonN_MC_Phi, MuonN_MC_Eta, MuonN_MC_Pt;
  extern Float_t MuonF_MC_E, MuonF_MC_Px, MuonF_MC_Py, MuonF_MC_Pz, MuonF_MC_Phi, MuonF_MC_Eta, MuonF_MC_Pt;
  extern Float_t MuonL_MC_E, MuonL_MC_Px, MuonL_MC_Py, MuonL_MC_Pz, MuonL_MC_Phi, MuonL_MC_Eta, MuonL_MC_Pt;
  extern Float_t MuonS_MC_E, MuonS_MC_Px, MuonS_MC_Py, MuonS_MC_Pz, MuonS_MC_Phi, MuonS_MC_Eta, MuonS_MC_Pt;
  extern Float_t Mmumu_Muons_MC;


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

int FillMM(TRootMuon* mymuon1, TRootMuon* mymuon2, TLorentzVector* correctedMuon1, TLorentzVector* correctedMuon2, bool doMC, bool doR9Rescaling, TClonesArray* mcParticles)
{

      // Fill muons stuff
      TRootMuon *leadingMuon;
      TRootMuon *subleadingMuon;
      TLorentzVector *correctedleadingMuon;
      TLorentzVector *correctedsubleadingMuon;
      if( (correctedMuon1->Pt()) > (correctedMuon2->Pt()) )      {
        leadingMuon = mymuon1;
        correctedleadingMuon = correctedMuon1;
        subleadingMuon = mymuon2;
        correctedsubleadingMuon = correctedMuon2;
      } else {
        leadingMuon = mymuon2;
        correctedleadingMuon = correctedMuon2;
        subleadingMuon = mymuon1;
        correctedsubleadingMuon = correctedMuon1;
      }
      MuonL_Pt = correctedleadingMuon->Pt();
      cerr<<endl<<" in FILLMMG : correctedleadingMuon->Pt() = "<<correctedleadingMuon->Pt()<<endl;
      MuonL_Eta = leadingMuon->Eta();
      MuonL_Phi = leadingMuon->Phi();
      MuonL_E = correctedleadingMuon->Energy();
      MuonL_Px = correctedleadingMuon->Px();
      MuonL_Py = correctedleadingMuon->Py();
      MuonL_Pz = correctedleadingMuon->Pz();
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
      MuonS_E = correctedsubleadingMuon->Energy();
      MuonS_Px = correctedsubleadingMuon->Px();
      MuonS_Py = correctedsubleadingMuon->Py();
      MuonS_Pz = correctedsubleadingMuon->Pz();
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



    TRootMuon *minusMuon;
    TRootMuon *plusMuon;
    TLorentzVector *correctedminusMuon;
    TLorentzVector *correctedplusMuon;

    if( mymuon1->charge()>0 ){

	plusMuon  = (TRootMuon*) mymuon1;
	minusMuon = (TRootMuon*) mymuon2;
	correctedplusMuon  = (TLorentzVector*) correctedMuon1;
	correctedminusMuon = (TLorentzVector*) correctedMuon2;
    }
    else
    {
	minusMuon = (TRootMuon*) mymuon1;
	plusMuon  = (TRootMuon*) mymuon2;
	correctedminusMuon = (TLorentzVector*) correctedMuon1;
	correctedplusMuon  = (TLorentzVector*) correctedMuon2;

    }

    //-----END of : Implementation of s with MZ_Surface -----//


    // FILLING MINITREE INFORMATION
    MuonM_Pt = correctedminusMuon->Pt();
    MuonM_Eta = minusMuon->Eta();
    MuonM_Phi = minusMuon->Phi();
    MuonM_E = correctedminusMuon->Energy();
    MuonM_Px = correctedminusMuon->Px();
    MuonM_Py = correctedminusMuon->Py();
    MuonM_Pz = correctedminusMuon->Pz();
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
    MuonP_E = correctedplusMuon->Energy();
    MuonP_Px = correctedplusMuon->Px();
    MuonP_Py = correctedplusMuon->Py();
    MuonP_Pz = correctedplusMuon->Pz();
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
   
        //MuonLMC->Clear();
        //MuonSMC->Clear();
	delete MuonLMC;
	MuonLMC = 0;
	delete MuonSMC;
	MuonSMC = 0;


    }// end doMC

	return 0;
}
