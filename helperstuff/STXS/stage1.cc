// stage1_reco:
// ggH_0J =0
// ggH_1J_PTH_0_60 =1
// ggH_1J_PTH_60_120 =2
// ggH_1J_PTH_120_200=3
// ggH_1J_PTH_200=4
// ggH_2J_PTH_0_60 =5
// ggH_2J_PTH_60_120 =6
// ggH_2J_PTH_120_200=7
// ggH_2J_PTH_200=8
// VBF_VBF = 9
// VBF_PTJET1_GT200_1J = 10
// VH_lep_0_150 =11
// VH_Had = 12
// ttH_Lep = 13
// ttH_Hap = 14
// VBF_1j = 15;
// VBF_2j = 16;
// VBF_PTJET1_GT200_2J = 17
// VBF_ptj_GT200_VH = 18
// VBF_rest_VH =19 
// ggH_VBF_2j=20 
// ggH_VBF_3j=21 
// VBF_VBF_3j=22 
// VH_Lep_GT150=23 

// stage1p1_reco:
// ggH_0J =0
// ggH_1J_PTH_0_60 =1
// ggH_1J_PTH_60_120 =2
// ggH_1J_PTH_120_200=3
// ggH_1J_PTH_200=4
// ggH_2J_PTH_0_60 =5
// ggH_2J_PTH_60_120 =6
// ggH_2J_PTH_120_200=7
// ggH_2J_PTH_200=8
// ggH_VBF_2j=9 
// ggH_VBF_3j=10 

// VBF_1j = 11;
// VBF_GT200_1J = 12;
// VBF_2j_mjj_120_350 = 13;
// VBF_2j_mjj_350_700_2j = 14;
// VBF_2j_mjj_GT700_2j = 15;
// VBF_2j_mjj_350_700_3j = 16;
// VBF_2j_mjj_GT700_3j = 17;
// VBF_2j = 18;
// VBF_GT200_2J = 19;

// VBF_GT200_VH = 20;
// VH_Had = 21; 
// VBF_rest_VH =22;

// VH_lep_0_75 =23;
// VH_lep_75_150 =24;
// VH_Lep_GT150=25;
// ttH_Lep = 26 ;
// ttH_Hap = 27;
#include "TString.h"
#include "TH1F.h"

float bins_hpt4[]={0,60,120,200};

int stage1_reco(int Njets, float pTj1, float mjj, float deta_jj, float H_pt, string &reco_catName,int category, float D1jet){
	static auto hpt_bin = make_unique<TH1F>("hpt_bin","",3, bins_hpt4);
	hpt_bin->SetDirectory(0);
	if(category>2 ){
		switch (category){
			case 3: reco_catName = "VH_Lep"; break; 
			case 4: reco_catName = "VH_Had"; break; 
			case 6: reco_catName = "TTH_Had"; break; 
			case 5: reco_catName = "TTH_Lep"; break; 
		}
		return  8+category; 
	}
	else if(category!=0){
		int vbfTopo=0;
		if (Njets<2) vbfTopo=0; 
		vbfTopo = mjj > 400.0 && deta_jj > 2.8;

		//if (pTj1>200 && H_pt<200)   {reco_catName = "VBF_ptj_GT200"; return 10;}

		//		if (pTj1>200 )   {reco_catName = "VBF_ptj_GT200"; return 10;}
		if (pTj1>200 )   {
			if(category==1)
			{reco_catName = "VBF_ptj_GT200_1J"; return 10;}
			else if (category==2)
			{reco_catName = "VBF_ptj_GT200_2J"; return 17;}
		}
		else if(vbfTopo) 	{reco_catName = "VBF_VBF"; return 9;}
		else {
			if(category==1){
				if(D1jet>0.7)
				{reco_catName = "VBF_1j"; return 15;}
				else
					goto lineuntag;
			}
			else if (category == 2)
			{reco_catName = "VBF_2j"; return 16;}
		}
	}
	else{
		if (Njets==0)        {reco_catName = "ggH_0j"; return 0;}
		else if (Njets==1)   {
lineuntag:
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_1j_0_60"; break;
				case 2: reco_catName = "ggH_1j_60_120"; break;
				case 3: reco_catName = "ggH_1j_120_200"; break;
				case 4: reco_catName = "ggH_1j_GT200"; break;
			}
			return binpt;
		} 
		else if (Njets >= 2){
			if ( 60 < mjj && mjj < 120 ) {reco_catName = "VH_Had";return 12;}
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_2j_0_60"; break;
				case 2: reco_catName = "ggH_2j_60_120"; break;
				case 3: reco_catName = "ggH_2j_120_200"; break;
				case 4: reco_catName = "ggH_2j_GT200"; break;
			}
			return binpt+4;
		}
	}
	throw std::exception();
}

int stage1_reco_sync(int Njets, float pTj1, float mjj, float deta_jj, float H_pt, string &reco_catName,int category, float D1jet,float pt_hjj){
	static auto hpt_bin = make_unique<TH1F>("hpt_bin","",3, bins_hpt4);
	hpt_bin->SetDirectory(0);
	int vbfTopo=0;
	if (Njets<2) vbfTopo=0; 
	vbfTopo = mjj > 400.0 && deta_jj > 2.8;
	//	if(category>2 && category!=4){
	if(category>4 ){
		switch (category){
			//	case 3: reco_catName = "VH_Lep"; break; 
			//	case 4: reco_catName = "VH_Had"; break; 
			case 6: reco_catName = "TTH_Had"; break; 
			case 5: reco_catName = "TTH_Lep"; break; 
		}
		return  8+category; 
	}
	else if(category==3){
		if(H_pt<150){reco_catName = "VH_Lep_0_150"; return 11;}
		else {reco_catName = "VH_Lep_GT150"; return 23;}
	}
	else if(category!=0){
		if (pTj1>200 )   {
			if(category==1)
			{reco_catName = "VBF_ptj_GT200_1J"; return 10;}
			else if (category==2)
			{reco_catName = "VBF_ptj_GT200_2J"; return 17;}
			else if (category==4)
			{reco_catName = "VBF_ptj_GT200_VH"; return 18;}
		}
		else if(vbfTopo) 	{
			if(category==4)
			{reco_catName = "VBF_rest_VH"; return 19;}
			else
			{ 
				if(pt_hjj < 25)
				{reco_catName = "VBF_VBF_2j"; return 9;}
				else
				{reco_catName = "VBF_VBF_3j"; return 22;}
			}
		}
		else if ( 60 < mjj && mjj < 120 && category==4){
			reco_catName = "VH_Had"; return 12;
		}
		else {
			if(category==1){
				if(D1jet>0.7)
				{reco_catName = "VBF_1j"; return 15;}
				else
					goto lineuntag;
			}
			else if (category==2)
			{reco_catName = "VBF_2j"; return 16;}
			else
			{reco_catName = "VBF_rest_VH"; return 19;}

		}
	}
	else{
		if (Njets==0)        {reco_catName = "ggH_0j"; return 0;}
		else if (Njets==1)   {
lineuntag:
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_1j_0_60"; break;
				case 2: reco_catName = "ggH_1j_60_120"; break;
				case 3: reco_catName = "ggH_1j_120_200"; break;
				case 4: reco_catName = "ggH_1j_GT200"; break;
			}
			return binpt;
		} 
		else if ( Njets>=2) {
			//			if ( 60 < mjj && mjj < 120 ) {reco_catName = "VH_Had";return 12;}
			if(vbfTopo) {
				if(pt_hjj<25)
				{reco_catName = "ggH_VBF_2j";return 20;}
				else
				{reco_catName = "ggH_VBF_3j";return 21;}
			}
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_2j_0_60"; break;
				case 2: reco_catName = "ggH_2j_60_120"; break;
				case 3: reco_catName = "ggH_2j_120_200"; break;
				case 4: reco_catName = "ggH_2j_GT200"; break;
			}
			return binpt+4;
		}
	}
	throw std::exception();
}
int stage1_reco_stage1(int Njets, float pTj1, float mjj, float deta_jj, float H_pt, string &reco_catName,int category){
	static auto hpt_bin = make_unique<TH1F>("hpt_bin","",3, bins_hpt4);
	if(category>2 && category!=4){
		switch (category){
			case 3: reco_catName = "VH_Lep"; break; 
				//			case 4: reco_catName = "VH_Had"; break; 
			case 6: reco_catName = "TTH_Had"; break; 
			case 5: reco_catName = "TTH_Lep"; break; 
		}
		return  8+category; 
	}
	else{ 
		int vbfTopo=0;
		if (Njets<2) vbfTopo=0; 
		vbfTopo = mjj > 400.0 && deta_jj > 2.8;

		if (pTj1>200 && H_pt<200)   {reco_catName = "VBF_ptj_GT200"; return 10;}
		//if (pTj1>200 )   {reco_catName = "VBF_ptj_GT200"; return 10;}
		else if(vbfTopo) 	{reco_catName = "VBF_VBF"; return 9;}
		if (Njets==0)        {reco_catName = "ggH_0j"; return 0;}
		else if (Njets==1)   {
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_1j_0_60"; break;
				case 2: reco_catName = "ggH_1j_60_120"; break;
				case 3: reco_catName = "ggH_1j_120_200"; break;
				case 4: reco_catName = "ggH_1j_GT200"; break;
			}
			return binpt;
		} 
		else if (Njets>=2) {
			if ( 60 < mjj && mjj < 120 ) {reco_catName = "VH_Had";return 12;}
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_2j_0_60"; break;
				case 2: reco_catName = "ggH_2j_60_120"; break;
				case 3: reco_catName = "ggH_2j_120_200"; break;
				case 4: reco_catName = "ggH_2j_GT200"; break;
			}
			return binpt+4;
		}
	}
	throw std::exception();
}
int stage1_reco_1p1_sync(int Njets, float pTj1, float mjj, float deta_jj, float H_pt, string &reco_catName,int category, float D1jet,float pt_hjj){
	static auto hpt_bin = make_unique<TH1F>("hpt_bin","",3, bins_hpt4);
	int vbfTopo=0;
	if (Njets<2) vbfTopo=0; 
	vbfTopo = mjj > 400.0 && deta_jj > 2.8;
	if(category>4 ){
		switch (category){
			case 6: reco_catName = "TTH_Had"; break; 
			case 5: reco_catName = "TTH_Lep"; break; 
		}
		return  21+category; 
	}
	else if(category==3){
		if(H_pt<150 && H_pt>75){reco_catName = "VH_Lep_75_150"; return 24;}
		else if(H_pt<75){reco_catName = "VH_Lep_0_75"; return 23;}
		else {reco_catName = "VH_Lep_GT150"; return 25;}
	}
	else if(category!=0){
		if (H_pt>200 )   {
			if(category==1)
			{reco_catName = "VBF_GT200_1J"; return 12;}
			else if (category==2)
			{reco_catName = "VBF_GT200_2J"; return 19;}
			else if (category==4)
			{reco_catName = "VBF_GT200_VH"; return 20;}
		}
		else {
			if (category ==1){
				if(D1jet>0.7)
				{reco_catName = "VBF_1j"; return 11;}
				else
					goto lineuntag; 
			}
			else if (category ==2){
				if (mjj > 120 && mjj < 350 )
				{reco_catName = "VBF_2j_mjj_120_350"; return 13;}
				else if (mjj > 350 && mjj < 700 ){
					if (pt_hjj<25) {reco_catName = "VBF_2j_mjj_350_700_2j"; return 14;}
					else {reco_catName = "VBF_2j_mjj_350_700_3j"; return 16;}
				}
				else if (mjj > 700 ){
					if (pt_hjj<25) {reco_catName = "VBF_2j_mjj_GT700_2j"; return 15;}
					else {reco_catName = "VBF_2j_mjj_GT700_3j"; return 17;}
				}
				else 
				{reco_catName = "VBF_2j"; return 18;}
			}
			else if (category == 4){
				if ( 60 < mjj && mjj < 120)
				{reco_catName = "VH_Had"; return 21;}
				else
				{reco_catName = "VBF_rest_VH"; return 22;}
			}
		}
	}
	else {
		if (Njets==0)        {reco_catName = "ggH_0j"; return 0;}
		else if (Njets==1)   {
lineuntag:
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_1j_0_60"; break;
				case 2: reco_catName = "ggH_1j_60_120"; break;
				case 3: reco_catName = "ggH_1j_120_200"; break;
				case 4: reco_catName = "ggH_1j_GT200"; break;
			}
			return binpt;
		} 
		else if ( Njets>=2) {
			if(vbfTopo) {
				if(pt_hjj<25)
				{reco_catName = "ggH_VBF_2j";return 9;}
				else
				{reco_catName = "ggH_VBF_3j";return 10;}
			}
			int binpt = hpt_bin->FindBin(H_pt);
			switch(binpt){
				case 1: reco_catName = "ggH_2j_0_60"; break;
				case 2: reco_catName = "ggH_2j_60_120"; break;
				case 3: reco_catName = "ggH_2j_120_200"; break;
				case 4: reco_catName = "ggH_2j_GT200"; break;
			}
			return binpt+4;
		}
	}
	throw std::exception();
}
int convertCat (int htxs_stage1_red_cat){
	switch (htxs_stage1_red_cat){
		case 101: return 16; //ggH_VBF_2j
		case 102: return 18; //ggH_VBF_3j
		case 103: return 0; //ggH_0j 
		case 104: return 1;// ggH_1j_0_60 
		case 105: return 2;// ggH_1j_60_120
		case 106: return 3; // ggH_1j_120_200 
		case 107: return 4; //ggH_1j_200 
		case 108: return 5; //ggH_2j_0_60 
		case 109: return 6; //ggH_2j_60_120 
		case 110: return 7; //ggH_2j_120_200
		case 111: return 8;  //ggH_2j_200
		case 201: return 9; //VBF_VBF_2j
		case 202: return 19; //VBF_VBF_3j
		case 205: return 10; //VBF_jpt200
		case 301: return 11; //VH lep_0_150
		case 302: return 20; //VH lep_GT150
		case 203: return 12; //VH Had
		case 601: return 13; //ttH
		case 204: return 15; //VBF rest
		case 701: return 17; //BBH
		case -1: return -1; 
		case -2: return -2; 
		default : return 200;

	}
}

TString convertsys(TString inputname){
	if(inputname.Contains("ggH"))
		return "ggH";
	else if (inputname.Contains("Wminus")||inputname.Contains("Wplus"))
		return "WH";
	else if (inputname.Contains("ZH"))
		return "ZH";
	else if (inputname.Contains("VBF"))
		return "VBF";
	else if (inputname.Contains("ttH"))
		return "ttH";
	else if (inputname.Contains("bbH"))
		return "bbH";
	else if (inputname.Contains("ZZTo"))
		return "qqZZ";
	throw std::exception();
}
TString convertStage0Name(int category){
	TString categoryName_s;
	switch (category){
		case 0: categoryName_s = "ggH_untagged";break;
		case 1: categoryName_s = "VBF_1jet";break;
		case 2: categoryName_s = "VBF_2jets";break;
		case 3: categoryName_s = "VH_Lep";break;
		case 4: categoryName_s = "VH_Had";break;
		case 5: categoryName_s = "TTH_Lep";break;
		case 6: categoryName_s = "TTH_Had";break;

	}
	return categoryName_s;
}

int stage1_reco_sync(int Njets, float pTj1, float mjj, float deta_jj, float H_pt, int category, float D1jet, float pt_hjj){
	string reco_catName;
	return stage1_reco_sync(Njets, pTj1, mjj, deta_jj, H_pt, reco_catName, category, D1jet, pt_hjj);
}
