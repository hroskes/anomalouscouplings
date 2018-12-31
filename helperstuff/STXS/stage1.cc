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
// VBF_PTJET1_GT200 = 10
// VH_lep =11
// VH2JET = 12
// ttH_Lep = 13
// ttH_Hap = 14
// VBF_1j = 15;
// VBF_2j = 16;
// VBF_REST = ?

#include "TH1F.h"
#include "TString.h"

#include <memory>
#include <stdexcept>

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
		if (pTj1>200 )   {reco_catName = "VBF_ptj_GT200"; return 10;}
		else if(vbfTopo) 	{reco_catName = "VBF_VBF"; return 9;}
		else {
			if(category==1){
				if(D1jet>0.7)
				{reco_catName = "VBF_1j"; return 15;}
				else
					goto lineuntag;
			}
			else if (category==2)
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
		else {
			throw invalid_argument(TString("Njets = ") += Njets);
		}
	}
	throw runtime_error("?????");
}

int stage1_reco_stage1(int Njets, float pTj1, float mjj, float deta_jj, float H_pt, /*string &reco_catName,*/int category){
	static auto hpt_bin = make_unique<TH1F>("hpt_bin","",3, bins_hpt4);
	hpt_bin->SetDirectory(0);
	string reco_catName;
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
		else {
			throw invalid_argument(TString("Njets = ") += Njets);
		}
	}
}
int convertCat (int htxs_stage1_red_cat){
	switch (htxs_stage1_red_cat){
		case 101: return 16; //ggH_VBF
		case 103: return 0; //ggH_0j 
		case 104: return 1;// ggH_1j_0_60 
		case 105: return 2;// ggH_1j_60_120
		case 106: return 3; // ggH_1j_120_200 
		case 107: return 4; //ggH_1j_200 
		case 108: return 5; //ggH_2j_0_60 
		case 109: return 6; //ggH_2j_60_120 
		case 110: return 7; //ggH_2j_120_200
		case 111: return 8;  //ggH_2j_200
		case 201: return 9; //VBF_VBF
		case 205: return 10; //VBF_jpt200
		case 301: return 11; //VH lep
		case 203: return 12; //VH Had
		case 601: return 13; //ttH
		case 204: return 15; //VBF rest
		case -1: return -1; 
		case -2: return -2; 

	}
	throw invalid_argument(TString("htxs_stage1_red_cat = ")+=htxs_stage1_red_cat);
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
	else
		throw invalid_argument("inputname = "+inputname);
}
