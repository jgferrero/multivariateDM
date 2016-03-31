#include "MVA.h"

float subcuenta(TString sample);

void cuenta(){
 
	float sig = subcuenta("ttDM0001scalar0010");  
	float x00 = subcuenta("00_Fakes"          );  
	float x01 = subcuenta("01_Data"           );  
	float x02 = subcuenta("02_WZTo3LNu"       );  
	float x03 = subcuenta("03_ZZ"             );  
	float x04 = subcuenta("04_TTTo2L2Nu"      );  
	float x05 = subcuenta("05_ST"             );  
	float x06 = subcuenta("06_WW"             );  
	float x07 = subcuenta("07_ZJets"          );  
	float x08 = subcuenta("08_WJets"          );  
	float x09 = subcuenta("09_TTV"            );  
	float x10 = subcuenta("10_HWW"            );  
	float x11 = subcuenta("11_Wg"             );  
	float x12 = subcuenta("12_Zg"             );  
	float x13 = subcuenta("13_VVV"            );  
	float x14 = subcuenta("14_HZ"             );  
}

float subcuenta(TString sample){

	TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/minitrees/TTDM/" + sample + ".root", "READ" );

	TTree *tree = (TTree*) file -> Get( "latino" );

	float mva_ttDM0001scalar0010;

	tree -> SetBranchAddress( "mva_ttDM0001scalar0010", &mva_ttDM0001scalar0010 );

	int nentries = tree -> GetEntries();

	float mycounter = 0; 

	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		tree -> GetEntry(ievt);

		//if (mva_ttDM0001scalar0010  < 0.5 ) continue; 

		mycounter = mycounter + 1; //eventW; 

	}

	cout << sample << "    --- >  " << mycounter << endl; 

	return mycounter;

}
