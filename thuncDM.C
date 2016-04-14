#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "TCanvas.h"

const int goldberg =   9;  // number of muR/muF variations 
const int diabelli = 100;  // number of PDF variations
const int nbbin    =  20;  // number of bins either "mll" or "mth" are divide in 


//////////////  functions  /////////////////////////
void GetThSyst(TString sample);

////////////////////////////////////////////////////


//////////////  main  //////////////////////
void thuncDM(){

	GetThSyst("02_WZTo3LNu" );
	GetThSyst("03_ZZ"       );
	GetThSyst("04_TTTo2L2Nu");
	GetThSyst("05_ST"       );
	GetThSyst("06_WW"       );
	GetThSyst("07_ZJets"    );
	GetThSyst("09_TTV"      );
	GetThSyst("10_HWW"      );
	GetThSyst("11_Wg"       );
	//GetThSyst("12_Zg"       );
	//GetThSyst("13_VVV"      );
	GetThSyst("14_HZ"       );

}
////////////////////////////////////////////




////////////////////////////////////////////
void GetThSyst(TString sample){


	TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/minitrees_thunc/TTDM/" + sample + ".root", "READ" ); 

	TH1F *h_PDF = (TH1F*)file -> Get("h_pdfsum");
	TH1F *h_QCD = (TH1F*)file -> Get("h_qcdsum");
	
	TTree *tree = (TTree*) file -> Get("latino");

	// some of the variables below you need to set the cuts

	float MVA; 	tree -> SetBranchAddress( "mva_ttDM0001scalar0010", &MVA);

	vector<float> *weight; 	weight = 0; 	tree -> SetBranchAddress( "LHEweight", &weight );

        /////////////////////////////////////////////////////////////////////////////////

	float z[goldberg];  // RF  
	float y[diabelli];  // PDF 

	// just initialize 

	for( int a = 0; a < goldberg; a++ ){

		z[a] = 0;
	
	}

	for( int a = 0; a < diabelli; a++ ){

		y[a] = 0;

	}
	////////////////////////////////////////////////

	// getting into the tree...

	int nentries = tree -> GetEntries();    

	for( int i = 0; i < nentries; i++) {

		tree -> GetEntry(i);

		
		// cuts 
		// --------------------------------------------------------------------------------------------
		if( MVA < 0.16 ) continue;
		// --------------------------------------------------------------------------------------------	


		for( int a = 0; a < goldberg; a++ ){

			 z[a] = z[a] + weight -> at( a );      // RF
	
		}


		for( int a = 0; a < diabelli; a++ ){

			 y[a] = y[a] + weight -> at( a + 9 );  // PDF       WE KNOW THAT THE FIRST PDF VARIATION OCCUPIES THE 10th ELEMENT
	
		}


	}



	// displaying some results... 
	

	// for QCD 
	// -------

	// pseudo-GEN

	float QCD_GEN_ct  = h_QCD -> GetBinContent(1); 

	float QCD_GEN_up  = h_QCD -> GetBinContent(9);

	float QCD_GEN_do  = h_QCD -> GetBinContent(5);

	float QCD_GEN_ratio_up = QCD_GEN_up / QCD_GEN_ct ;
	float QCD_GEN_ratio_do = QCD_GEN_do / QCD_GEN_ct ;

	// RECO

	float QCD_RECO_ct  = z[0]; 

	float QCD_RECO_up  = z[8]; 

	float QCD_RECO_do  = z[4]; 

	float QCD_RECO_ratio_up = QCD_RECO_up / QCD_RECO_ct ;
	float QCD_RECO_ratio_do = QCD_RECO_do / QCD_RECO_ct ;

      

	// for PDF
	// -------

	// pseudo-GEN

	float PDF_GEN_mean = h_PDF -> Integral() / diabelli; 

	float PDF_GEN_mean_sq = 0;

	for( int a = 0; a < diabelli; a++ ){

		PDF_GEN_mean_sq += (   pow( h_PDF -> GetBinContent(a+1), 2 ) / diabelli   );  	

	}

	float PDF_GEN_SD = sqrt(   PDF_GEN_mean_sq - pow( PDF_GEN_mean, 2 )   ); 

	float PDF_GEN_ratio = 1 + ( PDF_GEN_SD / PDF_GEN_mean );


	// RECO

	float PDF_RECO_mean    = 0; 
	float PDF_RECO_mean_sq = 0;

	for( int a = 0; a < diabelli; a++ ){

		PDF_RECO_mean    += (      y[a]     / diabelli  );

		PDF_RECO_mean_sq += (  pow(y[a], 2) / diabelli  );

	}

	float PDF_RECO_SD  = sqrt(   PDF_RECO_mean_sq - pow( PDF_RECO_mean, 2 )   );

	float PDF_RECO_ratio = 1 + ( PDF_RECO_SD / PDF_RECO_mean );



 	cout << ""                                                                                                        << endl;
 	cout << "sample: " << sample                                                                                      << endl;
	cout << "-----------------------------------------------------------"                                             << endl; 
	cout << "  >> QCD up   -- xs = " << QCD_GEN_ratio_up << "      -- acc = " << QCD_RECO_ratio_up / QCD_GEN_ratio_up << endl;
	cout << "  >> QCD down -- xs = " << QCD_GEN_ratio_do << "      -- acc = " << QCD_RECO_ratio_do / QCD_GEN_ratio_up << endl;
	cout << "-----------------------------------------------------------"                                             << endl; 
	cout << "  >> PDF      -- xs = " << PDF_GEN_ratio    << "      -- acc = " << PDF_RECO_ratio    / PDF_GEN_ratio    << endl;
 	cout << ""                                                                                                        << endl;
 	cout << ""                                                                                                        << endl;

}  
////////////////////////////////////////////






