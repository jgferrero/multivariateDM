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
void thunc2(TString sample);
float  GetRECO(TString sample);  // RECO computation
float   GetGEN(TString sample);  // GEN computation
////////////////////////////////////////////////////


//////////////  main  //////////////////////
void thunc(){

	//thunc2("WWTo2L2Nu"              );
	//thunc2("WZTo3LNu"               );
	thunc2("GluGluHToWWTo2L2Nu_M126");
	thunc2("VBFHToWWTo2L2Nu_M126"   );

}
////////////////////////////////////////////


////////////////////////////////////////////
void thunc2(TString sample){

	float gen  =  GetGEN(sample);
	float reco = GetRECO(sample);

	cout << sample << "         -- xs = " << gen << "    -- acc = " << reco/gen << endl;  	

}
////////////////////////////////////////////


////////////////////////////////////////////
float GetRECO(TString sample){


	TFile* file = new TFile( "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox_76x/22Jan_25ns_mAODv2_MC__l2loose__hadd__bSFL2pTEff__l2tight/latino_" + sample + ".root", "READ" ); 
        /////////////////////////////////////////////////////////////////////////////////


	// protection 1
	bool TheHistoExists = false; 
	TheHistoExists = file -> GetListOfKeys() -> Contains("mcWeightExplainedOrdered");
	if ( TheHistoExists == false ) return -1; 

	// protection 2
	TH1F *histo = (TH1F*)file -> Get("mcWeightExplainedOrdered");
	int nlabel = histo -> GetEntries(); 
	if ( nlabel == 0 ) return -1; 
        /////////////////////////////////////////////////////////////////////////////////


	TTree *tree = (TTree*) file -> Get("latino");

	// some of the variables below you need to set the cuts

	float mll       ; tree -> SetBranchAddress( "mll"       , &mll        );
	float mth       ; tree -> SetBranchAddress( "mth"       , &mth        );
	float metPfType1; tree -> SetBranchAddress( "metPfType1", &metPfType1 );
	float ptll      ; tree -> SetBranchAddress( "ptll"      , &ptll       );

	vector<float> *weight                   ; weight                    = 0;  tree -> SetBranchAddress( "std_vector_LHE_weight"    , &weight                    );
	vector<float> *std_vector_lepton_pt     ; std_vector_lepton_pt      = 0;  tree -> SetBranchAddress( "std_vector_lepton_pt"     , &std_vector_lepton_pt      );
	vector<float> *std_vector_lepton_flavour; std_vector_lepton_flavour = 0;  tree -> SetBranchAddress( "std_vector_lepton_flavour", &std_vector_lepton_flavour );
	vector<float> *std_vector_jet_pt        ; std_vector_jet_pt         = 0;  tree -> SetBranchAddress( "std_vector_jet_pt"        , &std_vector_jet_pt         );
	vector<float> *std_vector_jet_cmvav2    ; std_vector_jet_cmvav2     = 0;  tree -> SetBranchAddress( "std_vector_jet_cmvav2"    , &std_vector_jet_cmvav2     );
        /////////////////////////////////////////////////////////////////////////////////


	float z_total[goldberg]          ;   // RF  global
	float z      [goldberg][nbbin][2];   // RF  per variable and bin
	float y_total[diabelli]          ;   // PDF global
	float y      [diabelli][nbbin][2];   // PDF per variable and bin


	// just initialize 

	for( int a = 0; a < goldberg; a++ ){

		z_total[a] = 0;

		for( int bin = 0; bin < nbbin; bin++ ){

			z[a][bin][0] = 0;

			z[a][bin][1] = 0;

		} 

	}

	for( int a = 0; a < diabelli; a++ ){

		y_total[a] = 0;

		for( int bin = 0; bin < nbbin; bin++ ){

			y[a][bin][0] = 0;

			y[a][bin][1] = 0;

		} 

	}
	////////////////////////////////////////////////


	// getting into the tree...

	int nentries = tree -> GetEntries();    

	for( int i = 0; i < nentries; i++) {

		tree -> GetEntry(i);

		
		// cuts according to Andrea
		// https://github.com/latinos/PlotsConfigurations/blob/master/Configurations/ggH/cuts.py#L5-L30
		// --------------------------------------------------------------------------------------------
		if( mll < 12                           ) continue;
		if( std_vector_lepton_pt -> at(0) < 20 ) continue;  
		if( std_vector_lepton_pt -> at(1) < 10 ) continue; 
		if( std_vector_lepton_pt -> at(1) < 10 ) continue; 
		if( metPfType1 < 20                    ) continue;
		if( ptll < 30                          ) continue;
   
		if( std_vector_lepton_flavour->at(0) * std_vector_lepton_flavour->at(1) != -11*13       ) continue;
		if( ( abs(std_vector_lepton_flavour->at(1)) == 13 && std_vector_lepton_pt->at(1) < 13 ) ) continue;
		if( mth < 60                                                                            ) continue; 
		if( std_vector_jet_pt -> at(0) > 30                                                     ) continue; 
		if( ( std_vector_jet_pt->at(0) > 20 && std_vector_jet_cmvav2->at(0) > -0.715 )          ) continue; 
		if( ( std_vector_jet_pt->at(1) > 20 && std_vector_jet_cmvav2->at(1) > -0.715 )          ) continue; 
		if( ( std_vector_jet_pt->at(2) > 20 && std_vector_jet_cmvav2->at(2) > -0.715 )          ) continue; 
		if( ( std_vector_jet_pt->at(3) > 20 && std_vector_jet_cmvav2->at(3) > -0.715 )          ) continue; 
		if( ( std_vector_jet_pt->at(4) > 20 && std_vector_jet_cmvav2->at(4) > -0.715 )          ) continue; 
		if( ( std_vector_jet_pt->at(5) > 20 && std_vector_jet_cmvav2->at(5) > -0.715 )          ) continue;
		if( ( std_vector_jet_pt->at(6) > 20 && std_vector_jet_cmvav2->at(6) > -0.715 )          ) continue;
		if( ( std_vector_jet_pt->at(7) > 20 && std_vector_jet_cmvav2->at(7) > -0.715 )          ) continue;
		if( ( std_vector_jet_pt->at(8) > 20 && std_vector_jet_cmvav2->at(8) > -0.715 )          ) continue;
		if( ( std_vector_jet_pt->at(9) > 20 && std_vector_jet_cmvav2->at(9) > -0.715 )          ) continue;
		// --------------------------------------------------------------------------------------------	


		for( int a = 0; a < goldberg; a++ ){

			 z_total[a] = z_total[a] + weight -> at( a );      // RF
	
		}


		for( int a = 0; a < diabelli; a++ ){

			 y_total[a] = y_total[a] + weight -> at( a + 9 );  // PDF       WE KNOW THAT THE FIRST PDF VARIATION OCCUPIES THE 10th ELEMENT
	
		}



		for( int bin = 0; bin < nbbin; bin++ ){

			// mll
			// ---

			if( (mll/10 - bin) >= 0 && (mll/10 - bin) < 1 ){    //  how to assign a bin to a branch value

				for( int a = 0; a < goldberg; a++ ){   

					z[a][bin][0] = z[a][bin][0] + weight -> at( a  );      // RF

				}

				for( int a = 0; a < diabelli; a++ ){

					 y[a][bin][0] = y[a][bin][0] + weight -> at( a + 9 );  // PDF       WE KNOW THAT THE FIRST PDF VARIATION OCCUPIES THE 10th ELEMENT

				}

			}

			// mth
			// ---
			if( (mth/20 - bin) >= 0 && (mth/20 - bin) < 1 ){    //  how to assign a bin to a branch value

				for( int a = 0; a < goldberg; a++ ){

					 z[a][bin][1] = z[a][bin][1] + weight -> at( a );     // RF

				}

				for( int a = 0; a < diabelli; a++ ){

					 y[a][bin][1] = y[a][bin][1] + weight -> at( a + 9 );  // PDF       WE KNOW THAT THE FIRST PDF VARIATION OCCUPIES THE 10th ELEMENT

				}

			}

		}

	}

	// displaying some results... 
	
	// for QCD... 

	float RFup = abs(z_total[8] - z_total[0]) / z_total[0];
	float RFdn = abs(z_total[4] - z_total[0]) / z_total[0];

	float quotient = z_total[8] / z_total[0];


	// for PDF

	float m = 0; float m2 = 0;

	for( int a = 0; a < diabelli; a++ ){

		m  = m  + y_total[a];

		m2 = m2 + pow(y_total[a], 2);

	}

	float mean = m/diabelli; 
	float unc  = sqrt(m2/diabelli - pow(mean, 2));
	float recoup = mean + unc; 
	//float quotient= recoup/mean;


	return quotient; 


	////////////////////
 

	// plotter  to the end
	// -------------------

	/*TH1F* h_central= new TH1F( "h_central", " ", 20, 0, 400 );          
	TH1F* h_up     = new TH1F( "h_up"    , " ", 20, 0, 400 );           
	TH1F* h_dn     = new TH1F( "h_dn"    , " ", 20, 0, 400 );           

	for( int bin = 0; bin < nbbin; bin++ ){

		float r = 0; float r2 = 0; 

		for( int a = 0; a < diabelli; a++ ){

			r  = r  + y[a][bin][1];    			

			r2 = r2 + pow( y[a][bin][1], 2 );    		 


		}

		float central = r/diabelli;                                       float central_ratio = 1; 
		float up      = central + sqrt( r2/diabelli - pow(central, 2) );  float up_ratio = 1; if (central != 0) up_ratio = up/central; 
		float dn      = central - sqrt( r2/diabelli - pow(central, 2) );  float dn_ratio = 1; if (central != 0) dn_ratio = dn/central;

		h_central-> Fill( 20*bin, central_ratio );        
		h_up     -> Fill( 20*bin, up_ratio      );          
		h_dn     -> Fill( 20*bin, dn_ratio      );       


	}


	TCanvas* c = new TCanvas( "pdf", "", 600, 600);
	h_central-> SetLineColor(kBlack);  h_central -> SetLineWidth(3);
	h_up     -> SetLineColor(kRed);
	h_dn     -> SetLineColor(kRed);

	h_up -> SetTitle(sample + "     PDF unc");
	h_up -> SetXTitle("mth" );                             
	h_up -> SetYTitle("rel. unc." );                   
	h_up -> SetMaximum( 1.1);
	h_up -> SetMinimum( 0.9);

	h_up     -> Draw("hist");
	h_dn     -> Draw("hist, same");
	h_central-> Draw("hist, same");

	c -> SaveAs("../<plots-subdirectory>/" + sample + "_mth_ratio.png");            
	c -> SaveAs("../<plots-subdirectory>/" + sample + "_mth_ratio.eps");              */
 
}  
////////////////////////////////////////////





////////////////////////////////////////////
float GetGEN(TString sample){

	TFile* file = new TFile( "/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox_76x/22Jan_25ns_mAODv2_MC__l2loose__hadd__bSFL2pTEff__l2tight/latino_" + sample + ".root", "READ" ); 
        /////////////////////////////////////////////////////////////////////////////////


	// protection 1
	bool TheHistoExists = false; 
	TheHistoExists = file -> GetListOfKeys() -> Contains("mcWeightExplainedOrdered");
	if ( TheHistoExists == false ) return -1.; 

	// protection 2
	TH1F *histo = (TH1F*)file -> Get("mcWeightExplainedOrdered");
	int nlabel = histo -> GetEntries(); 
	if ( nlabel == 0 ) return -1.; 
        /////////////////////////////////////////////////////////////////////////////////


	TTree *tree = (TTree*) file -> Get("latino"  );	

	TTree *tree2= (TTree*) file -> Get("mcweight");  // this is the tree for GEN !!

	vector<float> *weight; weight = 0;  tree -> SetBranchAddress( "std_vector_LHE_weight", &weight );   // just it's neccesary to check one branch

	////////////////////////////////////////////////

	float z_total[goldberg];   // RF, global
	float y_total[diabelli];   // PDF, global


	// just initialize 

	for( int a = 0; a < goldberg; a++ ){

		z_total[a] = 0;


	}

	for( int a = 0; a < diabelli; a++ ){

		y_total[a] = 0;


	}

	////////////////////////////////////////////////


	// getting into the tree...

	int nentries = tree -> GetEntries();    

	for( int i = 0; i < nentries; i++) {

		tree -> GetEntry(i);

		for( int a = 0; a < goldberg; a++ ){

			 z_total[a] = z_total[a] + weight -> at( a );      // RF
	
		}

		for( int a = 0; a < diabelli; a++ ){

			 y_total[a] = y_total[a] + weight -> at( a + 9 );  // PDF       WE KNOW THAT THE FIRST PDF VARIATION OCCUPIES THE 10th ELEMENT
	
		}


	}
	////////////////////////////////////////////////

	// displaying some results... 	

	// for QCD...

	float RFup = abs(z_total[8] - z_total[0]) / z_total[0];
	float RFdn = abs(z_total[4] - z_total[0]) / z_total[0];

	float quotient = z_total[8] / z_total[0];


	// for PDF

	float m = 0; float m2 = 0;

	for( int a = 0; a < diabelli; a++ ){

		m  = m  + y_total[a];

		m2 = m2 + pow( y_total[a], 2 );

	}

	float mean = m/diabelli; 
	float unc  = sqrt(m2/diabelli - pow(mean, 2));
	float genup= mean + unc;
	//float quotient = genup/mean;

	return quotient;
	
}
////////////////////////////////////////////
