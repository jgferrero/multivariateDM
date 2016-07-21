#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TColor.h"  
#include "TLatex.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TInterpreter.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "THStack.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraphPainter.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "cosmetics.h"


const int n = 50;   
const double x[n] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49 }; 

void LoadLimits(TString MLP);

//----------

void LimitOptimisation() {

	LoadLimits("01"); 
	LoadLimits("02");  
	LoadLimits("03");  
	LoadLimits("04");  
	LoadLimits("05");  
	LoadLimits("06");  
	LoadLimits("07");    

}

void LoadLimits(TString MLP){	

	cout << "                       " << endl; 
	cout << "                       " << endl; 
	cout << "   >> >> >> MLP " << MLP << endl;  
	cout << "                       " << endl; 


	gInterpreter->ExecuteMacro("GoodStyle.C");

	// ----------------------

	TFile* file = new TFile( "ttDM0001scalar00500_" + MLP + ".root", "READ" );

	TTree* tree = (TTree*) file -> Get( "limit" );

	Double_t limit; 

	tree -> SetBranchAddress( "limit", &limit );


	Double_t y_exp[n]; 

	// initialisation
	for(int j = 0; j < n; j++ ) {
		y_exp[j] = 0;
	} 

	int j = 0; 

	for ( long ievt = 0; ievt < tree -> GetEntries(); ievt++ ) { 

		tree -> GetEntry(ievt); 
		
	 	if (ievt%6==2) y_exp[j++] = limit;

	}

	float expected = 1e8; 
	float cut   = 0;
 	
	for( int j = 0; j < n; j++ ){

		//cout << x[j]*0.02 << " -- " << y_exp[j] << endl; 

		cout << y_exp[j] << endl; 
	
		if(   ( y_exp[j] < expected )  &&  ( y_exp[j] != 0 )   ){

			expected = y_exp[j];    //cout << expected << endl; 

			cut = j; 

		}

	}

	//cout << "   cut at " << 0.02*cut << "    expected limit =  " << expected << endl; 
	//cout << 0.02*cut << " & " << expected << endl; 
	
}





