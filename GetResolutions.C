#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TVector2.h"
#include "TChain.h"
#include "TColor.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "Math/Minimizer.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooVoigtian.h"
#include "RooFitResult.h"

using namespace RooFit;


enum{ parall, transv, MET, nvariable };
enum{ pT, sumET, NVtx, nparameter};
enum{ dat, GJets, QCD, WJets, ZGJets_old, ZNuNuGJets40130_old, ZGTo2LG_old, WGJets_old, WGToLNuG_old, TTGJets_old, TGJets_old, nprocess_old}; 

enum{ data2016D, 
      GJets40100, GJets100200, GJets200400, GJets400600, GJets600Inf, 
      QCD200300, QCD300500, QCD500700, QCD7001000, QCD10001500, QCD15002000, QCD2000Inf,
      WJets100200, WJets200400, WJets400600, WJets600800, WJets8001200, WJets12002500, WJets2500Inf,
      ZGJets, ZNuNuGJets40130, ZGTo2LG, 
      WGJets, WGToLNuG,
      TTGJets, TGJets,
      nprocess }; 

const bool DoFillHistograms = true  ;  
const bool DoGlobalPlots    = false ; 
const bool DoFit            = false ; 

const bool DoNotEventLoop   = false; 

const bool RunOverAllData = true; 
const bool RunOverAllMC   = false; 

      int _TotalEntries[nprocess]; 

const int MaxEntriesData    =2203625;
const int MaxEntriesMC      = 100000;

//const int MaxEntriesMCStored =  20000; // 20000 (11vii)

const TString NTuplaDir = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/gamma/APPROVAL/";

//const float TheLuminosity = 4.32421;   // fb⁻¹
const float TheLuminosity = 1.0;    // fb⁻¹  (4.35)


TString allhistosWriteTo  = "allhistos_02Aug"; 
TString allhistosReadFrom = "allhistos_12Jul"; 


const bool GJetsMC = false;

const float a = 0.5346                ; 
const float b = 0.2166                ; 
const float c = 2 * sqrt( 2 * log(2) ); 

const int HowManyBkgs = 3;   // provisional  // ideally, HowManyBkgs = 5

float          xs     [nprocess]; 
TString sampleID      [nprocess];
TString processID     [nprocess];
TString processIDfancy[nprocess];
Color_t ProcessColor  [nprocess];
 
TString variableID     [nvariable];
TString variableIDfancy[nvariable];
TString sigma_variable [nvariable];

TString parameterID     [nparameter];
TString parameterIDfancy[nparameter];

float maxpT    = 300.; float minpT    =50.; const int nbinpT    = 25; float _uppT   ;
float maxsumET =3000.; float minsumET = 0.; const int nbinsumET = 30; float _upsumET;  
float maxNVtx  =  40.; float minNVtx  = 0.; const int nbinNVtx  =  8; float _upNVtx ;

float maxMET   = 200.; float minMET   =    0.; const int nbinMET   = 40;
float maxuPara = 200.; float minuPara = -200.; const int nbinuPara = 80; 
float maxuPerp = 200.; float minuPerp = -200.; const int nbinuPerp = 80; 

/*float maxpT    = 200.; float minpT    = 0.; const int nbinpT    =  5; float _uppT   ;
float maxsumET = 300.; float minsumET = 0.; const int nbinsumET =  5; float _upsumET;  
float maxNVtx  =  20.; float minNVtx  = 0.; const int nbinNVtx  =  5; float _upNVtx ;*/

const int aa = nvariable;   // just for next blocks
const int bb = nprocess ;   // just for next blocks

// ----- write -----------------------------
TH1F*    h_global_W     [aa][bb]           ;
TH1F*    h_resol_pT_W   [aa][bb][nbinpT   ];
TH1F*    h_resol_sumET_W[aa][bb][nbinsumET];
TH1F*    h_resol_NVtx_W [aa][bb][nbinNVtx ]; 
// ----- read ------------------------------
TH1F*    h_global     [aa][bb]             ;
THStack* s_global     [aa]                 ;
TH1F*    Ratio        [aa]                 ;
TH1F*    h_resol_pT   [aa][bb+1][nbinpT   ];
TH1F*    h_resol_sumET[aa][bb+1][nbinsumET];
TH1F*    h_resol_NVtx [aa][bb+1][nbinNVtx ]; 
// -----------------------------------------


float xpT   [nbinpT   ]; float ypT_A   [nbinpT   ]; float epT_A   [nbinpT   ]; float chi2pT_A   [nbinpT   ]; 
			 float ypT_B   [nbinpT   ]; float epT_B   [nbinpT   ]; float chi2pT_B   [nbinpT   ];

float xsumET[nbinsumET]; float ysumET_A[nbinsumET]; float esumET_A[nbinsumET]; float chi2sumET_A[nbinsumET];
			 float ysumET_B[nbinsumET]; float esumET_B[nbinsumET]; float chi2sumET_B[nbinsumET];

float xNVtx [nbinNVtx ]; float yNVtx_A [nbinNVtx ]; float eNVtx_A [nbinNVtx ]; float chi2NVtx_A [nbinNVtx ];
			 float yNVtx_B [nbinNVtx ]; float eNVtx_B [nbinNVtx ]; float chi2NVtx_B [nbinNVtx ];



void Assign          (                                                                                                          ); 
void  FillHistograms (                                                                                                          );
void GlobalPlots     (                                                                                                          );
void OrderFits       (                                                                                                          );
void  FillHistogram  ( int process                                                                                              );
void GetResolution   ( int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi, TString GJetsOrigin   );
void PlotResolution  ( int ivar, int parameter                                                                                  );
void SetAxis         ( TH1* hist, TString xtitle, TString ytitle, Float_t xoffset, Float_t yoffset                              ); 
void DrawLatex       ( Font_t tfont, Float_t x, Float_t y, Float_t tsize, Short_t align, const char* text                       );
double GetFWHM       ( double sigma, double gamma                                                                               );
double GetFWHMerror  ( double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg );


void GetResolutions(){

	Assign(); 

	FillHistograms(); 

	GlobalPlots();

	OrderFits();


	cout << "            " << endl; 
	cout << "--- VALE ---" << endl; 
	cout << "            " << endl; 
}


void Assign(){


	// maria n-tuplas (approval version) -> https://dalfonso.web.cern.ch/dalfonso/ICHEP/sampleLocation.txt

	// x-sections: 
	//		https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
	//		http://mvesterb.web.cern.ch/mvesterb/met/things/samples.dat

	int i = 0; 

	sampleID[i] = "SinglePhoton_Run2016D_PromptReco_v2";           processID[i] = "data2016D"      ; processIDfancy[i++] = "data 2016-D"  ; 

	sampleID[i] = "GJets_HT40to100"        ; xs[i] =   20730000.0; processID[i] = "GJets40100"     ; processIDfancy[i] = "#gamma + jets"  ; ProcessColor[i++] = kYellow;
	sampleID[i] = "GJets_HT100to200"       ; xs[i] =    9226000.0; processID[i] = "GJets100200"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;
	sampleID[i] = "GJets_HT200to400"       ; xs[i] =    2300000.0; processID[i] = "GJets200400"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;
	sampleID[i] = "GJets_HT400to600"       ; xs[i] =     277400.0; processID[i] = "GJets400600"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;
	sampleID[i] = "GJets_HT600toInf"       ; xs[i] =      93380.0; processID[i] = "GJets600Inf"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;
	sampleID[i] = "QCD_HT200to300_ext"     ; xs[i] = 1735000000.0; processID[i] = "QCD200300"      ; processIDfancy[i] = "QCD"            ; ProcessColor[i++] = kGreen ;
	sampleID[i] = "QCD_HT300to500"         ; xs[i] =  366800000.0; processID[i] = "QCD300500"      ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
	sampleID[i] = "QCD_HT500to700"         ; xs[i] =   29370000.0; processID[i] = "QCD500700"      ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
	sampleID[i] = "QCD_HT700to1000_ext"    ; xs[i] =    6524000.0; processID[i] = "QCD7001000"     ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
	sampleID[i] = "QCD_HT1000to1500"       ; xs[i] =    1064000.0; processID[i] = "QCD10002000"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
	sampleID[i] = "QCD_HT1500to2000"       ; xs[i] =     119900.0; processID[i] = "QCD15002000"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
	sampleID[i] = "QCD_HT2000toInf"        ; xs[i] =      25240.0; processID[i] = "QCD2000Inf"     ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
	sampleID[i] = "WJetsToLNu_HT100to200"  ; xs[i] =    1627450.0; processID[i] = "WJets100200"    ; processIDfancy[i] = "W + jets"       ; ProcessColor[i++] = kBlue  ;
	sampleID[i] = "WJetsToLNu_HT200to400"  ; xs[i] =     435237.0; processID[i] = "WJets200400"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
	sampleID[i] = "WJetsToLNu_HT400to600"  ; xs[i] =      59191.0; processID[i] = "WJets400600"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
	sampleID[i] = "WJetsToLNu_HT600to800"  ; xs[i] =      14581.0; processID[i] = "WJets600800"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
	sampleID[i] = "WJetsToLNu_HT800to1200" ; xs[i] =       6655.0; processID[i] = "WJets8001200"   ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
	sampleID[i] = "WJetsToLNu_HT1200to2500"; xs[i] =       1608.1; processID[i] = "WJets12002500"  ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
	sampleID[i] = "WJetsToLNu_HT2500toInf" ; xs[i] =         38.9; processID[i] = "WJets2500Inf"   ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
	sampleID[i] = "ZGJets"                 ; xs[i] =        190.3; processID[i] = "ZGJets"         ; processIDfancy[i] = "Z#gamma + jets" ; ProcessColor[i++] = kViolet;
	sampleID[i] = "ZNuNuGJets40130"        ; xs[i] =       2816.0; processID[i] = "ZNuNuGJets40130"; processIDfancy[i] = " "              ; ProcessColor[i++] = kViolet;
	sampleID[i] = "ZGTo2LG"                ; xs[i] =     117864.0; processID[i] = "ZGTo2LG"        ; processIDfancy[i] = " "              ; ProcessColor[i++] = kViolet;
	sampleID[i] = "WGJets"                 ; xs[i] =        663.7; processID[i] = "WGJets"         ; processIDfancy[i] = "W#gamma + jets" ; ProcessColor[i++] = kCyan  ;
	sampleID[i] = "WGToLNuG"               ; xs[i] =     585800.0; processID[i] = "WGToLNuG"       ; processIDfancy[i] = " "              ; ProcessColor[i++] = kCyan  ;
	sampleID[i] = "TTGJets"                ; xs[i] =       3697.0; processID[i] = "TTGJets"        ; processIDfancy[i] = "tt#gamma + jets"; ProcessColor[i++] = kPink  ; 
	sampleID[i] = "TGJets"                 ; xs[i] =       3697.0; processID[i] = "TGJets"         ; processIDfancy[i] = "t#gamma + jets" ; ProcessColor[i++] = kPink  ;

	// ------------------------------------------------------------------------------------------



	variableID[parall] = "parallel"  ; 
	variableID[transv] = "transverse"; 
	variableID[MET   ] = "MET"       ;

	variableIDfancy[parall] = "u_{||}+#gamma_{T} [GeV]"; 
	variableIDfancy[transv] = "u_{#perp}  [GeV]"       ;
	variableIDfancy[MET   ] = "E_{T}^{miss} [GeV]"     ;

	sigma_variable[parall] = "#sigma( u_{||} ) [GeV]"   ;
	sigma_variable[transv] = "#sigma( u_{#perp} )  [GeV]";

	parameterID[pT   ] = "pT"   ;
	parameterID[sumET] = "sumET";
	parameterID[NVtx ] = "NVtx" ;

	parameterIDfancy[pT   ] = "#gamma q_{T} [GeV]";
	parameterIDfancy[sumET] = "#Sigma E_{T} [TeV]";
	parameterIDfancy[NVtx ] = "number of vertices";

}

void FillHistograms(){

	if( DoFillHistograms == false ) { cout << "   >>   FillHistograms deactivated !!!   " << endl;  return; }

	for( int j = 0; j < nprocess; j++ ){

		cout << "next n-tupla: " << processID[j] << endl; 

		FillHistogram( j );

	}


	if( DoNotEventLoop == false ){

		cout << "the writing step begins..." << endl; 

		TFile* allhistos = new TFile("histograms/" + allhistosWriteTo +".root", "recreate");

		for( int i = 0; i < nvariable; i++ ){

			for( int j = 0; j < nprocess; j++ ){ 

				h_global_W[i][j] -> Write();

				if( i == MET ) continue;
							
				for( int k = 0; k < nbinpT   ; k++) 	h_resol_pT_W   [i][j][k] -> Write();
				for( int k = 0; k < nbinsumET; k++)	h_resol_sumET_W[i][j][k] -> Write();
				for( int k = 0; k < nbinNVtx ; k++)	h_resol_NVtx_W [i][j][k] -> Write();

			}

		}

		allhistos -> Close();

		cout << "                          " << endl; 
		cout << "                          " << endl;
		cout << "    histograms filled !!! " << endl;
		cout << "                          " << endl;
		cout << "                          " << endl;

	}

}



void FillHistogram( int process ){

	// ----------------------------------------------------------------------------------------------------------------------------------------------

	/*TChain* tree = new TChain("METtree");

	TString path = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/gamma/"; 

	if( process == dat ){

		tree -> Add( path + "4fbV6/SinglePhoton_Run2016B_PromptReco_v2/METtree.root" );

	}


	else if( process == GJets){

		tree -> Add( path + "GJets_HT40to100/METtree.root"  ); 
		tree -> Add( path + "GJets_HT100to200/METtree.root" ); 
		tree -> Add( path + "GJets_HT200to400/METtree.root" ); 
		tree -> Add( path + "GJets_HT400to600/METtree.root" ); 
		//tree -> Add( path + "GJets_HT600toInf/METtree.root" ); 

	}


	else if( process == QCD ){

		//tree -> Add( path + "QCD_HT200to300/METtree.root"   ); 
		tree -> Add( path + "QCD_HT300to500/METtree.root"   ); 
		tree -> Add( path + "QCD_HT500to700_ext/METtree.root"   ); 
		tree -> Add( path + "QCD_HT700to1000_ext/METtree.root"  ); 
		tree -> Add( path + "QCD_HT1000to1500/METtree.root" ); 
		tree -> Add( path + "QCD_HT1500to2000/METtree.root" ); 
		tree -> Add( path + "QCD_HT2000toInf/METtree.root"  ); 			

	}


	else if( process == WJets ){

		tree -> Add( path + "WJetsToLNu_HT100to200/METtree.root"  ); 
		tree -> Add( path + "WJetsToLNu_HT200to400/METtree.root"  ); 
		tree -> Add( path + "WJetsToLNu_HT400to600/METtree.root"  ); 
		tree -> Add( path + "WJetsToLNu_HT600to800/METtree.root"  ); 
		tree -> Add( path + "WJetsToLNu_HT800to1200_ext/METtree.root" ); 
		tree -> Add( path + "WJetsToLNu_HT1200to2500/METtree.root"); 
		tree -> Add( path + "WJetsToLNu_HT2500toInf/METtree.root" ); 	

	}


	else if( process == ZGJets ){

		tree -> Add( path + "ZGJets/METtree.root"  ); 

	}

	
	else if( process == ZNuNuGJets40130 ){

		tree -> Add( path + "ZNuNuGJets40130/METtree.root"  ); 

	}


	else if( process == ZGTo2LG ){

		tree -> Add( path + "ZGTo2LG/METtree.root"  ); 

	}


	else if( process == WGJets ){

		tree -> Add( path + "WGJets/METtree.root"  ); 

	}


	else if( process == WGToLNuG ){

		tree -> Add( path + "WGToLNuG/METtree.root"  ); 

	}


	else if( process == TTGJets ){

		//tree -> Add( path + "TTGJets/METtree.root"  ); 
		
	}


	else if( process == TGJets ){

		tree -> Add( path + "TGJets/METtree.root"  ); 
		
	}


	else {

		return; 

	}*/

	// ----------------------------------------------------------------------------------------------------------------------------------------------


	TFile* TheFile = new TFile( NTuplaDir + sampleID[process] + "/METtree.root", "read"); 

	TTree* tree = (TTree*) TheFile -> Get( "METtree" );

	int   ngamma                            ;
	int   gamma_idCutBased                  ;
	float gamma_r9                          ;	
	float gamma_pt                          ;
	float gamma_eta                         ;
	float gamma_phi                         ;

	float gamma_sigmaIetaIeta ; 

	int   nLepGood10                        ;  
        float lep_pt                            ;  

	int   HBHENoiseFilter                   ;
	int   HBHENoiseIsoFilter                ;
	int   CSCTightHalo2015Filter            ;
	int   EcalDeadCellTriggerPrimitiveFilter;
	int   goodVertices                      ;
	int   eeBadScFilter                     ;

	int HLT_Photon30 ;
	int HLT_Photon50 ;
	int HLT_Photon75 ;
	int HLT_Photon90 ;
	int HLT_Photon120;

	int HLT_Photon30_Prescale ;
	int HLT_Photon50_Prescale ;
	int HLT_Photon75_Prescale ;
	int HLT_Photon90_Prescale ;
	int HLT_Photon120_Prescale;


	float met_sumEt                         ;
	int   nVert                             ;
	float met_pt                            ;
	float met_phi                           ;
	//float uPara                             ;
	//float uPerp                             ;
	float puWeight                          ;
	float genWeight	                        ;


	tree -> SetBranchAddress( "ngamma"                                 , &ngamma                             );
	tree -> SetBranchAddress( "gamma_idCutBased"                       , &gamma_idCutBased                   );
	tree -> SetBranchAddress( "gamma_r9"                               , &gamma_r9                           );
	tree -> SetBranchAddress( "gamma_pt"                               , &gamma_pt                           );
	tree -> SetBranchAddress( "gamma_eta"                              , &gamma_eta                          );
	tree -> SetBranchAddress( "gamma_phi"                              , &gamma_phi                          );

	tree -> SetBranchAddress( "gamma_sigmaIetaIeta"                    , &gamma_sigmaIetaIeta                );

	tree -> SetBranchAddress( "nLepGood10"                             , &nLepGood10                         );
	tree -> SetBranchAddress( "lep_pt"                                 , &lep_pt                             );

	tree -> SetBranchAddress( "Flag_HBHENoiseFilter"                   , &HBHENoiseFilter                    );
	tree -> SetBranchAddress( "Flag_HBHENoiseIsoFilter"                , &HBHENoiseIsoFilter                 );
	tree -> SetBranchAddress( "Flag_CSCTightHalo2015Filter"            , &CSCTightHalo2015Filter             );
	tree -> SetBranchAddress( "Flag_EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter );
	tree -> SetBranchAddress( "Flag_goodVertices"                      , &goodVertices                       );
	tree -> SetBranchAddress( "Flag_eeBadScFilter"                     , &eeBadScFilter                      );

	tree -> SetBranchAddress( "HLT_Photon30"                           , &HLT_Photon30                       );
	tree -> SetBranchAddress( "HLT_Photon50"                           , &HLT_Photon50                       );
	tree -> SetBranchAddress( "HLT_Photon75"                           , &HLT_Photon75                       );
	tree -> SetBranchAddress( "HLT_Photon90"                           , &HLT_Photon90                       );
	tree -> SetBranchAddress( "HLT_Photon120"                          , &HLT_Photon120                      );

	if( process == data2016D ){

		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v_Prescale" , &HLT_Photon30_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale" , &HLT_Photon50_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale" , &HLT_Photon75_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale" , &HLT_Photon90_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale", &HLT_Photon120_Prescale   );

	}


	tree -> SetBranchAddress( "met_sumEt"                              , &met_sumEt                          );
	tree -> SetBranchAddress( "nVert"                                  , &nVert                              );
	tree -> SetBranchAddress( "met_pt"                                 , &met_pt                             );
	tree -> SetBranchAddress( "met_phi"                                , &met_phi                            );
	//tree -> SetBranchAddress( "met_uPara_zll"                          , &uPara                              );
	//tree -> SetBranchAddress( "met_uPerp_zll"                          , &uPerp                              );
	tree -> SetBranchAddress( "puWeight"                               , &puWeight                           );
	tree -> SetBranchAddress( "genWeight"                              , &genWeight                          );



	for( int i = 0; i < nvariable; i++ ){

		if( i != MET ) h_global_W[i][process] = new TH1F( "h_global_" + variableID[i] + "_" + processID[process], variableID[i], nbinuPara, minuPara, maxuPara );
		if( i == MET ) h_global_W[i][process] = new TH1F( "h_global_" + variableID[i] + "_" + processID[process], variableID[i], nbinMET  , minMET  , maxMET   );

		if( i == MET) continue;

		for( int k = 0; k < nbinpT; k++ ){

			h_resol_pT_W   [i][process][k] = new TH1F( Form("h_resol_pT_"    + variableID[i] + "_" + processID[process] + "_%d", k), "resolution pT"        , nbinuPara, minuPara, maxuPara );

		}

		for( int k = 0; k < nbinsumET; k++ ){

			h_resol_sumET_W[i][process][k] = new TH1F( Form("h_resol_sumET_" + variableID[i] + "_" + processID[process] + "_%d", k), "resolution sum ET"    , nbinuPara, minuPara, maxuPara );

		}

		for( int k = 0; k < nbinNVtx; k++ ){

			h_resol_NVtx_W [i][process][k] = new TH1F( Form("h_resol_NVtx_"  + variableID[i] + "_" + processID[process] + "_%d", k), "resolution num of vtx", nbinuPara, minuPara, maxuPara );

		}

	}


	int nentries = tree -> GetEntries();//cout << " nentries = " << nentries << endl;

	_TotalEntries[process] = nentries; 

	if ( DoNotEventLoop == true ) return; 

	if ( process == data2016D  &&  RunOverAllData == false ) nentries = MaxEntriesData; 
		
	if ( process != data2016D  &&  RunOverAllMC   == false ) nentries = MaxEntriesMC  ;

	float baseW = 1.0;	if ( process != data2016D ) baseW = xs[process]/nentries; 

 
	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		tree -> GetEntry(ievt);	

		cout << "\n ievt = " << ievt << endl; 
		
		//if(  process == data2016D  &&  ( ievt == 301922 || ievt == 743589 || ievt == 2493920 )  ) continue;

		//if( ievt%10000 == 0 )  cout << "  >>  evt " << ievt << endl;

			//cout << "ngamma = " << ngamma << endl;

		if ( ngamma < 1                       ) continue;

			//cout << "id = " <<  gamma_idCutBased << endl;

		if ( gamma_idCutBased != 3            ) continue; 

			//cout << "gamma_r9 = " << gamma_r9 << endl;

		if ( gamma_r9 < 0.9 || gamma_r9 > 1.0 ) continue;
 
			//cout << "gamma_pt = " << gamma_pt <<  endl;

		if ( gamma_pt < 50                    ) continue; 

			//cout << "gamma_eta = " << gamma_eta << endl;

		if ( gamma_eta > 1.5                  ) continue; 



		if ( nLepGood10 > 0                   ) continue; 

		// MET filters
		if ( HBHENoiseFilter                    != 1 ) continue; 
		if ( HBHENoiseIsoFilter                 != 1 ) continue; 
		if ( CSCTightHalo2015Filter             != 1 ) continue; 
		if ( EcalDeadCellTriggerPrimitiveFilter != 1 ) continue; 
		if ( goodVertices                       != 1 ) continue; 
		if ( eeBadScFilter                      != 1 ) continue; 	


		float eventW = 1.0 ;

		// triggers
		if( process == data2016D ){

			//cout << HLT_Photon30 << "  -- " << HLT_Photon50 << "  -- " << HLT_Photon75 << "  -- " << HLT_Photon90 << "  -- " << HLT_Photon120 << "  -- " << endl;

			//if( HLT_Photon30 == 0  &&  HLT_Photon50 == 0  &&  HLT_Photon75 == 0  &&  HLT_Photon90 == 0  &&  HLT_Photon120 == 0 ) continue;


			if( HLT_Photon120 != 0 ) eventW *= HLT_Photon120_Prescale;

			else{

				if( HLT_Photon90 != 0 )	eventW *= HLT_Photon90_Prescale;

				else{

					if( HLT_Photon75 != 0 )	eventW *= HLT_Photon75_Prescale;

					else{

						if( HLT_Photon50 != 0 )	eventW *= HLT_Photon50_Prescale;

						else{

							if( HLT_Photon30 != 0 )	eventW *= HLT_Photon30_Prescale;


						}

					}

				}
	
			}



		}


 
		if( process != data2016D ){

			eventW *= baseW    ;

			eventW *= puWeight ;

			eventW *= genWeight;

		}


		if( gamma_pt < minpT || gamma_pt >= maxpT || met_sumEt < minsumET || met_sumEt >= maxsumET || nVert < minNVtx || nVert >= maxNVtx) continue;

		//cout << ievt << " -- " << ngamma << " -- " << gamma_pt << " -- " << met_sumEt << " -- " << nVert << endl;

		int l = floor(    nbinpT    * (gamma_pt  - minpT   ) / (maxpT    - minpT   )    ); //cout << " l = " << l << endl;
		int m = floor(    nbinsumET * (met_sumEt - minsumET) / (maxsumET - minsumET)    ); //cout << " m = " << m << endl;
		int n = floor(    nbinNVtx  * (nVert     - minNVtx ) / (maxNVtx  - minNVtx )    ); //cout << " n = " << n << endl;



		TVector2 qT, ET, uT; 

		ET.SetMagPhi( met_pt  , met_phi   );

		qT.SetMagPhi( gamma_pt, gamma_phi );

		uT = -1* ( ET + qT ); 

		float uPara = (  uT.Px() * qT.Px() + uT.Py() * qT.Py()  ) / qT.Mod() + qT.Mod(); // cout <<  uPara << endl;
		float uPerp = (  uT.Px() * qT.Py() - uT.Py() * qT.Px()  ) / qT.Mod();


		h_global_W     [parall][process]    -> Fill( uPara , eventW ); 
		h_global_W     [transv][process]    -> Fill( uPerp , eventW );
		h_global_W     [MET   ][process]    -> Fill( met_pt, eventW );

		h_resol_pT_W   [parall][process][l] -> Fill( uPara , eventW ); 
		h_resol_pT_W   [transv][process][l] -> Fill( uPerp , eventW );	

		h_resol_sumET_W[parall][process][m] -> Fill( uPara , eventW ); 
		h_resol_sumET_W[transv][process][m] -> Fill( uPerp , eventW );

		h_resol_NVtx_W [parall][process][n] -> Fill( uPara , eventW );
		h_resol_NVtx_W [transv][process][n] -> Fill( uPerp , eventW );			


	}

}



void GlobalPlots(){

	if( DoGlobalPlots == false ) { cout << "   >>   GlobalPlots deactivated !!!   " << endl;  return; }

	TFile* allhistos = new TFile( "histograms/" + allhistosReadFrom + ".root", "read" );

	for( int i = 0; i < nvariable; i++ ){

		for( int j = 0; j < nprocess; j++ ){

			h_global[i][j] = (TH1F*) allhistos -> Get("h_global_" + variableID[i] + "_" + processID[j]);

		}

	
		h_global[i][dat] -> SetMarkerStyle(20);

		h_global[i][dat] -> SetMarkerColor(kBlack);

		h_global[i][dat] -> SetLineColor(kBlack);

	}


	float TheCurrentLuminosity = TheLuminosity; 

	//( RunOverAllData == true )   ?   TheCurrentLuminosity  =  TheLuminosity   :   TheCurrentLuminosity  =  TheLuminosity * MaxEntriesData/_TotalEntries[dat]  ; 


	float weight[nprocess]; 

	for( int j = 1; j < nprocess; j++ ){

		if( j == TTGJets ) continue;

		 weight[j] = xs[j]/_TotalEntries[j];

		for( int i = 0; i < nvariable; i++ ){ 		

			h_global[i][j] -> Scale( weight[j]*TheCurrentLuminosity ); 
 			//h_global[i][j] -> Scale( TheCurrentLuminosity ); 

			//h_global[i][j] -> SetFillColor( ProcessColor[j] );
			h_global[i][j] -> SetFillColorAlpha( ProcessColor[j], 0.7 );
		
		}		

	}

	
	for( int i = 0; i < nvariable; i++ ){

			h_global[i][TTGJets] -> Add( h_global[i][TGJets         ] );

			h_global[i][ZGJets ] -> Add( h_global[i][ZNuNuGJets40130] );
			h_global[i][ZGJets ] -> Add( h_global[i][ZGTo2LG        ] );

			h_global[i][WGJets ] -> Add( h_global[i][WGToLNuG       ] );

	}


	for( int i = 0; i < nvariable; i++ ){

		s_global[i]  = new THStack( variableID[i], variableID[i] );

		//s_global[i] -> Add( h_global[i][TTGJets        ] );
		s_global[i] -> Add( h_global[i][TGJets         ] );
		s_global[i] -> Add( h_global[i][ZGJets         ] );
		//s_global[i] -> Add( h_global[i][ZNuNuGJets40130] );
		//s_global[i] -> Add( h_global[i][ZGTo2LG        ] );
		s_global[i] -> Add( h_global[i][WJets          ] );
		s_global[i] -> Add( h_global[i][WGJets         ] );
		//s_global[i] -> Add( h_global[i][WGToLNuG       ] );
		s_global[i] -> Add( h_global[i][QCD            ] );
		s_global[i] -> Add( h_global[i][GJets          ] );

	}



	for( int i = 0; i < nvariable; i++ ){

		Ratio[i] = (TH1F*) h_global[i][dat] -> Clone( "ratio_" + variableID[i] );

	}



	for( int r = 0; r < nbinuPara; r++ ){

		TH1F* h_allMC = (TH1F*) h_global[parall][GJets] -> Clone( "h_allMC_" + variableID[parall] );

		for( int j = 2; j < nprocess; j++ ) h_allMC->Add( h_global[parall][j] ); 

		float TheRatio = h_global[parall][dat]->GetBinContent(r) / h_allMC->GetBinContent(r); 
		float TheError = h_global[parall][dat]->GetBinError(r)   / h_allMC->GetBinContent(r);

		Ratio[parall] -> SetBinContent(r, TheRatio);
		Ratio[parall] -> SetBinError  (r, TheError);

	}


	for( int r = 0; r < nbinuPerp; r++ ){

		TH1F* h_allMC = (TH1F*) h_global[transv][GJets] -> Clone( "h_allMC_" + variableID[transv] );

		for( int j = 2; j < nprocess; j++ ) h_allMC->Add( h_global[transv][j] ); 

		float TheRatio = h_global[transv][dat]->GetBinContent(r) / h_allMC->GetBinContent(r); 
		float TheError = h_global[transv][dat]->GetBinError(r)   / h_allMC->GetBinContent(r);

		Ratio[transv] -> SetBinError  (r, TheError);

	}
	
	for( int r = 0; r < nbinMET; r++ ){

		TH1F* h_allMC = (TH1F*) h_global[MET][GJets] -> Clone( "h_allMC_" + variableID[MET] );

		for( int j = 2; j < nprocess; j++ ) h_allMC->Add( h_global[MET][j] ); 

		float TheRatio = h_global[MET][dat]->GetBinContent(r) / h_allMC->GetBinContent(r);
		float TheError = h_global[MET][dat]->GetBinError(r)   / h_allMC->GetBinContent(r);

		Ratio[MET] -> SetBinContent(r, TheRatio);
		Ratio[MET] -> SetBinError  (r, TheError);

	}



	for( int i = 0; i < nvariable; i++ ){

		TCanvas* c = new TCanvas( "canvas_" + variableID[i], variableID[i], 550, 720 );  

		TPad* pad1 = new TPad("pad1", "pad1", 0.05, 0.0, 1.00, 0.3);
		TPad* pad2 = new TPad("pad2", "pad2", 0.05, 0.3, 1.00, 1.0);

		pad1 -> SetTopMargin   (0.08);
		pad1 -> SetBottomMargin(0.35);
		pad1 -> Draw();

		pad2 -> SetTopMargin   (0.08);
		pad2 -> SetBottomMargin(0.02);
		pad2 -> Draw();

		// ---------------------------------------------------------------

		pad2 -> cd();

		pad2 -> SetLogy(); 

		h_global[i][dat] -> SetTitle("");

		h_global[i][dat] -> SetStats(false);   // it has priority over the gStyle->SetOptStats option

 		SetAxis( h_global[i][dat], "", Form("Events / %1.0f GeV", (maxuPara-minuPara)/nbinuPara), 1.5, 1.8 );

		h_global[i][dat] -> SetMinimum(0.1);

		//h_global[i][dat] -> SetXTitle( variableIDfancy[i] + " [GeV]" );
		//h_global[i][dat] -> SetYTitle( Form("Events / %1.0f GeV", (maxuPara-minuPara)/nbinuPara) );

		//h_global[i][dat] -> GetXaxis() ->SetTitleOffset(1.2);
		//h_global[i][dat] -> GetYaxis() ->SetTitleOffset(1.5);

		h_global[i][dat] -> Draw("e"); 

		s_global[i] -> Draw("hist same");

		h_global[i][dat] -> Draw("e same");

		// ---------------------------------------------------------------
		
		TLegend* TheLegend = new TLegend( 0.68, 0.65, 0.88, 0.88 );

		//TheLegend -> SetHeader("processes");

		TheLegend -> SetBorderSize(0);

		TheLegend -> SetTextSize(0.030);

		for( int j = 0; j < nprocess; j++ ){ 

			if(  j == TTGJets  ||  j == ZNuNuGJets40130  ||  j == ZGTo2LG  ||  j == WGToLNuG ) continue;

			( j == data2016D )   ?   TheLegend -> AddEntry( h_global[i][j], processIDfancy[j], "p" )   :   TheLegend -> AddEntry( h_global[i][j], processIDfancy[j], "f" );

		}

		TheLegend -> Draw();

		// ---------------------------------------------------------------

		//DrawLatex( tfont, x, y, tsize, align, text, setndc );

		DrawLatex( 61, 0.100, 0.945, 0.050, 11, "CMS"                                             );

     	 	DrawLatex( 52, 0.205, 0.945, 0.030, 11, "Preliminary"                                     );

		DrawLatex( 42, 0.900, 0.945, 0.050, 31, Form("%.3f fb^{-1} (13TeV, 2016)", TheLuminosity) );


		/*TLatex HeaderLeft, HeaderRight;

		HeaderLeft .SetTextAlign(11);   // left-bottom
		HeaderRight.SetTextAlign(31);	// right-bottom

		HeaderLeft .SetTextSize(0.03);
		HeaderRight.SetTextSize(0.03);

		HeaderLeft .SetTextFont(42);
		HeaderRight.SetTextFont(42);

		HeaderLeft .SetNDC();
		HeaderRight.SetNDC();

		HeaderLeft .DrawLatex ( 0.1, 0.94,       "CMS Preliminary"                              );
		HeaderRight.DrawLatex ( 0.9, 0.94, Form( "%4.2f fb^{-1} (13TeV, 2016)", TheLuminosity ) );*/

		// ---------------------------------------------------------------

		////pad2 -> GetFrame() -> DrawClone();

  		pad2 -> RedrawAxis();

		// ---------------------------------------------------------------

		pad1 -> cd();

		pad1 -> SetLogy();		

		Ratio[i] -> SetTitle("");

		Ratio[i] -> SetStats(false);   // it has priority over the gStyle->SetOptStats option

 		SetAxis( Ratio[i], variableIDfancy[i], "data / MC ", 1.4, 0.75 );

		/*Ratio[i] -> SetXTitle( " "       );
		Ratio[i] -> SetYTitle( "data/MC" );

		Ratio[i] -> GetXaxis() ->SetTitleOffset(1.2);
		Ratio[i] -> GetYaxis() ->SetTitleOffset(1.5);

		Ratio[i] -> GetXaxis()->SetLabelFont(43); 
   		Ratio[i] -> GetXaxis()->SetTitleFont(43);   
   		Ratio[i] -> GetYaxis()->SetLabelFont(43);   
   		Ratio[i] -> GetYaxis()->SetTitleFont(43);  

		Ratio[i] -> GetXaxis()->SetLabelSize(18.); 
   		Ratio[i] -> GetXaxis()->SetTitleSize(18.);   
   		Ratio[i] -> GetYaxis()->SetLabelSize(18.);   
   		Ratio[i] -> GetYaxis()->SetTitleSize(18.);   */

      		Ratio[i] -> GetYaxis() -> SetRangeUser(1e-4, 1e4);

		Ratio[i] -> Draw("ep");

		pad1 -> RedrawAxis();

      		// ---------------------------------------------------------------

		c -> SaveAs( "global/" + variableID[i] + ".pdf" );
		c -> SaveAs( "global/" + variableID[i] + ".png" );

	}

	//allhistos -> Close();
 
}

void OrderFits(){

	if( DoFit == false ) { cout << "   >>   OrderFits deactivated !!!   " << endl;  return; };


	// load histograms

	TFile* allhistos = new TFile( "histograms/" + allhistosReadFrom + ".root", "read" );

	for( int i = 0; i < nvariable; i++ ){

		if( i == MET ) continue; 

		for( int j = 0; j < nprocess; j++ ){

			for( int k = 0; k < nbinpT   ; k++){

				h_resol_pT   [i][j][k] = (TH1F*) allhistos -> Get( Form("h_resol_pT_"    + variableID[i] + "_" + processID[j] + "_%d", k) );

/* remove some day */		if( j != dat )	h_resol_pT[i][j][k] -> Scale( xs[j]/_TotalEntries[j] );

			}

			for( int k = 0; k < nbinsumET; k++){	

				h_resol_sumET[i][j][k] = (TH1F*) allhistos -> Get( Form("h_resol_sumET_" + variableID[i] + "_" + processID[j] + "_%d", k) );

/* remove some day */		if( j != dat )	h_resol_sumET[i][j][k] -> Scale( xs[j]/_TotalEntries[j] );

			}

			for( int k = 0; k < nbinNVtx ; k++){	

				h_resol_NVtx [i][j][k] = (TH1F*) allhistos -> Get( Form("h_resol_NVtx_"  + variableID[i] + "_" + processID[j] + "_%d", k) );

/* remove some day */		if( j != dat )	h_resol_NVtx[i][j][k] -> Scale( xs[j]/_TotalEntries[j] );

			}

		}

	}



	// create all-background template for the fit procedure

	for( int k = 0; k < nbinpT; k++ ){

		h_resol_pT   [parall][nprocess][k] = new TH1F( Form("h_resol_pT_"  + variableID[parall] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara ); 
		h_resol_pT   [transv][nprocess][k] = new TH1F( Form("h_resol_pT_"  + variableID[transv] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara ); 

	}

	for( int k = 0; k < nbinsumET; k++ ){

		h_resol_sumET[parall][nprocess][k] = new TH1F( Form("h_resol_sumET_" + variableID[parall] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara ); 
		h_resol_sumET[transv][nprocess][k] = new TH1F( Form("h_resol_sumET_" + variableID[transv] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara );

	}

	for( int k = 0; k < nbinNVtx; k++ ){

		h_resol_NVtx [parall][nprocess][k] = new TH1F( Form("h_resol_NVtx_"  + variableID[parall] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara ); 
		h_resol_NVtx [transv][nprocess][k] = new TH1F( Form("h_resol_NVtx_"  + variableID[transv] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara );

	}


	for( int i = 0; i < nvariable; i++ ){

	if( i == MET ) continue; 

		for( int j = 2; j < nprocess ; j++ ){

			for( int k = 0; k < nbinpT; k++ ){

				h_resol_pT   [i][nprocess][k] -> Add( h_resol_pT   [i][j][k] ); 

			}

			for( int k = 0; k < nbinsumET; k++ ){

				h_resol_sumET[i][nprocess][k] -> Add( h_resol_sumET[i][j][k] );

			}

			for( int k = 0; k < nbinNVtx; k++ ){

				h_resol_NVtx [i][nprocess][k] -> Add( h_resol_NVtx [i][j][k] );

			}

		}

	}

 
	// fits & final plots 

	for( int i = 0; i < nvariable; i++){
	//for( int i = 0; i < 1; i++){

		if( i == MET ) continue; 


		// ----- resolution: photon pT

		for( int k = 0; k < nbinpT; k++ ){
		//for( int k = 0; k < 1; k++ ){

			xpT[k] = minpT + (1.0*k+0.5)*(maxpT-minpT)/nbinpT;  

			GetResolution( i, pT, k, ypT_A[k],    epT_A[k],    chi2pT_A[k], "GJetsVoigtian" );
			//GetResolution( i, pT, k, ypT_B[k],    epT_B[k],    chi2pT_B[k], "GJetsFromMC"   );

			//cout << xpT[k] << " -- " << ypT_A[k] << " -- " << epT_A[k]<< " -- " << chi2pT_A[k] << endl;

		}
		
		PlotResolution( i, pT     );  cout << "final plot for " << parameterID[pT]    << " done !! " << endl; 


		// ----- resolution: sumEt

		//for( int k = 0; k < nbinsumET; k++ ){

		//	xsumET[k] = minsumET + (1.0*k+0.5)*(maxsumET-minsumET)/nbinsumET;  

		//	xsumET[k] = xsumET[k]/1000;  

		//	GetResolution( i, sumET, k, ysumET_A[k], esumET_A[k], chi2sumET_A[k], "GJetsVoigtian" );
			//GetResolution( i, sumET, k, ysumET_B[k], esumET_B[k], chi2sumET_B[k], "GJetsFromMC"   );

		//}

		//PlotResolution( i, sumET  );  cout << "final plot for " << parameterID[sumET] << " done !! " << endl; 


		// ----- resolution: number of vertices

		for( int k = 0; k < nbinNVtx; k++ ){

			xNVtx[k] = minNVtx + (1.0*k+0.5)*(maxNVtx-minNVtx)/nbinNVtx;

			GetResolution( i, NVtx, k, yNVtx_A[k],  eNVtx_A[k],  chi2NVtx_A[k], "GJetsVoigtian"  );
			//GetResolution( i, NVtx, k, yNVtx_B[k],  eNVtx_B[k],  chi2NVtx_B[k], "GJetsFromMC"    );

		}

		PlotResolution( i, NVtx   );  cout << "final plot for "  << parameterID[NVtx]  << " done !! " << endl;


	}

	cout << "                 " << endl; 
	cout << "exiting OrderFits" << endl; 
	cout << "                 " << endl; 
}




void PlotResolution( int ivar, int parameter ){

	TCanvas* c = new TCanvas("resol_" + variableID[ivar] + "_" + parameterID[parameter], variableID[ivar] + "_" + parameterID[parameter], 600, 600);

	TMultiGraph*  TheMultiGraph = new TMultiGraph();   // if >1 TGraphErrors plotted on the same canvas

	TGraphErrors* TheGraph; 

	if( parameter == pT    ) TheGraph = new TGraphErrors(nbinpT   , xpT   , ypT_A   , 0, epT_A    );
	if( parameter == sumET ) TheGraph = new TGraphErrors(nbinsumET, xsumET, ysumET_A, 0, esumET_A );
	if( parameter == NVtx  ) TheGraph = new TGraphErrors(nbinNVtx , xNVtx , yNVtx_A , 0, eNVtx_A  );
 

	TheGraph -> SetTitle("");

	TheGraph -> GetXaxis() -> SetTitle( parameterIDfancy[parameter] );

	TheGraph -> GetYaxis() -> SetTitle( sigma_variable[ivar] );

	TheGraph -> GetYaxis() -> SetRangeUser( 0.0, 50.0 );

	TheGraph -> SetMarkerStyle( 21 );

	TheGraph -> SetMarkerColor( kBlue );
	
	TheGraph -> SetLineColor( kBlue );

	TheGraph -> Draw("AP");

	/*TLatex tex;
	tex.SetTextAlign(13);
	tex.SetTextSize(0.03);
	tex.SetNDC();
	tex.DrawLatex ( 0.1, 0.95, "CMS Preliminary                     2.318 fb^{-1} (13TeV)" );*/

	DrawLatex( 61, 0.100, 0.945, 0.050, 11, "CMS"                                             );
 	DrawLatex( 52, 0.205, 0.945, 0.030, 11, "Preliminary"                                     );
	DrawLatex( 42, 0.900, 0.945, 0.050, 31, Form("%.3f fb^{-1} (13TeV, 2016)", TheLuminosity) );


	c -> SaveAs( "resol/" + variableID[ivar] + "_" + parameterID[parameter] + ".pdf" ); 
	c -> SaveAs( "resol/" + variableID[ivar] + "_" + parameterID[parameter] + ".png" ); 

}



void GetResolution(int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi, TString GJetsOrigin ){

	int ivar = whichvar; 

	RooRealVar x           ( "variable"   , "variable"          ,  0  , -500, 500 );   // random variable; redefine its range
	RooRealVar Gauss_mean  ( "Gauss_mean" , "Gauss mean"        ,  0  ,  -10, 10  );  
	RooRealVar BW_gamma    ( "BW_gamma"   , "Breit-Wigner gamma",  2.3,    0, 100 );  
	RooRealVar Gauss_sigma ( "Gauss_sigma", "Gauss sigma"       , 10  ,    0, 100 );  


	TH1F* h1;
	TH1F* h2; 
	TH1F* h3;

	int min; int width; 

	if( parameter == pT    ){

		h1 = (TH1F*)  h_resol_pT[ivar][dat  ][ibin] -> Clone();
		h2 = (TH1F*)  h_resol_pT[ivar][QCD  ][ibin] -> Clone();
		h3 = (TH1F*)  h_resol_pT[ivar][GJets][ibin] -> Clone();

		min = minpT; 
		width = (maxpT - minpT)/nbinpT; 

	}

	if( parameter == sumET ){

		h1 = (TH1F*)  h_resol_sumET[ivar][dat  ][ibin] -> Clone();
		h2 = (TH1F*)  h_resol_sumET[ivar][QCD  ][ibin] -> Clone();
		h3 = (TH1F*)  h_resol_sumET[ivar][GJets][ibin] -> Clone();

		min = minsumET; 
		width = (maxsumET - minsumET)/nbinsumET; 

	}

	if( parameter == NVtx  ){

		h1 = (TH1F*)  h_resol_NVtx[ivar][dat  ][ibin] -> Clone();
		h2 = (TH1F*)  h_resol_NVtx[ivar][QCD  ][ibin] -> Clone();
		h3 = (TH1F*)  h_resol_NVtx[ivar][GJets][ibin] -> Clone();

		min = minNVtx; 
		width = (maxNVtx - minNVtx)/nbinNVtx; 

	}

	int low = min + ibin*width; cout << low << endl;
	int up  = low + width; 

	// convert TH1F to RooDataHist
	RooDataHist  data( "data" , "data"      , x,  h1 ); 
	RooDataHist  bckg( "bkg"  , "background", x,  h2 );
	RooDataHist gjets( "gjets", "gammajets" , x,  h3 );


	// convert  RooDataHist into template (PDF)
	RooHistPdf bkgpdf  ( "bkgpdf"  , "bkgpdf"  , x, bckg  );
	RooHistPdf gjetspdf( "gjetspdf", "gjetspdf", x, gjets );

	// pure PDF (do not import signal, but recreate it)
	RooVoigtian voigt ("", "", x, Gauss_mean, BW_gamma, Gauss_sigma ); 

	// fraction 
	RooRealVar fbkg( "fbkg", "bkg fraction", 0.5, 0., 1. );

	// model = (1-f)*voigt + f*bkg 
	RooAddPdf modelA( "modelA", "modelA", RooArgList( voigt   , bkgpdf ), fbkg); 
	RooAddPdf modelB( "modelB", "modelB", RooArgList( gjetspdf, bkgpdf ), fbkg); 


	RooFitResult* result; 

	if( GJetsOrigin == "GJetsVoigtian" ) result = modelA.fitTo( data, Extended(kFALSE), RooFit::Save(kTRUE) );   
	if( GJetsOrigin == "GJetsFromMC"   ) result = modelB.fitTo( data, Extended(kFALSE), RooFit::Save(kTRUE) );   

	// --------  -------  -------  -------  -------  -------  -------  -------  ------- 

	double sigma  = Gauss_sigma.getVal()  ;   double gamma  = BW_gamma.getVal()   ;
  	double esigma = Gauss_sigma.getError();   double egamma = BW_gamma.getError ();

  	double Vss = result -> correlation( Gauss_sigma, Gauss_sigma );  double Vsg = result -> correlation( Gauss_sigma, BW_gamma );
  	double Vgs = result -> correlation( BW_gamma   , Gauss_sigma );  double Vgg = result -> correlation( BW_gamma   , BW_gamma );
  	
 	double FWHM  = GetFWHM     (sigma, gamma                                    );
  	double eFWHM = GetFWHMerror(sigma, gamma, esigma, egamma, Vss, Vsg, Vgs, Vgg);

	resol  =  FWHM/2.3546;
	eresol = eFWHM/2.3546;
 
	// --------  -------  -------  -------  -------  -------  -------  -------  ------- 

	TCanvas* mycanvas = new TCanvas("c", "c", 600, 600);

	RooPlot* frame = x.frame( Title(" ") );

	frame -> GetXaxis() -> SetRangeUser( -200, 200 );

	frame -> GetXaxis() -> SetTitle(variableIDfancy[ivar]);

	data.plotOn( frame );

	modelA.plotOn( frame, Components(bkgpdf), LineColor(kRed)  , LineStyle(kDashed), FillColor(kRed)  , DrawOption("F") );

	modelA.plotOn( frame, Components(voigt) , LineColor(kGreen), LineStyle(kDashed), FillColor(kGreen), DrawOption("L") );

	modelA.plotOn( frame,                     LineColor(kBlue) , LineStyle(kDashed), FillColor(kBlue) , DrawOption("L") );

	frame -> Draw();

	float chi2 = frame -> chiSquare();
	chichi = chi2; 

	/*TLatex mylatex;
	mylatex.SetTextAlign(13);
	mylatex.SetTextSize(0.03);
	mylatex.SetNDC();
	mylatex.DrawLatex ( 0.1, 0.95, Form("CMS Preliminary  2.318 fb^{-1} (13TeV)   #chi^{2} = %5.2f", chi2) );*/

	DrawLatex( 61, 0.100, 0.945, 0.050, 11, "CMS"                                             );
 	DrawLatex( 52, 0.205, 0.945, 0.030, 11, "Preliminary"                                     );
	DrawLatex( 42, 0.900, 0.945, 0.050, 31, Form("%.3f fb^{-1} (13TeV, 2016)", TheLuminosity) );
	DrawLatex( 42, 0.840, 0.900, 0.030, 31, Form("#chi^{2} = %5.2f ", chi2)                   );

	mycanvas -> SaveAs( Form( "fit/" + variableID[ivar] + "_" + parameterID[parameter] + "_%dto%d.pdf", low, up ) );
	mycanvas -> SaveAs( Form( "fit/" + variableID[ivar] + "_" + parameterID[parameter] + "_%dto%d.png", low, up ) );
}

 
double GetFWHM( double sigma, double gamma ){

	double Gauss_FWHM = c * sigma;
  	double BW_FWHM    = 2 * gamma;

  	return a * BW_FWHM + sqrt ( b * pow( BW_FWHM, 2 ) + pow( Gauss_FWHM, 2 ) );

}


double GetFWHMerror( double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg ){

	double Gauss_FWHM = c * sigma;
  	double BW_FWHM    = 2 * gamma;
	double sqroot     = sqrt ( b * pow( BW_FWHM, 2 ) + pow( Gauss_FWHM, 2 ) );

	// partial derivatives FWHM(voigtian) w.r.t sigma and gamma
	double dFWHMdsigma = c * (      Gauss_FWHM / sqroot );
	double dFWHMdgamma = 2 * ( a + b * BW_FWHM / sqroot );

	// esigma * esigma * pow( Vss, 2 ) gives covariance(sigma, sigma) etc
	double p1 = dFWHMdsigma * dFWHMdsigma * esigma * esigma * pow( Vss, 2 );
	double p2 = dFWHMdsigma * dFWHMdgamma * esigma * egamma * pow( Vsg, 2 );
	double p3 = dFWHMdgamma * dFWHMdsigma * egamma * esigma * pow( Vgs, 2 );
	double p4 = dFWHMdgamma * dFWHMdgamma * egamma * egamma * pow( Vgg, 2 );

	return sqrt ( p1 + p2 + p3 + p4 );

}



void SetAxis( TH1* hist, TString xtitle, TString ytitle, Float_t xoffset, Float_t yoffset ){

	gPad -> cd();

	gPad -> Update();

	// See https://root.cern.ch/doc/master/classTAttText.html#T4
	Float_t padw = gPad -> XtoPixel( gPad -> GetX2() );
	Float_t padh = gPad -> YtoPixel( gPad -> GetY1() );

	Float_t size = (padw < padh) ? padw : padh;

	size = 20. / size;   // like this label size is always 20 pixels

	TAxis* xaxis = (TAxis*) hist -> GetXaxis();
	TAxis* yaxis = (TAxis*) hist -> GetYaxis();

	xaxis -> SetTitleOffset(xoffset);
	yaxis -> SetTitleOffset(yoffset);

	xaxis -> SetLabelSize(size);
	yaxis -> SetLabelSize(size);
	xaxis -> SetTitleSize(size);
	yaxis -> SetTitleSize(size);

	xaxis -> SetTitle(xtitle);
	yaxis -> SetTitle(ytitle);

	//yaxis -> CenterTitle();

	gPad -> GetFrame() -> DrawClone();
	gPad -> RedrawAxis();

}



void DrawLatex( Font_t tfont, Float_t x, Float_t y, Float_t tsize, Short_t align, const char* text ){

	TLatex* tl = new TLatex(x, y, text);

	tl -> SetNDC      ( true );
	tl -> SetTextAlign( align);
	tl -> SetTextFont ( tfont);
	tl -> SetTextSize ( tsize);

	tl -> Draw("same");

}



