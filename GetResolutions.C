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
enum{ dat, GJets, QCD, WJets, ZGJets, ZNuNuGJets40130, ZGTo2LG, WGJets, WGToLNuG, TTGJets, TGJets, nprocess}; 


const bool DoNotEventLoop = false; 

const bool RunOverAllData = true; 
const bool RunOverAllMC   = true; 

      int _TotalEntries[nprocess]; 

const int MaxEntriesData    =   10000;
const int MaxEntriesMC      =     100;

const int MaxEntriesMCStored =  20000; // 20000 (11vii)


const float TheLuminosity = 4.32421;   // fb⁻¹


TString allhistosWriteTo  = "allhistos_12Jul"; 
TString allhistosReadFrom = "allhistos_11Jul"; 


const bool SkipFit = true; 

const bool GJetsMC = false;

const float a = 0.5346                ; 
const float b = 0.2166                ; 
const float c = 2 * sqrt( 2 * log(2) ); 

const int HowManyBkgs = 3;   // provisional  // ideally, HowManyBkgs = 5

float          xs   [nprocess]; 
TString processID   [nprocess];
TColor  processcolor[nprocess];
 
TString variableID     [nvariable];
TString variableIDfancy[nvariable];
TString yTitle         [nvariable];

TString parameterID[nparameter];

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
TH1F*    h_global     [aa][bb]           ;
THStack* s_global     [aa]               ;
TH1F*    h_resol_pT   [aa][bb][nbinpT   ];
TH1F*    h_resol_sumET[aa][bb][nbinsumET];
TH1F*    h_resol_NVtx [aa][bb][nbinNVtx ]; 
// -----------------------------------------

double GetFWHM      ( double sigma, double gamma                                                                               );
double GetFWHMerror ( double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg );
void  FillHistograms(                                                                                                          );
void  FillHistogram ( int process                                                                                              );
void GetResolution  ( int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi                        );
void GetGlobalPlots (                                                                                                          );

void GetResolutions(){

	processID[dat            ] = "data"           ; 
	processID[GJets          ] = "GJets"          ;  //processcolor[GJets  ] = kYellow ;
	processID[QCD            ] = "QCD"            ;  //processcolor[QCD    ] = kSpring ;
 	processID[WJets          ] = "WJets"          ;  //processcolor[WJets  ] = kViolet ;
	processID[ZGJets         ] = "ZGJets"         ;  //processcolor[ZGJets ] = kMagenta;
	processID[ZNuNuGJets40130] = "ZNuNuGJets40130";  
	processID[ZGTo2LG        ] = "ZGTo2LG"        ;  
	processID[WGJets         ] = "WGJets"         ;  //processcolor[WG     ] = kCyan   ;
	processID[WGToLNuG       ] = "WGToLNuG"       ;  
	processID[TTGJets        ] = "TTGJets"        ;  //processcolor[TTGJets] = kRed    ;
	processID[TGJets         ] = "TGJets"         ;  


	// --- x-sections ---------------------------------------------------------------------------

	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
	// http://mvesterb.web.cern.ch/mvesterb/met/things/samples.dat

	// ------------------------------------------------------------------------------------------

	xs[GJets          ] =    20790      + 9238 + 2305 + 274.4 + 93.46; 
	xs[QCD            ] = 27990000      + 1712000 + 347700 + 32100 + 6831 + 1207 + 119.9 + 25.24;
	xs[WJets          ] =    61526.7   ;
	xs[ZGJets         ] =        0.1903;
	xs[ZNuNuGJets40130] =        2.816 ;
	xs[ZGTo2LG        ] =      117.864 ;
	xs[WGJets         ] =        0.6637; // 405.271;
	xs[WGToLNuG       ] =      585.800 ;
	xs[TTGJets        ] =        3.967 ;
	xs[TGJets         ] =        3.967 ;

	// ------------------------------------------------------------------------------------------

	variableID[parall] = "parallel"  ; 
	variableID[transv] = "transverse"; 
	variableID[MET   ] = "MET"       ;

	variableIDfancy[parall] = "u_{||}+#gamma_{T}"; 
	variableIDfancy[transv] = "u_{#perp}"        ;
	variableIDfancy[MET   ] = "E_{T}^{miss}"     ;

	yTitle[parall] = "#sigma(u_{||}) [GeV]"   ;
	yTitle[transv] = "#sigma(u_{#perp}) [GeV]";

	parameterID[pT   ] = "pT"   ;
	parameterID[sumET] = "sumET";
	parameterID[NVtx ] = "NVtx" ;

	FillHistograms(); 

	cout << "                          " << endl; 
	cout << "                          " << endl;
	cout << "    histograms filled !!! " << endl;
	cout << "                          " << endl;
	cout << "                          " << endl;

	//GetGlobalPlots();

	if( SkipFit == true ) return;

	//for( int j = 2; j < HowManyBkgs+2 ; j++ ){
	for( int j = 2; j < nprocess ; j++ ){

		for( int i = 0; i < nbinpT; i++ ){

			h_resol_pT[parall][j][i] -> Scale( xs[j]/h_resol_pT[parall][j][i]->Integral() );
			h_resol_pT[transv][j][i] -> Scale( xs[j]/h_resol_pT[transv][j][i]->Integral() );

			if( j > 2 )  h_resol_pT[parall][QCD][i] -> Add( h_resol_pT[parall][j][i] );

		}

		for( int i = 0; i < nbinsumET; i++ ){

			h_resol_sumET[parall][j][i] -> Scale( xs[j]/h_resol_sumET[parall][j][i]->Integral() );
			h_resol_sumET[transv][j][i] -> Scale( xs[j]/h_resol_sumET[transv][j][i]->Integral() );

			if( j > 2 )  h_resol_sumET[parall][QCD][i] -> Add( h_resol_sumET[parall][j][i] );

		}

		for( int i = 0; i < nbinNVtx; i++ ){

			h_resol_NVtx[parall][j][i] -> Scale( xs[j]/h_resol_NVtx[parall][j][i]->Integral() );
			h_resol_NVtx[transv][j][i] -> Scale( xs[j]/h_resol_NVtx[transv][j][i]->Integral() );

			if( j > 2 )  h_resol_NVtx[parall][QCD][i] -> Add( h_resol_NVtx[parall][j][i] );

		}

	}

	float xpT   [nbinpT   ]; float ypT   [nbinpT   ]; float epT   [nbinpT   ]; float chi2pT   [nbinpT   ];
	float xsumET[nbinsumET]; float ysumET[nbinsumET]; float esumET[nbinsumET]; float chi2sumET[nbinsumET];
	float xNVtx [nbinNVtx ]; float yNVtx [nbinNVtx ]; float eNVtx [nbinNVtx ]; float chi2NVtx [nbinNVtx ];

	for( int k = 0; k < 1; k++){

		for( int i = 0; i < nbinpT; i++ ){

			xpT[i] = minpT + (1.0*i+0.5)*(maxpT-minpT)/nbinpT;  //cout << xpT[i] << endl;

			GetResolution( k, 1, i, ypT[i],    epT[i],    chi2pT[i]    );

		}

		for( int i = 0; i < nbinsumET; i++ ){

			xsumET[i] = minsumET + (1.0*i+0.5)*(maxsumET-minsumET)/nbinsumET;  
			xsumET[i] = xsumET[i]/1000;  

			GetResolution( k, 2, i, ysumET[i], esumET[i], chi2sumET[i] );

		}

		for( int i = 0; i < nbinNVtx; i++ ){

			xNVtx[i] = minNVtx + (1.0*i+0.5)*(maxNVtx-minNVtx)/nbinNVtx;

			GetResolution( k, 3, i, yNVtx[i],  eNVtx[i],  chi2NVtx[i]  );

		}



		TGraphErrors* graphpT    = new TGraphErrors(nbinpT   , xpT   , ypT   , 0, epT    );
		TGraphErrors* graphsumET = new TGraphErrors(nbinsumET, xsumET, ysumET, 0, esumET );
		TGraphErrors* graphNVtx  = new TGraphErrors(nbinNVtx , xNVtx , yNVtx , 0, eNVtx  );

		graphpT    -> SetTitle("");
		graphsumET -> SetTitle("");
		graphNVtx  -> SetTitle("");

		graphpT    -> GetXaxis()->SetTitle("#gamma q_{T} [GeV]");
		graphsumET -> GetXaxis()->SetTitle("#Sigma E_{T} [TeV]");
		graphNVtx  -> GetXaxis()->SetTitle("number of vertices");

		graphpT    -> GetYaxis()->SetTitle(yTitle[k]);
		graphsumET -> GetYaxis()->SetTitle(yTitle[k]);
		graphNVtx  -> GetYaxis()->SetTitle(yTitle[k]);

		graphpT    -> GetYaxis()->SetRangeUser(0.0, 40.0);

		graphpT    -> SetMarkerStyle(21);
		graphsumET -> SetMarkerStyle(21);
		graphNVtx  -> SetMarkerStyle(21);

		graphpT    -> SetMarkerColor(kGreen+3);
		graphsumET -> SetMarkerColor(kGreen+3);
		graphNVtx  -> SetMarkerColor(kGreen+3);

		graphpT    -> SetLineColor(kGreen+3);
		graphsumET -> SetLineColor(kGreen+3);
		graphNVtx  -> SetLineColor(kGreen+3);

		TLatex tex;
		tex.SetTextAlign(13);
		tex.SetTextSize(0.03);
		tex.SetNDC();

		TCanvas* cpT = new TCanvas("cpT", "cpT", 600, 600);

		graphpT -> Draw("AP");

		tex.DrawLatex ( 0.1, 0.95, "CMS Preliminary                     2.318 fb^{-1} (13TeV)" );

		cpT -> SaveAs("resol/" + variableID[k] + "_pT.pdf"); 
		cpT -> SaveAs("resol/" + variableID[k] + "_pT.png"); 


		TCanvas* csumET = new TCanvas("csumET", "csumET", 600, 600);

		graphsumET -> Draw("AP");

		tex.DrawLatex ( 0.1, 0.95, "CMS Preliminary                     2.318 fb^{-1} (13TeV)" );

		csumET -> SaveAs("resol/" + variableID[k] + "_sumET.pdf");
		csumET -> SaveAs("resol/" + variableID[k] + "_sumET.png");


		TCanvas* cNVtx = new TCanvas("cNVtx", "cNVtx", 600, 600);

		graphNVtx -> Draw("AP");

		tex.DrawLatex ( 0.1, 0.95, "CMS Preliminary                     2.318 fb^{-1} (13TeV)" );

		cNVtx -> SaveAs("resol/" + variableID[k] + "_NVtx.pdf");
		cNVtx -> SaveAs("resol/" + variableID[k] + "_NVtx.png");



	}

}



void FillHistograms(){

	for( int j = 0; j < nprocess; j++ ){
	//for( int j = 0; j < HowManyBkgs+2; j++ ){

		cout << "process: " << processID[j] << endl; 

		FillHistogram( j );

	}

	cout << "the writing step begins..." << endl; 

	if( DoNotEventLoop == false ){

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

	}

}



void FillHistogram( int process ){

	TChain* tree = new TChain("METtree");

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

	}



	int   ngamma                            ;
	int   gamma_idCutBased                  ;
	float gamma_r9                          ;	
	float gamma_pt                          ;
	float gamma_eta                         ;
	float gamma_phi                         ;

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

	if ( process == dat  &&  RunOverAllData == false ) nentries = MaxEntriesData; 
		
	if ( process != dat  &&  RunOverAllMC   == false ) nentries = MaxEntriesMC  ;
	  
	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		if( ievt%10000 == 0 )  cout << "  >> 10k more... " << endl;

		tree -> GetEntry(ievt);	//cout << " ievt = " << ievt << endl;

		if ( ngamma < 1                       ) continue;
		if ( gamma_idCutBased != 3            ) continue; 
		if ( gamma_r9 < 0.9 || gamma_r9 > 1.0 ) continue; 
		if ( gamma_pt < 50                    ) continue; 
		if ( gamma_eta > 1.5                  ) continue; 

		if ( nLepGood10 > 0                   ) continue; 
		
		// MET filters
		if ( HBHENoiseFilter                    != 1 ) continue; 
		if ( HBHENoiseIsoFilter                 != 1 ) continue; 
		if ( CSCTightHalo2015Filter             != 1 ) continue; 
		if ( EcalDeadCellTriggerPrimitiveFilter != 1 ) continue; 
		if ( goodVertices                       != 1 ) continue; 
		if ( eeBadScFilter                      != 1 ) continue; 		

		// triggers
		if( process == dat ){

			//cout << HLT_Photon30 << "  -- " << HLT_Photon50 << "  -- " << HLT_Photon75 << "  -- " << HLT_Photon90 << "  -- " << HLT_Photon120 << "  -- " << endl;

			if( HLT_Photon30 == 0  &&  HLT_Photon50 == 0  &&  HLT_Photon75 == 0  &&  HLT_Photon90 == 0  &&  HLT_Photon120 == 0 ) continue;

		}



		float eventW = 1.0 ;
 
		if( process != dat ){

			eventW *= puWeight ;

			eventW *= genWeight;

		}


		if( gamma_pt < minpT || gamma_pt >= maxpT || met_sumEt < minsumET || met_sumEt >= maxsumET || nVert < minNVtx || nVert >= maxNVtx) continue;

		//cout << ievt << " -- " << ngamma << " -- " << gamma_pt << " -- " << met_sumEt << " -- " << nVert << endl;

		int l = floor(    nbinpT    * (gamma_pt  - minpT   ) / (maxpT    - minpT   )    ); //cout << l << endl;
		int m = floor(    nbinsumET * (met_sumEt - minsumET) / (maxsumET - minsumET)    ); //cout << m << endl;
		int n = floor(    nbinNVtx  * (nVert     - minNVtx ) / (maxNVtx  - minNVtx )    ); //cout << n << endl;


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

	//TCanvas* c = new TCanvas( "mycanvas", "mycanvas", 600, 600 );

	//h_resol[parall][process] -> Draw();

}



void GetGlobalPlots(){

	//if( RunOverAllData == false ) { cout << "   >>   RunOverAllData is deactivated !!!   " << endl;  return; }

	TFile* allhistos = new TFile( "histograms/" + allhistosReadFrom + ".root", "read" );

	for( int i = 0; i < nvariable; i++ ){

		for( int j = 0; j < nprocess; j++ ){

			h_global[i][j] = (TH1F*) allhistos -> Get("h_global_" + variableID[i] + "_" + processID[j]);

		}

	
		h_global[i][dat] -> SetMarkerStyle(20);

		h_global[i][dat] -> SetMarkerColor(kBlack);

		h_global[i][dat] -> SetLineColor(kBlack);

	}


	float TheCurrentLuminosity; 

	( RunOverAllData == true )   ?   TheCurrentLuminosity  =  TheLuminosity   :   TheCurrentLuminosity  =  TheLuminosity * MaxEntriesData/_TotalEntries[dat]  ; 


	float weight[nprocess]; 

	for( int j = 1; j < nprocess; j++ ){

		( RunOverAllMC == true )   ?   weight[j] = xs[j]/_TotalEntries[j]   :    weight[j] = xs[j]/MaxEntriesMCStored  ; 		

		h_global[parall][j] -> Scale( weight[j]*TheCurrentLuminosity );  
		h_global[transv][j] -> Scale( weight[j]*TheCurrentLuminosity ); 

		h_global[parall][j] -> SetFillColor( j+1 );
		h_global[transv][j] -> SetFillColor( j+1 );	

		//h_global[parall][j] -> SetMarkerStyle( 21 );
		//h_global[transv][j] -> SetMarkerStyle( 21 );

		//h_global[parall][j] -> SetMarkerColor( j+1 );
		//h_global[transv][j] -> SetMarkerColor( j+1 );


	}


	s_global[parall]  = new THStack( variableID[parall], variableID[parall] );
	s_global[transv]  = new THStack( variableID[transv], variableID[transv] );


	for( int i = 0; i < nvariable; i++ ){

		//s_global[i] -> Add( h_global[i][TTGJets        ] );
		s_global[i] -> Add( h_global[i][TGJets         ] );
		s_global[i] -> Add( h_global[i][ZGJets         ] );
		s_global[i] -> Add( h_global[i][ZNuNuGJets40130] );
		s_global[i] -> Add( h_global[i][ZGTo2LG        ] );
		s_global[i] -> Add( h_global[i][WJets          ] );
		s_global[i] -> Add( h_global[i][WGJets         ] );
		s_global[i] -> Add( h_global[i][WGToLNuG       ] );
		s_global[i] -> Add( h_global[i][QCD            ] );
		s_global[i] -> Add( h_global[i][GJets          ] );

	}

	
	TCanvas* c[nvariable]; 

	for( int i = 0; i < nvariable; i++ ){

		c[i] = new TCanvas( "canvas_" + variableID[i], "one canvas", 600, 600 );

		c[i] -> SetLogy(); 

		h_global[i][dat] -> SetTitle("");

		h_global[i][dat] -> SetStats(false);   // it has priority over the gStyle->SetOptStats option

		h_global[i][dat] -> SetXTitle( variableIDfancy[i] + " [GeV]" );
		h_global[i][dat] -> SetYTitle( Form("Events / %1.0f GeV", (maxuPara-minuPara)/nbinuPara) );

		h_global[i][dat] -> GetXaxis() ->SetTitleOffset(1.2);
		h_global[i][dat] -> GetYaxis() ->SetTitleOffset(1.5);

		h_global[i][dat] -> Draw("E1"); 

		s_global[i] -> Draw("hist same");

		h_global[i][dat] -> Draw("E1 same");

		// ---------------------------------------------------------------
		
		TLegend* TheLegend = new TLegend( 0.65, 0.65, 0.88, 0.88 );

		//TheLegend -> SetHeader("processes");

		TheLegend -> SetBorderSize(0);

		TheLegend -> SetTextSize(0.025);

		for( int j = 0; j < nprocess; j++ ){ 

			( j == dat )   ?   TheLegend -> AddEntry( h_global[i][j], processID[j], "p" )   :   TheLegend -> AddEntry( h_global[i][j], processID[j], "f" );

		}

		TheLegend -> Draw();

		// ---------------------------------------------------------------

		TLatex HeaderLeft, HeaderRight;

		HeaderLeft .SetTextAlign(11);   // left-bottom
		HeaderRight.SetTextAlign(31);	// right-bottom

		HeaderLeft .SetTextSize(0.03);
		HeaderRight.SetTextSize(0.03);

		HeaderLeft .SetTextFont(42);
		HeaderRight.SetTextFont(42);

		HeaderLeft .SetNDC();
		HeaderRight.SetNDC();

		HeaderLeft .DrawLatex ( 0.1, 0.92,       "CMS Preliminary"                              );
		HeaderRight.DrawLatex ( 0.9, 0.92, Form( "%4.2f fb^{-1} (13TeV, 2016)", TheLuminosity ) );

		// ---------------------------------------------------------------

		c[i] -> SaveAs( "global/" + variableID[i] + ".pdf" );
		c[i] -> SaveAs( "global/" + variableID[i] + ".png" );

	}

	//allhistos -> Close();
 
}



void GetResolution(int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi){

	int k = whichvar; 


	RooRealVar x           ( "variable"   , "variable"          ,  0  , -500, 500 );   // random variable; redefine its range
	RooRealVar Gauss_mean  ( "Gauss_mean" , "Gauss mean"        ,  0  ,  -10, 10  );  
	RooRealVar BW_gamma    ( "BW_gamma"   , "Breit-Wigner gamma",  2.3,    0, 100 );  
	RooRealVar Gauss_sigma ( "Gauss_sigma", "Gauss sigma"       , 10  ,    0, 100 );  


	TH1F* h1;
	TH1F* h2; 
	TH1F* h3;

	int min; int width; 

	if( parameter == 1 ){

		h1 = (TH1F*)  h_resol_pT[k][dat  ][ibin] -> Clone();
		h2 = (TH1F*)  h_resol_pT[k][QCD  ][ibin] -> Clone();
		h3 = (TH1F*)  h_resol_pT[k][GJets][ibin] -> Clone();

		min = minpT; 
		width = (maxpT - minpT)/nbinpT; 

	}

	if( parameter == 2 ){

		h1 = (TH1F*)  h_resol_sumET[k][dat  ][ibin] -> Clone();
		h2 = (TH1F*)  h_resol_sumET[k][QCD  ][ibin] -> Clone();
		h3 = (TH1F*)  h_resol_sumET[k][GJets][ibin] -> Clone();

		min = minsumET; 
		width = (maxsumET - minsumET)/nbinsumET; 

	}

	if( parameter == 3 ){

		h1 = (TH1F*)  h_resol_NVtx[k][dat  ][ibin] -> Clone();
		h2 = (TH1F*)  h_resol_NVtx[k][QCD  ][ibin] -> Clone();
		h3 = (TH1F*)  h_resol_NVtx[k][GJets][ibin] -> Clone();

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

	if( GJetsMC == false ) result = modelA.fitTo( data, Extended(kFALSE), RooFit::Save(kTRUE) );   
	if( GJetsMC == true  ) result = modelB.fitTo( data, Extended(kFALSE), RooFit::Save(kTRUE) );   

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

	frame -> GetXaxis() -> SetTitle(variableIDfancy[k]);

	data.plotOn( frame );

	modelA.plotOn( frame, Components(bkgpdf), LineColor(kRed)  , LineStyle(kDashed), FillColor(kRed)  , DrawOption("F") );

	modelA.plotOn( frame, Components(voigt) , LineColor(kGreen), LineStyle(kDashed), FillColor(kGreen), DrawOption("L") );

	modelA.plotOn( frame,                     LineColor(kBlue) , LineStyle(kDashed), FillColor(kBlue) , DrawOption("L") );

	frame -> Draw();

	float chi2 = frame -> chiSquare();
	chichi = chi2; 

	TLatex mylatex;
	mylatex.SetTextAlign(13);
	mylatex.SetTextSize(0.03);
	mylatex.SetNDC();
	mylatex.DrawLatex ( 0.1, 0.95, Form("CMS Preliminary  2.318 fb^{-1} (13TeV)   #chi^{2} = %5.2f", chi2) );

	mycanvas -> SaveAs( Form( "fit/" + variableID[k] + "_" + parameterID[parameter-1] + "_%dto%d.pdf", low, up ) );
	mycanvas -> SaveAs( Form( "fit/" + variableID[k] + "_" + parameterID[parameter-1] + "_%dto%d.png", low, up ) );
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





