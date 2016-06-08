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


const float a = 0.5346                ; 
const float b = 0.2166                ; 
const float c = 2 * sqrt( 2 * log(2) ); 

const float maxpT = 70; 
const int width = 5; 
float _up = 0; 

enum{ parall, transv, nvariable };
enum{ dat, QCD, TTGJets, WG, WJets, ZGJets, nprocess}; 

TString processID [nprocess] ; 
TString variableID[nvariable];
float xs[nprocess]; 

TH1F* h_resol[nvariable][nprocess];



double GetFWHM     ( double sigma, double gamma                                                                               );
double GetFWHMerror( double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg );
void   GetHistograms(); 
void GetHistogram( int process );
void GetResolution();

void GetResolutions(){

	//int nbin = maxpT/width; 
 
	//float output[2][nbin];

	//for( int i = 0; i < nbin; i++ ){

	//	_up = ( i + 1 ) * width;

		GetResolution();

	//}

}

void GetResolution(){

	processID[dat    ] = "data"    ; 
	processID[QCD    ] = "QCD"     ; 
	processID[TTGJets] = "TTGJets" ; 
	processID[WG     ] = "WG"      ;
 	processID[WJets  ] = "WJets"   ; 
	processID[ZGJets ] = "ZGJets"  ; 

	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
	xs[QCD    ] = 27990000 + 1712000 + 347700 + 32100 + 6831 + 1207 + 119.9 + 25.24;
	xs[TTGJets] =     3.967;
	xs[WG     ] =   405.271;
	xs[WJets  ] = 61526.7  ;
	//xs[ZGJets ] = ;

	variableID[parall] = "parallel"  ; 
	variableID[transv] = "transverse"; 

	RooRealVar x           ( "variable"   , "variable"          ,  0  , -500, 500 );   // random variable; redefine its range
	RooRealVar Gauss_mean  ( "Gauss_mean" , "Gauss mean"        ,  0  ,  -10, 10  );  
	RooRealVar BW_gamma    ( "BW_gamma"   , "Breit-Wigner gamma",  2.3,    0, 100 );  
	RooRealVar Gauss_sigma ( "Gauss_sigma", "Gauss sigma"       , 10  ,    0, 100 );  

	GetHistograms();

	for( int j = 1; j < 3 ; j++ ){

		h_resol[parall][j] -> Scale( xs[j]/h_resol[parall][j]->Integral() );

		if( j > 1 )  h_resol[parall][QCD] -> Add( h_resol[parall][j] );
	}
 
	// convert TH1F to RooDataHist
	RooDataHist data ( "data", "data"      , x,  h_resol[parall][dat] ); 
	RooDataHist bckg ( "bkg" , "background", x,  h_resol[parall][QCD] ); 

	// convert  RooDataHist into template (PDF)
	RooHistPdf bkgpdf ( "bkgpdf", "bkgpdf", x, bckg );

	// pure PDF (do not import signal, but recreate it)
	RooVoigtian voigt ("", "", x, Gauss_mean, BW_gamma, Gauss_sigma ); 
	
	// fraction 
	RooRealVar fbkg( "fbkg", "bkg fraction", 0.5, 0., 1. );

	// model = (1-f)*voigt + f*bkg 
	RooAddPdf model( "model", "model", RooArgList( voigt, bkgpdf ), fbkg);

	RooFitResult* result = model.fitTo( data, Extended(kFALSE), RooFit::Save(kTRUE) );   

	// --------  -------  -------  -------  -------  -------  -------  -------  ------- 

	double sigma  = Gauss_sigma.getVal()  ;   double gamma  = BW_gamma.getVal()   ;
  	double esigma = Gauss_sigma.getError();   double egamma = BW_gamma.getError ();

  	double Vss = result -> correlation( Gauss_sigma, Gauss_sigma );  double Vsg = result -> correlation( Gauss_sigma, BW_gamma );
  	double Vgs = result -> correlation( BW_gamma   , Gauss_sigma );  double Vgg = result -> correlation( BW_gamma   , BW_gamma );
  	
 	double FWHM  = GetFWHM     (sigma, gamma                                    );
  	double eFWHM = GetFWHMerror(sigma, gamma, esigma, egamma, Vss, Vsg, Vgs, Vgg);

  	//cout << "  -- sigma " <<  sigma << endl;
  	//cout << "  -- gamma " <<  gamma << endl;

	cout << "                                     " << endl; 
	cout << "                                     " << endl; 
	cout << "                                     " << endl; 
	cout << _up << " GeV -> FWHM = "  << FWHM << " +/- " << eFWHM << endl;   
	cout << "                                     " << endl; 
	cout << "                                     " << endl; 
	cout << "                                     " << endl; 

	// --------  -------  -------  -------  -------  -------  -------  -------  ------- 

	RooPlot* frame = x.frame( Title("mytitle") );

        frame -> GetXaxis() -> SetRangeUser( -200, 200 );

	data.plotOn( frame );

        model.plotOn( frame, Components(bkgpdf), LineColor(kRed)  , LineStyle(kDashed), FillColor(kRed)    , DrawOption("F") );

        model.plotOn( frame, Components(voigt) , LineColor(kGreen), LineStyle(kDashed), FillColor(kGreen+1), DrawOption("L") );

	frame -> Draw();
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

void GetHistograms(){

	//for( int j = 0; j < nprocess; j++ ){
	for( int j = 0; j < 3; j++ ){

		GetHistogram( j );

	}

}

void GetHistogram( int process ){

	//TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/latino_" + processID[j] + ".root", "READ" );

	//TTree *tree = (TTree*) file -> Get( "latino" );

	TChain* tree = new TChain("latino");

	TString path = "~/eos/cms/store/group/phys_higgs/cmshww/amassiro/TEST/latino_"; 

	if( process == dat ){

		tree -> Add( path + "SinglePhoton_00__part0.root" );
		tree -> Add( path + "SinglePhoton_00__part1.root" );
		tree -> Add( path + "SinglePhoton_00__part2.root" );
		tree -> Add( path + "SinglePhoton_00__part3.root" );
		tree -> Add( path + "SinglePhoton_00__part4.root" );

		tree -> Add( path + "SinglePhoton_01__part0.root" );
		tree -> Add( path + "SinglePhoton_01__part1.root" );
		tree -> Add( path + "SinglePhoton_01__part2.root" );
		tree -> Add( path + "SinglePhoton_01__part3.root" );
		tree -> Add( path + "SinglePhoton_01__part4.root" );

		tree -> Add( path + "SinglePhoton_02__part0.root" );
		tree -> Add( path + "SinglePhoton_02__part1.root" );
		tree -> Add( path + "SinglePhoton_02__part2.root" );
		tree -> Add( path + "SinglePhoton_02__part3.root" );
		tree -> Add( path + "SinglePhoton_02__part4.root" );	

		tree -> Add( path + "SinglePhoton_03__part0.root" );
		tree -> Add( path + "SinglePhoton_03__part1.root" );
		tree -> Add( path + "SinglePhoton_03__part2.root" );
		tree -> Add( path + "SinglePhoton_03__part3.root" );
		tree -> Add( path + "SinglePhoton_03__part4.root" );

		tree -> Add( path + "SinglePhoton_04__part0.root" );
		tree -> Add( path + "SinglePhoton_04__part1.root" );
		tree -> Add( path + "SinglePhoton_04__part2.root" );
		tree -> Add( path + "SinglePhoton_04__part3.root" );
		tree -> Add( path + "SinglePhoton_04__part4.root" );

		tree -> Add( path + "SinglePhoton_05__part0.root" );
		tree -> Add( path + "SinglePhoton_05__part1.root" );
		tree -> Add( path + "SinglePhoton_05__part2.root" );
		tree -> Add( path + "SinglePhoton_05__part3.root" );
		tree -> Add( path + "SinglePhoton_05__part4.root" );

		tree -> Add( path + "SinglePhoton_06__part0.root" );
		tree -> Add( path + "SinglePhoton_06__part1.root" );
		tree -> Add( path + "SinglePhoton_06__part2.root" );

	}

	else if( process == QCD ){

		tree -> Add( path + "QCD100To200__part0.root"  ); 
		tree -> Add( path + "QCD100To200__part1.root"  ); 
		tree -> Add( path + "QCD100To200__part2.root"  ); 
		tree -> Add( path + "QCD100To200__part3.root"  ); 
		tree -> Add( path + "QCD100To200__part4.root"  ); 

		tree -> Add( path + "QCD200To300__part0.root"  ); 
		tree -> Add( path + "QCD200To300__part1.root"  ); 
		tree -> Add( path + "QCD200To300__part2.root"  ); 
		tree -> Add( path + "QCD200To300__part3.root"  ); 
		tree -> Add( path + "QCD200To300__part4.root"  ); 	

		tree -> Add( path + "QCD700To1000__part0.root"  ); 
		tree -> Add( path + "QCD700To1000__part1.root"  ); 
		tree -> Add( path + "QCD700To1000__part2.root"  ); 
		tree -> Add( path + "QCD700To1000__part3.root"  ); 
		tree -> Add( path + "QCD700To1000__part4.root"  ); 

		tree -> Add( path + "QCD1000To1500__part0.root" ); 
		tree -> Add( path + "QCD1000To1500__part1.root" ); 
		tree -> Add( path + "QCD1000To1500__part2.root" ); 

		tree -> Add( path + "QCD1500To2000__part0.root" ); 
		tree -> Add( path + "QCD1500To2000__part1.root" ); 
		tree -> Add( path + "QCD1500To2000__part2.root" ); 

		tree -> Add( path + "QCD2000ToInf__part0.root"  ); 
		tree -> Add( path + "QCD2000ToInf__part1.root"  ); 

	}

	else if( process == TTGJets ){

		tree -> Add( path + "TTGJets_00__part0.root"  ); 
		tree -> Add( path + "TTGJets_00__part1.root"  ); 
		tree -> Add( path + "TTGJets_00__part2.root"  ); 
		tree -> Add( path + "TTGJets_00__part3.root"  ); 
		tree -> Add( path + "TTGJets_00__part4.root"  ); 

		tree -> Add( path + "TTGJets_01.root"         ); 

	}

	else if( process == WG ){

		tree -> Add( path + "WGToLNuG.root"  ); 

	}

	else if( process == WJets ){

		tree -> Add( path + "WJetsToLNu_00__part0.root"  ); 
		tree -> Add( path + "WJetsToLNu_00__part1.root"  ); 
		tree -> Add( path + "WJetsToLNu_00__part2.root"  ); 
		tree -> Add( path + "WJetsToLNu_00__part3.root"  ); 
		tree -> Add( path + "WJetsToLNu_00__part4.root"  ); 

		tree -> Add( path + "WJetsToLNu_01__part0.root"  ); 
		tree -> Add( path + "WJetsToLNu_01__part1.root"  ); 
		tree -> Add( path + "WJetsToLNu_01__part2.root"  ); 
		tree -> Add( path + "WJetsToLNu_01__part3.root"  ); 
		tree -> Add( path + "WJetsToLNu_01__part4.root"  ); 

		tree -> Add( path + "WJetsToLNu_02__part0.root"  ); 
		tree -> Add( path + "WJetsToLNu_02__part1.root"  ); 
		tree -> Add( path + "WJetsToLNu_02__part2.root"  ); 
		tree -> Add( path + "WJetsToLNu_02__part3.root"  ); 
		tree -> Add( path + "WJetsToLNu_02__part4.root"  ); 

		tree -> Add( path + "WJetsToLNu_03__part0.root"  ); 
		tree -> Add( path + "WJetsToLNu_03__part1.root"  ); 
		tree -> Add( path + "WJetsToLNu_03__part2.root"  ); 
		tree -> Add( path + "WJetsToLNu_03__part3.root"  ); 
		tree -> Add( path + "WJetsToLNu_03__part4.root"  ); 

		tree -> Add( path + "WJetsToLNu_04__part0.root"  ); 
		tree -> Add( path + "WJetsToLNu_04__part1.root"  ); 
		tree -> Add( path + "WJetsToLNu_04__part2.root"  ); 
		tree -> Add( path + "WJetsToLNu_04__part3.root"  ); 
		tree -> Add( path + "WJetsToLNu_04__part4.root"  ); 

	}


	else if( process == ZGJets ){

		tree -> Add( path + "ZNuNuGJets.root"  ); 


	}

	else {

		return; 

	}

	
	float baseW        ; 
	float metPfType1   ;
	float metPfType1Phi; 
	
	vector<float> *std_vector_photon_pt;  std_vector_photon_pt  = 0;
	vector<float> *std_vector_photon_phi; std_vector_photon_phi = 0;


	tree -> SetBranchAddress( "baseW"                , &baseW                 );
	tree -> SetBranchAddress( "metPfType1"           , &metPfType1            );
	tree -> SetBranchAddress( "metPfType1Phi"        , &metPfType1Phi         );
	tree -> SetBranchAddress( "std_vector_photon_pt" , &std_vector_photon_pt  );
	tree -> SetBranchAddress( "std_vector_photon_phi", &std_vector_photon_phi );


	for( int i = 0; i < nvariable; i++ ){

			h_resol[i][process] = new TH1F( "h_resol_" + variableID[i] + "_" + processID[process], "resolution", 400, -200, 200 );

	}


	int nentries = tree -> GetEntries();

	if ( nentries > 3000 ) nentries = 3000; 


	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		tree -> GetEntry(ievt);	

		if ( std_vector_photon_pt->at(0) != 0 ) {

			//if ( std_vector_photon_pt->at(0) < ( _up - width ) || std_vector_photon_pt->at(0) > _up ) continue; 	

			TVector2 qT, ET, uT; 

			ET.SetMagPhi( metPfType1, metPfType1Phi );

			qT.SetMagPhi( std_vector_photon_pt->at(0), std_vector_photon_phi->at(0) );

			uT = -1* ( ET + qT ); 

			TVector2 u_parallel = uT.Proj(qT); //cout << qT.DeltaPhi(u_parallel) << endl; 
			
			float u_parall; ( qT.Mod() > u_parallel.Mod() )   ?   u_parall = (qT+u_parallel).Mod()   :   u_parall = -1* (qT+u_parallel).Mod();   // cout << u_parall << endl;

			float u_transv = sqrt( pow( uT.Mod(), 2) - pow( u_parallel.Mod(), 2) );

				h_resol[parall][process] -> Fill( u_parall ); 
				h_resol[transv][process] -> Fill( u_transv );			

		}

	}


	//cout << _up << endl; 

	//TCanvas* c = new TCanvas( "mycanvas", "mycanvas", 600, 600 );

	//h_resol[parall][process] -> Draw();

}
