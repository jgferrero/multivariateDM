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

enum{ parall, transv, nvariable };
enum{ pT, sumET, NVtx, nparameter};
enum{ dat, GJets, QCD, TTGJets, WG, WJets, ZGJets, nprocess}; 
const int HowManyBkgs = 2;   // provisional

float          xs[nprocess]; 
TString processID[nprocess];
 
TString variableID     [nvariable];
TString variableIDfancy[nvariable];
TString yTitle         [nvariable];

TString parameterID[nparameter];

float maxpT    = 800.; float minpT    = 0.; const int nbinpT    = 20; float _uppT   ;
float maxsumET =3000.; float minsumET = 0.; const int nbinsumET = 30; float _upsumET;  
float maxNVtx  =  40.; float minNVtx  = 0.; const int nbinNVtx  =  8; float _upNVtx ;

/*float maxpT    = 200.; float minpT    = 0.; const int nbinpT    =  5; float _uppT   ;
float maxsumET = 300.; float minsumET = 0.; const int nbinsumET =  5; float _upsumET;  
float maxNVtx  =  20.; float minNVtx  = 0.; const int nbinNVtx  =  5; float _upNVtx ;*/

const int aa = nvariable;   // just for next block
const int bb = nprocess ;   // just for next block

TH1F* h_resol_pT   [aa][bb][nbinpT   ];
TH1F* h_resol_sumET[aa][bb][nbinsumET];
TH1F* h_resol_NVtx [aa][bb][nbinNVtx ]; 

const int GJetsMC = false;

double GetFWHM      ( double sigma, double gamma                                                                               );
double GetFWHMerror ( double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg );
void   GetHistograms(                                                                                                          );
void   GetHistogram ( int process                                                                                              );
void GetResolution  (int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi                         );

void GetResolutions(){

	processID[dat    ] = "data"    ; 
	processID[GJets  ] = "GJets"   ; 
	processID[QCD    ] = "QCD"     ; 
	processID[TTGJets] = "TTGJets" ; 
	processID[WG     ] = "WG"      ;
 	processID[WJets  ] = "WJets"   ; 
	processID[ZGJets ] = "ZGJets"  ; 

	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
	xs[GJets  ] = 20790 + 9238 + 2305 + 274.4 + 93.46; 
	xs[QCD    ] = 27990000 + 1712000 + 347700 + 32100 + 6831 + 1207 + 119.9 + 25.24;
	xs[TTGJets] =     3.967;
	xs[WG     ] =   405.271;
	xs[WJets  ] = 61526.7  ;
	//xs[ZGJets ] = ? ;

	variableID[parall] = "parallel"  ; 
	variableID[transv] = "transverse"; 

	variableIDfancy[parall] = "u_{||}+q_{T}"; 
	variableIDfancy[transv] = "u_{#perp}"   ; 

	yTitle[parall] = "#sigma(u_{||}) [GeV]"   ;
	yTitle[transv] = "#sigma(u_{#perp}) [GeV]";

	parameterID[pT   ] = "pT"   ;
	parameterID[sumET] = "sumET";
	parameterID[NVtx ] = "NVtx" ;

	GetHistograms();

	cout << "    histograms loaded !!! " << endl;

	for( int j = 2; j < HowManyBkgs+2 ; j++ ){

		for( int i = 0; i < nbinpT; i++ ){

			h_resol_pT[parall][j][i] -> Scale( xs[j]/h_resol_pT[parall][j][i]->Integral() );
			h_resol_pT[transv][j][i] -> Scale( xs[j]/h_resol_pT[transv][j][i]->Integral() );

			//if( j > 2 )  h_resol[parall][QCD][] -> Add( h_resol[parall][j][] );

		}

		for( int i = 0; i < nbinsumET; i++ ){

			h_resol_sumET[parall][j][i] -> Scale( xs[j]/h_resol_sumET[parall][j][i]->Integral() );
			h_resol_sumET[transv][j][i] -> Scale( xs[j]/h_resol_sumET[transv][j][i]->Integral() );

			//if( j > 2 )  h_resol[parall][QCD][] -> Add( h_resol[parall][j][] );

		}

		for( int i = 0; i < nbinNVtx; i++ ){

			h_resol_NVtx[parall][j][i] -> Scale( xs[j]/h_resol_NVtx[parall][j][i]->Integral() );
			h_resol_NVtx[transv][j][i] -> Scale( xs[j]/h_resol_NVtx[transv][j][i]->Integral() );

			//if( j > 2 )  h_resol[parall][QCD][] -> Add( h_resol[parall][j][] );

		}

	}



	float outputpT   [4][nbinpT   ]; float xpT   [nbinpT   ]; float ypT   [nbinpT   ]; float epT   [nbinpT   ]; float chi2pT   [nbinpT   ];
	float outputsumET[4][nbinsumET]; float xsumET[nbinsumET]; float ysumET[nbinsumET]; float esumET[nbinsumET]; float chi2sumET[nbinsumET];
	float outputNVtx [4][nbinNVtx ]; float xNVtx [nbinNVtx ]; float yNVtx [nbinNVtx ]; float eNVtx [nbinNVtx ]; float chi2NVtx [nbinNVtx ];

	for( int k = 0; k < nvariable; k++){

		for( int i = 0; i < nbinpT; i++ ){

			xpT[i] = minpT + (1.0*i+0.5)*(maxpT-minpT)/nbinpT;  cout << xpT[i] << endl;

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
	RooAddPdf modelA( "modelA", "modelB", RooArgList( voigt   , bkgpdf ), fbkg); 
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



void GetHistograms(){

	//for( int j = 0; j < nprocess; j++ ){
	for( int j = 0; j < HowManyBkgs+2; j++ ){

		GetHistogram( j );

		//cout << j << " success  " << endl;

	}

}


void GetHistogram( int process ){

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

	else if( process == GJets){

		tree -> Add( path + "GJets40To100__part1.root"  ); 
		tree -> Add( path + "GJets40To100__part2.root"  ); 

		tree -> Add( path + "GJets100To200__part0.root" );
		tree -> Add( path + "GJets100To200__part1.root" );
		tree -> Add( path + "GJets100To200__part2.root" );
 
		tree -> Add( path + "GJets200To400__part0.root" ); 
		tree -> Add( path + "GJets200To400__part1.root" ); 
		tree -> Add( path + "GJets200To400__part2.root" ); 
		tree -> Add( path + "GJets200To400__part3.root" ); 
		tree -> Add( path + "GJets200To400__part4.root" ); 

		tree -> Add( path + "GJets400To600__part0.root" ); 
		tree -> Add( path + "GJets400To600__part1.root" ); 

		tree -> Add( path + "GJets600ToInf__part0.root" ); 
		tree -> Add( path + "GJets600ToInf__part1.root" ); 

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


	float metPfType1     ;
	float metPfType1Phi  ;
	float metPfType1SumEt;
	float nvtx           ;
	
	vector<float> *std_vector_photon_pt;  std_vector_photon_pt  = 0;
	vector<float> *std_vector_photon_phi; std_vector_photon_phi = 0;


	tree -> SetBranchAddress( "metPfType1"           , &metPfType1            );
	tree -> SetBranchAddress( "metPfType1Phi"        , &metPfType1Phi         );
	tree -> SetBranchAddress( "std_vector_photon_pt" , &std_vector_photon_pt  );
	tree -> SetBranchAddress( "std_vector_photon_phi", &std_vector_photon_phi );
	tree -> SetBranchAddress( "metPfType1SumEt"      , &metPfType1SumEt       );
	tree -> SetBranchAddress( "nvtx"                 , &nvtx                  );


	for( int i = 0; i < nvariable; i++ ){

		for( int k = 0; k < nbinpT; k++ ){

			h_resol_pT   [i][process][k] = new TH1F( Form("h_resol_pT_"    + variableID[i] + "_" + processID[process] + "_%d", k), "resolution pT"        , 80, -200, 200 );

		}

		for( int k = 0; k < nbinsumET; k++ ){

			h_resol_sumET[i][process][k] = new TH1F( Form("h_resol_sumET_" + variableID[i] + "_" + processID[process] + "_%d", k), "resolution sum ET"    , 80, -200, 200 );

		}

		for( int k = 0; k < nbinNVtx; k++ ){

			h_resol_NVtx [i][process][k] = new TH1F( Form("h_resol_NVtx_"  + variableID[i] + "_" + processID[process] + "_%d", k), "resolution num of vtx", 80, -200, 200 );

		}

	}

	int nentries = tree -> GetEntries(); 

	if ( process == dat                    ) nentries = 200000; 
	if ( process != dat && nentries > 50000) nentries =  50000; 

	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		tree -> GetEntry(ievt);	

		float photonpT = std_vector_photon_pt->at(0);

		if ( photonpT > 0 ) {

			if( photonpT >= maxpT || metPfType1SumEt >= maxsumET || nvtx >= maxNVtx) continue;

			int l = floor(    nbinpT    * (photonpT        - minpT   ) / (maxpT    - minpT   )    ); //cout << l << endl;
			int m = floor(    nbinsumET * (metPfType1SumEt - minsumET) / (maxsumET - minsumET)    ); //cout << m << endl;
			int n = floor(    nbinNVtx  * (nvtx            - minNVtx ) / (maxNVtx  - minNVtx )    ); //cout << n << endl;

			TVector2 qT, ET, uT; 

			ET.SetMagPhi( metPfType1, metPfType1Phi );

			qT.SetMagPhi( std_vector_photon_pt->at(0), std_vector_photon_phi->at(0) );

			uT = -1* ( ET + qT ); 

			TVector2 u_parallel = uT.Proj(qT); //cout << qT.DeltaPhi(u_parallel) << endl; 
			
			float u_parall; ( qT.Mod() > u_parallel.Mod() )   ?   u_parall = (qT+u_parallel).Mod()                                    :   u_parall = -1* (qT+u_parallel).Mod();   // cout << u_parall << endl;

			float u_transv; ( qT.Mod() > u_parallel.Mod() )   ?   u_transv = sqrt( pow( uT.Mod(), 2) - pow( u_parallel.Mod(), 2) )    :   u_transv = -1* sqrt( pow( uT.Mod(), 2) - pow( u_parallel.Mod(), 2) );

			h_resol_pT   [parall][process][l] -> Fill( u_parall ); 
			h_resol_pT   [transv][process][l] -> Fill( u_transv );	

			h_resol_sumET[parall][process][m] -> Fill( u_parall ); 
			h_resol_sumET[transv][process][m] -> Fill( u_transv );

			h_resol_NVtx [parall][process][n] -> Fill( u_parall );
			h_resol_NVtx [transv][process][n] -> Fill( u_transv );	


		}


	}

	//TCanvas* c = new TCanvas( "mycanvas", "mycanvas", 600, 600 );

	//h_resol[parall][process] -> Draw();

}

