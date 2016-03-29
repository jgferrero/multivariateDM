#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "TBranch.h"

#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "TPython.h"   // linacre
#include "StopTreeLooper.h"

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


enum {tt, WW, bk, dt};

enum {data, bkg};

enum {dist3, dist2, dist1, dist, prox};
//enum {dist, prox};

enum {njet, nbjet};


const int        nDM = 8    ;   // DM samples
const float        L = 0.553;   // luminosity
const int       nbin = 40   ;   // MVA output histogram
const float      inf = -0.5 ;   // MVA output histogram
const float      sup = 1.5  ;   // MVA output histogram
const int        nDD = 4    ;   // number datadriven estimated bkgs plus two: the rest of the bkgs (merged) & data
const int        nSB = 2    ;   // number of sidebands
const int       nBox = 2    ;   // number of boxes per CR (equal to the number of datadriven estimated bkgs)
const int  bandwidth = 2    ;   // side band width
const int        nCf = 6    ;   // number of boxes configurations


TString DM      [nDM]                       ;   // DM samples ids
TH1F* MVAdata   [nDM]                       ;
TH1F* MVAsig    [nDM]                       ;
TH1F* MVAbkg    [nDM]                       ; 
TPad* pad       [nDM]                       ;   // MVA output
TH1F* MVAtt     [nDM]                       ;
TH1F* MVAWW     [nDM]                       ;
float y_sig     [nDM]                       ;   // yields
float y_bkg     [nDM]                       ;   // yields
float frontier  [nDM][nSB+1]                ;   // sidebands borders (frontier[nSB] = optimised cut)
TH2F* njnb      [nDM][nSB]  [nDD]           ;   // CR 2D-histograms
TCanvas* c2     [nDM]                       ;   // control regions 
TString estimated_bkg       [nDD]           ;   // estimated bkgs samples ids
TString  SBid        [nSB]                  ;   // SB ids    
//int wall     [4][nDM][nSB]       [nBox]   ;   // box boundaries (by hand): [njet-down, njet-up, nbjet-down, nbjet-up] 
int wall     [4]                 [nBox][nCf];   // box boundaries (by hand): [njet-down, njet-up, nbjet-down, nbjet-up] 
float population[nDM][nSB]  [nDD][nBox][nCf];   // populations inside the boxes at the CRs  
float SF        [nDM][nSB]  [nDD-2]    [nCf];   // scale factors 
float SF_e      [nDM][nSB]  [nDD-2]    [nCf];   // scale factors errors
float SF_cov    [nDM][nSB]             [nCf];   // scale factors covariance (case just two factors) 
TH1F* CRhisto[2][nDM][nSB]  [2]             ;   // CR 1D-histograms
TString    id               [2]             ;   // CR 1D-histograms
TString varID[2]                            ;   // CR 1D-histograms
TCanvas* c3     [nDM]                       ;   // CR 1D-histograms

// functions
// ---------
void assign(); 

TH2F* fillhistogram(TString sample); 

TH1F* getMVAoutput(TString DM, TString sample);

int GetExtremumScore(TH1F* MVAsig, TH1F* MVAbkg); 

TH2F* GetCRhistogram(int i, int SB, int bkg);

TH1F* GetCRsubhistogram(int i, int SB, int variable, int set);

float GetCRpopulation(int i, int SB, int DD, int box, int cfg);

void assign_top();

double sqr(double x);
double  rd(double x);
double  th(double x);


void assign() {

	int iDM = 0;

	DM[iDM++] = "ttDM1scalar20"    ;
	DM[iDM++] = "ttDM1scalar50"    ;
	DM[iDM++] = "ttDM1scalar500"   ;
	DM[iDM++] = "ttDM10scalar10"   ; 
	DM[iDM++] = "ttDM50scalar50"   ;
	DM[iDM++] = "ttDM50scalar200"  ;
	DM[iDM++] = "ttDM50scalar300"  ;
	DM[iDM++] = "ttDM150scalar200" ;  

	estimated_bkg[tt] = "TTJets"   ;
	estimated_bkg[WW] = "WWTo2L2Nu";
	estimated_bkg[bk] = "minorbkg" ;
	estimated_bkg[dt] = "data"     ;

	id[data] = "data"  ;
	id[bkg]  = "allbkg";

	varID[njet]  = "number of jets"  ;
	varID[nbjet] = "number of b-jets";

	int iSB = 0;
	//SBid[iSB++] = "sideband trans-3"  ;
	//SBid[iSB++] = "sideband trans-2"  ;
	//SBid[iSB++] = "sideband trans-1"  ;
	SBid[iSB++] = "sideband distal"   ;
	SBid[iSB++] = "sideband proximal" ;


int i = 0; 
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ 
//                                       box WW                                                                                           box tt
// ------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------- 
//     njet down              njet up                nbjet down             nbjet up                    njet down              njet up                nbjet down             nbjet up                 
// ------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------- 
wall[0][0][i]=0;        wall[1][0][i]= 1;        wall[2][0][i]=0;         wall[3][0][i]=1;        wall[0][1][i]=2;       wall[1][1][i]=4;        wall[2][1][i]=1;        wall[3][1][i++]=3;

wall[0][0][i]=1;        wall[1][0][i]=10;        wall[2][0][i]=0;         wall[3][0][i]=1;        wall[0][1][i]=2;       wall[1][1][i]=4;        wall[2][1][i]=1;        wall[3][1][i++]=3;       

wall[0][0][i]=0;        wall[1][0][i]=10;        wall[2][0][i]=0;         wall[3][0][i]=1;        wall[0][1][i]=2;       wall[1][1][i]=4;        wall[2][1][i]=1;        wall[3][1][i++]=3;        

wall[0][0][i]=0;        wall[1][0][i]= 1;        wall[2][0][i]=0;         wall[3][0][i]=1;        wall[0][1][i]=1;       wall[1][1][i]=5;        wall[2][1][i]=1;        wall[3][1][i++]=3;

wall[0][0][i]=1;        wall[1][0][i]=10;        wall[2][0][i]=0;         wall[3][0][i]=1;        wall[0][1][i]=1;       wall[1][1][i]=5;        wall[2][1][i]=1;        wall[3][1][i++]=3;    
   
wall[0][0][i]=0;        wall[1][0][i]=10;        wall[2][0][i]=0;         wall[3][0][i]=1;        wall[0][1][i]=1;       wall[1][1][i]=5;        wall[2][1][i]=1;        wall[3][1][i++]=3;       
// ------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------- 
}

TH1F* getMVAoutput(TString DM, TString sample) {

	TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/application_outputs/OPT/h_Cuts_" + DM + "_" + sample + ".root", "READ" );

	TH1F* MVAoutput = (TH1F*) file -> Get("hCuts");

	return MVAoutput;
}




