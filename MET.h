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

enum{ Zee, Zmumu, Gamma, nkanal }; 
enum{ parall, transv, MET, nvariable };
enum{ pT, sumET, NVtx, nparameter};


enum{ DoubleEG2016B  , DoubleEG2016C  , DoubleEG2016D  , 
      DY_ee, TT_ee, WW_ee, WZTo2L2Q_ee, WZTo3LNu_ee, ZZTo4L_ee, ZZTo2L2Q_ee, ZZTo2L2Nu_ee, 

      DoubleMuon2016B, DoubleMuon2016C, DoubleMuon2016D, 
      DY_mm, TT_mm, WW_mm, WZTo2L2Q_mm, WZTo3LNu_mm, ZZTo4L_mm, ZZTo2L2Q_mm, ZZTo2L2Nu_mm, 

      SinglePhoton2016B, SinglePhoton2016C, SinglePhoton2016D, 
      GJets40100, GJets100200, GJets200400, GJets400600, GJets600Inf, 
      QCD200300, QCD300500, QCD500700, QCD7001000, QCD10001500, QCD15002000, QCD2000Inf,
      WJets100200, WJets200400, WJets400600, WJets600800, WJets8001200, WJets12002500, WJets2500Inf,
      ZGJets, ZNuNuGJets40130, ZGTo2LG, 
      WGJets, WGToLNuG,
      TTGJets, TGJets,

      nprocess };   // APPROVAL



TString allhistosWriteTo  = "11Aug_ALL"; 
TString allhistosReadFrom = "11Aug_ALL"; 

const float TheLuminosity = 12.895;   //  5.9202 + 2.6453 + 4.3296 fb⁻¹  

const bool DoFillHistograms = false;  
const bool DoGlobalPlots    = false; 
const bool DoFit            = true ; 

const bool SkipEventLoop    = false; 

const bool RunOverAllData = true; 
const bool RunOverAllMC   = true; 

      int _TotalEntries[nprocess]; 

const int MaxEntriesData    = 100;
const int MaxEntriesMC      = 10000;

//const int MaxEntriesMCStored =  20000; // 20000 (11vii)

TString NTuplaDir[nkanal]; 
TString kanalID  [nkanal]; 


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
bool    isData        [nprocess];
int     kanal         [nprocess];
 

TString xx             [nvariable];
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
TH1F*    h_global     [aa][bb]           ;
TH1F*    h_data       [aa]               ;
TH1F*    h_mc         [aa]               ;
THStack* s_global     [aa]               ;
TH1F*    Ratio        [aa]               ;
TH1F*    Unc          [aa]               ;
TH1F*    h_resol_pT   [aa][bb][nbinpT   ];
TH1F*    h_resol_pT_fit   [aa][3 ][nbinpT   ][nkanal];
TH1F*    h_resol_sumET[aa][bb][nbinsumET];
TH1F*    h_resol_sumET_fit[aa][3 ][nbinsumET][nkanal];
TH1F*    h_resol_NVtx [aa][bb][nbinNVtx ]; 
TH1F*    h_resol_NVtx_fit [aa][3 ][nbinNVtx ][nkanal]; 
// -----------------------------------------


float xpT   [nbinpT   ]; float ypT_A_Zee   [nbinpT   ]; float epT_A_Zee   [nbinpT   ]; float chi2pT_A_Zee   [nbinpT   ]; 
			 float ypT_A_Zmm   [nbinpT   ]; float epT_A_Zmm   [nbinpT   ]; float chi2pT_A_Zmm   [nbinpT   ]; 
			 float ypT_A_lum   [nbinpT   ]; float epT_A_lum   [nbinpT   ]; float chi2pT_A_lum   [nbinpT   ]; 
			 float ypT_B       [nbinpT   ]; float epT_B       [nbinpT   ]; float chi2pT_B       [nbinpT   ];

float xsumET[nbinsumET]; float ysumET_A_Zee[nbinsumET]; float esumET_A_Zee[nbinsumET]; float chi2sumET_A_Zee[nbinsumET];
			 float ysumET_A_Zmm[nbinsumET]; float esumET_A_Zmm[nbinsumET]; float chi2sumET_A_Zmm[nbinsumET];
			 float ysumET_A_lum[nbinsumET]; float esumET_A_lum[nbinsumET]; float chi2sumET_A_lum[nbinsumET];
			 float ysumET_B    [nbinsumET]; float esumET_B    [nbinsumET]; float chi2sumET_B    [nbinsumET];

float xNVtx [nbinNVtx ]; float yNVtx_A_Zee [nbinNVtx ]; float eNVtx_A_Zee [nbinNVtx ]; float chi2NVtx_A_Zee [nbinNVtx ];
			 float yNVtx_A_Zmm [nbinNVtx ]; float eNVtx_A_Zmm [nbinNVtx ]; float chi2NVtx_A_Zmm [nbinNVtx ];
			 float yNVtx_A_lum [nbinNVtx ]; float eNVtx_A_lum [nbinNVtx ]; float chi2NVtx_A_lum [nbinNVtx ];
			 float yNVtx_B     [nbinNVtx ]; float eNVtx_B     [nbinNVtx ]; float chi2NVtx_B     [nbinNVtx ];



void Assign          (                                                                                                          ); 
void  FillHistograms (                                                                                                          );
void GlobalPlots     ( int ch                                                                                                   );
void OrderFits       (                                                                                                          );
void  FillHistogram  ( int process                                                                                              );
void GetResolution   ( int ch, int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi, TString GJetsOrigin );
void PlotResolution  ( int ivar, int parameter                                                                                  );
void SetAxis         ( TH1* hist, TString xtitle, TString ytitle, Float_t xoffset, Float_t yoffset                              ); 
void DrawLatex       ( Font_t tfont, Float_t x, Float_t y, Float_t tsize, Short_t align, const char* text                       );
double GetFWHM       ( double sigma, double gamma                                                                               );
double GetFWHMerror  ( double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg );




void Assign(){

	NTuplaDir[Zee  ] = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/diELE/APPROVAL/" ;
	NTuplaDir[Zmumu] = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/diMUON/APPROVAL/";
	NTuplaDir[Gamma] = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/gamma/APPROVAL/" ;

	kanalID[Zee  ] = "Zee"  ;
	kanalID[Zmumu] = "Zmm"  ;
	kanalID[Gamma] = "Gamma";


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	int i = 0; 


isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "data";  sampleID[i] = "DoubleEG_Run2016B_PromptReco_v2"  ; processID[i++] = "DoubleEG2016B";                                       
isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "    ";  sampleID[i] = "DoubleEG_Run2016C_PromptReco_v2"  ; processID[i++] = "DoubleEG2016C";
isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "    ";  sampleID[i] = "DoubleEG_Run2016D_PromptReco_v2"  ; processID[i++] = "DoubleEG2016D";

isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "DY"      ; ProcessColor[i] = kYellow; sampleID[i] = "DYJetsToLL_M50" ; processID[i] = "DYJetsToLL_M50" ; xs[i++] = 6024000.0;     
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "tt"      ; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_DiLepton"; processID[i] = "TTJets_DiLepton"; xs[i++] =   87310.0;         
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "dibosons"; ProcessColor[i] = kBlue  ; sampleID[i] = "WWTo2L2Nu"      ; processID[i] = "WWTo2L2Nu"      ; xs[i++] =   12178.0;        
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo2L2Q"       ; processID[i] = "WZTo2L2Q"       ; xs[i++] =    5595.0;           
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo3LNu"       ; processID[i] = "WZTo3LNu"       ; xs[i++] =    4429.0;         
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo4L"         ; processID[i] = "ZZTo4L"         ; xs[i++] =    1256.0;          
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Q"       ; processID[i] = "ZZTo2L2Q"       ; xs[i++] =    3280.0;          
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Nu"      ; processID[i] = "ZZTo2L2Nu"      ; xs[i++] =     564.0; 



isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "data";  sampleID[i] = "DoubleMuon_Run2016B_PromptReco_v2"; processID[i++] = "DoubleMuon2016B";
isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "    ";  sampleID[i] = "DoubleMuon_Run2016C_PromptReco_v2"; processID[i++] = "DoubleMuon2016C";
isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "    ";  sampleID[i] = "DoubleMuon_Run2016D_PromptReco_v2"; processID[i++] = "DoubleMuon2016D";

isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "DY"      ; ProcessColor[i] = kYellow; sampleID[i] = "DYJetsToLL_M50" ; processID[i] = "DYJetsToLL_M50" ; xs[i++] = 6024000.0;     
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "tt"      ; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_DiLepton"; processID[i] = "TTJets_DiLepton"; xs[i++] =   87310.0;         
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "dibosons"; ProcessColor[i] = kBlue  ; sampleID[i] = "WWTo2L2Nu"      ; processID[i] = "WWTo2L2Nu"      ; xs[i++] =   12178.0;        
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo2L2Q"       ; processID[i] = "WZTo2L2Q"       ; xs[i++] =    5595.0;           
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo3LNu"       ; processID[i] = "WZTo3LNu"       ; xs[i++] =    4429.0;         
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo4L"         ; processID[i] = "ZZTo4L"         ; xs[i++] =    1256.0;          
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Q"       ; processID[i] = "ZZTo2L2Q"       ; xs[i++] =    3280.0;          
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Nu"      ; processID[i] = "ZZTo2L2Nu"      ; xs[i++] =     564.0; 



isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016B_PromptReco_v2";           processID[i] = "SinglePhoton2016B"      ; processIDfancy[i++] = "data";
isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016C_PromptReco_v2";           processID[i] = "SinglePhoton2016C"      ; processIDfancy[i++] = "    ";
isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016D_PromptReco_v2";           processID[i] = "SinglePhoton2016D"      ; processIDfancy[i++] = "    "; 

isData[i] = 0; kanal[i] = 2; sampleID[i] = "GJets_HT40to100"        ; xs[i] =   20730000.0; processID[i] = "GJets40100"     ; processIDfancy[i] = "#gamma + jets"  ; ProcessColor[i++] = kYellow;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "GJets_HT100to200"       ; xs[i] =    9226000.0; processID[i] = "GJets100200"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "GJets_HT200to400"       ; xs[i] =    2300000.0; processID[i] = "GJets200400"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "GJets_HT400to600"       ; xs[i] =     277400.0; processID[i] = "GJets400600"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "GJets_HT600toInf"       ; xs[i] =      93380.0; processID[i] = "GJets600Inf"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kYellow;

isData[i] = 0; kanal[i] = 2; sampleID[i] = "QCD_HT200to300_ext"     ; xs[i] = 1735000000.0; processID[i] = "QCD200300"      ; processIDfancy[i] = "QCD"            ; ProcessColor[i++] = kGreen ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "QCD_HT300to500"         ; xs[i] =  366800000.0; processID[i] = "QCD300500"      ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "QCD_HT500to700"         ; xs[i] =   29370000.0; processID[i] = "QCD500700"      ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "QCD_HT700to1000_ext"    ; xs[i] =    6524000.0; processID[i] = "QCD7001000"     ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "QCD_HT1000to1500"       ; xs[i] =    1064000.0; processID[i] = "QCD10002000"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "QCD_HT1500to2000"       ; xs[i] =     119900.0; processID[i] = "QCD15002000"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "QCD_HT2000toInf"        ; xs[i] =      25240.0; processID[i] = "QCD2000Inf"     ; processIDfancy[i] = " "              ; ProcessColor[i++] = kGreen ;

isData[i] = 0; kanal[i] = 2; sampleID[i] = "WJetsToLNu_HT100to200"  ; xs[i] =    1627450.0; processID[i] = "WJets100200"    ; processIDfancy[i] = "W + jets"       ; ProcessColor[i++] = kBlue  ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WJetsToLNu_HT200to400"  ; xs[i] =     435237.0; processID[i] = "WJets200400"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WJetsToLNu_HT400to600"  ; xs[i] =      59191.0; processID[i] = "WJets400600"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WJetsToLNu_HT600to800"  ; xs[i] =      14581.0; processID[i] = "WJets600800"    ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WJetsToLNu_HT800to1200" ; xs[i] =       6655.0; processID[i] = "WJets8001200"   ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WJetsToLNu_HT1200to2500"; xs[i] =       1608.1; processID[i] = "WJets12002500"  ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WJetsToLNu_HT2500toInf" ; xs[i] =         38.9; processID[i] = "WJets2500Inf"   ; processIDfancy[i] = " "              ; ProcessColor[i++] = kBlue  ;

isData[i] = 0; kanal[i] = 2; sampleID[i] = "ZGJets"                 ; xs[i] =        190.3; processID[i] = "ZGJets"         ; processIDfancy[i] = "Z#gamma + jets" ; ProcessColor[i++] = kViolet;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "ZNuNuGJets40130"        ; xs[i] =       2816.0; processID[i] = "ZNuNuGJets40130"; processIDfancy[i] = " "              ; ProcessColor[i++] = kViolet;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "ZGTo2LG"                ; xs[i] =     117864.0; processID[i] = "ZGTo2LG"        ; processIDfancy[i] = " "              ; ProcessColor[i++] = kViolet;

isData[i] = 0; kanal[i] = 2; sampleID[i] = "WGJets"                 ; xs[i] =        663.7; processID[i] = "WGJets"         ; processIDfancy[i] = "W#gamma + jets" ; ProcessColor[i++] = kCyan  ;
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WGToLNuG"               ; xs[i] =     585800.0; processID[i] = "WGToLNuG"       ; processIDfancy[i] = " "              ; ProcessColor[i++] = kCyan  ;

isData[i] = 0; kanal[i] = 2; sampleID[i] = "TTGJets"                ; xs[i] =       3697.0; processID[i] = "TTGJets"        ; processIDfancy[i] = "tt#gamma + jets"; ProcessColor[i++] = kPink  ; 
isData[i] = 0; kanal[i] = 2; sampleID[i] = "TGJets"                 ; xs[i] =       3697.0; processID[i] = "TGJets"         ; processIDfancy[i] = "               "; ProcessColor[i  ] = kPink  ;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	variableID[parall] = "parallel"  ; 
	variableID[transv] = "transverse"; 
	variableID[MET   ] = "MET"       ;

	variableIDfancy[parall] = "u_{||}+V_{T} [GeV]"; 
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


