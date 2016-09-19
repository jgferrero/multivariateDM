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
 
enum{ parall_G, transv_G, MET, VpT, nvert, nvariable_G };
enum{ parall_R, transv_R, scale,           nvariable_R }; 

enum{ pT, sumET, NVtx, nparameter};

enum{ nominal, JECu, JECd, UnEu, UnEd, nsystematic }; 


enum{ DoubleEG2016B, DoubleEG2016C, DoubleEG2016D, DoubleEG2016E, DoubleEG2016F, 
      DY_ee, TT_ee, TTsemiT_ee, TTsemiTbar_ee, 
      WW_ee, WZTo2L2Q_ee, WZTo3LNu_ee, ZZTo4L_ee, ZZTo2L2Q_ee, ZZTo2L2Nu_ee, WWW_ee, WWZ_ee, WZZ_ee, ZZZ_ee, 
      SingleTop1_ee, SingleTop2_ee, SingleTop3_ee, SingleTop4_ee, SingleTop5_ee,

      DoubleMuon2016B, DoubleMuon2016C, DoubleMuon2016D, DoubleMuon2016E, DoubleMuon2016F, 
      DY_mm, TT_mm, TTsemiT_mm, TTsemiTbar_mm,
      WW_mm, WZTo2L2Q_mm, WZTo3LNu_mm, ZZTo4L_mm, ZZTo2L2Q_mm, ZZTo2L2Nu_mm, WWW_mm, WWZ_mm, WZZ_mm, ZZZ_mm,
      SingleTop1_mm, SingleTop2_mm, SingleTop3_mm, SingleTop4_mm, SingleTop5_mm,

      SinglePhoton2016B, SinglePhoton2016C, SinglePhoton2016D, SinglePhoton2016E, SinglePhoton2016F,
      GJets40100, GJets100200, GJets200400, GJets400600, GJets600Inf, 
      QCD200300, QCD300500, QCD500700, QCD7001000, QCD10001500, QCD15002000, QCD2000Inf,
      WJets100200, WJets200400, WJets400600, WJets600800, WJets8001200, WJets12002500, WJets2500Inf,
      ZGJets, ZNuNuGJets40130, ZGTo2LG, 
      WGJets, WGToLNuG,
      TTGJets, TGJets,

      nprocess };   // APPROVAL ( Maria should be aware of those ones commented, 10.ix )


TString allhistosWriteTo  = "sync"; //""160910_further-processes"; //20Aug_qT-gt-50
TString allhistosReadFrom = "160914_Zee_yes-pu"; // "17Aug_new-xs_pu-rw" 

const float TheUnifiedLuminosity = 20.110; 
      float TheLuminosity[nkanal];

const bool DoFillHistograms = 1;  
const bool DoGlobalPlots    = 0; 
const bool DoFit            = 0; 

const bool SkipEventLoop    = 0; 

const bool RunOverAllData = 1; 
const bool RunOverAllMC   = 1; 

      int _TotalEntries[nprocess]; 

const int MaxEntriesData    = 100;
const int MaxEntriesMC      = 1000;

const int alpha = DoubleEG2016B;
const int omega = DoubleEG2016B;

const bool PrintSync = true; 


TString NTuplaDir   [nkanal]; 
TString kanalID     [nkanal]; 
TString kanalIDfancy[nkanal]; 


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
 

TString variableID_G     [nvariable_G];
TString variableIDfancy_G[nvariable_G];

TString variableID_R     [nvariable_R];
TString variableIDfancy_R[nvariable_R];
TString sigma_variable   [nvariable_R];

TString parameterID     [nparameter];
TString parameterIDfancy[nparameter];

TString systematicID[nsystematic];

float maxpT    = 300.; float minpT    = 0.; const int nbinpT    = 30; float _uppT   ;
float maxsumET =3000.; float minsumET = 0.; const int nbinsumET = 30; float _upsumET;  
float maxNVtx  =  40.; float minNVtx  = 0.; const int nbinNVtx  =  8; float _upNVtx ;

float maxuPara = 200.; float minuPara = -200.; const int nbinuPara = 80; 
float maxuPerp = 200.; float minuPerp = -200.; const int nbinuPerp = 80; 
float maxscale =   5.; float minscale =   -3.; const int nbinscale =400; 
float maxMET   = 200.; float minMET   =    0.; const int nbinMET   = 40;
float maxnvert =  50.; float minnvert =    0.; const int nbinnvert = 50; 


/*float maxpT    = 200.; float minpT    = 0.; const int nbinpT    =  5; float _uppT   ;
float maxsumET = 300.; float minsumET = 0.; const int nbinsumET =  5; float _upsumET;  
float maxNVtx  =  20.; float minNVtx  = 0.; const int nbinNVtx  =  5; float _upNVtx ;*/

const int aa = nvariable_G;   // just for next blocks
const int bb = nvariable_R;   // just for next blocks
const int zz = nprocess   ;   // just for next blocks

// ----- write ----------------------------------------------------
TH1F*    h_global_W       [aa][zz]                   [nsystematic];
TH1F*    h_resol_pT_W     [bb][zz][nbinpT   ]        [nsystematic];
TH1F*    h_resol_sumET_W  [bb][zz][nbinsumET]        [nsystematic];
TH1F*    h_resol_NVtx_W   [bb][zz][nbinNVtx ]        [nsystematic]; 
// ----- read -----------------------------------------------------
TH1F*    h_global         [aa][zz]                   [nsystematic];
TH1F*    h_data           [aa]                       [nsystematic];
TH1F*    h_mc             [aa]                       [nsystematic];
THStack* s_global         [aa]                                    ;
TH1F*    Ratio            [aa]                                    ;
TH1F*    Unc              [aa]                                    ;
TH1F*    h_resol_pT       [bb][zz][nbinpT   ]        [nsystematic];
TH1F*    h_resol_pT_fit   [bb][3 ][nbinpT   ][nkanal][nsystematic];
TH1F*    h_resol_sumET    [bb][zz][nbinsumET]        [nsystematic];
TH1F*    h_resol_sumET_fit[bb][3 ][nbinsumET][nkanal][nsystematic];
TH1F*    h_resol_NVtx     [bb][zz][nbinNVtx ]        [nsystematic]; 
TH1F*    h_resol_NVtx_fit [bb][3 ][nbinNVtx ][nkanal][nsystematic];
TH1F*    h_SR[3][nparameter][nkanal];   // 3 histos{ up pad, ratio, ratio-unc } x 3 param x 3 kanales
TString  histoID[3]; 
// ----------------------------------------------------------------
		 


void Assign          (                                                                                                          ); 
void  FillHistograms (                                                                                                          );
void GlobalPlots     ( int ch                                                                                                   );
void OrderFits       (                                                                                                          );
void  FillHistogram  ( int process                                                                                              );
void GetResolution   ( int ch, int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi, TString GJetsOrigin, int syst );
void PlotResolution  ( int ivar, int param                                                                                      );
void SetAxis         ( TH1* hist, TString xtitle, TString ytitle, Float_t xoffset, Float_t yoffset                              ); 
void DrawLatex       ( Font_t tfont, Float_t x, Float_t y, Float_t tsize, Short_t align, const char* text                       );
double GetFWHM       ( double sigma, double gamma                                                                               );
double GetFWHMerror  ( double sigma, double gamma, double esigma, double egamma, double Vss, double Vsg, double Vgs, double Vgg );
void BuildRatio      ( int ivar, int nbin                                                                                       );
void GetPUreweighting( int ch                                                                                                   );
void GetRatioValues  ( float x, float ex, float y, float ey, float y_JECu, float y_JECd, float y_UnEu, float y_UnEd, float& z, float& ez_data, float& ez_mc );
void MoveOverflows   ( TH1* hist                                                                                                );

void Assign(){

	TheLuminosity[Zee  ] = 5.885 + 2.646 + 4.329;// + 4.049 ;//+ 3.159;  // fb⁻¹  
 	TheLuminosity[Zmumu] = 5.920 + 2.645 + 4.330 + 4.050 + 3.160;  // fb⁻¹  
	TheLuminosity[Gamma] = 5.9   + 2.65  + 4.35  + 4.050 + 3.16 ;  // fb⁻¹     

	NTuplaDir[Zee  ] = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/diELE/APPROVAL/" ;
	NTuplaDir[Zmumu] = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/diMUON/APPROVAL/";
	NTuplaDir[Gamma] = "~/eos/cms/store/group/phys_jetmet/dalfonso/ICHEP/gamma/APPROVAL/" ;

	kanalID[Zee  ] = "Zee"  ;
	kanalID[Zmumu] = "Zmm"  ;
	kanalID[Gamma] = "Gamma";

	kanalIDfancy[Zee  ] = "Z #rightarrow ee"    ;
	kanalIDfancy[Zmumu] = "Z #rightarrow #mu#mu";
	kanalIDfancy[Gamma] = "#gamma + jets"       ;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// https://dalfonso.web.cern.ch/dalfonso/ICHEP/sampleLocation.txt
// https://github.com/GuillelmoGomezCeballos/MitAnalysisRunII/blob/master/data/xs.dat
// http://mvesterb.web.cern.ch/mvesterb/met/things/samples.dat

	int i = 0; 

isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "data";  sampleID[i] = "DoubleEG_Run2016B_PromptReco_v2"  ; processID[i++] = "DoubleEG2016B";    
isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "    ";  sampleID[i] = "DoubleEG_Run2016C_PromptReco_v2"  ; processID[i++] = "DoubleEG2016C";
isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "    ";  sampleID[i] = "DoubleEG_Run2016D_PromptReco_v2"  ; processID[i++] = "DoubleEG2016D";
isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "    ";  sampleID[i] = "DoubleEG_Run2016E_PromptReco_v2"  ; processID[i++] = "DoubleEG2016E";
isData[i] = 1; kanal[i] = 0; processIDfancy[i] = "    ";  sampleID[i] = "DoubleEG_Run2016F_PromptReco_v1"  ; processID[i++] = "DoubleEG2016F";

isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "DY"; ProcessColor[i] = kYellow; sampleID[i] = "DYJetsToLL_M50"             ; processID[i] = "DYJetsToLL_M50"     ; xs[i++] = 6025200. ;  
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "tt"; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_DiLepton"            ; processID[i] = "TTJets_DiLepton"    ; xs[i++] =   88287.7;  
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "  "; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_SingleLeptonFromT"   ; processID[i] = "TTJets_SemiFromT"   ; xs[i++] =  176575.4;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "  "; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_SingleLeptonFromTbar"; processID[i] = "TTJets_SemiFromTbar"; xs[i++] =  176575.4;

isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "dibosons"; ProcessColor[i] = kBlue  ; sampleID[i] = "WWTo2L2Nu"      ; processID[i] = "WWTo2L2Nu"      ; xs[i++] =   12177.6 ;        
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo2L2Q"       ; processID[i] = "WZTo2L2Q"       ; xs[i++] =    6204.9 ;  // 5595.0*1.109
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo3LNu"       ; processID[i] = "WZTo3LNu"       ; xs[i++] =    4911.8 ;  // 4429.0*1.109
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo4L"         ; processID[i] = "ZZTo4L"         ; xs[i++] =    1256.  ;          
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Q"       ; processID[i] = "ZZTo2L2Q"       ; xs[i++] =    3542.0 ;  // 3220.0*1.1  
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Nu"      ; processID[i] = "ZZTo2L2Nu"      ; xs[i++] =     564.  ; 
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WWW"            ; processID[i] = "WWW"            ; xs[i++] =     208.6 ;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WWZ"            ; processID[i] = "WWZ"            ; xs[i++] =     165.1 ;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZZ"            ; processID[i] = "WZZ"            ; xs[i++] =      55.7 ;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZZ"            ; processID[i] = "ZZZ"            ; xs[i++] =      14.0 ;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "T_tWch"         ; processID[i] = "T_tWch"         ; xs[i++] =   35600.  ;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TBar_tWch"      ; processID[i] = "TBar_tWch"      ; xs[i++] =   35600.  ;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TToLep_tch"     ; processID[i] = "TToLep_tch"     ; xs[i++] =   44070.48;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TBarToLep_tch"  ; processID[i] = "TBarToLep_tch"  ; xs[i++] =   26227.8 ;
isData[i] = 0; kanal[i] = 0; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TToLep_sch"     ; processID[i] = "TToLep_sch"     ; xs[i++] =    3680.64;


isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "data";  sampleID[i] = "DoubleMuon_Run2016B_PromptReco_v2"; processID[i++] = "DoubleMuon2016B";
isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "    ";  sampleID[i] = "DoubleMuon_Run2016C_PromptReco_v2"; processID[i++] = "DoubleMuon2016C";
isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "    ";  sampleID[i] = "DoubleMuon_Run2016D_PromptReco_v2"; processID[i++] = "DoubleMuon2016D";
isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "    ";  sampleID[i] = "DoubleMuon_Run2016E_PromptReco_v2"; processID[i++] = "DoubleMuon2016E";
isData[i] = 1; kanal[i] = 1; processIDfancy[i] = "    ";  sampleID[i] = "DoubleMuon_Run2016F_PromptReco_v1"; processID[i++] = "DoubleMuon2016F";

isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "DY"; ProcessColor[i] = kYellow; sampleID[i] = "DYJetsToLL_M50"             ; processID[i] = "DYJetsToLL_M50"     ; xs[i++] = 6025200. ;  
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "tt"; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_DiLepton"            ; processID[i] = "TTJets_DiLepton"    ; xs[i++] =   88287.7;  
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "  "; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_SingleLeptonFromT"   ; processID[i] = "TTJets_SemiFromT"   ; xs[i++] =  176575.4;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "  "; ProcessColor[i] = kRed   ; sampleID[i] = "TTJets_SingleLeptonFromTbar"; processID[i] = "TTJets_SemiFromTbar"; xs[i++] =  176575.4;

isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "dibosons"; ProcessColor[i] = kBlue  ; sampleID[i] = "WWTo2L2Nu"      ; processID[i] = "WWTo2L2Nu"      ; xs[i++] =   12177.6 ;        
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo2L2Q"       ; processID[i] = "WZTo2L2Q"       ; xs[i++] =    6204.9 ;  // 5595.0*1.109
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZTo3LNu"       ; processID[i] = "WZTo3LNu"       ; xs[i++] =    4911.8 ;  // 4429.0*1.109
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo4L"         ; processID[i] = "ZZTo4L"         ; xs[i++] =    1256.  ;          
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Q"       ; processID[i] = "ZZTo2L2Q"       ; xs[i++] =    3542.0 ;  // 3220.0*1.1  
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZTo2L2Nu"      ; processID[i] = "ZZTo2L2Nu"      ; xs[i++] =     564.  ; 
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WWW"            ; processID[i] = "WWW"            ; xs[i++] =     208.6 ;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WWZ"            ; processID[i] = "WWZ"            ; xs[i++] =     165.1 ;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "WZZ"            ; processID[i] = "WZZ"            ; xs[i++] =      55.7 ;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "ZZZ"            ; processID[i] = "ZZZ"            ; xs[i++] =      14.0 ;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "T_tWch"         ; processID[i] = "T_tWch"         ; xs[i++] =   35600.  ;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TBar_tWch"      ; processID[i] = "TBar_tWch"      ; xs[i++] =   35600.  ;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TToLep_tch"     ; processID[i] = "TToLep_tch"     ; xs[i++] =   44070.48;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TBarToLep_tch"  ; processID[i] = "TBarToLep_tch"  ; xs[i++] =   26227.8 ;
isData[i] = 0; kanal[i] = 1; processIDfancy[i] = "        "; ProcessColor[i] = kBlue  ; sampleID[i] = "TToLep_sch"     ; processID[i] = "TToLep_sch"     ; xs[i++] =    3680.64;


isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016B_PromptReco_v2";           processID[i] = "SinglePhoton2016B"      ; processIDfancy[i++] = "data";
isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016C_PromptReco_v2";           processID[i] = "SinglePhoton2016C"      ; processIDfancy[i++] = "    ";
isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016D_PromptReco_v2";           processID[i] = "SinglePhoton2016D"      ; processIDfancy[i++] = "    "; 
isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016E_PromptReco_v2";           processID[i] = "SinglePhoton2016E"      ; processIDfancy[i++] = "    "; 
isData[i] = 1; kanal[i] = 2; sampleID[i] = "SinglePhoton_Run2016F_PromptReco_v1";           processID[i] = "SinglePhoton2016F"      ; processIDfancy[i++] = "    ";  // not yet

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

isData[i] = 0; kanal[i] = 2; sampleID[i] = "ZGJets"         ; xs[i] =    190.3; processID[i] = "ZGJets"         ; processIDfancy[i] = "Z#gamma + jets" ; ProcessColor[i++] = kViolet; 
isData[i] = 0; kanal[i] = 2; sampleID[i] = "ZNuNuGJets40130"; xs[i] =   3928.4; processID[i] = "ZNuNuGJets40130"; processIDfancy[i] = " "              ; ProcessColor[i++] = kViolet; // 2806.0*1.4
isData[i] = 0; kanal[i] = 2; sampleID[i] = "ZGTo2LG"        ; xs[i] = 123900. ; processID[i] = "ZGTo2LG"        ; processIDfancy[i] = " "              ; ProcessColor[i++] = kViolet; 

isData[i] = 0; kanal[i] = 2; sampleID[i] = "WGJets"         ; xs[i] =   929.18; processID[i] = "WGJets"         ; processIDfancy[i] = "W#gamma + jets" ; ProcessColor[i++] = kCyan  ; // 663.7*1.4
isData[i] = 0; kanal[i] = 2; sampleID[i] = "WGToLNuG"       ; xs[i] = 585960. ; processID[i] = "WGToLNuG"       ; processIDfancy[i] = " "              ; ProcessColor[i++] = kCyan  ; // 514000.*1.14

isData[i] = 0; kanal[i] = 2; sampleID[i] = "TTGJets"        ; xs[i] =   4620.8; processID[i] = "TTGJets"        ; processIDfancy[i] = "tt#gamma + jets"; ProcessColor[i++] = kPink  ; // 1444.*3.2 
isData[i] = 0; kanal[i] = 2; sampleID[i] = "TGJets"         ; xs[i] =   2967. ; processID[i] = "TGJets"         ; processIDfancy[i] = "               "; ProcessColor[i  ] = kPink  ;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	variableID_G[parall_G] = "parallel"  ; 
	variableID_G[transv_G] = "transverse"; 
	variableID_G[MET     ] = "MET"       ;
	variableID_G[VpT     ] = "VpT"       ;
	variableID_G[nvert   ] = "nvert"     ;

	variableIDfancy_G[parall_G] = "u_{||} + q_{T} [GeV]"   ; 
	variableIDfancy_G[transv_G] = "u_{#perp}  [GeV]"       ;
	variableIDfancy_G[MET     ] = "E_{T}^{miss} [GeV]"     ;
	variableIDfancy_G[VpT     ] = "q_{T} [GeV]"            ;
	variableIDfancy_G[nvert   ] = "Number of Vertices"     ;


	variableID_R[parall_R] = "parallel"  ; 
	variableID_R[transv_R] = "transverse"; 
	variableID_R[scale   ] = "scale"     ;

	variableIDfancy_R[parall_R] = "u_{||} + q_{T} [GeV]"   ; 
	variableIDfancy_R[transv_R] = "u_{#perp}  [GeV]"       ;
	variableIDfancy_R[scale   ] = "- u_{||} / q_{T}"       ;

	sigma_variable[parall_R] = "#sigma( u_{||} ) [GeV]"    ;
	sigma_variable[transv_R] = "#sigma( u_{#perp} )  [GeV]";
	sigma_variable[scale   ] = "- #LT u_{||} / q_{T} #GT"  ;



	parameterID[pT   ] = "pT"   ;
	parameterID[sumET] = "sumET";
	parameterID[NVtx ] = "NVtx" ;

	parameterIDfancy[pT   ] = "q_{T} [GeV]";
	parameterIDfancy[sumET] = "#Sigma E_{T} [TeV]";
	parameterIDfancy[NVtx ] = "Number of Vertices";


	histoID[0] = "up"            ;
	histoID[1] = "ratio"         ;
	histoID[2] = "ratio-unc"     ;


	
	systematicID[nominal] = "nominal"; 
	systematicID[JECu   ] = "JEC-up" ; 
	systematicID[JECd   ] = "JEC-do" ; 
	systematicID[UnEu   ] = "UnE-up" ; 
	systematicID[UnEd   ] = "UnE-do" ; 

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



void BuildRatio( int ivar, int nbin ) {

	float dtSum_last = 0; 

	for( int r = -1; r < nbin+1; r++ ){

		float dtSum = h_data[ivar][nominal] -> GetBinContent(r+1);
		float dtErr = h_data[ivar][nominal] -> GetBinError(r+1)  ;
		float mcSum = h_mc  [ivar][nominal] -> GetBinContent(r+1);
		float mcErr = h_mc  [ivar][nominal] -> GetBinError(r+1)  ;

		float dtSum_JECu = h_data[ivar][JECu] -> GetBinContent(r+1);
		float dtErr_JECu = h_data[ivar][JECu] -> GetBinError(r+1)  ;
		float mcSum_JECu = h_mc  [ivar][JECu] -> GetBinContent(r+1);
		float mcErr_JECu = h_mc  [ivar][JECu] -> GetBinError(r+1)  ;

		float dtSum_JECd = h_data[ivar][JECd] -> GetBinContent(r+1);
		float dtErr_JECd = h_data[ivar][JECd] -> GetBinError(r+1)  ;
		float mcSum_JECd = h_mc  [ivar][JECd] -> GetBinContent(r+1);
		float mcErr_JECd = h_mc  [ivar][JECd] -> GetBinError(r+1)  ;

		float mcSum_UnEu = h_mc  [ivar][UnEu] -> GetBinContent(r+1);
		float mcSum_UnEd = h_mc  [ivar][UnEd] -> GetBinContent(r+1);


		float TheRatio 	= 999; 
		float TheError  = 999;
		float TheUncErr = 999; 

		if( mcSum > 0. ) {

			TheRatio =          dtSum / mcSum;
			TheError =          dtErr / mcSum;
			TheUncErr= TheRatio * sqrt(  pow( mcErr/mcSum, 2 )  +  pow( ((mcSum_JECu-mcSum_JECd)/2)/mcSum, 2 )  +  pow( ((mcSum_UnEu-mcSum_UnEd)/2)/mcSum, 2 )  );  

		}

		Ratio[ivar] -> SetBinContent(r+1, TheRatio);  
		Ratio[ivar] -> SetBinError  (r+1, TheError);

		Unc[ivar] -> SetBinContent(r+1, 1.0);
		Unc[ivar] -> SetBinError  (r+1, TheUncErr);

		dtSum_last = dtSum; 

		//if( ivar == nvert )  cout << r+1 << " -- " << dtSum << " -- " << mcSum << " -- "  << endl;


	}

	//MoveOverflows( Ratio[ivar] );
	//MoveOverflows( Unc  [ivar] );	 

}



void GetPUreweighting( int ch ) {

	TFile* myfile = new TFile("pileup/PU-reweighting_160914_Zee-no-pu.root", "update");

	TH1F* h_pu = (TH1F*) h_data[nvert][nominal] -> Clone( "h_pu-reweighting_" + kanalID[ch] );

	float dt_total = h_data[nvert][nominal]->Integral(); 
	float mc_total = h_mc  [nvert][nominal]->Integral();

	for( int i = 1; i < h_pu->GetNbinsX()+1; i++ ) { 

		float dt = h_data[nvert][nominal]->GetBinContent(i); 
		float mc = h_mc  [nvert][nominal]->GetBinContent(i);  

		float sf;  ( mc > 0. )  ?  sf = (dt/dt_total)/(mc/mc_total)  :  sf = 1.0; 	

		h_pu -> SetBinContent( i, sf ); 

	}

	h_pu -> Write();

	myfile->Close();
	
} 



void MoveOverflows(TH1* hist)
{

  Float_t xmin = hist->GetBinLowEdge(1); 
  Float_t xmax = hist->GetBinLowEdge(hist->GetNbinsX()+1);  	//if( xmax == 50.)  xmax = 45.;

  int nentries = hist->GetEntries();
  int nbins    = hist->GetNbinsX();
  
  TAxis* xaxis = (TAxis*)hist->GetXaxis();


  // Underflow
  //----------------------------------------------------------------------------
  if (xmin != -999)
    {
      Int_t   firstBin = -1;
      Float_t firstVal = 0;
      Float_t firstErr = 0;
      
      for (Int_t i=0; i<=nbins+1; i++)
	{
	  if (xaxis->GetBinLowEdge(i) < xmin)
	    {
	      firstVal += hist->GetBinContent(i);
	      firstErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      hist->SetBinContent(i, 0);
	      hist->SetBinError  (i, 0);
	    }
	  else if (firstBin == -1)
	    {
	      firstVal += hist->GetBinContent(i);
	      firstErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      firstBin = i;
	    }
	}

      firstErr = sqrt(firstErr);
  
      hist->SetBinContent(firstBin, firstVal);
      hist->SetBinError  (firstBin, firstErr);
    }


  // Overflow
  //----------------------------------------------------------------------------
  if (xmax != -999)
    {
      Int_t   lastBin = -1;
      Float_t lastVal = 0;
      Float_t lastErr = 0;
      
      for (Int_t i=nbins+1; i>=0; i--)
	{
	  Float_t lowEdge = xaxis->GetBinLowEdge(i);
      
	  if (lowEdge >= xmax)
	    {
	      lastVal += hist->GetBinContent(i);
	      lastErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      hist->SetBinContent(i, 0);
	      hist->SetBinError  (i, 0);
	    }
	  else if (lastBin == -1)
	    {
	      lastVal += hist->GetBinContent(i);
	      lastErr += (hist->GetBinError(i)*hist->GetBinError(i));
	      lastBin = i;
	    }
	}

      lastErr = sqrt(lastErr);
  
      hist->SetBinContent(lastBin, lastVal);
      hist->SetBinError  (lastBin, lastErr);
    }

  hist->SetEntries(nentries);
}


