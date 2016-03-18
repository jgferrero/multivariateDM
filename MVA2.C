// root -l -x -q optim2.C+ 

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TBranch.h"
#include "TVector3.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"
//#include "/home/jgarciaf/root/tmva/test/TMVAGui.C"

void optim22(TString DM);
void optim23(TString DM, TString sample);

// main
void MVA2(){

	optim22("ttDM1scalar20"   );
	optim22("ttDM1scalar50"   );
	optim22("ttDM1scalar500"  );
	optim22("ttDM10scalar10"  ); 
	optim22("ttDM50scalar50"  );
	optim22("ttDM50scalar200" );
	optim22("ttDM50scalar300" );
	optim22("ttDM50scalar300" );
	optim22("ttDM150scalar200");

}


void optim22(TString DM){

	optim23( DM, DM                   );
	optim23( DM, "TTJets"             );
	optim23( DM, "DYJetsToLL_M-10to50"); 
	optim23( DM, "DYJetsToLL_M-50"    );
	optim23( DM, "ST_tW_top"          );
	optim23( DM, "ST_tW_antitop"      );
	optim23( DM, "WWTo2L2Nu"          );
	optim23( DM, "WZTo3LNu"           );  
	optim23( DM, "ZZTo4L"             ); 
	optim23( DM, "WJetsToLNu"         );

	optim23( DM, "DoubleEG"       );
	optim23( DM, "DoubleMuon"     );
	optim23( DM, "MuonEG"         );
	optim23( DM, "SingleElectron" );
	optim23( DM, "SingleMuon"     );

}


void optim23(TString DM, TString sample){


/////////
// READER
/////////           
// analogously to the Factory, the communication between the user application and the MVA methods
// is interfaced by the TMVA Reader , which is created by the user:
// TMVA::Reader* reader = new TMVA::Reader( "<options>" );

TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );   
/////////////////////////////////


//////////////////
// INPUT VARIABLES
//////////////////
// the  user  registers  the  names  of  the  input  variables  with  the  Reader.   
// they  are  required  to  be the same (and in the same order) as the names used for training 
// (this requirement is not actually mandatory, but enforced to ensure the consistency between training and application). 
// together with the name is given the address of a local variable, which carries the updated input values during the event loop
float channel, MET, mll, njet, nbjet15, lep1Pt, lep2Pt, jetpt1, jetpt2;
float  dphilep1lep2, dphilep1jet1, dphilep1jet2, dphilep1MET, dphilep2jet1, dphilep2jet2, dphilep2MET, dphijet1jet2, dphijet1MET, dphijet2MET, dphilep1lep2MET, dphijet1jet2MET;

reader -> AddVariable( "channel"  , &channel   ); 
reader -> AddVariable( "MET"      , &MET       ); 
reader -> AddVariable( "mll"      , &mll       ); 
reader -> AddVariable( "njet"     , &njet      ); 
reader -> AddVariable( "nbjet15"  , &nbjet15   ); 
reader -> AddVariable( "lep1Pt"   , &lep1Pt    ); 
reader -> AddVariable( "lep2Pt"   , &lep2Pt    ); 
reader -> AddVariable( "jetpt1"   , &jetpt1    ); 
reader -> AddVariable( "jetpt2"   , &jetpt2    ); 

reader -> AddVariable ("dphilep1lep2"   , &dphilep1lep2  );
reader -> AddVariable ("dphilep1jet1"   , &dphilep1jet1  );
reader -> AddVariable ("dphilep1jet2"   , &dphilep1jet2  );
reader -> AddVariable ("dphilep1MET"    , &dphilep1MET   );
reader -> AddVariable ("dphilep2jet1"   , &dphilep2jet1  );
reader -> AddVariable ("dphilep2jet2"   , &dphilep2jet2  );
reader -> AddVariable ("dphilep2MET"    , &dphilep2MET   );
reader -> AddVariable ("dphijet1jet2"   , &dphijet1jet2  );
reader -> AddVariable ("dphijet1MET"    , &dphijet1MET   );
reader -> AddVariable ("dphijet2MET"    , &dphijet2MET   );
reader -> AddVariable ("dphilep1lep2MET", &dphilep1lep2MET);
//reader -> AddVariable ("dphijet1jet2MET", &dphijet1jet2MET);

/////////////////////////////////


//////////////////////
// BOOKING MVA METHODS 
//////////////////////
// the selected MVA methods are booked with the Reader 
// using the weight files from the preceding training job:
// reader->BookMVA( "<YourMethodName>", "<path/JobName_MethodName.weights.xml>" );
// the first argument is a user defined name to ditinguish between methods 
// (it does not need to be the same as for training, although this could be a useful choice

//reader -> BookMVA( "Cuts", "training_test_outputs/weightsOPT/" + DM + "_Cuts.weights.xml" ); 
//reader -> BookMVA( "Cuts", "training_test_outputs/weightsOPT/" + DM + "_LD.weights.xml" ); 
reader -> BookMVA( "Cuts", "training_test_outputs/weightsOPT/" + DM + "_MLP.weights.xml" ); 
/////////////////////////////////


//////////////////////////
// REQUESTING MVA RESPONSE
//////////////////////////
// output histogram
//Int_t nbin = 100;
//TH1F *hMLP;
//hMLP  = new TH1F( "hMLP", "MLP",  nbin, -0.3, 2.2 );   //  EXTREMELY IMPORTANT: ANN OUTBUT BETWEEN 0 & 1 !!!

Int_t nbin = 40;
TH1F *hCuts;
hCuts  = new TH1F( "hCuts", "Cuts",  nbin, -0.5, 1.5 );


TFile *input(0);
input = TFile::Open( "minitrees/noNJetCut/bonsai_" + sample + ".root", "UPDATE" );
//input = TFile::Open( "minitrees/bonsai_" + sample + ".root", "UPDATE" );

TTree *theTree = (TTree*)input -> Get("latino");

cout << " ... opening file: " << sample << endl;
cout << " ... processing: " << theTree -> GetEntries() << " events" << endl;

// output branch 
float MVAresponse; 
TBranch *MVABranch = theTree -> Branch( "MVA"+DM, &MVAresponse, "MVAresponse/F" );
//theTree -> Branch( "MVA"+DM, &MVAresponse, "MVAresponse/F" );

//load "auxliar variables":                                                                          !?
//float <my_aux_var>;                                                                                !?
//float weight                                                                                       !?
//theTree->SetBranchAddress( "TChannel",  &<my_aux_var> ); //     CORTES                             !?
//theTree->SetBranchAddress( "Weight"  ,  &weight       ); //     PESADO funcion Fillhistogram       !?

float var01, var02, var03, var04, var05, var06, var07, var08, var09, var10, var11, var12, var13, var14, var15, var16, var17, var18, var19, var20;

theTree -> SetBranchAddress( "channel"  , &var01 );
theTree -> SetBranchAddress( "MET"      , &var02 );
theTree -> SetBranchAddress( "mll"      , &var03 );
theTree -> SetBranchAddress( "njet"     , &var04 );
theTree -> SetBranchAddress( "nbjet15"  , &var05 );
theTree -> SetBranchAddress( "lep1Pt"   , &var06 );
theTree -> SetBranchAddress( "lep2Pt"   , &var07 );
theTree -> SetBranchAddress( "jetpt1"   , &var08 );
theTree -> SetBranchAddress( "jetpt2"   , &var09 );

theTree -> SetBranchAddress( "dphilep1lep2", &var10 );
theTree -> SetBranchAddress( "dphilep1jet1", &var11 );
theTree -> SetBranchAddress( "dphilep1jet2", &var12 );
theTree -> SetBranchAddress( "dphilep1MET" , &var13 );
theTree -> SetBranchAddress( "dphilep2jet1", &var14 );
theTree -> SetBranchAddress( "dphilep2jet2", &var15 );
theTree -> SetBranchAddress( "dphilep2MET" , &var16 );
theTree -> SetBranchAddress( "dphijet1jet2", &var17 );
theTree -> SetBranchAddress( "dphijet1MET" , &var18 );
theTree -> SetBranchAddress( "dphijet2MET" , &var19 );
theTree -> SetBranchAddress( "dphilep1lep2MET" , &var20 );
//theTree -> SetBranchAddress( "dphijet1jet2MET" , &var21 );

float baseW;
theTree -> SetBranchAddress( "baseW"  , &baseW   );


//watch
//TStopwatch sw;   
//sw.Start();

//loop
int ntotal = 0; if ( theTree -> GetEntries() > 60000 ) ntotal = 60000; else ntotal = theTree -> GetEntries();

for ( Long64_t ievt = 0; ievt < theTree -> GetEntries(); ievt++ ) {

	if ( ievt % 500000 == 0 ) cout << "--- ... Processing event: " << ievt << endl;

	theTree -> GetEntry(ievt);

	channel   = var01;
	MET       = var02;
	mll       = var03;
	njet      = var04;
	nbjet15   = var05;
	lep1Pt    = var06;
	lep2Pt    = var07;
	jetpt1    = var08;
	jetpt2    = var09;
	dphilep1lep2 = var10;
	dphilep1jet1 = var11;
	dphilep1jet2 = var12;
	dphilep1MET  = var13;
	dphilep2jet1 = var14;
	dphilep2jet2 = var15;
	dphilep2MET  = var16;
	dphijet1jet2 = var17;
	dphijet1MET  = var18;
	dphijet2MET  = var19;
	dphilep1lep2MET = var20;
	//dphijet1jet2MET = var21;


	//if ( MET > 120 && var2_1 > 20  &&  var2_2 > 20 && eta1 < 2.4  &&  eta2 < 2.4  &&  flv1*flv2 < 0  &&  njet > 1 && ( channel > 1  ||  ( mll < 76 || mll > 106 ) )  ) {

// the rectangular cut classifier is special since it returns a binary answer for a given set of input variables and cuts. 
// the user must specify the desired signal efficiency to define the working pointaccording to which the Reader will choose the cuts:
// Bool_t passed = reader->EvaluateMVA( "Cuts", signalEfficiency );

			//float MVAresponse = reader -> EvaluateMVA("Cuts", 0.5);
			MVAresponse = reader -> EvaluateMVA("Cuts");
			hCuts -> Fill(MVAresponse, baseW);                   // histo
			MVABranch -> Fill();                                 // branch

}    
                                            

//watch
//sw.Stop();
//cout << "--- End of event loop: " << endl; ; sw.Print();

/////////////////////////////////

//to add the MVA-branch to the input file
theTree -> Write("", TObject::kOverwrite);

input -> Close();





//new file for the MVA-histogram
//TFile *target  =  TFile::Open ("application_outputs/CR/h_MLP_"+sample+".root","UPDATE" );
TFile *target  =  TFile::Open ("application_outputs/OPT/h_Cuts_" + DM + "_" + sample + ".root","UPDATE" );   

//hMLP -> Write();
hCuts -> Write();                                                                                            // comment when just it's intended to compute CR steps

cout << "--- Created root file: \"h_Cuts_" + DM + "_" + sample + ".root\" containing the MVA output histograms" << endl;
	
target -> Close();





delete reader;

std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;

} 




