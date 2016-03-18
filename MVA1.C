//root -l -x optim1.C+ 

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"
//#include "TMVA/TMVAGui.h"


void optim12(TString DM); 


//////////////
//// main ///////////////////////////////////////////
//////////////

void MVA1(){

	optim12("ttDM1scalar20"   );
	optim12("ttDM1scalar50"   );
	optim12("ttDM1scalar500"  );
	optim12("ttDM10scalar10"  ); 
	optim12("ttDM50scalar50"  );
	optim12("ttDM50scalar200" );
	optim12("ttDM50scalar300" );
	optim12("ttDM50scalar300" );
	optim12("ttDM150scalar200");

        // proposed id
	/*optim12 ("ttDM001scalar0010");
	optim12 ("ttDM001scalar0020");
	optim12 ("ttDM001scalar0050");
	optim12 ("ttDM001scalar0100");
	optim12 ("ttDM001scalar0200");
	optim12 ("ttDM001scalar0300");
	optim12 ("ttDM001scalar0500");
	optim12 ("ttDM001scalar1000");
	optim12 ("ttDM010scalar0010");
	optim12 ("ttDM010scalar0050");
	optim12 ("ttDM010scalar0100");
	optim12 ("ttDM050scalar0050");
	optim12 ("ttDM050scalar0200");
	optim12 ("ttDM050scalar0300");
	optim12 ("ttDM150scalar0200");
	//ttDM150scalar0500 not available
	//ttDM150scalar1000 not available
	optim12 ("ttDM500scalar0500");
	optim12 ("ttDM001pseudo0010");
	optim12 ("ttDM001pseudo0020");
	optim12 ("ttDM001pseudo0050");
	optim12 ("ttDM001pseudo0100");
	optim12 ("ttDM001pseudo0200");
	optim12 ("ttDM001pseudo0300");
	optim12 ("ttDM001pseudo0500");
	//ttDM001pseudo1000 not available
	optim12 ("ttDM010pseudo0010");
	optim12 ("ttDM010pseudo0050");
	optim12 ("ttDM010pseudo0100");
	optim12 ("ttDM050pseudo0050");
	optim12 ("ttDM050pseudo0200");
	optim12 ("ttDM050pseudo0300");
	optim12 ("ttDM150pseudo0200");
	optim12 ("ttDM150pseudo0500");
	//ttDM150pseudo1000 not available
	optim12 ("ttDM500pseudo0500");*/

}

/////////////////////////////////////////////////////

/////////////////
//// optim12 ////////////////////////////////////////
/////////////////

void optim12(TString DM){

// (blocks from fig 7 sec 3.1 TMVA Users Guide)  


                           /////////////////////////////////
/////////////////////////////  block 0: ROOT TARGET FILE  //////////////////////////////////////////////////////////////////////// 
                           /////////////////////////////////   

TString outfileName( "training_test_outputs/trainOPT/" + DM + ".root" );   

TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

(TMVA::gConfig().GetIONames()).fWeightFileDir = "training_test_outputs/weightsOPT";  // weights output file
   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  



                           ////////////////////////
/////////////////////////////  block 1: FACTORY  ///////////////////////////////////////////////////////////////////////////////// 
                           ////////////////////////  

// TMVA::Factory* factory = new TMVA::Factory( "<JobName>", outputFile, "<options>" );
// list of transformations to test; formatting example: Transformations=I;D;P;U;G,D, for identity, decorrelation, PCA, Uniform and  Gaussianisation followed by decorrelation transformations.
// the list of transformations contains a default set of data preprocessing steps for test and visualisation purposes only

TMVA::Factory *factory = new TMVA::Factory( DM, outputFile,    
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  



                           //////////////////////////////////////////////////////
/////////////////////////////  block 2: INITIALIALISE TRAINING AND TEST TREES  /////////////////////////////////////////////////// 
                           //////////////////////////////////////////////////////

// the input data sets used for training and testing of the multivariate methods need to be handed to the Factory
// data trees can be provided specifially for the purpose of either training or testing or for both purposes.
// in the latter case the factory then splits the tree into one part for training, the other for testing (see also sec 3.1.4)
   
TString path = "minitrees/bonsai_";

// sig
//----
TString sigFile    = path + DM + ".root";

TFile *inputSFile  = TFile::Open(sigFile, "read");
TTree *signal      = (TTree*)inputSFile -> Get("latino");
cout << "--- TMVAnalysis/Accessing Signal: " << sigFile << endl;

// bkgs
//-----
TString bkg1File   = path + "TTJets.root"             ;
TString bkg2File   = path + "DYJetsToLL_M-10to50.root"; 
TString bkg3File   = path + "DYJetsToLL_M-50.root"    ;
TString bkg4File   = path + "ST_tW_top.root"          ;
TString bkg5File   = path + "ST_tW_antitop.root"      ;
TString bkg6File   = path + "WWTo2L2Nu.root"          ;
TString bkg7File   = path + "WZTo3LNu.root"           ;  
TString bkg8File   = path + "ZZTo4L.root"             ; 
TString bkg9File   = path + "WJetsToLNu.root"         ;

TFile *inputB1File = TFile::Open(bkg1File, "read");
TTree *background1 = (TTree*)inputB1File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 1: " << bkg1File << endl;

TFile *inputB2File = TFile::Open(bkg2File, "read");
TTree *background2 = (TTree*)inputB2File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 2 : " << bkg2File << endl;

TFile *inputB3File = TFile::Open(bkg3File, "read");
TTree *background3 = (TTree*)inputB3File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 3: " << bkg3File << endl;

TFile *inputB4File = TFile::Open(bkg4File, "read");
TTree *background4 = (TTree*)inputB4File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 4: " << bkg4File << endl;
 
TFile *inputB5File = TFile::Open(bkg5File, "read");
TTree *background5 = (TTree*)inputB5File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 5: " << bkg5File << endl;

TFile *inputB6File = TFile::Open(bkg6File, "read");
TTree *background6 = (TTree*)inputB6File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 6: " << bkg6File << endl;

TFile *inputB7File = TFile::Open(bkg7File, "read");
TTree *background7 = (TTree*)inputB7File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 7: " << bkg7File << endl;

TFile *inputB8File = TFile::Open(bkg8File, "read");
TTree *background8 = (TTree*)inputB8File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 8: " << bkg8File << endl;

TFile *inputB9File = TFile::Open(bkg9File, "read");
TTree *background9 = (TTree*)inputB9File -> Get("latino");
cout << "--- TMVAnalysis/Accessing Bkg 9: " << bkg9File << endl;
 
// set the event weights (these weights are applied in addition to individual event weights that can be specified)

//float lumi = 1.269;
float lumi = 0.553;

// sig
Double_t sigWeight  = lumi;

// bkg
//Double_t denominator = SUM_{i=1}^n  xs_i/ Nevents_i 
//Double_t bkg1Weight  = xs_1 / Nevents_bkg1 / denominator;
//Double_t bkg2Weight  = xs_2 / Nevents_bkg2 / denominator;
// ... 
//Double_t bkg1Weight  = xs_n / Nevents_bkgn / denominator;

Double_t bkg1Weight = lumi;
Double_t bkg2Weight = lumi;
Double_t bkg3Weight = lumi;
Double_t bkg4Weight = lumi;
Double_t bkg5Weight = lumi;
Double_t bkg6Weight = lumi;
Double_t bkg7Weight = lumi;
Double_t bkg8Weight = lumi;
Double_t bkg9Weight = lumi;

// register the trees
factory->AddSignalTree    (signal,      sigWeight );
factory->AddBackgroundTree(background1, bkg1Weight);
factory->AddBackgroundTree(background2, bkg2Weight);
factory->AddBackgroundTree(background3, bkg3Weight);
factory->AddBackgroundTree(background4, bkg4Weight);
factory->AddBackgroundTree(background5, bkg5Weight);
factory->AddBackgroundTree(background6, bkg6Weight);
factory->AddBackgroundTree(background7, bkg7Weight);
factory->AddBackgroundTree(background8, bkg8Weight);
factory->AddBackgroundTree(background9, bkg9Weight);

// rather  than  having  only  global  weighting  factors  for  individual  input  trees  which  allow  to  scale
// them  to  the  same  luminosity,  individual  event  weights  can  be  applied  as  well.
// these  weights should be available event-by-event, i.e.  as a column or a function of columns of the input data sets.
// to specify the weights to be used for the training use the command:
//	factory->SetWeightExpression( "<YourWeightExpression>" );
// or if you have different expressions (variables) used as weights in the signal and background trees:
//	factory->SetSignalWeightExpression( "<YourSignalWeightExpression>" );
//	factory->SetBackgroundWeightExpression( "<YourBackgroundWeightExpression>" );
factory->SetWeightExpression( "baseW" ); // be sure it works properly!!   
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                           /////////////////////////////
/////////////////////////////  block 3: ADDVARIABLES  //////////////////////////////////////////////////////////////////////////// 
                           /////////////////////////////  

// specify input variables
// must be TFormula-compliant functions of branches in the training trees are registered. 
// it takes the variable name (string), which must have a correspondence in the input ROOT tree or input text file, and optionally a number type ('F' (default) and 'I'). (sec 3.1.3 TMVA Users Guide) 
// [be careful with the order: it must be respected at the application step (optim2.C)]
// factory->AddVariable( "<YourVar1>+<YourVar2>", "Pretty Title", "Unit", 'F' ); see Code Example 13, TMVA Users Guide

factory->AddVariable ("channel"  , "", "", 'F');
factory->AddVariable ("MET"      , "", "", 'F');
factory->AddVariable ("mll"      , "", "", 'F');
factory->AddVariable ("njet"     , "", "", 'F');
factory->AddVariable ("nbjet15"  , "", "", 'F');
factory->AddVariable ("lep1Pt"   , "", "", 'F');
factory->AddVariable ("lep2Pt"   , "", "", 'F');
factory->AddVariable ("jetpt1"   , "", "", 'F');
factory->AddVariable ("jetpt2"   , "", "", 'F');

factory->AddVariable ("dphilep1lep2"   , "", "", 'F');
factory->AddVariable ("dphilep1jet1"   , "", "", 'F');
factory->AddVariable ("dphilep1jet2"   , "", "", 'F');
factory->AddVariable ("dphilep1MET"    , "", "", 'F');
factory->AddVariable ("dphilep2jet1"   , "", "", 'F');
factory->AddVariable ("dphilep2jet2"   , "", "", 'F');
factory->AddVariable ("dphilep2MET"    , "", "", 'F');
factory->AddVariable ("dphijet1jet2"   , "", "", 'F');
factory->AddVariable ("dphijet1MET"    , "", "", 'F');
factory->AddVariable ("dphijet2MET"    , "", "", 'F');
factory->AddVariable ("dphilep1lep2MET", "", "", 'F');
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  


                           /////////////////////////////////////////
/////////////////////////////  BLOCK 4: PRESEL CUTS & PREPARATION //////////////////////////////////////////////////////////////// 
                           /////////////////////////////////////////

// it is possible to apply selection requirements (cuts) upon the input events.  
// these requirements can depend on any variable present in the input data sets, 
// i.e., they are not restricted to the variables used by the methods
// [they must be added to the application step (optim2.C) too]

//TCut mycut = "pfType1Met>120 && pt1>20.0 && pt2>20.0 && std_vector_lepton_eta[0]<2.4 && std_vector_lepton_eta[1]<2.4 && std_vector_lepton_flavour[0]*std_vector_lepton_flavour[1]<0 && njet>1 && mll>20 && ( channel>1 || (mll<76 || mll>106) )";   // latinos channel assignment: channel==0 mm, channel==1 ee, channel==2 em, channel==3 me
   
// prepare the training and test trees.
// if nTrainSignal=0 and nTest Signal=0,  the  signal  sample  is  split in  half  for training  and  testing.   
// the same rules apply to background.   
// since zero is default, not specifying anything corresponds to splitting the samples in two halves.
// the option SplitMode defines how the training and test samples are selected from the source trees.
// with SplitMode=Random,  events are selected randomly.  
// with SplitMode=Alternate,  events are chosen  in  alternating  turns  for  the  training  and  test  samples  
// as  they  occur  in  the  source  trees until the desired numbers of training and test events are selected. 
   
//factory->PrepareTrainingAndTestTree(mycut,     
//				":nTrain_Signal=0:nTest_Signal=0:nTrain_Background=2000:nTest_Background=2000:SplitMode=Alternate:!V");

//no cuts
factory->PrepareTrainingAndTestTree("",     
				":nTrain_Signal=0:nTest_Signal=0:nTrain_Background=2000:nTest_Background=2000:SplitMode=Alternate:!V");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  


                           //////////////////////////
/////////////////////////////  BLOCK 5: BOOK MVA  //////////////////////////////////////////////////////////////////////////////// 
                           //////////////////////////
 
// selected MVA methods are booked through a type identifer and a user-defined unique name, 
// and configuration options are specified via an option string. 
// for "rectangular cut optimisation" see Option Table 10, pg 61 TMVA Users Guide

//factory->BookMethod( TMVA::Types::kCuts, "Cuts", "FitMethod=GA:EffMethod=EffSel:VarProp[0]=FSmart:VarProp[1]=FSmart:VarProp[2]=FSmart:VarProp[3]=FSmart" );
//factory->BookMethod( TMVA::Types::kLD, "LD" );
//factory->BookMethod( TMVA::Types::kMLP, "MLP",  "H:!V:NeuronType=linear :VarTransform=N:NCycles=600:HiddenLayers=1:TestRate=5:!UseRegulator" );
factory->BookMethod( TMVA::Types::kMLP, "MLP",    "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=600:HiddenLayers=25,10:TestRate=5:!UseRegulator" );
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  


                           /////////////////////////////////////////////////////
/////////////////////////////  BLOCKS 6, 7 & 8: TRAIN, TEST & EVALUATE MVAS  ///////////////////////////////////////////////////// 
                           /////////////////////////////////////////////////////

// Train MVAs using the set of training events
factory->TrainAllMethods();

// ---- Evaluate all MVAs using the set of test events
factory->TestAllMethods();

// ----- Evaluate and compare performance of all configured MVAs
factory->EvaluateAllMethods();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Save the output
outputFile->Close();

cout << "==> Wrote root file: " << outputFile->GetName() << endl;
cout << "==> TMVAClassification is done!" << endl;

delete factory;


} // optim12 end
  
