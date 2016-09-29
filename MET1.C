#include "MET.h"

void MET1(){

	if( DoFillHistograms == false ) { cout << "   >>   FillHistograms deactivated !!!   " << endl;  return; }

	Assign(); 

	TH1::SetDefaultSumw2(); 

	FillHistograms(); 

}


void FillHistograms(){

	//for( int j = 0; j < nprocess; j++ ){
	for( int j = alpha; j < omega+1; j++ ){

		//if ( j == DoubleEG2016B     || j == DoubleEG2016C     || j == DoubleEG2016D     ) continue; 
		//if ( j == DoubleMuon2016B   || j == DoubleMuon2016C   || j == DoubleMuon2016D   ) continue; 
		//if ( j == SinglePhoton2016B || j == SinglePhoton2016C || j == SinglePhoton2016D ) continue; 
		//if ( j != DoubleEG2016E && j != DoubleEG2016F && j != DoubleMuon2016E && j != DoubleMuon2016F && j != SinglePhoton2016E && j != SinglePhoton2016F ) continue; 

		cout << "next n-tupla: " << processID[j] << endl; 

		FillHistogram( j );

		if( SkipEventLoop == false ){

			TFile* allhistos = new TFile("histograms/" + allhistosWriteTo +".root", "update");

			for( int s = 0; s < nsystematic; s++ ){ 

				for( int i = 0; i < nvariable_G; i++ ){

					h_global_W[i][j][s] -> Write();  

				}


				for( int i = 0; i < nvariable_R; i++ ){		
						
					for( int k = 0; k < nbinpT   ; k++) 	h_resol_pT_W   [i][j][k][s] -> Write();
					for( int k = 0; k < nbinsumET; k++)	h_resol_sumET_W[i][j][k][s] -> Write();
					for( int k = 0; k < nbinNVtx ; k++)	h_resol_NVtx_W [i][j][k][s] -> Write();

				}

			}

			allhistos -> Close();

		}

	}   // nprocess

	cout << "                          " << endl; 
	cout << "                          " << endl;
	cout << "    histograms filled !!! " << endl;
	cout << "                          " << endl;
	cout << "                          " << endl;

}



void FillHistogram( int process ){


	ofstream sync;   // for the sync with Leonora

	sync.open( "sync/160926_" + processID[process] + "_weights.txt" );

	//sync << "RUN       \t LUMI     \t EVT   \t Zll pT  \t lep-1 pT \t lep-2 pT \t MET       \n \n \n"; 

	sync << "baseW \t puW \t genW \t eventW (total)       \n \n \n"; 

	int TheKanal = kanal[process]; 


	TFile* kfactorFile = new TFile( "kfactor.root", "read");

	TF1*  kfactor_Gjets = (TF1*)  kfactorFile -> Get( "gj-40"    );
	TH1F* kfactor_ZnnG  = (TH1F*) kfactorFile -> Get( "znng-130" );  
	TH1F* kfactor_ZllG  = (TH1F*) kfactorFile -> Get( "zllg-130" );  


	TFile* f_pu = new TFile( "pileup/PU-reweighting_160926_definitivo_runs-E-F.root", "read");

	TH1F* h_pu[nkanal];

	h_pu[Zee  ] = (TH1F*) f_pu -> Get( "h_pu-reweighting_" + kanalID[Zee  ] );
	h_pu[Zmumu] = (TH1F*) f_pu -> Get( "h_pu-reweighting_" + kanalID[Zmumu] );
	h_pu[Gamma] = (TH1F*) f_pu -> Get( "h_pu-reweighting_" + kanalID[Gamma] );
 

	TFile* TheFile = new TFile( NTuplaDir[TheKanal] + sampleID[process] + "/METtree.root", "read"); 

	TH1F* TheSumGenWeight; float TheDenominator; 

	if( !isData[process] ) { TH1F* TheSumGenWeight = (TH1F*) TheFile -> Get( "SumGenWeights" ); TheDenominator = TheSumGenWeight -> Integral(); }  // cout << "total nentries from the histo = " << TheDenominator << endl;

	TTree* tree = (TTree*) TheFile -> Get( "METtree" );

	tree -> SetBranchStatus("*", 0);  // disable all branches


	UInt_t    run ; 
	UInt_t    lumi;
	ULong64_t evt ; 

	int nLepGood20                          ;
	float zll_pt                            ;
	float zll_phi                           ; 
	float zll_mass                          ; 

	int   ngamma                            ;
	int   gamma_idCutBased[3]               ;
	float gamma_hOverE[3]			;
	float gamma_sigmaIetaIeta[3]		;
	float gamma_r9[3]                       ;	
	float gamma_pt[3]                       ;
	float gamma_eta[3]                      ;
	float gamma_phi[3]                      ;

	int   nLepGood10                        ;  
	float lep_pt[4]                         ;  
	float lep_eta[4]                        ;  

	int   HBHENoiseFilter                   ;
	int   HBHENoiseIsoFilter                ;
	int   CSCTightHalo2015Filter            ;
	int   EcalDeadCellTriggerPrimitiveFilter;
	int   goodVertices                      ;
	int   eeBadScFilter                     ;
	int   globalTightHalo2016Filter         ;
	float badMuonFilter                     ;
	float badChargedHadronFilter            ;

	int HLT_Mu17_Mu8                        ;
	int HLT_Mu17_TkMu8                      ;
	int HLT_Ele23_Ele12                     ; 

	int HLT_DoubleMu	;
	int HLT_DoubleEG	;
	int HLT_SingleMu	;
	int HLT_SingleEle	;

	int HLT_Photon30                        ;
	int HLT_Photon50                        ;
	int HLT_Photon75                        ;
	int HLT_Photon90                        ;
	int HLT_Photon120                       ;

	int HLT_Photon30_Prescale               ;
	int HLT_Photon50_Prescale               ;
	int HLT_Photon75_Prescale               ;
	int HLT_Photon90_Prescale               ;
	int HLT_Photon120_Prescale              ;

	float met_sumEt                         ;
	int   nVert                             ;
	float met_pt                            ;
	float met_phi                           ;
	float met_jecUp_pt                      ;
	float met_jecUp_phi                     ; 
	float met_jecDown_pt                    ;
	float met_jecDown_phi                   ; 
	float met_UnEup_pt                      ;
	float met_UnEup_phi                     ;
	float met_UnEdown_pt                    ;
	float met_UnEdown_phi                   ;

	float puWeight                          ;
	float genWeight	                        ;

	int   nGenPart                          ; 
	int   GenPart_pdgId[40]                 ;
	float GenPart_pt[40]                    ;



	tree -> SetBranchStatus( "run" , 1 );	tree -> SetBranchAddress( "run" , &run  );
	tree -> SetBranchStatus( "lumi", 1 );	tree -> SetBranchAddress( "lumi", &lumi );
	tree -> SetBranchStatus( "evt" , 1 );	tree -> SetBranchAddress( "evt" , &evt  );


	if( TheKanal == Zee || TheKanal == Zmumu ){

		tree -> SetBranchStatus("nLepGood20", 1);	tree -> SetBranchAddress( "nLepGood20"                             , &nLepGood20                         );
		tree -> SetBranchStatus("zll_pt"    , 1);	tree -> SetBranchAddress( "zll_pt"                                 , &zll_pt                             );
		tree -> SetBranchStatus("zll_phi"   , 1);	tree -> SetBranchAddress( "zll_phi"                                , &zll_phi                            );
		tree -> SetBranchStatus("zll_mass"  , 1);	tree -> SetBranchAddress( "zll_mass"                               , &zll_mass                           );

	}

	if( TheKanal == Gamma ){

		tree -> SetBranchStatus("ngamma"             , 1);	tree -> SetBranchAddress( "ngamma"                                 , &ngamma                             );
		tree -> SetBranchStatus("gamma_idCutBased"   , 1);	tree -> SetBranchAddress( "gamma_idCutBased"                       , &gamma_idCutBased                   );
		tree -> SetBranchStatus("gamma_hOverE"       , 1);	tree -> SetBranchAddress( "gamma_hOverE"                           , &gamma_hOverE                       );
		tree -> SetBranchStatus("gamma_sigmaIetaIeta", 1);	tree -> SetBranchAddress( "gamma_sigmaIetaIeta"                    , &gamma_sigmaIetaIeta                );
		tree -> SetBranchStatus("gamma_r9"           , 1);	tree -> SetBranchAddress( "gamma_r9"                               , &gamma_r9                           );
		tree -> SetBranchStatus("gamma_pt"           , 1);	tree -> SetBranchAddress( "gamma_pt"                               , &gamma_pt                           );
		tree -> SetBranchStatus("gamma_eta"          , 1);	tree -> SetBranchAddress( "gamma_eta"                              , &gamma_eta                          );
		tree -> SetBranchStatus("gamma_phi"          , 1);	tree -> SetBranchAddress( "gamma_phi"                              , &gamma_phi                          );

		tree -> SetBranchStatus("nLepGood10"         , 1);	tree -> SetBranchAddress( "nLepGood10"                             , &nLepGood10                         );

	}

	tree -> SetBranchStatus("lep_pt"                                 , 1);		tree -> SetBranchAddress( "lep_pt"                                 , &lep_pt                             );
	tree -> SetBranchStatus("lep_eta"                                , 1);		tree -> SetBranchAddress( "lep_eta"                                , &lep_eta                            );

	tree -> SetBranchStatus("Flag_HBHENoiseFilter"                   , 1);		tree -> SetBranchAddress( "Flag_HBHENoiseFilter"                   , &HBHENoiseFilter                    );
	tree -> SetBranchStatus("Flag_HBHENoiseIsoFilter"                , 1);		tree -> SetBranchAddress( "Flag_HBHENoiseIsoFilter"                , &HBHENoiseIsoFilter                 );
	tree -> SetBranchStatus("Flag_CSCTightHalo2015Filter"            , 1);		tree -> SetBranchAddress( "Flag_CSCTightHalo2015Filter"            , &CSCTightHalo2015Filter             );
	tree -> SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter", 1);		tree -> SetBranchAddress( "Flag_EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter );
	tree -> SetBranchStatus("Flag_goodVertices"                      , 1);		tree -> SetBranchAddress( "Flag_goodVertices"                      , &goodVertices                       );
	tree -> SetBranchStatus("Flag_eeBadScFilter"                     , 1);		tree -> SetBranchAddress( "Flag_eeBadScFilter"                     , &eeBadScFilter                      );
	tree -> SetBranchStatus("Flag_globalTightHalo2016Filter"         , 1);		tree -> SetBranchAddress( "Flag_globalTightHalo2016Filter"         , &globalTightHalo2016Filter          );
	tree -> SetBranchStatus("Flag_badMuonFilter"                     , 1);		tree -> SetBranchAddress( "Flag_badMuonFilter"                     , &badMuonFilter                      );
	tree -> SetBranchStatus("Flag_badChargedHadronFilter"            , 1);		tree -> SetBranchAddress( "Flag_badChargedHadronFilter"            , &badChargedHadronFilter             );

	tree -> SetBranchStatus("HLT_Photon30" , 1);		tree -> SetBranchAddress( "HLT_Photon30"                           , &HLT_Photon30                       );
	tree -> SetBranchStatus("HLT_Photon50" , 1);		tree -> SetBranchAddress( "HLT_Photon50"                           , &HLT_Photon50                       );
	tree -> SetBranchStatus("HLT_Photon75" , 1);		tree -> SetBranchAddress( "HLT_Photon75"                           , &HLT_Photon75                       );
	tree -> SetBranchStatus("HLT_Photon90" , 1);		tree -> SetBranchAddress( "HLT_Photon90"                           , &HLT_Photon90                       );
	tree -> SetBranchStatus("HLT_Photon120", 1);		tree -> SetBranchAddress( "HLT_Photon120"                          , &HLT_Photon120                      );

	if( isData[process] ){

		tree -> SetBranchStatus("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"      , 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"      , &HLT_Mu17_Mu8    );
		tree -> SetBranchStatus("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"    , 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"    , &HLT_Mu17_TkMu8  );
		tree -> SetBranchStatus("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Ele23_Ele12 );

		tree -> SetBranchStatus("HLT_DoubleMu" , 1);	tree -> SetBranchAddress( "HLT_DoubleMu" , &HLT_DoubleMu  );
		tree -> SetBranchStatus("HLT_DoubleEG" , 1);	tree -> SetBranchAddress( "HLT_DoubleEG" , &HLT_DoubleEG  );
		tree -> SetBranchStatus("HLT_SingleMu" , 1);	tree -> SetBranchAddress( "HLT_SingleMu" , &HLT_SingleMu  );
		tree -> SetBranchStatus("HLT_SingleEle", 1);	tree -> SetBranchAddress( "HLT_SingleEle", &HLT_SingleEle );


		tree -> SetBranchStatus("HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v_Prescale" , 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon30_Prescale    );
		tree -> SetBranchStatus("HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale" , 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon50_Prescale    );
		tree -> SetBranchStatus("HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale" , 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon75_Prescale    );
		tree -> SetBranchStatus("HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale" , 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon90_Prescale    );
		tree -> SetBranchStatus("HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale", 1);	tree -> SetBranchAddress( "HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale"  , &HLT_Photon120_Prescale   );

	}

	tree -> SetBranchStatus("met_sumEt"      , 1);		tree -> SetBranchAddress( "met_sumEt"                              , &met_sumEt                          );
	tree -> SetBranchStatus("nVert"          , 1);		tree -> SetBranchAddress( "nVert"                                  , &nVert                              );
	tree -> SetBranchStatus("met_pt"         , 1);		tree -> SetBranchAddress( "met_pt"                                 , &met_pt                             );
	tree -> SetBranchStatus("met_phi"        , 1);		tree -> SetBranchAddress( "met_phi"                                , &met_phi                            );
	tree -> SetBranchStatus("met_jecUp_pt"   , 1);		tree -> SetBranchAddress( "met_jecUp_pt"                           , &met_jecUp_pt                       );
	tree -> SetBranchStatus("met_jecUp_phi"  , 1);		tree -> SetBranchAddress( "met_jecUp_phi"                          , &met_jecUp_phi                      );
	tree -> SetBranchStatus("met_jecDown_pt" , 1);		tree -> SetBranchAddress( "met_jecDown_pt"                         , &met_jecDown_pt                     );
	tree -> SetBranchStatus("met_jecDown_phi", 1);		tree -> SetBranchAddress( "met_jecDown_phi"                        , &met_jecDown_phi                    );

	tree -> SetBranchStatus("met_shifted_UnclusteredEnUp_pt"   , 1); tree -> SetBranchAddress( "met_shifted_UnclusteredEnUp_pt"   , &met_UnEup_pt   );
	tree -> SetBranchStatus("met_shifted_UnclusteredEnUp_phi"  , 1); tree -> SetBranchAddress( "met_shifted_UnclusteredEnUp_phi"  , &met_UnEup_phi  );
	tree -> SetBranchStatus("met_shifted_UnclusteredEnDown_pt" , 1); tree -> SetBranchAddress( "met_shifted_UnclusteredEnDown_pt" , &met_UnEdown_pt );
	tree -> SetBranchStatus("met_shifted_UnclusteredEnDown_phi", 1); tree -> SetBranchAddress( "met_shifted_UnclusteredEnDown_phi", &met_UnEdown_phi);

	
	if(  !isData[process] ){

		tree -> SetBranchStatus("puWeight" , 1);	tree -> SetBranchAddress( "puWeight"                               , &puWeight                           );
		tree -> SetBranchStatus("genWeight", 1);	tree -> SetBranchAddress( "genWeight"                              , &genWeight                          );

	}

	if( process == GJets40100 || process == GJets100200 || process == GJets200400|| process == GJets400600 || process == GJets600Inf || process == ZGJets || process == ZGTo2LG ){

			tree -> SetBranchStatus("nGenPart"     , 1);	tree -> SetBranchAddress( "nGenPart"     , &nGenPart      );
			tree -> SetBranchStatus("GenPart_pdgId", 1);	tree -> SetBranchAddress( "GenPart_pdgId", &GenPart_pdgId );
			tree -> SetBranchStatus("GenPart_pt"   , 1);	tree -> SetBranchAddress( "GenPart_pt"   , &GenPart_pt    );

	}
	

	for( int s = 0; s < nsystematic; s++ ){
	for( int i = 0; i < nvariable_G; i++ ){

		if( i == parall_R ) h_global_W[i][process][s] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process] + "_" + systematicID[s], variableID_G[i], nbinuPara, minuPara, maxuPara );
		if( i == transv_R ) h_global_W[i][process][s] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process] + "_" + systematicID[s], variableID_G[i], nbinuPerp, minuPerp, maxuPerp );
		if( i == MET      ) h_global_W[i][process][s] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process] + "_" + systematicID[s], variableID_G[i], nbinMET  , minMET  , maxMET   );
		if( i == VpT      ) h_global_W[i][process][s] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process] + "_" + systematicID[s], variableID_G[i], nbinpT   , minpT   , maxpT    );
		if( i == nvert    ) h_global_W[i][process][s] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process] + "_" + systematicID[s], variableID_G[i], nbinnvert, minnvert, maxnvert );
		
	}  // variable_G	



		for( int k = 0; k < nbinpT; k++ ){

			h_resol_pT_W   [parall_R][process][k][s] = new TH1F( Form("h_resol_pT_" + variableID_R[parall_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution pT"        , nbinuPara, minuPara, maxuPara );

			h_resol_pT_W   [transv_R][process][k][s] = new TH1F( Form("h_resol_pT_" + variableID_R[transv_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution pT"        , nbinuPerp, minuPerp, maxuPerp );

			h_resol_pT_W   [scale   ][process][k][s] = new TH1F( Form("h_resol_pT_" + variableID_R[scale   ] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution pT"        , nbinscale, minscale, maxscale );

		}

		for( int k = 0; k < nbinsumET; k++ ){

			h_resol_sumET_W[parall_R][process][k][s] = new TH1F( Form("h_resol_sumET_" + variableID_R[parall_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution sum ET"    , nbinuPara, minuPara, maxuPara );

			h_resol_sumET_W[transv_R][process][k][s] = new TH1F( Form("h_resol_sumET_" + variableID_R[transv_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution sum ET"    , nbinuPerp, minuPerp, maxuPerp );

			h_resol_sumET_W[scale   ][process][k][s] = new TH1F( Form("h_resol_sumET_" + variableID_R[scale   ] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution sum ET"    , nbinscale, minscale, maxscale );

		}

		for( int k = 0; k < nbinNVtx; k++ ){

			h_resol_NVtx_W [parall_R][process][k][s] = new TH1F( Form("h_resol_NVtx_"  + variableID_R[parall_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution num of vtx", nbinuPara, minuPara, maxuPara );

			h_resol_NVtx_W [transv_R][process][k][s] = new TH1F( Form("h_resol_NVtx_"  + variableID_R[transv_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution num of vtx", nbinuPerp, minuPerp, maxuPerp );

			h_resol_NVtx_W [scale   ][process][k][s] = new TH1F( Form("h_resol_NVtx_"  + variableID_R[scale   ] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d_" + systematicID[s], k), "resolution num of vtx", nbinscale, minscale, maxscale );

		}

	}  // systematic 


	int nentries = tree -> GetEntries(); //cout << " nentries = " << nentries << endl;

	_TotalEntries[process] = nentries; 

	if (  isData[process]  &&  RunOverAllData == false ) nentries = MaxEntriesData; 
		
	if ( !isData[process]  &&  RunOverAllMC   == false ) nentries = MaxEntriesMC  ;

	float baseW = 1.0;   if ( !isData[process] ) baseW = xs[process]/TheDenominator;  //cout << xs[process] << " -- " << nentries << " -- " << TheDenominator << " -- " << baseW << endl; 

	if ( SkipEventLoop == true ) return; 
 
	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {



		if( ievt%100000 == 0 )  cout << processID[process] << "   >>  evt " << ievt << endl;

		tree -> GetEntry(ievt);	

		float eventW = 1.0 ;



		if ( TheKanal == Zee || TheKanal == Zmumu ){

			if ( nLepGood20 != 2                     ) continue;

			if (  abs(lep_eta[0]) > 2.4             ||  abs(lep_eta[1]) > 2.4             ) continue; 	 

			if (  abs( abs(lep_eta[0])-1.5 ) < 0.1  ||  abs( abs(lep_eta[1])-1.5 ) < 0.1  ) continue;


			if ( zll_mass < 81. || zll_mass > 101.   ) continue; 

			if ( lep_pt[0] < 20. || lep_pt[1] < 20.  ) continue; 	//cout << lep_pt[0] << " -- " << lep_pt[1] << endl; 

			//if ( zll_pt < 50.                        ) continue;    //cout << zll_pt << endl;



			if( isData[process] ){   // triggers

				if( HLT_DoubleMu == 0 && HLT_DoubleEG == 0 && HLT_SingleMu == 0 && HLT_SingleEle == 0 ) continue;

			}
	
		}



		if ( TheKanal == Gamma ){ 

			if ( ngamma != 1                            ) continue;

			//if ( gamma_idCutBased[0] != 3                                      ) continue; 
			if ( gamma_hOverE[0] > 0.05                                          ) continue; 
			if ( gamma_sigmaIetaIeta[0] < 0.005 || gamma_sigmaIetaIeta[0] > 0.01 ) continue; 

			if ( gamma_r9[0] < 0.9 || gamma_r9[0] > 1.0 ) continue;
			if ( gamma_pt[0] < 50                       ) continue; 
			if ( abs(gamma_eta[0]) > 1.4442             ) continue; 

			if ( nLepGood10 != 0                        ) continue; 

			if( isData[process] ){   // triggers


				if( HLT_Photon30 == 0  &&  HLT_Photon50 == 0  &&  HLT_Photon75 == 0  &&  HLT_Photon90 == 0  &&  HLT_Photon120 == 0 ) continue;


				if( HLT_Photon120 == 1 ) eventW *= HLT_Photon120_Prescale;

				else{

					if( HLT_Photon90 == 1 )	eventW *= HLT_Photon90_Prescale;

					else{

						if( HLT_Photon75 == 1 )	eventW *= HLT_Photon75_Prescale;

						else{

							if( HLT_Photon50 == 1 )	eventW *= HLT_Photon50_Prescale;

							else{

								if( HLT_Photon30 == 1 )	eventW *= HLT_Photon30_Prescale;


							}

						}

					}
	
				}

			}

		}



		// MET filters
		// -----------------------------------------------------
		if ( HBHENoiseFilter                    != 1 ) continue; 
		if ( HBHENoiseIsoFilter                 != 1 ) continue; 
		//if ( CSCTightHalo2015Filter             != 1 ) continue; 
		if ( EcalDeadCellTriggerPrimitiveFilter != 1 ) continue; 
		if ( goodVertices                       != 1 ) continue; 
		if ( eeBadScFilter                      != 1 ) continue; 	
		if ( globalTightHalo2016Filter          != 1 ) continue;
		if ( badMuonFilter                      != 1 ) continue;
		if ( badChargedHadronFilter             != 1 ) continue;
		// -----------------------------------------------------


		// ----- sync -------------------------			



		//if(  ( process == DoubleEG2016B  &&  run == 273158 )  ||  ( process == DoubleMuon2016B  &&  run == 273158 )  ){



		//if(  process == DY_mm  &&  ( lumi > 153860 && lumi <= 155000 ) ){

		//	sync << Form( "%12u * %12u * %12llu * %12f * %12f * %12f * %12f *\n", run, lumi, evt, zll_pt, lep_pt[0], lep_pt[1], met_pt );

		//}




		//} 

		//if(  process == DoubleEG2016C    &&  run == 275658  )
		//if(  process == DoubleEG2016D    &&  run == 276315  ) 
		

		// ------------------------------------

 
		if( !isData[process] ){

			eventW *= baseW    ; 	

			//cout << nVert << " -- " << h_pu[TheKanal] -> GetBinContent( h_pu[TheKanal]->FindBin(nVert) ) << endl;	

			//int nVert_rw;  ( nVert >= 38 )  ?  nVert_rw = 38  :  nVert_rw = nVert; 

			float puW = h_pu[TheKanal] -> GetBinContent( h_pu[TheKanal]->FindBin(nVert) );

			eventW *= puW;   		

			eventW *= genWeight;// /abs(genWeight);

		}



		// ----- k-factor ----------------------------------------
 
		if( process == GJets40100 || process == GJets100200 || process == GJets200400|| process == GJets400600 || process == GJets600Inf || process == ZGJets || process == ZGTo2LG ){

			int genindex = -1; 

      			for( int j = 0; j < nGenPart; j++){

	  			if( abs( GenPart_pdgId[j] ) == 22 ){

	      				genindex = j;

	      				break;
	    
				}

			}

			float gamma_pt_GEN = GenPart_pt[genindex];

			if( process == GJets40100 || process == GJets100200 || process == GJets200400|| process == GJets400600 || process == GJets600Inf ){

				eventW *= kfactor_Gjets -> Eval( gamma_pt_GEN );

			}

			else if( process == ZGJets ){

				float TheFactor = kfactor_ZnnG -> GetBinContent( kfactor_ZnnG->FindBin(gamma_pt_GEN) ); 

				if( TheFactor > 1. )  eventW *= TheFactor; 

			}

			else if( process == ZGTo2LG ){

				float TheFactor = kfactor_ZllG -> GetBinContent( kfactor_ZllG->FindBin(gamma_pt_GEN) ); 

				if( TheFactor > 1. )  eventW *= TheFactor; 

			}

			
		}

		// -------------------------------------------------------

		float boson_pt  = -999.;
		float boson_phi = -999.;	
	
		if ( TheKanal == Zee || TheKanal == Zmumu ) { boson_pt = zll_pt     ;  boson_phi = zll_phi     ; }   

		if ( TheKanal == Gamma )                    { boson_pt = gamma_pt[0];  boson_phi = gamma_phi[0]; }


		TVector2 qT, ET, uT; 

		float _met_pt, _met_phi; 

		for( int s = 0; s < nsystematic; s++ ){

			if( s == nominal ) { _met_pt = met_pt        ; _met_phi = met_phi;         }
			if( s == JECu    ) { _met_pt = met_jecUp_pt  ; _met_phi = met_jecUp_phi;   }
			if( s == JECd    ) { _met_pt = met_jecDown_pt; _met_phi = met_jecDown_phi; }
			if( s == UnEu    ) { _met_pt = met_UnEup_pt  ; _met_phi = met_UnEup_phi;   }
			if( s == UnEd    ) { _met_pt = met_UnEdown_pt; _met_phi = met_UnEdown_phi; }


			ET.SetMagPhi( _met_pt  , _met_phi   );

			qT.SetMagPhi( boson_pt, boson_phi );

			uT = -1* ( ET + qT ); 


			float uPara = (  uT.Px() * qT.Px() + uT.Py() * qT.Py()  ) / qT.Mod();

			float Scale = -1. * uPara/qT.Mod();   // cout << abs(uPara) << " -- " << qT.Mod() <<  " -- " << Scale << endl; 
	 
			uPara += qT.Mod();   // cout <<  uPara << endl;

			float uPerp = (  uT.Px() * qT.Py() - uT.Py() * qT.Px()  ) / qT.Mod();


			if(run==273158 && s==nominal) sync << Form( "*%12llu*%12f*%12f*\n", evt, uPara, uPerp );


			h_global_W     [parall_G][process][s]    -> Fill( uPara   , eventW );  
			h_global_W     [transv_G][process][s]    -> Fill( uPerp   , eventW );
			h_global_W     [MET     ][process][s]    -> Fill( _met_pt , eventW );
			h_global_W     [VpT     ][process][s]    -> Fill( boson_pt, eventW );
			h_global_W     [nvert   ][process][s]    -> Fill( nVert   , eventW );


			if( boson_pt < minpT || boson_pt >= maxpT || met_sumEt < minsumET || met_sumEt >= maxsumET || nVert < minNVtx || nVert >= maxNVtx) continue;

			int l = floor(    nbinpT    * (boson_pt    - minpT   ) / (maxpT    - minpT   )    ); //cout << " l = " << l << endl;
			int m = floor(    nbinsumET * (met_sumEt   - minsumET) / (maxsumET - minsumET)    ); //cout << " m = " << m << endl;
			int n = floor(    nbinNVtx  * (nVert       - minNVtx ) / (maxNVtx  - minNVtx )    ); //cout << " n = " << n << endl;

			h_resol_pT_W   [parall_R][process][l][s] -> Fill( uPara , eventW ); 
			h_resol_pT_W   [transv_R][process][l][s] -> Fill( uPerp , eventW );	
			h_resol_pT_W   [scale   ][process][l][s] -> Fill( Scale , eventW );	

			h_resol_sumET_W[parall_R][process][m][s] -> Fill( uPara , eventW ); 
			h_resol_sumET_W[transv_R][process][m][s] -> Fill( uPerp , eventW );
			h_resol_sumET_W[scale   ][process][m][s] -> Fill( Scale , eventW );

			h_resol_NVtx_W [parall_R][process][n][s] -> Fill( uPara , eventW );
			h_resol_NVtx_W [transv_R][process][n][s] -> Fill( uPerp , eventW );
			h_resol_NVtx_W [scale   ][process][n][s] -> Fill( Scale , eventW );

		}   // systematic

	}   // ievt

	sync.close();

}
