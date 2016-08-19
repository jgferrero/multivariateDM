#include "MET.h"

void MET1(){

	if( DoFillHistograms == false ) { cout << "   >>   FillHistograms deactivated !!!   " << endl;  return; }

	Assign(); 

	FillHistograms(); 

}


void FillHistograms(){

	for( int j = 0; j < nprocess; j++ ){
	//for( int j = alpha; j < omega+1; j++ ){

		cout << "next n-tupla: " << processID[j] << endl; 

		FillHistogram( j );

	}


	if( SkipEventLoop == false ){

		cout << "the writing step begins..." << endl; 

		TFile* allhistos = new TFile("histograms/" + allhistosWriteTo +".root", "update");

		for( int j = 0; j < nprocess; j++ ){
		//for( int j = alpha; j < omega+1; j++ ){ 

			for( int i = 0; i < nvariable_G; i++ ){

				h_global_W[i][j] -> Write();  

			}


			for( int i = 0; i < nvariable_R; i++ ){		
						
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

	int TheKanal = kanal[process]; 


	TFile* kfactorFile = new TFile( "kfactor.root", "read");

	TF1*  kfactor_Gjets = (TF1*)  kfactorFile -> Get( "gj-40"    );
	TH1F* kfactor_ZnnG  = (TH1F*) kfactorFile -> Get( "znng-130" );  
	TH1F* kfactor_ZllG  = (TH1F*) kfactorFile -> Get( "zllg-130" );  


	TFile* f_pu = new TFile( "pileup/PU-reweighting_15Aug_old-xs.root", "read");

	TH1F* h_pu[nkanal];

	h_pu[Zee  ] = (TH1F*) f_pu -> Get( "h_pu-reweighting_" + kanalID[Zee  ] );
	h_pu[Zmumu] = (TH1F*) f_pu -> Get( "h_pu-reweighting_" + kanalID[Zmumu] );
	h_pu[Gamma] = (TH1F*) f_pu -> Get( "h_pu-reweighting_" + kanalID[Gamma] );
 

	TFile* TheFile = new TFile( NTuplaDir[TheKanal] + sampleID[process] + "/METtree.root", "read"); 

	TH1F* TheCount = (TH1F*) TheFile -> Get( "Count" );  float TheDenominator = TheCount -> GetEntries(); // cout << "total nentries from the histo = " << TheDenominator << endl;

	TTree* tree = (TTree*) TheFile -> Get( "METtree" );


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

	float puWeight                          ;
	float genWeight	                        ;

	int   nGenPart                          ; 
	int   GenPart_pdgId[40]                 ;
	float GenPart_pt[40]                    ;


	if( TheKanal == Zee || TheKanal == Zmumu ){

		tree -> SetBranchAddress( "nLepGood20"                             , &nLepGood20                         );
		tree -> SetBranchAddress( "zll_pt"                                 , &zll_pt                             );
		tree -> SetBranchAddress( "zll_phi"                                , &zll_phi                            );
		tree -> SetBranchAddress( "zll_mass"                               , &zll_mass                           );

	}

	if( TheKanal == Gamma ){

		tree -> SetBranchAddress( "ngamma"                                 , &ngamma                             );
		tree -> SetBranchAddress( "gamma_idCutBased"                       , &gamma_idCutBased                   );
		tree -> SetBranchAddress( "gamma_hOverE"                           , &gamma_hOverE                       );
		tree -> SetBranchAddress( "gamma_sigmaIetaIeta"                    , &gamma_sigmaIetaIeta                );
		tree -> SetBranchAddress( "gamma_r9"                               , &gamma_r9                           );
		tree -> SetBranchAddress( "gamma_pt"                               , &gamma_pt                           );
		tree -> SetBranchAddress( "gamma_eta"                              , &gamma_eta                          );
		tree -> SetBranchAddress( "gamma_phi"                              , &gamma_phi                          );

		tree -> SetBranchAddress( "nLepGood10"                             , &nLepGood10                         );
	}

	tree -> SetBranchAddress( "lep_pt"                                 , &lep_pt                             );


	tree -> SetBranchAddress( "Flag_HBHENoiseFilter"                   , &HBHENoiseFilter                    );
	tree -> SetBranchAddress( "Flag_HBHENoiseIsoFilter"                , &HBHENoiseIsoFilter                 );
	tree -> SetBranchAddress( "Flag_CSCTightHalo2015Filter"            , &CSCTightHalo2015Filter             );
	tree -> SetBranchAddress( "Flag_EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter );
	tree -> SetBranchAddress( "Flag_goodVertices"                      , &goodVertices                       );
	tree -> SetBranchAddress( "Flag_eeBadScFilter"                     , &eeBadScFilter                      );
	tree -> SetBranchAddress( "Flag_globalTightHalo2016Filter"         , &globalTightHalo2016Filter          );
	tree -> SetBranchAddress( "Flag_badMuonFilter"                     , &badMuonFilter                      );
	tree -> SetBranchAddress( "Flag_badChargedHadronFilter"            , &badChargedHadronFilter             );

	tree -> SetBranchAddress( "HLT_Photon30"                           , &HLT_Photon30                       );
	tree -> SetBranchAddress( "HLT_Photon50"                           , &HLT_Photon50                       );
	tree -> SetBranchAddress( "HLT_Photon75"                           , &HLT_Photon75                       );
	tree -> SetBranchAddress( "HLT_Photon90"                           , &HLT_Photon90                       );
	tree -> SetBranchAddress( "HLT_Photon120"                          , &HLT_Photon120                      );

	if( isData[process] ){

		tree -> SetBranchAddress( "HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"      , &HLT_Mu17_Mu8    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"    , &HLT_Mu17_TkMu8  );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Ele23_Ele12 );

		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon30_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon30_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon50_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon50_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon75_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon75_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon90_R9Id90_HE10_IsoM_v_Prescale"   , &HLT_Photon90_Prescale    );
		tree -> SetBranchAddress( "HLT_BIT_HLT_Photon120_R9Id90_HE10_IsoM_v_Prescale"  , &HLT_Photon120_Prescale   );

	}

	tree -> SetBranchAddress( "met_sumEt"                              , &met_sumEt                          );
	tree -> SetBranchAddress( "nVert"                                  , &nVert                              );
	tree -> SetBranchAddress( "met_pt"                                 , &met_pt                             );
	tree -> SetBranchAddress( "met_phi"                                , &met_phi                            );

	if(  !isData[process] ){

		tree -> SetBranchAddress( "puWeight"                               , &puWeight                           );
		tree -> SetBranchAddress( "genWeight"                              , &genWeight                          );

	}

	if( process == GJets40100 || process == GJets100200 || process == GJets200400|| process == GJets400600 || process == GJets600Inf || process == ZGJets || process == ZGTo2LG ){

			tree -> SetBranchAddress( "nGenPart"     , &nGenPart      );
			tree -> SetBranchAddress( "GenPart_pdgId", &GenPart_pdgId );
			tree -> SetBranchAddress( "GenPart_pt"   , &GenPart_pt    );

	}
	


	for( int i = 0; i < nvariable_G; i++ ){

		if( i == parall_R ) h_global_W[i][process] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process], variableID_G[i], nbinuPara, minuPara, maxuPara );
		if( i == transv_R ) h_global_W[i][process] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process], variableID_G[i], nbinuPerp, minuPerp, maxuPerp );
		if( i == MET      ) h_global_W[i][process] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process], variableID_G[i], nbinMET  , minMET  , maxMET   );
		if( i == VpT      ) h_global_W[i][process] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process], variableID_G[i], nbinpT   , minpT   , maxpT    );
		if( i == nvert    ) h_global_W[i][process] = new TH1F( "h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[process], variableID_G[i], nbinnvert, minnvert, maxnvert );

	}	

	//for( int i = 0; i < nvariable_R; i++ ){

		for( int k = 0; k < nbinpT; k++ ){

			h_resol_pT_W   [parall_R][process][k] = new TH1F( Form("h_resol_pT_" + variableID_R[parall_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution pT"        , nbinuPara, minuPara, maxuPara );

			h_resol_pT_W   [transv_R][process][k] = new TH1F( Form("h_resol_pT_" + variableID_R[transv_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution pT"        , nbinuPerp, minuPerp, maxuPerp );

			h_resol_pT_W   [scale   ][process][k] = new TH1F( Form("h_resol_pT_" + variableID_R[scale   ] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution pT"        , nbinscale, minscale, maxscale );

		}

		for( int k = 0; k < nbinsumET; k++ ){

			h_resol_sumET_W[parall_R][process][k] = new TH1F( Form("h_resol_sumET_" + variableID_R[parall_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution sum ET"    , nbinuPara, minuPara, maxuPara );

			h_resol_sumET_W[transv_R][process][k] = new TH1F( Form("h_resol_sumET_" + variableID_R[transv_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution sum ET"    , nbinuPerp, minuPerp, maxuPerp );

			h_resol_sumET_W[scale   ][process][k] = new TH1F( Form("h_resol_sumET_" + variableID_R[scale   ] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution sum ET"    , nbinscale, minscale, maxscale );

		}

		for( int k = 0; k < nbinNVtx; k++ ){

			h_resol_NVtx_W [parall_R][process][k] = new TH1F( Form("h_resol_NVtx_"  + variableID_R[parall_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution num of vtx", nbinuPara, minuPara, maxuPara );

			h_resol_NVtx_W [transv_R][process][k] = new TH1F( Form("h_resol_NVtx_"  + variableID_R[transv_R] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution num of vtx", nbinuPerp, minuPerp, maxuPerp );

			h_resol_NVtx_W [scale   ][process][k] = new TH1F( Form("h_resol_NVtx_"  + variableID_R[scale   ] + "_" + kanalID[TheKanal] + "_"  + processID[process] + "_%d", k), "resolution num of vtx", nbinscale, minscale, maxscale );

		}

	//}


	int nentries = tree -> GetEntries();//cout << " nentries = " << nentries << endl;

	_TotalEntries[process] = nentries; 

	if (  isData[process]  &&  RunOverAllData == false ) nentries = MaxEntriesData; 
		
	if ( !isData[process]  &&  RunOverAllMC   == false ) nentries = MaxEntriesMC  ;

	float baseW = 1.0;   if ( !isData[process] ) baseW = xs[process]/TheDenominator;  //cout << xs[process] << " -- " << nentries << " -- " << TheDenominator << " -- " << baseW << endl; 

	if ( SkipEventLoop == true ) return; 
 
	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {



		if( ievt%10000 == 0 )  cout << processID[process] << "   >>  evt " << ievt << endl;

		tree -> GetEntry(ievt);	

		float eventW = 1.0 ;



		if ( TheKanal == Zee || TheKanal == Zmumu ){

			if ( nLepGood20 < 2                    ) continue; 

			if ( zll_mass < 81. || zll_mass > 101. ) continue; 

			if ( lep_pt[0] < 25 || lep_pt[1] < 20  ) continue; 	//cout << lep_pt[0] << " -- " << lep_pt[1] << endl; 


			if( isData[process] ){   // triggers

				if( HLT_Mu17_Mu8 == 0 && HLT_Mu17_TkMu8 == 0 && HLT_Ele23_Ele12 == 0 ) continue;

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
		if ( HBHENoiseFilter                    != 1 ) continue; 
		if ( HBHENoiseIsoFilter                 != 1 ) continue; 
		//if ( CSCTightHalo2015Filter             != 1 ) continue; 
		if ( EcalDeadCellTriggerPrimitiveFilter != 1 ) continue; 
		if ( goodVertices                       != 1 ) continue; 
		if ( eeBadScFilter                      != 1 ) continue; 	
		if ( globalTightHalo2016Filter          != 1 ) continue;
		if ( badMuonFilter                      != 1 ) continue;
		if ( badChargedHadronFilter             != 1 ) continue;


 
		if( !isData[process] ){

			eventW *= baseW    ; 	

			//cout << nVert << " -- " << h_pu[TheKanal] -> GetBinContent( h_pu[TheKanal]->FindBin(nVert) ) << endl;	

			eventW *= puWeight * h_pu[TheKanal] -> GetBinContent( h_pu[TheKanal]->FindBin(nVert) );  		

			eventW *= genWeight/abs(genWeight);

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

		float boson_pt, boson_phi;	
	
		if ( TheKanal == Zee || TheKanal == Zmumu ) boson_pt = zll_pt     ;  boson_phi = zll_phi     ;   

		if ( TheKanal == Gamma )                    boson_pt = gamma_pt[0];  boson_phi = gamma_phi[0];


		TVector2 qT, ET, uT; 

		ET.SetMagPhi( met_pt  , met_phi   );

		qT.SetMagPhi( boson_pt, boson_phi );

		uT = -1* ( ET + qT ); 


		float uPara = (  uT.Px() * qT.Px() + uT.Py() * qT.Py()  ) / qT.Mod();

		float Scale = abs(uPara)/qT.Mod();   // cout << abs(uPara) << " -- " << qT.Mod() <<  " -- " << Scale << endl; 
 
		uPara += qT.Mod();   // cout <<  uPara << endl;

		float uPerp = (  uT.Px() * qT.Py() - uT.Py() * qT.Px()  ) / qT.Mod();



		h_global_W     [parall_G][process]    -> Fill( uPara   , eventW );  
		h_global_W     [transv_G][process]    -> Fill( uPerp   , eventW );
		h_global_W     [MET     ][process]    -> Fill( met_pt  , eventW );
		h_global_W     [VpT     ][process]    -> Fill( boson_pt, eventW );
		h_global_W     [nvert   ][process]    -> Fill( nVert   , eventW );


		if( boson_pt < minpT || boson_pt >= maxpT || met_sumEt < minsumET || met_sumEt >= maxsumET || nVert < minNVtx || nVert >= maxNVtx) continue;

		int l = floor(    nbinpT    * (boson_pt    - minpT   ) / (maxpT    - minpT   )    ); //cout << " l = " << l << endl;
		int m = floor(    nbinsumET * (met_sumEt   - minsumET) / (maxsumET - minsumET)    ); //cout << " m = " << m << endl;
		int n = floor(    nbinNVtx  * (nVert       - minNVtx ) / (maxNVtx  - minNVtx )    ); //cout << " n = " << n << endl;

		h_resol_pT_W   [parall_R][process][l] -> Fill( uPara , eventW ); 
		h_resol_pT_W   [transv_R][process][l] -> Fill( uPerp , eventW );	
		h_resol_pT_W   [scale   ][process][l] -> Fill( Scale , eventW );	

		h_resol_sumET_W[parall_R][process][m] -> Fill( uPara , eventW ); 
		h_resol_sumET_W[transv_R][process][m] -> Fill( uPerp , eventW );
		h_resol_sumET_W[scale   ][process][m] -> Fill( Scale , eventW );

		h_resol_NVtx_W [parall_R][process][n] -> Fill( uPara , eventW );
		h_resol_NVtx_W [transv_R][process][n] -> Fill( uPerp , eventW );
		h_resol_NVtx_W [scale   ][process][n] -> Fill( Scale , eventW );


	}

}
