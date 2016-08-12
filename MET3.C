#include "MET.h"

void MET3(){

	if( DoFit == false ) { cout << "   >>   OrderFits deactivated !!!   " << endl;  return; };

	Assign(); 


	// load histograms

	TFile* allhistos = new TFile( "histograms/" + allhistosReadFrom + ".root", "read" );

	for( int i = 0; i < nvariable; i++ ){

		if( i == MET ) continue; 

		for( int j = 0; j < nprocess; j++ ){

			int TheKanal = kanal[j];

			for( int k = 0; k < nbinpT   ; k++){

				h_resol_pT   [i][j][k] = (TH1F*) allhistos -> Get( Form("h_resol_pT_"    + variableID[i] + "_" + kanalID[TheKanal] + "_" + processID[j] + "_%d", k) );
			}

			for( int k = 0; k < nbinsumET; k++){	

				h_resol_sumET[i][j][k] = (TH1F*) allhistos -> Get( Form("h_resol_sumET_" + variableID[i] + "_" + kanalID[TheKanal] + "_"+ processID[j] + "_%d", k) );

			}

			for( int k = 0; k < nbinNVtx ; k++){	

				h_resol_NVtx [i][j][k] = (TH1F*) allhistos -> Get( Form("h_resol_NVtx_"  + variableID[i] + "_" + kanalID[TheKanal] + "_"+ processID[j] + "_%d", k) );

			}

		}

	}




	// create all-background template for the fit procedure

	for( int k = 0; k < nbinpT; k++ ){

		for( int m = 0; m < nkanal; m++ ){

			h_resol_pT_fit   [parall][0][k][m] = new TH1F( Form("h_resol_pT_"  + variableID[parall] + "_" + kanalID[m] + "_data_%d",   k), "data   ", nbinuPara, minuPara, maxuPara );
 			h_resol_pT_fit   [parall][1][k][m] = new TH1F( Form("h_resol_pT_"  + variableID[parall] + "_" + kanalID[m] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara );
 			h_resol_pT_fit   [parall][2][k][m] = new TH1F( Form("h_resol_pT_"  + variableID[parall] + "_" + kanalID[m] + "_signal_%d", k), "signal ", nbinuPara, minuPara, maxuPara );
			h_resol_pT_fit   [transv][0][k][m] = new TH1F( Form("h_resol_pT_"  + variableID[transv] + "_" + kanalID[m] + "_data_%d",   k), "data   ", nbinuPerp, minuPerp, maxuPerp ); 
			h_resol_pT_fit   [transv][1][k][m] = new TH1F( Form("h_resol_pT_"  + variableID[transv] + "_" + kanalID[m] + "_allbkg_%d", k), "all bkg", nbinuPerp, minuPerp, maxuPerp ); 
			h_resol_pT_fit   [transv][2][k][m] = new TH1F( Form("h_resol_pT_"  + variableID[transv] + "_" + kanalID[m] + "_signal_%d", k), "signal ", nbinuPerp, minuPerp, maxuPerp ); 
		}

	}

	for( int k = 0; k < nbinsumET; k++ ){

		for( int m = 0; m < nkanal; m++ ){

			h_resol_sumET_fit[parall][0][k][m] = new TH1F( Form("h_resol_sumET_" + variableID[parall] + "_" + kanalID[m] + "_data_%d",   k), "data   ", nbinuPara, minuPara, maxuPara );
	 		h_resol_sumET_fit[parall][1][k][m] = new TH1F( Form("h_resol_sumET_" + variableID[parall] + "_" + kanalID[m] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara );
 			h_resol_sumET_fit[parall][2][k][m] = new TH1F( Form("h_resol_sumET_" + variableID[parall] + "_" + kanalID[m] + "_signal_%d", k), "signal ", nbinuPara, minuPara, maxuPara );
			h_resol_sumET_fit[transv][0][k][m] = new TH1F( Form("h_resol_sumET_" + variableID[transv] + "_" + kanalID[m] + "_data_%d",   k), "data   ", nbinuPerp, minuPerp, maxuPerp );
			h_resol_sumET_fit[transv][1][k][m] = new TH1F( Form("h_resol_sumET_" + variableID[transv] + "_" + kanalID[m] + "_allbkg_%d", k), "all bkg", nbinuPerp, minuPerp, maxuPerp );
			h_resol_sumET_fit[transv][2][k][m] = new TH1F( Form("h_resol_sumET_" + variableID[transv] + "_" + kanalID[m] + "_signal_%d", k), "signal ", nbinuPerp, minuPerp, maxuPerp );

		}

	}

	for( int k = 0; k < nbinNVtx; k++ ){

		for( int m = 0; m < nkanal; m++ ){

			h_resol_NVtx_fit [parall][0][k][m] = new TH1F( Form("h_resol_NVtx_"  + variableID[parall] + "_" + kanalID[m] + "_data_%d",   k), "data   ", nbinuPara, minuPara, maxuPara );
	 		h_resol_NVtx_fit [parall][1][k][m] = new TH1F( Form("h_resol_NVtx_"  + variableID[parall] + "_" + kanalID[m] + "_allbkg_%d", k), "all bkg", nbinuPara, minuPara, maxuPara );
 			h_resol_NVtx_fit [parall][2][k][m] = new TH1F( Form("h_resol_NVtx_"  + variableID[parall] + "_" + kanalID[m] + "_signal_%d", k), "signal ", nbinuPara, minuPara, maxuPara );
			h_resol_NVtx_fit [transv][0][k][m] = new TH1F( Form("h_resol_NVtx_"  + variableID[transv] + "_" + kanalID[m] + "_data_%d",   k), "data   ", nbinuPerp, minuPerp, maxuPerp );
			h_resol_NVtx_fit [transv][1][k][m] = new TH1F( Form("h_resol_NVtx_"  + variableID[transv] + "_" + kanalID[m] + "_allbkg_%d", k), "all bkg", nbinuPerp, minuPerp, maxuPerp );
			h_resol_NVtx_fit [transv][2][k][m] = new TH1F( Form("h_resol_NVtx_"  + variableID[transv] + "_" + kanalID[m] + "_signal_%d", k), "signal ", nbinuPerp, minuPerp, maxuPerp );

		}

	}

	for( int i = 0; i < nvariable; i++ ){

	if( i == MET ) continue; 

		for( int j = 2; j < nprocess ; j++ ){

			for( int k = 0; k < nbinpT; k++ ){

				h_resol_pT_fit   [i][0][k][Zee] -> Add( h_resol_pT   [i][DoubleEG2016B][k] ); 
				h_resol_pT_fit   [i][0][k][Zee] -> Add( h_resol_pT   [i][DoubleEG2016C][k] ); 
				h_resol_pT_fit   [i][0][k][Zee] -> Add( h_resol_pT   [i][DoubleEG2016D][k] ); 
				h_resol_pT_fit   [i][2][k][Zee] -> Add( h_resol_pT   [i][DY_ee        ][k] ); 
				h_resol_pT_fit   [i][1][k][Zee] -> Add( h_resol_pT   [i][TT_ee        ][k] ); 
				h_resol_pT_fit   [i][1][k][Zee] -> Add( h_resol_pT   [i][WW_ee        ][k] ); 
				h_resol_pT_fit   [i][1][k][Zee] -> Add( h_resol_pT   [i][WZTo2L2Q_ee  ][k] ); 
				h_resol_pT_fit   [i][1][k][Zee] -> Add( h_resol_pT   [i][WZTo3LNu_ee  ][k] ); 
				h_resol_pT_fit   [i][1][k][Zee] -> Add( h_resol_pT   [i][ZZTo4L_ee    ][k] ); 
				h_resol_pT_fit   [i][1][k][Zee] -> Add( h_resol_pT   [i][ZZTo2L2Q_ee  ][k] ); 
				h_resol_pT_fit   [i][1][k][Zee] -> Add( h_resol_pT   [i][ZZTo2L2Nu_ee ][k] ); 

				h_resol_pT_fit   [i][0][k][Zmumu] -> Add( h_resol_pT   [i][DoubleMuon2016B][k] ); 
				h_resol_pT_fit   [i][0][k][Zmumu] -> Add( h_resol_pT   [i][DoubleMuon2016C][k] ); 
				h_resol_pT_fit   [i][0][k][Zmumu] -> Add( h_resol_pT   [i][DoubleMuon2016D][k] ); 
				h_resol_pT_fit   [i][2][k][Zmumu] -> Add( h_resol_pT   [i][DY_mm          ][k] ); 
				h_resol_pT_fit   [i][1][k][Zmumu] -> Add( h_resol_pT   [i][TT_mm          ][k] ); 
				h_resol_pT_fit   [i][1][k][Zmumu] -> Add( h_resol_pT   [i][WW_mm          ][k] ); 
				h_resol_pT_fit   [i][1][k][Zmumu] -> Add( h_resol_pT   [i][WZTo2L2Q_mm    ][k] ); 
				h_resol_pT_fit   [i][1][k][Zmumu] -> Add( h_resol_pT   [i][WZTo3LNu_mm    ][k] ); 
				h_resol_pT_fit   [i][1][k][Zmumu] -> Add( h_resol_pT   [i][ZZTo4L_mm      ][k] ); 
				h_resol_pT_fit   [i][1][k][Zmumu] -> Add( h_resol_pT   [i][ZZTo2L2Q_mm    ][k] ); 
				h_resol_pT_fit   [i][1][k][Zmumu] -> Add( h_resol_pT   [i][ZZTo2L2Nu_mm   ][k] );

				h_resol_pT_fit   [i][0][k][Gamma] -> Add( h_resol_pT   [i][SinglePhoton2016B][k] ); 
				h_resol_pT_fit   [i][0][k][Gamma] -> Add( h_resol_pT   [i][SinglePhoton2016C][k] ); 
				h_resol_pT_fit   [i][0][k][Gamma] -> Add( h_resol_pT   [i][SinglePhoton2016D][k] ); 
				h_resol_pT_fit   [i][2][k][Gamma] -> Add( h_resol_pT   [i][GJets40100 ][k] ); 
				h_resol_pT_fit   [i][2][k][Gamma] -> Add( h_resol_pT   [i][GJets100200][k] ); 
				h_resol_pT_fit   [i][2][k][Gamma] -> Add( h_resol_pT   [i][GJets200400][k] ); 
				h_resol_pT_fit   [i][2][k][Gamma] -> Add( h_resol_pT   [i][GJets400600][k] ); 
				h_resol_pT_fit   [i][2][k][Gamma] -> Add( h_resol_pT   [i][GJets600Inf][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][QCD200300      ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][QCD300500      ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][QCD500700      ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][QCD7001000     ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][QCD10001500    ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][QCD15002000    ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][QCD2000Inf     ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WJets100200    ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WJets200400    ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WJets400600    ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WJets600800    ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WJets8001200   ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WJets12002500  ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WJets2500Inf   ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][ZGJets         ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][ZNuNuGJets40130][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][ZGTo2LG        ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WGJets         ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][WGToLNuG       ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][TTGJets        ][k] ); 
				h_resol_pT_fit   [i][1][k][Gamma] -> Add( h_resol_pT   [i][TGJets         ][k] ); 

			}

			for( int k = 0; k < nbinsumET; k++ ){

				h_resol_sumET_fit   [i][0][k][Zee] -> Add( h_resol_sumET   [i][DoubleEG2016B][k] ); 
				h_resol_sumET_fit   [i][0][k][Zee] -> Add( h_resol_sumET   [i][DoubleEG2016C][k] ); 
				h_resol_sumET_fit   [i][0][k][Zee] -> Add( h_resol_sumET   [i][DoubleEG2016D][k] ); 
				h_resol_sumET_fit   [i][2][k][Zee] -> Add( h_resol_sumET   [i][DY_ee        ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zee] -> Add( h_resol_sumET   [i][TT_ee        ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zee] -> Add( h_resol_sumET   [i][WW_ee        ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zee] -> Add( h_resol_sumET   [i][WZTo2L2Q_ee  ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zee] -> Add( h_resol_sumET   [i][WZTo3LNu_ee  ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zee] -> Add( h_resol_sumET   [i][ZZTo4L_ee    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zee] -> Add( h_resol_sumET   [i][ZZTo2L2Q_ee  ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zee] -> Add( h_resol_sumET   [i][ZZTo2L2Nu_ee ][k] ); 

				h_resol_sumET_fit   [i][0][k][Zmumu] -> Add( h_resol_sumET   [i][DoubleMuon2016B][k] ); 
				h_resol_sumET_fit   [i][0][k][Zmumu] -> Add( h_resol_sumET   [i][DoubleMuon2016C][k] ); 
				h_resol_sumET_fit   [i][0][k][Zmumu] -> Add( h_resol_sumET   [i][DoubleMuon2016D][k] ); 
				h_resol_sumET_fit   [i][2][k][Zmumu] -> Add( h_resol_sumET   [i][DY_mm          ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zmumu] -> Add( h_resol_sumET   [i][TT_mm          ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zmumu] -> Add( h_resol_sumET   [i][WW_mm          ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zmumu] -> Add( h_resol_sumET   [i][WZTo2L2Q_mm    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zmumu] -> Add( h_resol_sumET   [i][WZTo3LNu_mm    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zmumu] -> Add( h_resol_sumET   [i][ZZTo4L_mm      ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zmumu] -> Add( h_resol_sumET   [i][ZZTo2L2Q_mm    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Zmumu] -> Add( h_resol_sumET   [i][ZZTo2L2Nu_mm   ][k] );

				h_resol_sumET_fit   [i][0][k][Gamma] -> Add( h_resol_sumET   [i][SinglePhoton2016B][k] ); 
				h_resol_sumET_fit   [i][0][k][Gamma] -> Add( h_resol_sumET   [i][SinglePhoton2016C][k] ); 
				h_resol_sumET_fit   [i][0][k][Gamma] -> Add( h_resol_sumET   [i][SinglePhoton2016D][k] ); 
				h_resol_sumET_fit   [i][2][k][Gamma] -> Add( h_resol_sumET   [i][GJets40100 ][k] ); 
				h_resol_sumET_fit   [i][2][k][Gamma] -> Add( h_resol_sumET   [i][GJets100200][k] ); 
				h_resol_sumET_fit   [i][2][k][Gamma] -> Add( h_resol_sumET   [i][GJets200400][k] ); 
				h_resol_sumET_fit   [i][2][k][Gamma] -> Add( h_resol_sumET   [i][GJets400600][k] ); 
				h_resol_sumET_fit   [i][2][k][Gamma] -> Add( h_resol_sumET   [i][GJets600Inf][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][QCD200300      ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][QCD300500      ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][QCD500700      ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][QCD7001000     ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][QCD10001500    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][QCD15002000    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][QCD2000Inf     ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WJets100200    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WJets200400    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WJets400600    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WJets600800    ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WJets8001200   ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WJets12002500  ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WJets2500Inf   ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][ZGJets         ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][ZNuNuGJets40130][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][ZGTo2LG        ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WGJets         ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][WGToLNuG       ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][TTGJets        ][k] ); 
				h_resol_sumET_fit   [i][1][k][Gamma] -> Add( h_resol_sumET   [i][TGJets         ][k] ); 

			}

			for( int k = 0; k < nbinNVtx; k++ ){

				h_resol_NVtx_fit   [i][0][k][Zee] -> Add( h_resol_NVtx   [i][DoubleEG2016B][k] ); 
				h_resol_NVtx_fit   [i][0][k][Zee] -> Add( h_resol_NVtx   [i][DoubleEG2016C][k] ); 
				h_resol_NVtx_fit   [i][0][k][Zee] -> Add( h_resol_NVtx   [i][DoubleEG2016D][k] ); 
				h_resol_NVtx_fit   [i][2][k][Zee] -> Add( h_resol_NVtx   [i][DY_ee        ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zee] -> Add( h_resol_NVtx   [i][TT_ee        ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zee] -> Add( h_resol_NVtx   [i][WW_ee        ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zee] -> Add( h_resol_NVtx   [i][WZTo2L2Q_ee  ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zee] -> Add( h_resol_NVtx   [i][WZTo3LNu_ee  ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zee] -> Add( h_resol_NVtx   [i][ZZTo4L_ee    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zee] -> Add( h_resol_NVtx   [i][ZZTo2L2Q_ee  ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zee] -> Add( h_resol_NVtx   [i][ZZTo2L2Nu_ee ][k] ); 

				h_resol_NVtx_fit   [i][0][k][Zmumu] -> Add( h_resol_NVtx   [i][DoubleMuon2016B][k] ); 
				h_resol_NVtx_fit   [i][0][k][Zmumu] -> Add( h_resol_NVtx   [i][DoubleMuon2016C][k] ); 
				h_resol_NVtx_fit   [i][0][k][Zmumu] -> Add( h_resol_NVtx   [i][DoubleMuon2016D][k] ); 
				h_resol_NVtx_fit   [i][2][k][Zmumu] -> Add( h_resol_NVtx   [i][DY_mm          ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zmumu] -> Add( h_resol_NVtx   [i][TT_mm          ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zmumu] -> Add( h_resol_NVtx   [i][WW_mm          ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zmumu] -> Add( h_resol_NVtx   [i][WZTo2L2Q_mm    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zmumu] -> Add( h_resol_NVtx   [i][WZTo3LNu_mm    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zmumu] -> Add( h_resol_NVtx   [i][ZZTo4L_mm      ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zmumu] -> Add( h_resol_NVtx   [i][ZZTo2L2Q_mm    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Zmumu] -> Add( h_resol_NVtx   [i][ZZTo2L2Nu_mm   ][k] );

				h_resol_NVtx_fit   [i][0][k][Gamma] -> Add( h_resol_NVtx   [i][SinglePhoton2016B][k] ); 
				h_resol_NVtx_fit   [i][0][k][Gamma] -> Add( h_resol_NVtx   [i][SinglePhoton2016C][k] ); 
				h_resol_NVtx_fit   [i][0][k][Gamma] -> Add( h_resol_NVtx   [i][SinglePhoton2016D][k] ); 
				h_resol_NVtx_fit   [i][2][k][Gamma] -> Add( h_resol_NVtx   [i][GJets40100 ][k] ); 
				h_resol_NVtx_fit   [i][2][k][Gamma] -> Add( h_resol_NVtx   [i][GJets100200][k] ); 
				h_resol_NVtx_fit   [i][2][k][Gamma] -> Add( h_resol_NVtx   [i][GJets200400][k] ); 
				h_resol_NVtx_fit   [i][2][k][Gamma] -> Add( h_resol_NVtx   [i][GJets400600][k] ); 
				h_resol_NVtx_fit   [i][2][k][Gamma] -> Add( h_resol_NVtx   [i][GJets600Inf][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][QCD200300      ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][QCD300500      ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][QCD500700      ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][QCD7001000     ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][QCD10001500    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][QCD15002000    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][QCD2000Inf     ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WJets100200    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WJets200400    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WJets400600    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WJets600800    ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WJets8001200   ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WJets12002500  ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WJets2500Inf   ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][ZGJets         ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][ZNuNuGJets40130][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][ZGTo2LG        ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WGJets         ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][WGToLNuG       ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][TTGJets        ][k] ); 
				h_resol_NVtx_fit   [i][1][k][Gamma] -> Add( h_resol_NVtx   [i][TGJets         ][k] ); 

			}

		}

	}

 
	// fits & final plots 

	for( int i = 0; i < nvariable; i++){

		if( i == MET ) continue; 


		// ----- resolution: photon pT

		for( int k = 0; k < nbinpT; k++ ){

			xpT[k] = minpT + (1.0*k+0.5)*(maxpT-minpT)/nbinpT;  

			GetResolution( Zee  , i, pT, k, ypT_A_Zee[k],    epT_A_Zee[k],    chi2pT_A_Zee[k], "GJetsVoigtian" );
			GetResolution( Zmumu, i, pT, k, ypT_A_Zmm[k],    epT_A_Zmm[k],    chi2pT_A_Zmm[k], "GJetsVoigtian" );
			GetResolution( Gamma, i, pT, k, ypT_A_lum[k],    epT_A_lum[k],    chi2pT_A_lum[k], "GJetsVoigtian" );

			//GetResolution( i, pT, k, ypT_B[k],  epT_B[k],  chi2pT_B[k], "GJetsFromMC"    );

		}
		
		PlotResolution( i, pT     );  cout << "final plot for " << parameterID[pT]    << " done !! " << endl; 


		// ----- resolution: sumEt

		//for( int k = 0; k < nbinsumET; k++ ){

		//	xsumET[k] = minsumET + (1.0*k+0.5)*(maxsumET-minsumET)/nbinsumET;  

		//	xsumET[k] = xsumET[k]/1000;  

		//	GetResolution( i, sumET, k, ysumET_A[k], esumET_A[k], chi2sumET_A[k], "GJetsVoigtian" );
			//GetResolution( i, sumET, k, ysumET_B[k], esumET_B[k], chi2sumET_B[k], "GJetsFromMC"   );

		//}

		//PlotResolution( i, sumET  );  cout << "final plot for " << parameterID[sumET] << " done !! " << endl; 


		// ----- resolution: number of vertices

		for( int k = 0; k < nbinNVtx; k++ ){

			xNVtx[k] = minNVtx + (1.0*k+0.5)*(maxNVtx-minNVtx)/nbinNVtx;

			GetResolution( Zee  , i, NVtx, k, yNVtx_A_Zee[k],    eNVtx_A_Zee[k],    chi2NVtx_A_Zee[k], "GJetsVoigtian" );
			GetResolution( Zmumu, i, NVtx, k, yNVtx_A_Zmm[k],    eNVtx_A_Zmm[k],    chi2NVtx_A_Zmm[k], "GJetsVoigtian" );
			GetResolution( Gamma, i, NVtx, k, yNVtx_A_lum[k],    eNVtx_A_lum[k],    chi2NVtx_A_lum[k], "GJetsVoigtian" );

			//GetResolution( i, NVtx, k, yNVtx_B[k],  eNVtx_B[k],  chi2NVtx_B[k], "GJetsFromMC"    );

		}

		PlotResolution( i, NVtx   );  cout << "final plot for "  << parameterID[NVtx]  << " done !! " << endl;


	}

	cout << "                 " << endl; 
	cout << "exiting OrderFits" << endl; 
	cout << "                 " << endl; 
}




void PlotResolution( int ivar, int parameter ){

	TCanvas* c = new TCanvas("resol_" + variableID[ivar] + "_" + parameterID[parameter], variableID[ivar] + "_" + parameterID[parameter], 600, 600);

	TMultiGraph*  TheMultiGraph = new TMultiGraph();   // if >1 TGraphErrors plotted on the same canvas

	TGraphErrors* TheGraph[3]; 

	if( parameter == pT    ){
 
		TheGraph[Zee  ] = new TGraphErrors(nbinpT   , xpT   , ypT_A_Zee   , 0, epT_A_Zee    );
		TheGraph[Zmumu] = new TGraphErrors(nbinpT   , xpT   , ypT_A_Zmm   , 0, epT_A_Zmm    );
		TheGraph[Gamma] = new TGraphErrors(nbinpT   , xpT   , ypT_A_lum   , 0, epT_A_lum    );

	}

	if( parameter == sumET ){

		 TheGraph[Zee  ] = new TGraphErrors(nbinsumET, xsumET, ysumET_A_Zee, 0, esumET_A_Zee );
		 TheGraph[Zmumu] = new TGraphErrors(nbinsumET, xsumET, ysumET_A_Zmm, 0, esumET_A_Zmm );
		 TheGraph[Gamma] = new TGraphErrors(nbinsumET, xsumET, ysumET_A_lum, 0, esumET_A_lum );

	}


	if( parameter == NVtx  ){

		TheGraph[Zee  ] = new TGraphErrors(nbinNVtx , xNVtx , yNVtx_A_Zee , 0, eNVtx_A_Zee  );
		TheGraph[Zmumu] = new TGraphErrors(nbinNVtx , xNVtx , yNVtx_A_Zmm , 0, eNVtx_A_Zmm  );
		TheGraph[Gamma] = new TGraphErrors(nbinNVtx , xNVtx , yNVtx_A_lum , 0, eNVtx_A_lum  );

	}
 

	TheGraph[0] -> SetTitle("Z #rightarrow ee");
	TheGraph[1] -> SetTitle("Z #rightarrow #mu#mu");
	TheGraph[2] -> SetTitle("#gamma");

	TheGraph[0] -> SetMarkerStyle( 22 );   TheGraph[0] -> SetMarkerSize( 1.2*TheGraph[0]->GetMarkerSize() );
	TheGraph[1] -> SetMarkerStyle( 23 );   TheGraph[1] -> SetMarkerSize( 1.2*TheGraph[0]->GetMarkerSize() );
	TheGraph[2] -> SetMarkerStyle( 21 );

	TheGraph[0] -> SetMarkerColor( kBlue );
	TheGraph[1] -> SetMarkerColor( kRed  );
	TheGraph[2] -> SetMarkerColor( kGreen);
	
	TheGraph[0] -> SetLineColor( kBlue  );
	TheGraph[1] -> SetLineColor( kRed   );
	TheGraph[2] -> SetLineColor( kGreen );

	TheGraph[0] -> SetFillStyle(0);
	TheGraph[1] -> SetFillStyle(0);
	TheGraph[2] -> SetFillStyle(0);

	TheMultiGraph -> Add( TheGraph[0] );
	TheMultiGraph -> Add( TheGraph[1] );
	TheMultiGraph -> Add( TheGraph[2] );

	TheMultiGraph -> Draw("AP");

	TheMultiGraph -> GetXaxis() -> SetTitle( parameterIDfancy[parameter] );
	TheMultiGraph -> GetYaxis() -> SetTitle( sigma_variable[ivar] );
	TheMultiGraph -> GetYaxis() -> SetRangeUser( 0.0, 50.0 );

	c-> BuildLegend(); 

	//TheGraph[0] -> Draw("AP"      );
	//TheGraph[1] -> Draw("AP, same");
	//TheGraph[2] -> Draw("AP, same");

	//TLatex tex;
	//tex.SetTextAlign(13);
	//tex.SetTextSize(0.03);
	//tex.SetNDC();
	//tex.DrawLatex ( 0.1, 0.95, "CMS Preliminary                     2.318 fb^{-1} (13TeV)" );

	DrawLatex( 61, 0.100, 0.945, 0.050, 11, "CMS"                                             );
 	DrawLatex( 52, 0.205, 0.945, 0.030, 11, "Preliminary"                                     );
	DrawLatex( 42, 0.900, 0.945, 0.050, 31, Form("%.3f fb^{-1} (13TeV, 2016)", TheLuminosity) );


	c -> SaveAs( "resol/" + variableID[ivar] + "_" + parameterID[parameter] + ".pdf" ); 
	c -> SaveAs( "resol/" + variableID[ivar] + "_" + parameterID[parameter] + ".png" ); 

}


void GetResolution( int ch, int whichvar, int parameter, int ibin, float& resol, float& eresol, float& chichi, TString GJetsOrigin ){

	int ivar = whichvar; 

	RooRealVar x           ( "variable"   , "variable"          ,  0  , -500, 500 );   // random variable; redefine its range
	RooRealVar Gauss_mean  ( "Gauss_mean" , "Gauss mean"        ,  0  ,  -10, 10  );  
	RooRealVar BW_gamma    ( "BW_gamma"   , "Breit-Wigner gamma",  2.3,    0, 100 );  
	RooRealVar Gauss_sigma ( "Gauss_sigma", "Gauss sigma"       , 10  ,    0, 100 );  


	TH1F* h1;
	TH1F* h2; 
	TH1F* h3;

	int min; int width; 

	if( parameter == pT    ){

		h1 = (TH1F*)  h_resol_pT_fit[ivar][0][ibin][ch] -> Clone();
		h2 = (TH1F*)  h_resol_pT_fit[ivar][1][ibin][ch] -> Clone();
		h3 = (TH1F*)  h_resol_pT_fit[ivar][2][ibin][ch] -> Clone();

		min = minpT; 
		width = (maxpT - minpT)/nbinpT; 

	}

	if( parameter == sumET ){

		h1 = (TH1F*)  h_resol_sumET_fit[ivar][0][ibin][ch] -> Clone();
		h2 = (TH1F*)  h_resol_sumET_fit[ivar][1][ibin][ch] -> Clone();
		h3 = (TH1F*)  h_resol_sumET_fit[ivar][2][ibin][ch] -> Clone();

		min = minsumET; 
		width = (maxsumET - minsumET)/nbinsumET; 

	}

	if( parameter == NVtx  ){

		h1 = (TH1F*)  h_resol_NVtx_fit[ivar][0][ibin][ch] -> Clone();
		h2 = (TH1F*)  h_resol_NVtx_fit[ivar][1][ibin][ch] -> Clone();
		h3 = (TH1F*)  h_resol_NVtx_fit[ivar][2][ibin][ch] -> Clone();

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

	if( GJetsOrigin == "GJetsVoigtian" ) result = modelA.fitTo( data, Extended(kFALSE), RooFit::Save(kTRUE) );   
	if( GJetsOrigin == "GJetsFromMC"   ) result = modelB.fitTo( data, Extended(kFALSE), RooFit::Save(kTRUE) );   

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

	frame -> GetXaxis() -> SetTitle(variableIDfancy[ivar]);

	data.plotOn( frame );

	modelA.plotOn( frame, Components(bkgpdf), LineColor(kRed)  , LineStyle(kDashed), FillColor(kRed)  , DrawOption("F") );

	modelA.plotOn( frame, Components(voigt) , LineColor(kGreen), LineStyle(kDashed), FillColor(kGreen), DrawOption("L") );

	modelA.plotOn( frame,                     LineColor(kBlue) , LineStyle(kDashed), FillColor(kBlue) , DrawOption("L") );

	frame -> Draw();

	float chi2 = frame -> chiSquare();
	chichi = chi2; 

	/*TLatex mylatex;
	mylatex.SetTextAlign(13);
	mylatex.SetTextSize(0.03);
	mylatex.SetNDC();
	mylatex.DrawLatex ( 0.1, 0.95, Form("CMS Preliminary  2.318 fb^{-1} (13TeV)   #chi^{2} = %5.2f", chi2) );*/

	DrawLatex( 61, 0.100, 0.945, 0.050, 11, "CMS"                                             );
 	DrawLatex( 52, 0.205, 0.945, 0.030, 11, "Preliminary"                                     );
	DrawLatex( 42, 0.900, 0.945, 0.050, 31, Form("%.3f fb^{-1} (13TeV, 2016)", TheLuminosity) );
	DrawLatex( 42, 0.840, 0.900, 0.030, 31, Form("#chi^{2} = %5.2f ", chi2)                   );

	mycanvas -> SaveAs( Form( "fit/" + kanalID[ch] + "_" + variableID[ivar] + "_" + parameterID[parameter] +  "_%dto%d.pdf", low, up ) );
	mycanvas -> SaveAs( Form( "fit/" + kanalID[ch] + "_" + variableID[ivar] + "_" + parameterID[parameter] +  "_%dto%d.png", low, up ) );
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





