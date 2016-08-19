#include "MET.h"

void MET2(){

	if( DoGlobalPlots == false ) { cout << "   >>   GlobalPlots deactivated !!!   " << endl;  return; }


	Assign();



	TFile* allhistos = new TFile( "histograms/" + allhistosReadFrom + ".root", "read" );

	for( int i = 0; i < nvariable_G; i++ ){

		for( int j = 0; j < nprocess; j++ ){
		//for( int j = alpha; j < omega+1; j++ ){

			int TheKanal = kanal[j]; 

			h_global[i][j] = (TH1F*) allhistos -> Get("h_global_" + variableID_G[i] + "_" + kanalID[TheKanal] + "_" + processID[j]);

		}

	}


	for( int j = 0; j < nprocess; j++ ){
	//for( int j = alpha; j < omega+1; j++ ){

		if( isData[j] ) continue; 

		for( int i = 0; i < nvariable_G; i++ ){ 	

 			h_global[i][j] -> Scale( TheLuminosity ); 

			h_global[i][j] -> SetFillColorAlpha( ProcessColor[j], 0.7 );
		
		}		

	}

	for( int i = 0; i < nvariable_G; i++ ){

			h_global[i][WW_ee]->Add(h_global[i][WZTo2L2Q_ee ]);
			h_global[i][WW_ee]->Add(h_global[i][WZTo3LNu_ee ]);
			h_global[i][WW_ee]->Add(h_global[i][ZZTo4L_ee   ]);
			h_global[i][WW_ee]->Add(h_global[i][ZZTo2L2Q_ee ]);
			h_global[i][WW_ee]->Add(h_global[i][ZZTo2L2Nu_ee]);

	
			h_global[i][WW_mm]->Add(h_global[i][WZTo2L2Q_mm ]);
			h_global[i][WW_mm]->Add(h_global[i][WZTo3LNu_mm ]);
			h_global[i][WW_mm]->Add(h_global[i][ZZTo4L_mm   ]);
			h_global[i][WW_mm]->Add(h_global[i][ZZTo2L2Q_mm ]);
			h_global[i][WW_mm]->Add(h_global[i][ZZTo2L2Nu_mm]);


			h_global[i][GJets40100]->Add(h_global[i][GJets100200]);
			h_global[i][GJets40100]->Add(h_global[i][GJets200400]);
			h_global[i][GJets40100]->Add(h_global[i][GJets400600]);
			h_global[i][GJets40100]->Add(h_global[i][GJets600Inf]);  

			h_global[i][QCD200300]->Add(h_global[i][QCD300500  ]);
			h_global[i][QCD200300]->Add(h_global[i][QCD500700  ]);
			h_global[i][QCD200300]->Add(h_global[i][QCD7001000 ]);
			h_global[i][QCD200300]->Add(h_global[i][QCD10001500]);
			h_global[i][QCD200300]->Add(h_global[i][QCD15002000]);
			h_global[i][QCD200300]->Add(h_global[i][QCD2000Inf ]);

			h_global[i][WJets100200]->Add(h_global[i][WJets200400  ]);
			h_global[i][WJets100200]->Add(h_global[i][WJets400600  ]);
			h_global[i][WJets100200]->Add(h_global[i][WJets600800  ]);
			h_global[i][WJets100200]->Add(h_global[i][WJets8001200 ]);
			h_global[i][WJets100200]->Add(h_global[i][WJets12002500]);
			h_global[i][WJets100200]->Add(h_global[i][WJets2500Inf ]);

			h_global[i][ZGJets ] -> Add( h_global[i][ZNuNuGJets40130] );
			h_global[i][ZGJets ] -> Add( h_global[i][ZGTo2LG        ] );

			h_global[i][WGJets ] -> Add( h_global[i][WGToLNuG       ] );

			h_global[i][TTGJets] -> Add( h_global[i][TGJets         ] );  

	}


	GlobalPlots( Zee   ); 
	GlobalPlots( Zmumu ); 
	GlobalPlots( Gamma ); 

}


void GlobalPlots( int ch ){  

	int runB, runC, runD; 

	for( int i = 0; i < nvariable_G; i++ ){

		if ( ch == Zee   ) { runB =DoubleEG2016B    ; runC =DoubleEG2016C    ; runD =DoubleEG2016D    ; }
		if ( ch == Zmumu ) { runB =DoubleMuon2016B  ; runC =DoubleMuon2016C  ; runD =DoubleMuon2016D  ; }
		if ( ch == Gamma ) { runB =SinglePhoton2016B; runC =SinglePhoton2016C; runD =SinglePhoton2016D; }

		h_data[i] = (TH1F*) h_global[i][runB] -> Clone( "h_data_" + variableID_G[i] );

		h_data[i] -> Add( h_global[i][runC] );
		h_data[i] -> Add( h_global[i][runD] );

		h_data[i] -> SetMarkerStyle(20);

		h_data[i] -> SetMarkerColor(kBlack);

		h_data[i] -> SetLineColor(kBlack);

	}


	//-- h_mc 
 
	for( int i = 0; i < nvariable_G; i++ ){

		if( ch == Zee ){

			h_mc[i] = (TH1F*) h_global[i][DY_ee] -> Clone( "h_mc_" + variableID_G[i] );

			h_mc[i] -> Add( h_global[i][TT_ee] );
			h_mc[i] -> Add( h_global[i][WW_ee] );
		
		}

		if( ch == Zmumu ){

			h_mc[i] = (TH1F*) h_global[i][DY_mm] -> Clone( "h_mc_" + variableID_G[i] );

			h_mc[i] -> Add( h_global[i][TT_mm] );
			h_mc[i] -> Add( h_global[i][WW_mm] );
	
		}


		if( ch == Gamma ){

			h_mc[i] = (TH1F*) h_global[i][GJets40100] -> Clone( "h_mc_" + variableID_G[i] );

			h_mc[i] -> Add( h_global[i][QCD200300  ] );
			h_mc[i] -> Add( h_global[i][WJets100200] );
			h_mc[i] -> Add( h_global[i][WGJets     ] );
			h_mc[i] -> Add( h_global[i][ZGJets     ] );
			h_mc[i] -> Add( h_global[i][TTGJets    ] );			

		}

	}	



	//-- PU-reweighting

	//GetPUreweighting( ch );



	//-- stack
 
	for( int i = 0; i < nvariable_G; i++ ){

		s_global[i]  = new THStack( variableID_G[i], variableID_G[i] );

		if(  ch == Zee ){

			s_global[i] -> Add( h_global[i][TT_ee] );
			s_global[i] -> Add( h_global[i][WW_ee] );
			s_global[i] -> Add( h_global[i][DY_ee] );

		}

		if(  ch == Zmumu ){

			s_global[i] -> Add( h_global[i][TT_mm] );
			s_global[i] -> Add( h_global[i][WW_mm] );
			s_global[i] -> Add( h_global[i][DY_mm] );

		}

		if( ch == Gamma ){	

			s_global[i] -> Add( h_global[i][TTGJets        ] );
			s_global[i] -> Add( h_global[i][ZGJets         ] );
			s_global[i] -> Add( h_global[i][WJets100200    ] );
			s_global[i] -> Add( h_global[i][WGJets         ] );
			s_global[i] -> Add( h_global[i][QCD200300      ] );
			s_global[i] -> Add( h_global[i][GJets40100     ] );

		}

	}


	//-- ratio

	for( int i = 0; i < nvariable_G; i++ ){

		Ratio[i] = (TH1F*) h_data[i] -> Clone( "ratio_" + variableID_G[i] );
		Unc  [i] = (TH1F*) h_mc  [i] -> Clone( "ratio_" + variableID_G[i] ); 

		Unc[i]->SetLineColor(11);
		Unc[i]->SetFillColor(11);

	}

	BuildRatio( parall_G, nbinuPara ); 
	BuildRatio( transv_G, nbinuPerp ); 
	BuildRatio( MET     , nbinMET   ); 
	BuildRatio( VpT     , nbinpT    ); 
	BuildRatio( nvert   , nbinnvert ); 

	

	for( int i = 0; i < nvariable_G; i++ ){

		TCanvas* c = new TCanvas( "canvas_" + variableID_G[i], variableID_G[i], 550, 720 );  

		TPad* pad1 = new TPad("pad1", "pad1", 0.05, 0.0, 1.00, 0.3);
		TPad* pad2 = new TPad("pad2", "pad2", 0.05, 0.3, 1.00, 1.0);

		pad1 -> SetTopMargin   (0.08);
		pad1 -> SetBottomMargin(0.35);
		pad1 -> Draw();

		pad2 -> SetTopMargin   (0.08);
		pad2 -> SetBottomMargin(0.02);
		pad2 -> Draw();

		// ---------------------------------------------------------------

		pad2 -> cd();

		pad2 -> SetLogy(); 

		h_data[i] -> SetTitle("");

		h_data[i] -> SetStats(false);   // it has priority over the gStyle->SetOptStats option

 		SetAxis( h_data[i], "", Form("Events / %1.0f GeV", (maxuPara-minuPara)/nbinuPara), 1.5, 1.8 );

		h_data[i] -> SetMinimum(0.1);


		h_data[i] -> Draw("e"); 

		s_global[i] -> Draw("hist same");

		h_data[i] -> Draw("e same");

		// ---------------------------------------------------------------
		
		TLegend* TheLegend = new TLegend( 0.68, 0.65, 0.88, 0.88 );

		TheLegend -> SetBorderSize(0);

		TheLegend -> SetTextSize(0.030);

		for( int j = 0; j < nprocess; j++ ){ 


			if(  ch == Zee && ( j == DoubleEG2016B || j ==  DY_ee || j == TT_ee || j == WW_ee ) ) {

			( isData[j] )   ?   TheLegend -> AddEntry( h_data[i], processIDfancy[j], "p" )   :   TheLegend -> AddEntry( h_global[i][j], processIDfancy[j], "f" );

			}

			if(  ch == Zmumu && ( j == DoubleMuon2016B || j ==  DY_mm || j == TT_mm || j == WW_mm ) ) {

			( isData[j] )   ?   TheLegend -> AddEntry( h_data[i], processIDfancy[j], "p" )   :   TheLegend -> AddEntry( h_global[i][j], processIDfancy[j], "f" );

			}

			if(  ch == Gamma && ( j == SinglePhoton2016B || j ==  GJets40100 || j == QCD200300 || j == WJets100200 || j == WJets100200 || j == ZGJets || j == WGJets || j == TTGJets ) ) {

			( isData[j] )   ?   TheLegend -> AddEntry( h_data[i], processIDfancy[j], "p" )   :   TheLegend -> AddEntry( h_global[i][j], processIDfancy[j], "f" );

			}

			else continue; 

		}


		TheLegend -> Draw();

		// ---------------------------------------------------------------

		DrawLatex( 61, 0.100, 0.945, 0.050, 11, "CMS"                                             );

     	 	DrawLatex( 52, 0.205, 0.945, 0.030, 11, "Preliminary"                                     );

		DrawLatex( 42, 0.900, 0.945, 0.050, 31, Form("%.3f fb^{-1} (13TeV, 2016)", TheLuminosity) );

		// ---------------------------------------------------------------
		
  		pad2 -> RedrawAxis();

		// ---------------------------------------------------------------

		pad1 -> cd();

		//pad1 -> SetLogy();		

		Ratio[i] -> SetTitle("");

		Ratio[i] -> SetStats(false);   // it has priority over the gStyle->SetOptStats option

 		SetAxis( Ratio[i], variableIDfancy_G[i], "data / MC ", 1.4, 0.75 );

      		//Ratio[i] -> GetYaxis() -> SetRangeUser(1e-1, 1e1);
      		Ratio[i] -> GetYaxis() -> SetRangeUser(0., 2.);

		Ratio[i] -> Draw("ep");
		Unc[i]   -> Draw("e2, same");
		Ratio[i] -> Draw("ep, same");

		pad1 -> RedrawAxis();

      		// ---------------------------------------------------------------

		c -> SaveAs( "global/" + kanalID[ch] + "_" + variableID_G[i] + ".pdf" );
		c -> SaveAs( "global/" + kanalID[ch] + "_" + variableID_G[i] + ".png" );

	}
 
}



