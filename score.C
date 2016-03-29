
#include "MVA.h"


void score(){


 	gInterpreter->ExecuteMacro("PaperStyle.C");

	assign(); 

	TCanvas* c = new TCanvas( "MVA", "MVA", 900, 900 );
	c -> Divide( 3, 3 );


	for ( int i = 0; i < nDM; i++) {
	//for ( int i = 2; i < 3; i++) {



		MVAdata[i] = getMVAoutput( DM[i], "data"      );		
		MVAsig [i] = getMVAoutput( DM[i], DM[i]       );   MVAsig [i] -> Scale(L);
		MVAbkg [i] = getMVAoutput( DM[i], "bkg"       );   MVAbkg [i] -> Scale(L); 
		MVAtt  [i] = getMVAoutput( DM[i], "TTJets"    );   MVAtt  [i] -> Scale(L);   float total_tt =  MVAtt [i] -> Integral();  //cout << total_tt << endl; 
		MVAWW  [i] = getMVAoutput( DM[i], "WWTo2L2Nu" );   MVAWW  [i] -> Scale(L);   float total_WW =  MVAWW [i] -> Integral();  //cout << total_WW << endl; 


		cout << "  sample: " << DM[i]        << endl; 
		cout << "------------------------"   << endl;

		int k2 = GetExtremumScore( MVAsig[i], MVAbkg[i] );             	
		frontier[i][nSB] = inf + ( sup - inf ) * k2 / nbin;

		//cout << frontier[i][cut] << endl;

		y_sig[i] = MVAsig[i] -> Integral(k2, nbin+1);
		y_bkg[i] = MVAbkg[i] -> Integral(k2, nbin+1);
 
		//cout << "------------------------"   << endl;
		//cout << "  sig yield = " << y_sig[i] << endl; 
		//cout << "  bkg yield = " << y_bkg[i] << endl; 
		//cout << "------------------------"   << endl;
		//cout << "                        "   << endl;
		//cout << "                        "   << endl;



		c2[i] = new TCanvas( "CRs_" + DM[i], "control regions", 1200, 600 );
		c2[i] -> Divide( nSB, 1 );

		c3[i] = new TCanvas( "CRs_" + DM[i] + " (splitted)", "control regions", nSB*400, 2*400 );
		c3[i] -> Divide( nSB, 2 );

		TPad*    padCR[nSB]   ;
		TLegend* legCR[nSB]   ;
		TPad*   padCR2[nSB][2];

		for ( int j = 0; j < nSB; j++ ) {   // sideband

			frontier[i][j]  = inf + ( sup - inf ) * ( k2 - ( nSB - j ) * bandwidth ) / nbin;

		}

		for ( int j = 0; j < nSB; j++ ) {   // sideband


			njnb[i][j][dt] = GetCRhistogram( i, j, dt );
			njnb[i][j][WW] = GetCRhistogram( i, j, WW );
			njnb[i][j][tt] = GetCRhistogram( i, j, tt );

			njnb[i][j][tt] -> SetMarkerColor(kRed) ;
			njnb[i][j][WW] -> SetMarkerColor(1)    ;

			njnb[i][j][tt] -> SetLineColor(kRed)   ;
			njnb[i][j][WW] -> SetLineColor(1)      ;

			njnb[i][j][tt] -> SetLineWidth(3)      ;
			njnb[i][j][WW] -> SetLineWidth(3)      ;
	
			njnb[i][j][tt] -> SetMarkerStyle(20)   ;
			njnb[i][j][WW] -> SetMarkerStyle(33)   ;

			njnb[i][j][tt] -> SetTitle( DM[i] + " - " + SBid[j] );
			njnb[i][j][tt] -> SetXTitle("number of jets"  );
			njnb[i][j][tt] -> SetYTitle("number of b-jets");
				
			padCR[j] = (TPad*) c2[i] -> GetPad(j+1);
		
			padCR[j] -> cd();

			njnb[i][j][tt] -> Draw("box"      );
			njnb[i][j][WW] -> Draw("box, same");

			legCR[j] = new TLegend( 0.20,0.70,0.50,0.88 );
			legCR[j] -> AddEntry( njnb[i][j][tt], " tt", "f" );	
			legCR[j] -> AddEntry( njnb[i][j][WW], " WW", "f" );
			legCR[j] -> Draw();



			for ( int v = 0; v < 2; v++ ) {   // variable

				for ( int s = 0; s < 2; s++ ) {   // sample

					CRhisto[v][i][j][s] = GetCRsubhistogram( i, j, v, s);  

				}

				CRhisto[v][i][j][bkg]  -> SetFillColor( kGray );

				CRhisto[v][i][j][bkg]  -> SetFillStyle( 1001 );

				CRhisto[v][i][j][data] -> SetMarkerStyle( kFullCircle );

				CRhisto[v][i][j][data] -> SetMarkerSize( 1.0 );

				CRhisto[v][i][j][bkg]  -> SetTitle( varID[v] + " on " + SBid[j] );

				
				padCR2[j][v] = (TPad*) c3[i] -> GetPad( j + nSB*v + 1 ); 

				padCR2[j][v] -> cd();
		
				padCR2[j][v] -> SetLogy(); 

				float ymax =  CRhisto[v][i][j][bkg] -> GetMaximum();
		
				CRhisto[v][i][j][bkg]  -> Draw("hist"       );
				CRhisto[v][i][j][data] -> Draw("ep,   same");

				CRhisto[v][i][j][bkg] -> SetMaximum(2e1 * ymax);
				CRhisto[v][i][j][bkg] -> SetMinimum(1e-4);

			}   // v


			/////////////////////////////////////////////////////////////

			//int cfg = 2; 

			for (int l = 0; l < nCf; l++ ) {   // box configuration
				
				for ( int k = 0; k < nBox; k++ ) {   // box color

					population[i][j][dt][k][l] = GetCRpopulation( i, j, dt, k, l); //cout << population[i][j][dt][k] << endl;
					population[i][j][tt][k][l] = GetCRpopulation( i, j, tt, k, l); //cout << population[i][j][tt][k] << endl;
					population[i][j][WW][k][l] = GetCRpopulation( i, j, WW, k, l); //cout << population[i][j][WW][k] << endl;
					population[i][j][bk][k][l] = GetCRpopulation( i, j, WW, k, l); //cout << population[i][j][WW][k] << endl;

				}   // end box color (k)

				float n1 = population[i][j][dt][0][l] - population[i][j][bk][0][l];       float n1_e2 = population[i][j][dt][0][l] + population[i][j][bk][0][l];
				float n2 = population[i][j][dt][1][l] - population[i][j][bk][1][l];       float n2_e2 = population[i][j][dt][1][l] + population[i][j][bk][1][l];

				float a11 = population[i][j][WW][0][l];  float a12 = population[i][j][tt][0][l];
				float a21 = population[i][j][WW][1][l];  float a22 = population[i][j][tt][1][l];

				float det = a11 * a22  -  a12 * a21; 

				// inverse 
				float b11 =   a22/det;    float b12 = - a12/det;
				float b21 = - a21/det;    float b22 =   a11/det;

				SF[i][j][WW][l] = b11 * n1  +  b12 * n2;   
				SF[i][j][tt][l] = b21 * n1  +  b22 * n2;


				SF_e  [i][j][WW][l] = sqrt(  b11 * b11 * n1_e2   +   b12 * b12 * n2_e2  );      
				SF_e  [i][j][tt][l] = sqrt(  b21 * b21 * n1_e2   +   b22 * b22 * n2_e2  );
			 	SF_cov[i][j][l]     =     (  b11 * b21 * n1_e2   +   b12 * b22 * n2_e2  ) / ( SF_e  [i][j][WW][l] * SF_e  [i][j][tt][l] );



				cout << "                                    "                                                                                                           << endl;
				cout <<  SBid[j] << "  configuration " << l                                                                                                              << endl;
				cout << "------------------------------------"                                                                                                           << endl;
cout<<"  WW: "<<setprecision(3)<<SF[i][j][WW][l]<<" +/- "<<setprecision(3)<< SF_e[i][j][WW][l]<<"    tt: "<<setprecision(3)<<SF[i][j][tt][l]<<" +/- "<<setprecision(3)<<SF_e[i][j][tt][l]<< endl;
cout<<"  WW: "<<SF[i][j][WW][l]*total_WW<<" +/- "<< SF_e[i][j][WW][l]*total_WW <<"  tt: "<<SF[i][j][tt][l]*total_tt<<" +/- "<< SF_e[i][j][tt][l]*total_tt << endl;
				cout << "                                    "                                                                                                           << endl;
				cout << "  with correlation = " <<  SF_cov[i][j][l]                                                                                                      << endl;
				cout << "                                    "                                                                                                           << endl;
				cout << "                                    "                                                                                                           << endl;
                                                                                                                                    
			}   // end box configuration (l)
			
		}   // end sideband (j)

		//for (int l = 0; l < nCf; l++ ) {

			
		//cout << setprecision(3) <<SF[i][4][tt][l] << " +/- " << setprecision(3) << SF_e[i][4][tt][l] << "   " 
                //     << setprecision(3) <<SF[i][3][tt][l] << " +/- " << setprecision(3) << SF_e[i][3][tt][l] << "   " 
                //     << setprecision(3) <<SF[i][2][tt][l] << " +/- " << setprecision(3) << SF_e[i][2][tt][l] << "   " 
                //     << setprecision(3) <<SF[i][1][tt][l] << " +/- " << setprecision(3) << SF_e[i][1][tt][l] << "   " 
                //     << setprecision(3) <<SF[i][0][tt][l] << " +/- " << setprecision(3) << SF_e[i][0][tt][l] << endl;

		//cout << setprecision(3) <<SF[i][4][WW][l] << " +/- " << setprecision(3) << SF_e[i][4][WW][l] << "   " 
                //     << setprecision(3) <<SF[i][3][WW][l] << " +/- " << setprecision(3) << SF_e[i][3][WW][l] << "   " 
                //     << setprecision(3) <<SF[i][2][WW][l] << " +/- " << setprecision(3) << SF_e[i][2][WW][l] << "   " 
                //     << setprecision(3) <<SF[i][1][WW][l] << " +/- " << setprecision(3) << SF_e[i][1][WW][l] << "   " 
                //     << setprecision(3) <<SF[i][0][WW][l] << " +/- " << setprecision(3) << SF_e[i][0][WW][l] << endl;
		
		//}


		c2[i] -> SaveAs( "plots/CRplots/" + DM[i] + ".png" );
		c3[i] -> SaveAs( "plots/CRplots/splitted_" + DM[i] + ".png");

		//////////////////////////////////////////////////////////////////////

	}   // i

	//c ->GetFrame()->DrawClone();
	//c -> SaveAs("MVA.png");
}







TH1F* getMVAoutput(TString DM, TString sample) {

	TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/application_outputs/OPT/h_Cuts_" + DM + "_" + sample + ".root", "READ" );

	TH1F* MVAoutput = (TH1F*) file -> Get("hCuts");

	return MVAoutput;
}



int GetExtremumScore(TH1F* MVAsig, TH1F* MVAbkg){

	float S0 = MVAsig -> Integral(0, nbin+1);   
	float B0 = MVAbkg -> Integral(0, nbin+1);  

	float S[nbin+2];
	float B[nbin+2];

	float eff[nbin+2];
	float rej[nbin+2];

	float ratio = 0;
	int k2 = 0;

	for ( int k = 0; k < nbin + 2; k++ ) {

		if ( k != nbin+1 ) {
	
			S[k] = MVAsig -> Integral(k, nbin+1);

			B[k] = MVAbkg -> Integral(k, nbin+1);

		}
						
		else {

			S[k] = 0; B[k] = 0; 

		}
		
		eff[k] = S[k]/S0;
		 
		rej[k] = 1 - B[k]/B0;

		//cout << k << " -- " << S[k] << " -- " << B[k] << " -- " << S[k]/sqrt(S[k]+B[k]) << endl; 

		if ( S[k] != 0  &&  B[k] != 0 ) {
		
						
			//if ( eff[k] * S[k]/(S[k]+B[k])  >  ratio ) {

			//	ratio = eff[k] * S[k]/(S[k]+B[k]);

			//	k2 = k;
	
			//}


			//if ( S[k]/sqrt(S[k]+B[k] ) > ratio ) {

			//	ratio = S[k]/sqrt(S[k]+B[k]);

			//	k2 = k;
	
			//}


			if ( S[k]/sqrt(B[k]) > ratio ) {

				ratio = S[k]/sqrt(B[k]);

				k2 = k;

			}


			//if ( S[k]/B[k] > ratio ) {

			//	ratio = S[k]/B[k];

			//	k2 = k;

			//}

			
			//if ( sqrt(S[k]+B[k]) - sqrt(B[k]) > ratio ) {

			//	ratio = sqrt(S[k]+B[k]) - sqrt(B[k]);

			//	k2 = k; 

			//}

		}


	}

	//SB = new TGraph(n, cut, score);	

	float maxR = ratio;
	
	float cut = inf + ( sup - inf ) * k2 / nbin;

	//cout << "  the extremum on the score is " << maxR << " for the cut at " << cut << endl;
	//cout << "  signal efficiency = " << eff[k2] << "// bkg rejection = " << rej[k2]   << endl; 

	return k2;
}




TH2F* GetCRhistogram(int i, int SB, int DD) {

	TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/minitrees/noNJetCut/bonsai_" + estimated_bkg[DD] + ".root", "READ" );
	//TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/minitrees/bonsai_" + estimated_bkg[DD] + ".root", "READ" );

	TTree *tree = (TTree*) file -> Get("latino");

	float MVA, njet, nbjet15, baseW;
 
	tree -> SetBranchAddress( "MVA" + DM[i], &MVA     );
	tree -> SetBranchAddress( "njet"       , &njet    );
	tree -> SetBranchAddress( "nbjet15"    , &nbjet15 );
	tree -> SetBranchAddress( "baseW"      , &baseW   );

	int nentries = tree -> GetEntries();

	TH2F* njnb = new TH2F( "njnb", "# b-jets vs # jets", 10, 0, 10, 10, 0, 10 );

	//cout << SB << " - " << SBid[SB] << " - " << frontier[i][SB] << " - " << frontier[i][SB+1] << endl;

	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		tree -> GetEntry(ievt);

		//if ( DD != dt ) baseW = baseW * L; 

		if( MVA > frontier[i][SB] && MVA < frontier[i][SB+1] )  njnb -> Fill( njet, nbjet15, baseW ); 

	}

	return njnb;

}

TH1F* GetCRsubhistogram(int i, int SB, int variable, int set) {

	TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/minitrees/noNJetCut/bonsai_" + id[set] + ".root", "READ" );
	//TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/minitrees/bonsai_" + id[set] + ".root", "READ" );

	TTree *tree = (TTree*) file -> Get("latino");

	float MVA, njet, nbjet15, baseW;
 
	tree -> SetBranchAddress( "MVA" + DM[i], &MVA     );
	tree -> SetBranchAddress( "njet"       , &njet    );
	tree -> SetBranchAddress( "nbjet15"    , &nbjet15 );
	tree -> SetBranchAddress( "baseW"      , &baseW   );

	int nentries = tree -> GetEntries();

	TH1F* CRhisto = new TH1F( "control histogram", "", 10, 0, 10 );


	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		tree -> GetEntry(ievt);

		//if ( DD != dt ) baseW = baseW * L; 

		if( MVA > frontier[i][SB] && MVA < frontier[i][SB+1] ) { 

			if ( variable == 0 ) CRhisto -> Fill( njet   , baseW ); 
			if ( variable == 1 ) CRhisto -> Fill( nbjet15, baseW );
		}

	}

	return CRhisto;

}

float GetCRpopulation(int i, int SB, int DD, int box, int cfg) {

	TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/minitrees/noNJetCut/bonsai_" + estimated_bkg[DD] + ".root", "READ" );
	//TFile* file = new TFile( "/gpfs/csic_projects/cms/jgarciaf/MVA/minitrees/bonsai_" + estimated_bkg[DD] + ".root", "READ" );

	TTree *tree = (TTree*) file -> Get("latino");

	float MVA, njet, nbjet15, baseW;
 
	tree -> SetBranchAddress( "MVA" + DM[i], &MVA     );
	tree -> SetBranchAddress( "njet"       , &njet    );
	tree -> SetBranchAddress( "nbjet15"    , &nbjet15 );
	tree -> SetBranchAddress( "baseW"      , &baseW   );

	int nentries = tree -> GetEntries();

	float counter = 0; 

	for ( Long64_t ievt = 0; ievt < nentries; ievt++ ) {

		tree -> GetEntry(ievt);

		if( MVA > frontier[i][SB] && MVA < frontier[i][SB+1] ) {

			//if( njet >= wall[0][i][SB][box] && njet < wall[1][i][SB][box] && nbjet15 >= wall[2][i][SB][box] && nbjet15 < wall[3][i][SB][box] ) counter = counter + baseW;
		if( njet >= wall[0][box][cfg] && njet < wall[1][box][cfg] && nbjet15 >= wall[2][box][cfg] && nbjet15 < wall[3][box][cfg] ) counter = counter + baseW;		
			//if( njet >= 0 && njet < 2 && nbjet15 >= 0 && nbjet15 < 2 ) counter = counter + baseW;
		} 

	}
	
	if ( DD != dt) counter = counter * L;

	return counter;

}

