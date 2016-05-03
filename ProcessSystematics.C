const int nsyst    = 22;  // includes nominal, as many as yields files per ttDM point
const int nsamples = 12;  // number of processes considered (includes data & signal)

enum{ nominal, LepMupTup, LepMupTdo, LepElepTup, LepElepTdo, METup, METdo, JESMaxup, JESMaxdo, Btagup, Btagdo, Idisoup, Idisodo, Triggerup, Triggerdo, xsQCDup, acQCDup, xsQCDdo, acQCDdo, xsPDF, acPDF, stat };
enum{          LepMupT,              LepElepT,               MET,          JESMax,             Btag,           Idiso,            Trigger,                         QCD,                              PDF       };
enum{ data, sig, fakes, WZ, ZZ, TT, ST, WW, Zjets, TTV, Wg, HZ };

TString systID[nsyst]; 
TString sampleID[nsamples]; 

//-- just for (tex) tables of yields to the AN-16-011
float S[12][13];  // [ttDM scalar points][process+1] 
float P[ 8][13];  // [ttDM pseudo points][process+1] 

void ProcessSystematics2(int l, TString ttDM);  // first argument neccesary but unimportant


// just the matrix I'll intend to build and manage with this script...  
//
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//		nom	LepMupTup	LepMupTdo	LepElepTup	LepElepTdo	METup		METdo		JESMaxup	JESMaxdo  	... 
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//data     	
//signal
//00_Fakes
//02_WZTo3LNu
//03_ZZ
//04_TTTo2L2Nu
//05_ST
//06_WW
//07_ZJets
//09_TTV
//11_Wg
//14_HZ
//--------------------------------------------------------------------------------------------------------------------------------------------------------------


void ProcessSystematics(){

	// the assignment of TStrings... 

	sampleID[data ] = "data" ; 
	sampleID[sig  ] = "sig"  ; 
	sampleID[fakes] = "fakes"; 
	sampleID[WZ   ] = "WZ"   ; 
	sampleID[ZZ   ] = "ZZ"   ; 
	sampleID[TT   ] = "TT"   ; 
	sampleID[ST   ] = "ST"   ; 
	sampleID[WW   ] = "WW"   ; 
	sampleID[Zjets] = "Zjets"; 
	sampleID[TTV  ] = "TTV"  ; 
	sampleID[Wg   ] = "Wg"   ; 
	sampleID[HZ   ] = "HZ"   ; 

	systID[nominal   ] = "nominal"   ; 
	systID[LepMupTup ] = "LepMupTup" ;  
	systID[LepMupTdo ] = "LepMupTdo" ;  
	systID[LepElepTup] = "LepElepTup";  
	systID[LepElepTdo] = "LepElepTdo";  
	systID[METup     ] = "METup"     ;  
	systID[METdo     ] = "METdo"     ;  
	systID[JESMaxup  ] = "JESMaxup"  ;  
	systID[JESMaxdo  ] = "JESMaxdo"  ;  
	systID[Btagup    ] = "Btagup"    ;  
	systID[Btagdo    ] = "Btagdo"    ;  
	systID[Idisoup   ] = "Idisoup"   ;  
	systID[Idisodo   ] = "Idisodo"   ;  
	systID[Triggerup ] = "Triggerup" ;  
	systID[Triggerdo ] = "Triggerdo" ;  
	systID[xsQCDup   ] = "xsQCDup"   ;  
	systID[acQCDup   ] = "acQCDup"   ;  
	systID[xsQCDdo   ] = "xsQCDdo"   ;  
	systID[acQCDdo   ] = "acQCDdo"   ;  
	systID[xsPDF     ] = "xsPDF"     ;  
	systID[acPDF     ] = "acPDF"     ;  
	systID[stat      ] = "stat"      ;  



	// repeating for evry ttDM point...  

	ProcessSystematics2( 0, "ttDM0001scalar0010");
	//ProcessSystematics2( 1, "ttDM0001scalar0050");
	//ProcessSystematics2( 2, "ttDM0001scalar0100");
	//ProcessSystematics2( 3, "ttDM0001scalar0200");
	//ProcessSystematics2( 4, "ttDM0001scalar0300");
	//ProcessSystematics2( 5, "ttDM0001scalar0500");
	////ProcessSystematics2( 6, "ttDM0010scalar0010");
	//ProcessSystematics2( 7, "ttDM0010scalar0050");
	//ProcessSystematics2( 8, "ttDM0010scalar0100");
	////ProcessSystematics2( 9, "ttDM0050scalar0050");
	//ProcessSystematics2(10, "ttDM0050scalar0200");
	///ProcessSystematics2(11, "ttDM0050scalar0300");
	////ProcessSystematics2("ttDM0150scalar0200");
	////ProcessSystematics2("ttDM0500scalar0500");

	//ProcessSystematics2( 0, "ttDM0001pseudo0010");
	//ProcessSystematics2( 1, "ttDM0001pseudo0020");
	//ProcessSystematics2( 2, "ttDM0001pseudo0050");
	//ProcessSystematics2( 3, "ttDM0001pseudo0100");
	//ProcessSystematics2( 4, "ttDM0001pseudo0200");
	//ProcessSystematics2( 5, "ttDM0001pseudo0500");
	////ProcessSystematics2("ttDM0010pseudo0100");
	//ProcessSystematics2( 6, "ttDM0050pseudo0200");
	//ProcessSystematics2( 7, "ttDM0150pseudo0200");
	////ProcessSystematics2("ttDM0150pseudo0500");


	//-- (tex) tables of yields to the AN-16-011 
	//-- output goes to subfolder tmva/tex
	//-- scalar separated from pseudo

	/*ofstream sYields;

	sYields.open("/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/tmva/tex/mva_sYields.txt");

	sYields << Form("WW               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",

					  S[0][WW], S[1][WW], S[2][WW], S[3][WW], S[4][WW], S[5][WW], S[6][WW], S[7][WW], S[8][WW], S[9][WW], S[10][WW], S[11][WW]                           );

	sYields << Form("WZ               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",

					  S[0][WZ], S[1][WZ], S[2][WZ], S[3][WZ], S[4][WZ], S[5][WZ], S[6][WZ], S[7][WZ], S[8][WZ], S[9][WZ], S[10][WZ], S[11][WZ]                           );

	sYields << Form("ZZ               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",

					  S[0][ZZ], S[1][ZZ], S[2][ZZ], S[3][ZZ], S[4][ZZ], S[5][ZZ], S[6][ZZ], S[7][ZZ], S[8][ZZ], S[9][ZZ], S[10][ZZ], S[11][ZZ]                           );
          
	sYields << Form("W$\\gamma$        & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",  

					  S[0][Wg], S[1][Wg], S[2][Wg], S[3][Wg], S[4][Wg], S[5][Wg], S[6][Wg], S[7][Wg], S[8][Wg], S[9][Wg], S[10][Wg], S[11][Wg]                           );     

	sYields << Form("Z+jets           & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",

				S[0][Zjets], S[1][Zjets], S[2][Zjets], S[3][Zjets], S[4][Zjets], S[5][Zjets], S[6][Zjets], S[7][Zjets], S[8][Zjets], S[9][Zjets], S[10][Zjets], S[11][Zjets] );
              
	sYields << Form("$\\ttbar$V        & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",

					  S[0][TTV], S[1][TTV], S[2][TTV], S[3][TTV], S[4][TTV], S[5][TTV], S[6][TTV], S[7][TTV], S[8][TTV], S[9][TTV], S[10][TTV], S[11][TTV]               );
                
	sYields << Form("$\\ttbar$         & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n", 

					  S[0][TT], S[1][TT], S[2][TT], S[3][TT], S[4][TT], S[5][TT], S[6][TT], S[7][TT], S[8][TT], S[9][TT], S[10][TT], S[11][TT]                           );
               
	sYields << Form("tW               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",

					  S[0][ST], S[1][ST], S[2][ST], S[3][ST], S[4][ST], S[5][ST], S[6][ST], S[7][ST], S[8][ST], S[9][ST], S[10][ST], S[11][ST]                           );   

	sYields << Form("non-prompt       & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \\hline \n",

			        S[0][fakes], S[1][fakes], S[2][fakes], S[3][fakes], S[4][fakes], S[5][fakes], S[6][fakes], S[7][fakes], S[8][fakes], S[9][fakes], S[10][fakes], S[11][fakes] );     

	sYields << Form("total background & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\          \n", 

					  S[0][12], S[1][12], S[2][12], S[3][12], S[4][12], S[5][12], S[6][12], S[7][12], S[8][12], S[9][12], S[10][12], S[11][12]                           );

 	sYields << Form("data             & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \n",  

			 S[0][data], S[1][data], S[2][data], S[3][data], S[4][data], S[5][data], S[6][data], S[7][data], S[8][data], S[9][data], S[10][data], S[11][data]                    );  

	sYields << Form("signal           & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f \\\\ \\hline  \n",   

					  S[0][sig], S[1][sig], S[2][sig], S[3][sig], S[4][sig], S[5][sig], S[6][sig], S[7][sig], S[8][sig], S[9][sig], S[10][sig], S[11][sig]               );     

	sYields.close();*/


	/*ofstream pYields;

	pYields.open("/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/tmva/tex/mva_pYields.txt");

	pYields << Form("WW               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n",

					  P[0][WW], P[1][WW], P[2][WW], P[3][WW], P[4][WW], P[5][WW], P[6][WW], P[7][WW]                         );

	pYields << Form("WZ               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n",

					  P[0][WZ], P[1][WZ], P[2][WZ], P[3][WZ], P[4][WZ], P[5][WZ], P[6][WZ], P[7][WZ]                         );

	pYields << Form("ZZ               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n",

					  P[0][ZZ], P[1][ZZ], P[2][ZZ], P[3][ZZ], P[4][ZZ], P[5][ZZ], P[6][ZZ], P[7][ZZ]                         );
          
	pYields << Form("W$\\gamma$        & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n",  

					  P[0][Wg], P[1][Wg], P[2][Wg], P[3][Wg], P[4][Wg], P[5][Wg], P[6][Wg], P[7][Wg]                         );     

	pYields << Form("Z+jets           & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n",

				P[0][Zjets], P[1][Zjets], P[2][Zjets], P[3][Zjets], P[4][Zjets], P[5][Zjets], P[6][Zjets], P[7][Zjets]           );
              
	pYields << Form("$\\ttbar$V        & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n",

					  P[0][TTV], P[1][TTV], P[2][TTV], P[3][TTV], P[4][TTV], P[5][TTV], P[6][TTV], P[7][TTV]                 );
                
	pYields << Form("$\\ttbar$         & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n", 

					  P[0][TT], P[1][TT], P[2][TT], P[3][TT], P[4][TT], P[5][TT], P[6][TT], P[7][TT]                         );
               
	pYields << Form("tW               & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\ \n",

					  P[0][ST], P[1][ST], P[2][ST], P[3][ST], P[4][ST], P[5][ST], P[6][ST], P[7][ST]                         );   

	pYields << Form("non-prompt       & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                          \\\\ \\hline \n",

			        P[0][fakes], P[1][fakes], P[2][fakes], P[3][fakes], P[4][fakes], P[5][fakes], P[6][fakes], P[7][fakes]           );     

	pYields << Form("total background & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                 \\\\\n", 

					  P[0][12], P[1][12], P[2][12], P[3][12], P[4][12], P[5][12], P[6][12], P[7][12]                         );

 	pYields << Form("data             & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                                  \\\\ \n",  

			 P[0][data], P[1][data], P[2][data], P[3][data], P[4][data], P[5][data], P[6][data], P[7][data]                          );  

	pYields << Form("signal           & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f & %7.2f                          \\\\ \\hline  \n",   

					  P[0][sig], P[1][sig], P[2][sig], P[3][sig], P[4][sig], P[5][sig], P[6][sig], P[7][sig]                 );     

	pYields.close();*/

}



void ProcessSystematics2(int l, TString ttDM){

	cout << "                       " << endl;
	cout << "                       " << endl;
	cout << "                       " << endl;
	cout << " ===================== " << endl;
	cout << " == " << ttDM            << endl;
	cout << " ===================== " << endl;
	cout << "                       " << endl;

	float yield[nsyst][nsamples]; 


	//-- loading the yields from *.dat files... 

	for( int i = 0; i < nsyst; i++ ){  

		ifstream syst; 

	  	syst.open("/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS_old/tmva/yields/" + ttDM + "_" + systID[i] + ".dat"); 

		int j = 0; 

		while( 1 ){      

			if (!syst.good()) break;

			if (j > nsamples) break;
		
			syst >> yield[i][j];   

			j++;

		}

		syst.close();

	}



	if ( yield[nominal][fakes] < 0 )  yield[nominal][fakes] = 0; 				// protection 

	yield[acQCDup][ZZ] = 1.00; yield[acQCDdo][ZZ] = 1.00; yield[acPDF][ZZ] = 1.00;	        // protection 
	yield[acQCDup][ST] = 1.00; yield[acQCDdo][ST] = 1.00; yield[acPDF][ST] = 1.00;	        // protection 



	//-- for the weighting you'll need... 

	float totalBkg = yield[nominal][fakes] + yield[nominal][WZ]    + yield[nominal][ZZ]  + yield[nominal][TT] + yield[nominal][ST]    
                       + yield[nominal][WW]    + yield[nominal][Zjets] + yield[nominal][TTV] + yield[nominal][Wg] + yield[nominal][HZ];


	float otherBkg = totalBkg - yield[nominal][TT] - yield[nominal][WW] - yield[nominal][Zjets];   // bkgs directly from MC 


	//-- some (10) systematics will need a bit of postprocessing before being written in the datacards...

	float mean_yield[10][nsamples];

	float total_LepMupT  = 0;
	float total_LepElepT = 0;
	float total_MET      = 0;
	float total_JESMax   = 0;
	float total_Btag     = 0;
	float total_Idiso    = 0;
	float total_Trigger  = 0;
	float total_acQCD    = 0; 
	float total_acPDF    = 0; 

	
	//-- loop over processes

	for(int j = 0; j < nsamples; j++){

		//-- just printing the preavious readout on the screen...

		//cout << yield[nominal][j] << " -- " << yield[LepMupTup][j] << " -- " << yield[LepMupTdo][j] << " -- " << yield[LepElepTup][j] << " -- " << yield[LepElepTdo][j] << " -- " 
                     //                               << yield[METup]    [j] << " -- " << yield[METdo]    [j] << " -- " << yield[JESMaxup]  [j] << " -- " << yield[JESMaxdo]  [j] << " -- " 
                     //<< yield[Btagup] [j] << " -- " << yield[Btagdo]   [j] << " -- " << yield[Idisoup]  [j] << " -- " << yield[Idisodo]   [j] << " -- " << yield[Triggerup] [j] << " -- "  << yield[Triggerdo][j] 
		     //<< yield[xsQCDup][j] << " -- " << yield[acQCDup]  [j] << " -- " << yield[xsQCDdo]  [j] << " -- " << yield[acQCDdo]   [j] << " -- " << yield[xsPDF    ] [j] << " -- "  << yield[acPDF    ][j] << endl;
 

		float weight = yield[nominal][j]/totalBkg;   


		//-- some tricky stuff to get reasonable numbers... 

		mean_yield[LepMupT ][j] = (  abs( yield[LepMupTup ][j] - yield[LepMupTdo ][j] )  ) / ( yield[LepMupTup ][j] + yield[LepMupTdo ][j] );//* weight;
		mean_yield[LepElepT][j] = (  abs( yield[LepElepTup][j] - yield[LepElepTdo][j] )  ) / ( yield[LepElepTup][j] + yield[LepElepTdo][j] );//* weight;
		mean_yield[MET     ][j] = (  abs( yield[METup     ][j] - yield[METdo     ][j] )  ) / ( yield[METup     ][j] + yield[METdo     ][j] );//* weight;
		mean_yield[JESMax  ][j] = (  abs( yield[JESMaxup  ][j] - yield[JESMaxdo  ][j] )  ) / ( yield[JESMaxup  ][j] + yield[JESMaxdo  ][j] );//* weight;
		mean_yield[Btag    ][j] = (  abs( yield[Btagup    ][j] - yield[nominal][j] ) + abs( yield[Btagdo    ][j] - yield[nominal][j] )  ) / (  2* yield[nominal][j]  );//* weight;
		mean_yield[Idiso   ][j] = (  abs( yield[Idisoup   ][j] - yield[nominal][j] ) + abs( yield[Idisodo   ][j] - yield[nominal][j] )  ) / (  2* yield[nominal][j]  );//* weight; 
		mean_yield[Trigger ][j] = (  abs( yield[Triggerup ][j] - yield[nominal][j] ) + abs( yield[Triggerdo ][j] - yield[nominal][j] )  ) / (  2* yield[nominal][j]  );//* weight; 
		mean_yield[QCD     ][j] = (  abs( yield[acQCDup   ][j] - 1                 ) + abs( yield[acQCDdo   ][j] - 1                 )  ) / (  2                     );//* weight;   
		mean_yield[PDF     ][j] =    abs( yield[acPDF     ][j] - 1                 )                                                                                  ;//* weight;   
		mean_yield[9       ][j] =    yield[stat][j]/yield[nominal][j]; 



		if( j > 1 ){

			if( j == ST || j == WZ || j == ZZ || j == TTV ||j == Wg ||j == HZ ){   // MC bkgs

		 	total_LepMupT  += pow(mean_yield[LepMupT ][j], 2);  //cout << total_LepMupT << endl;
		 	total_LepElepT += pow(mean_yield[LepElepT][j], 2);  //cout << total_LepMupT << endl;
			total_MET      += pow(mean_yield[MET     ][j], 2);  //cout << total_LepMupT << endl;
		 	total_JESMax   += pow(mean_yield[JESMax  ][j], 2);  //cout << total_LepMupT << endl;
			total_Btag     += pow(mean_yield[Btag    ][j], 2);
			total_Idiso    += pow(mean_yield[Idiso   ][j], 2);
			total_Trigger  += pow(mean_yield[Trigger ][j], 2); 
			total_acQCD    += pow(mean_yield[QCD     ][j], 2);
			total_acPDF    += pow(mean_yield[PDF     ][j], 2);  //cout << total_acPDF << endl;

			}

		}
	
	}

	float total_LeppT = sqrt( total_LepMupT + total_LepElepT );   					      // combining both muon and electron LES for MC-bkgs
	float tt_LeppT    = sqrt(  pow( mean_yield[LepMupT][TT], 2) + pow( mean_yield[LepElepT][TT], 2)  );   // combining both muon and electron LES for TT (will go for signal)

	total_LepMupT  = sqrt(total_LepMupT );    //cout << total_LepMupT  << endl; 
	total_LepElepT = sqrt(total_LepElepT);    //cout << total_LepElepT << endl;
	total_MET      = sqrt(total_MET     );    //cout << total_MET      << endl;
	total_JESMax   = sqrt(total_JESMax  );    //cout << total_JESMax   << endl;
	total_Btag     = sqrt(total_Btag    );    //cout << total_Btag     << endl;
	total_Idiso    = sqrt(total_Idiso   );    //cout << total_Idiso    << endl;
	total_Trigger  = sqrt(total_Trigger );    //cout << total_Trigger  << endl;
	total_acQCD    = sqrt(total_acQCD   );    //cout << total_acQCD    << endl;
	total_acPDF    = sqrt(total_acPDF   );    //cout << total_acPDF    << endl;

	float total_lumi  = 0.027 * otherBkg / totalBkg; 
	float total_xsMC  = 0.020 * otherBkg / totalBkg;
	float total_boxes = sqrt(   pow( 0.07 * yield[nominal][TT] / totalBkg,  2 ) + pow( 0.32 * yield[nominal][ST] / totalBkg,  2 )   );
	float total_Zjets = 0.041 * yield[nominal][Zjets]/ totalBkg;
	float total_fakes = 0.30  * yield[nominal][fakes]/ totalBkg;
	

	//-- first systematics table on the note... 

	/*ofstream sysTable;

	sysTable.open(Form("/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/tmva/tex/%s_sysTable.txt", ttDM.Data()));

	sysTable << Form("lepton trigger        & %3.1f & %3.1f \\\\  \n", 100*total_Trigger, 100*mean_yield[Trigger][TT]);
	sysTable << Form("lepton ID and iso     & %3.1f & %3.1f \\\\  \n", 100*total_Idiso  , 100*mean_yield[Idiso  ][TT]);
	sysTable << Form("lepton energy scale   & %3.1f & %3.1f \\\\  \n", 100*total_LeppT  , 100*tt_LeppT               );
	sysTable << Form("jet energy scale      & %3.1f & %3.1f \\\\  \n", 100*total_JESMax , 100*mean_yield[JESMax ][TT]);
	sysTable << Form("b-tagging             & %3.1f & %3.1f \\\\  \n", 100*total_Btag   , 100*mean_yield[Btag   ][TT]);
	sysTable << Form("$\\MET$                & %3.1f & %3.1f \\\\ \n", 100*total_MET    , 100*mean_yield[MET    ][TT]);
	sysTable << Form("$\\ttbar$ and WW       & %3.1f & %3s \\\\   \n", 100*total_boxes  , " - "                      );
	sysTable << Form("Z+jets                & %3.1f & %3s \\\\    \n", 100*total_Zjets  , " - "                      );
	sysTable << Form("non-prompt            & %3.1f & %3s \\\\    \n", 100*total_fakes  , " - "                      );
	sysTable << Form("xs MC                 & %3.1f & %3s \\\\    \n", 100*total_xsMC   , " - "                      );
	sysTable << Form("QCD                   & %3.1f & %3.1f \\\\  \n", 100*total_acQCD  , 100*mean_yield[acQCD  ][TT]);
	sysTable << Form("PDF                   & %3.1f & %3.1f \\\\  \n", 100*total_acPDF  , 100*mean_yield[acPDF  ][TT]);
	sysTable << Form("luminosity            & %3.1f & %3.1f \\\\  \n", 100*total_lumi   , 100*0.027                  );

	sysTable.close();*/



	//-- first datacard used (compressed version) -> DEPRECATED ... 

	/*ofstream datacard;

	datacard.open(Form("/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/tmva/datacards/%s.txt", ttDM.Data()));

	datacard << "imax 1   number of channels\n";
	datacard << "jmax 1   number of backgrounds\n";
	datacard << "kmax 4   number of nuisance parameters\n";
	datacard << "------------\n";
	datacard << "\n";
	datacard << "bin 1\n";
	datacard << Form("observation %f\n", yield[nominal][data]);
	datacard << "------------\n";
	datacard << "\n";
	datacard << "bin             1           1\n";
	datacard << "process         DM          bkg\n";
	datacard << "process         0           1\n";
	datacard << Form("rate            %f    %f\n", yield[nominal][sig], totalBkg);
	datacard << "------------\n";
	datacard << Form("LepMupT    lnN	%f	%f \n", 1+mean_yield[LepMupT ][sig], 1+total_LepMupT );
	datacard << Form("LepElepT   lnN	%f	%f \n", 1+mean_yield[LepElepT][sig], 1+total_LepElepT);
	datacard << Form("MET        lnN	%f	%f \n", 1+mean_yield[MET     ][sig], 1+total_MET     );
	datacard << Form("JESMax     lnN	%f	%f \n", 1+mean_yield[JESMax  ][sig], 1+total_JESMax  );

	datacard.close();*/


	//-- protection against Combine module when the expected signal rate is too small... 
	
	TString mysuffix = "rescale001"; 

		if( yield[nominal][sig] < 0.1 ){

			yield[nominal][sig] = yield[nominal][sig] * 100; 

			mysuffix = "rescale100";

			cout << ttDM << " has been rescaled" << endl;
		
		} 

	//-- the actual datacard is constructed from here on... better with a loop? 

	ofstream datacard;

	datacard.open(Form("/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/tmva/datacards/%s_%s.txt", ttDM.Data(), mysuffix.Data() ));

	datacard << "imax 1   number of channels                                                                                                                                                                  \n" ;
	datacard << "jmax 10  number of backgrounds                                                                                                                                                               \n" ;
	datacard << "kmax 15  number of nuisance parameters                                                                                                                                                       \n" ;
	datacard << "------------                                                                                                                                                                                 \n" ;
	datacard << "                                                                                                                                                                                             \n" ;
	datacard << "bin 1                                                                                                                                                                                        \n" ;
	datacard << Form("observation %5.0f                                                                                                                                                 \n", yield[nominal][data]);
	datacard << "------------                                                                                                                                                                                 \n" ;
	datacard << "                                                                                                                                                                                             \n" ;
	datacard << "bin                1                1                1                1                1                1         	 1                1                1                1                1    \n" ;
	datacard << "process            DM               TT               ST               WW               DY               fakes            WZ               ZZ               TTV              Wg               HZ  \n" ;
	datacard << "process            0                1                2                3                4                5                6                7                8                9               10   \n" ;
	datacard << Form("rate          	%7.2f          %7.2f          %7.2f          %7.2f          %7.2f          %7.2f          %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n", 

                                                yield[nominal][sig], yield[nominal][TT], yield[nominal][ST] , yield[nominal][WW], yield[nominal][Zjets], yield[nominal][fakes],

		                                yield[nominal][WZ] , yield[nominal][ZZ], yield[nominal][TTV], yield[nominal][Wg], yield[nominal][HZ]                                                                 );

	datacard << "------------                                                                                                                                                                                 \n" ;


	datacard << Form("LepMupT    lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n", 

		    1+mean_yield[LepMupT][TT]  , "-"                       , 1+mean_yield[LepMupT][ST]  , "-"                       , "-"                       , "-"  ,

		    1+mean_yield[LepMupT][WZ]  , 1+mean_yield[LepMupT][ZZ] , 1+mean_yield[LepMupT][TT] , 1+mean_yield[LepMupT][TT] , 1+mean_yield[LepMupT][HZ]                                                      );

	datacard << Form("LepElepT   lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

		    1+mean_yield[LepElepT][TT] ,  "-"                      , 1+mean_yield[LepElepT][ST] , "-"                       , "-"                       , "-"  , 

		    1+mean_yield[LepElepT][WZ] , 1+mean_yield[LepElepT][ZZ], 1+mean_yield[LepElepT][TT], 1+mean_yield[LepElepT][TT], 1+mean_yield[LepElepT][HZ]                                                     );

	datacard << Form("MET        lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

	            1+mean_yield[MET][TT]      , "-"                       , 1+mean_yield[MET][ST]      , "-"                       , "-"                       , "-"  ,                         
  
		    1+mean_yield[MET][WZ]      , 1+mean_yield[MET][ZZ]     , 1+mean_yield[MET][TT]     , 1+mean_yield[MET][TT]     , 1+mean_yield[MET][HZ]                                                          );

	datacard << Form("JESMax     lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

	            1+mean_yield[JESMax][TT]   , "-"                       , 1+mean_yield[JESMax][ST]   , "-"                       , "-"                       , "-"  ,
  
		    1+mean_yield[JESMax][WZ]   , 1+mean_yield[JESMax][ZZ]  , 1+mean_yield[JESMax][TT]  , 1+mean_yield[JESMax][TT]  , 1+mean_yield[JESMax][HZ]                                                       );

	datacard << Form("Btag       lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

	            1+mean_yield[Btag][TT]     , "-"                       , 1+mean_yield[Btag][ST]     , "-"                       , "-"                       , "-"  ,
  
		    1+mean_yield[Btag][WZ]     , 1+mean_yield[Btag][ZZ]    , 1+mean_yield[Btag][TT]    , 1+mean_yield[Btag][TT]    , 1+mean_yield[Btag][HZ]                                                         );

	datacard << Form("Idiso      lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

	            1+mean_yield[Idiso][TT]    , "-"                       , 1+mean_yield[Idiso][ST]    , "-"                       , "-"                       , "-"  ,
  
		    1+mean_yield[Idiso][WZ]    , 1+mean_yield[Idiso][ZZ]   , 1+mean_yield[Idiso][TT]   , 1+mean_yield[Idiso][TT]   , 1+mean_yield[Idiso][HZ]                                                        );

	datacard << Form("Trigger    lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

	            1+mean_yield[Trigger][TT]  , "-"                       , 1+mean_yield[Trigger][ST]  , "-"                       , "-"                       , "-"  ,
  
		    1+mean_yield[Trigger][WZ]  , 1+mean_yield[Trigger][ZZ] , 1+mean_yield[Trigger][TT] , 1+mean_yield[Trigger][TT] , 1+mean_yield[Trigger][HZ]                                                      );

	datacard << Form("QCDacc     lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

	            1+mean_yield[QCD  ][TT]    , "-"                         , 1+mean_yield[QCD  ][ST]    , "-"                        , "-"                       , "-"  ,
    
		    1+mean_yield[QCD  ][WZ]    ,1+ mean_yield[QCD  ][ZZ]   , 1+mean_yield[QCD  ][TT]   , 1+mean_yield[QCD  ][TT]    , 1+mean_yield[QCD  ][TT]                                                       );
                                                   

	datacard << Form("PDFacc     lnN	%7.2f          %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",

	            1+mean_yield[PDF][TT]      , "-"                       , 1+mean_yield[PDF][ST]       , "-"                      , "-"                        , "-"  ,
  
		    1+mean_yield[PDF][WZ]      , 1.00                      , 1+mean_yield[PDF][TT]      , 1.00                     , 1+mean_yield[PDF][HZ]                                                               );

	datacard << Form("tt_WW      lnN	%5s            %7.2f          %5s            %7.2f          %5s            %5s            %5s            %5s            %5s            %5s            %5s         \n",

	            "-"                        , 1.07                      , "-"                         , 1.32                     , "-"                       , "-"  ,
  
		    "-"                        , "-"                       , "-"                         , "-"                      , "-"                                                                            );

	datacard << Form("DY         lnN	%5s            %5s            %5s            %5s            %7.2f          %5s            %5s            %5s            %5s            %5s            %5s         \n",

	            "-"                        , "-"                       , "-"                         , "-"                      , 1.041                     , "-"  ,
  
		    "-"                        , "-"                       , "-"                         , "-"                      , "-"                                                                            );

	datacard << Form("fakes      lnN	%5s            %5s            %5s            %5s            %5s            %7.2f            %5s            %5s            %5s            %5s            %5s         \n",

	            "-"                        , "-"                       , "-"                         , "-"                      , "-"                       , 1.30 ,
            
                    "-"                        , "-"                       , "-"                         , "-"                      , "-"                                                                            );

	datacard << Form("xsMC       lnN	%5s           %5s            %7.2f          %5s            %5s            %5s            %7.2f          %7.2f          %7.2f          %7.2f          %7.2f       \n",
 
                    "-"                        , "-"                       , 1+0.20                      , "-"                      , "-"                       , "-"  ,
  
		    1+0.20                     , 1+0.20                    , 1+0.20                      , 1+0.20                   , 1+0.20                                                                         );

	datacard << Form("stat       lnN	%7.3f          %5s            %7.3f          %5s            %5s            %5s            %7.3f          %7.3f          %7.3f          %7.3f          %7.3f       \n",

	            1+mean_yield[9  ][TT]    , "-"                         , 1+mean_yield[9  ][ST]    , "-"                        , "-"                       , "-"  ,
    
		    1+mean_yield[9  ][WZ]    ,1+ mean_yield[9  ][ZZ]   , 1+mean_yield[9  ][TT]   , 1+mean_yield[9  ][TT]    , 1+mean_yield[9  ][HZ]                                                       );

	datacard << Form("Lumi       lnN	%7.3f          %5s            %7.3f          %5s            %5s            %5s            %7.3f          %7.3f          %7.3f          %7.3f          %7.3f       \n",

	            1+0.027                    , "-"                       , 1+0.027                     , "-"                      , "-"                       , "-"  ,
  
		    1+0.027                    , 1+0.027                   , 1+0.027                     , 1+0.027                  , 1+0.027                                                                        );

	datacard.close();


	//-- txt tables for preliminary "monitoring" 
	//-- they fall into tmva/txt subfolder

	/*ofstream table;

	table.open(Form("/gpfs/csic_projects/cms/jgarciaf/CMSSW_7_6_3/src/AnalysisCMS/tmva/txt/%s.txt", ttDM.Data()));

	table << Form("%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\n",
	              "", "nom", "Muup", "Mudo", "Eleup", "Eledo", "METup", "METdo", "JESup", "JESdo", "Btagup", "Btagdo", "Idisoup", "Idisodo", "Trgup", "Trgdo", "acQCDup", "acQCDdo", "acPDF" );

	for(int j = 0; j < nsamples; j++){

	table << Form("%7s\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n",
                      sampleID[j].Data(), yield[nominal][j], yield[LepMupTup][j], yield[LepMupTdo][j], yield[LepElepTup][j], yield[LepElepTdo][j], 
                      yield[METup][j], yield[METdo][j], yield[JESMaxup][j], yield[JESMaxdo][j], yield[Btagup][j], yield[Btagdo][j], 
                      yield[Idisoup][j], yield[Idisodo][j], yield[Triggerup][j], yield[Triggerdo][j],                                    
		      yield[acQCDup][j]*yield[nominal][j], yield[acQCDdo][j]*yield[nominal][j], yield[acPDF][j]*yield[nominal][j]                                                  );

	}

	table << "\n"; 
 	table << "\n"; 
	table << "\n";

	//table << Form("%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\n",
	//              "", "nom", "Muup", "Mudo", "Eleup", "Eledo", "METup", "METdo", "JESup", "JESdo", "Btagup", "Btagdo", "Idisoup", "Idisodo", "Trgup", "Trgdo", "xsQCDup", "xsQCDdo", "acQCDup", "acQCDdo", "xsPDF", "acPDF" );

	//for(int j = 0; j < nsamples; j++){

	//table << Form("%7s\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n",
        //              sampleID[j].Data(), yield[nominal][j], (yield[LepMupTup][j]-yield[nominal][j])/yield[nominal][j],  (yield[LepMupTdo][j]-yield[nominal][j])/yield[nominal][j], 
        //                                                     (yield[LepElepTup][j]-yield[nominal][j])/yield[nominal][j], (yield[LepElepTdo][j]-yield[nominal][j])/yield[nominal][j], 
        //                                                     (yield[METup][j]-yield[nominal][j])/yield[nominal][j],      (yield[METdo][j]-yield[nominal][j])/yield[nominal][j], 
        //                                                     (yield[JESMaxup][j]-yield[nominal][j])/yield[nominal][j],   (yield[JESMaxdo][j]-yield[nominal][j])/yield[nominal][j], 
        //                                                     (yield[Btagup][j]-yield[nominal][j])/yield[nominal][j],     (yield[Btagdo][j]-yield[nominal][j])/yield[nominal][j], 
        //                                                     (yield[Idisoup][j]-yield[nominal][j])/yield[nominal][j],    (yield[Idisodo][j]-yield[nominal][j])/yield[nominal][j], 
        //                                                     (yield[Triggerup][j]-yield[nominal][j])/yield[nominal][j],  (yield[Triggerdo][j]-yield[nominal][j])/yield[nominal][j]      );

	//}


	table << Form("%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\t%7s\n",
	              "", "nom", "Mu", "Ele", "MET", "JES", "Btagup", "Btagdo", "Idisoup", "Idisodo", "Trgup", "Trgdo", "acQCDup", "acQCDdo", "acPDF" );

	for(int j = 0; j < nsamples; j++){

	table << Form("%7s\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n",
                      sampleID[j].Data(), yield[nominal][j], (yield[LepMupTup][j]-yield[LepMupTdo][j])/(yield[LepMupTdo][j]+yield[LepMupTdo][j]), 
							     (yield[LepElepTup][j]-yield[LepElepTdo][j])/(yield[LepElepTdo][j]+yield[LepElepTdo][j]), 
                                                             (yield[METup][j]-yield[METdo][j])/(yield[METdo][j]+yield[METdo][j]), 
                                                             (yield[JESMaxup][j]-yield[JESMaxdo][j])/(yield[JESMaxdo][j]+yield[JESMaxdo][j]), 
                                                             (yield[Btagup][j]-yield[nominal][j])/yield[nominal][j],     (yield[Btagdo][j]-yield[nominal][j])/yield[nominal][j], 
                                                             (yield[Idisoup][j]-yield[nominal][j])/yield[nominal][j],    (yield[Idisodo][j]-yield[nominal][j])/yield[nominal][j], 
                                                             (yield[Triggerup][j]-yield[nominal][j])/yield[nominal][j],  (yield[Triggerdo][j]-yield[nominal][j])/yield[nominal][j],      
							     yield[acQCDup][j], yield[acQCDdo][j], yield[acPDF][j]	                                                              );

	}


	table.close();*/



	//-- some mess for the (tex) tables of yields... 

	//for(int j = 0; j < nsamples; j++){
	
		//S[l][j] = yield[nominal][j];
		//P[l][j] = yield[nominal][j];

	//}

	//S[l][12] = yield[nominal][fakes] + yield[nominal][WZ]    + yield[nominal][ZZ]  + yield[nominal][TT] + yield[nominal][ST] 
        //         + yield[nominal][WW]    + yield[nominal][Zjets] + yield[nominal][TTV] + yield[nominal][Wg] + yield[nominal][HZ];

	//P[l][12] = yield[nominal][fakes] + yield[nominal][WZ]    + yield[nominal][ZZ]  + yield[nominal][TT] + yield[nominal][ST] 
        //         + yield[nominal][WW]    + yield[nominal][Zjets] + yield[nominal][TTV] + yield[nominal][Wg] + yield[nominal][HZ];

}


