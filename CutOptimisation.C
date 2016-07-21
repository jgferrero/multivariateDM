const int ncut = 50; 
enum{ data, sig, fakes, WZ, ZZ, TT, ST, WW, Zjets, TTV, Wg, HZ, nprocess }; 
const int nMLP = 7; 

TString MLP[nMLP]; 

void CutOptimisation(){

	int iMLP = 0; 

	MLP[iMLP++] = "01/"; 
	MLP[iMLP++] = "02/";
	MLP[iMLP++] = "03/";
	MLP[iMLP++] = "04/";
	MLP[iMLP++] = "05/";
	MLP[iMLP++] = "06/";
	MLP[iMLP++] = "07/";
	

	for( int k = 0; k < nMLP; k++ ){

	float cut = 0; 

	for( int i = 0; i < ncut; i++ ){

		cut = i * 1./ncut; 

		float yield[nprocess];

		ifstream syst; 

	  	syst.open(Form("/gpfs/csic_projects/cms/jgarciaf/CMSSW_8_0_5/src/AnalysisCMS/tmva/yields/" + MLP[k] + "ttDM0001scalar00500_%3.2f.dat", cut)); 
	
		int j = 0; 

		while( 1 ){      

			if (!syst.good()) break;

			if (j > nprocess) break;
		
			syst >> yield[j];     

			j++;

		}

		//cout << yield[data] << endl;

		syst.close();


		ofstream datacard;

		datacard.open(Form("/gpfs/csic_projects/cms/jgarciaf/CMSSW_8_0_5/src/AnalysisCMS/tmva/datacards/" + MLP[k] + "ttDM0001scalar00500_%3.2f.txt", cut));

		datacard << "imax 1   number of channels                                                                                                                                                                  \n" ;
		datacard << "jmax 10  number of backgrounds                                                                                                                                                               \n" ;
		datacard << "kmax 15  number of nuisance parameters                                                                                                                                                       \n" ;
		datacard << "------------                                                                                                                                                                                 \n" ;
		datacard << "                                                                                                                                                                                             \n" ;
		datacard << "bin 1                                                                                                                                                                                        \n" ;
		datacard << Form("observation %5.0f                                                                                                                                                          \n", yield[data]);
		datacard << "------------                                                                                                                                                                                 \n" ;
		datacard << "                                                                                                                                                                                             \n" ;
		datacard << "bin        \t     \t   1   \t   1   \t   1   \t   1   \t   1   \t   1   \t   1   \t   1   \t   1   \t   1   \t   1 \n" ;
		datacard << "process    \t     \t  DM   \t  TT   \t  ST   \t  WW   \t  DY   \t fakes \t  WZ   \t  ZZ   \t  TTV  \t  Wg   \t  HZ \n" ;
		datacard << "process    \t     \t   0   \t   1   \t   2   \t   3   \t   4   \t   5   \t   6   \t   7   \t   8   \t   9   \t  10 \n" ;
		datacard << Form("rate       \t     \t%7.3f \t%7.3f \t%7.3f \t%7.3f \t%7.3f \t%7.3f \t%7.3f \t%7.3f \t%7.3f \t%7.3f \t%7.3f \n", 

		                                        100*yield[sig], yield[TT], yield[ST] , yield[WW], yield[Zjets], yield[fakes],

				                        yield[WZ] , yield[ZZ], yield[TTV], yield[Wg], yield[HZ]                                                                 );

		datacard << "------------                                                                                                                                                                                 \n" ;

//datacard << "LepMupT    lnN	  1.001              -              1.007              -                -                -              1.008            1.000            1.001            1.001            1.011       \n" ;
//datacard << "LepElepT   lnN	  1.001              -              1.003              -                -                -              1.007            1.011            1.001            1.001            1.011       \n" ;
//datacard << "MET        lnN	  1.017              -              1.014              -                -                -              1.004            1.032            1.017            1.017            1.022       \n" ;
//datacard << "JESMax     lnN	  1.056              -              1.054              -                -                -              1.064            1.034            1.056            1.056            1.049       \n" ;
//datacard << "Btag       lnN	  1.034              -              1.040              -                -                -              1.027            1.017            1.034            1.034            1.026       \n" ;
//datacard << "Idiso      lnN	  1.022              -              1.022              -                -                -              1.030            1.030            1.022            1.022            1.021       \n" ;
//datacard << "Trigger    lnN	  1.000              -              1.000              -                -                -              1.001            1.001            1.000            1.000            1.001       \n" ;
//datacard << "QCDacc     lnN	  1.044              -              1.000              -                -                -              1.032            1.000            1.044            1.044            1.044       \n" ;
//datacard << "PDFacc     lnN	  1.000              -              1.000              -                -                -              1.001            1.000            1.000            1.000            1.003       \n" ;
//datacard << "tt_WW      lnN	    -              1.070              -              1.320              -                -                -                -                -                -                -         \n" ;
//datacard << "DY         lnN	    -                -                -                -              1.041              -                -                -                -                -                -         \n" ;
//datacard << "fakes      lnN	    -                -                -                -                -              1.300                -                -                -                -                -       \n" ; 
//datacard << "xsMC       lnN	    -               -              1.300              -                -                -              1.300            1.300            1.300            1.300            1.300        \n" ;
//datacard << "stat       lnN	  1.007              -              1.053              -                -                -              1.078            1.054            1.007            1.007            1.110       \n" ;
//datacard << "Lumi       lnN	  1.027              -              1.027              -                -                -              1.027            1.027            1.027            1.027            1.027       \n" ;


datacard << "LepMupT    lnN	  1.009              -              1.007              -                -                -              1.012            1.011            1.009            1.009            1.000       \n" ;       
datacard << "LepElepT   lnN	  1.003              -              1.003              -                -                -              1.001            1.005            1.003            1.003            1.000       \n" ;       
datacard << "MET        lnN	  1.053              -              1.037              -                -                -              1.040            1.027            1.053            1.053            1.018       \n" ;      
datacard << "JESMax     lnN	  1.019              -              1.022              -                -                -              1.009            1.035            1.019            1.019            1.074       \n" ;       
datacard << "Btag       lnN	  1.002              -              1.008              -                -                -              1.006            1.002            1.002            1.002            1.024       \n" ;       
datacard << "Idiso      lnN	  1.023              -              1.025              -                -                -              1.025            1.028            1.023            1.023            1.029       \n" ;       
datacard << "Trigger    lnN	  1.000              -              1.000              -                -                -              1.000            1.000            1.000            1.000            1.000       \n" ;       
datacard << "QCDacc     lnN	  1.006              -              1.000              -                -                -              1.024            1.000            1.006            1.006            1.006       \n" ;       
datacard << "PDFacc     lnN	  1.000              -              1.000              -                -                -              1.001            1.000            1.000            1.000            1.000       \n" ;       
datacard << "tt_WW      lnN	    -              1.070              -              1.320              -                -                -                -                -                -                -         \n" ;      
datacard << "DY         lnN	    -                -                -                -              1.041              -                -                -                -                -                -         \n" ;      
datacard << "fakes      lnN	    -                -                -                -                -              1.300                -                -                -                -                -       \n" ;       
datacard << "xsMC       lnN	    -               -              1.300              -                -                -              1.300            1.300            1.300            1.300            1.300        \n" ;       
datacard << "stat       lnN	  1.010              -              1.043              -                -                -              1.040            1.021            1.010            1.010            1.148       \n" ;       
datacard << "Lumi       lnN	  1.027              -              1.027              -                -                -              1.027            1.027            1.027            1.027            1.027       \n" ; 

		datacard.close();



	}  // i

	}  // k
}
