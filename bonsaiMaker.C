/******************************************************************************************************
  root -l -x -q add.C+ 
******************************************************************************************************/

#include "MVA.h"

enum {Loose, Tight};

struct Lepton {
	int index;
	int type; 
	int flavour;
  	TLorentzVector v;
};

Lepton  Lepton1;
Lepton  Lepton2;

const int ELECTRON_FLAVOUR = 11;
const int MUON_FLAVOUR     = 13;

const float ELECTRON_MASS =  0.000511;  // [GeV]
const float MUON_MASS     =  0.106;     // [GeV]

const float csvv2ivf_looseWP  = 0.605;

void add2(TString sample);

////////
/// main
////////

void bonsaiMaker() {

   	/*add2("ttDM1scalar20"      );
   	add2("ttDM1scalar50"      );
	add2("ttDM1scalar500"     );
	add2("ttDM10scalar10"     ); 
	add2("ttDM50scalar50"     );
	add2("ttDM50scalar200"    );
	add2("ttDM50scalar300"    );
	add2("ttDM50scalar300"    );
	add2("ttDM150scalar200"   );*/

	//add2("TTTo2L2Nu"           );
	//add2("DYJetsToLL_M-10to50"); 
	//add2("DYJetsToLL_M-50"    );
	//add2("ST_tW_top"          );
	//add2("ST_tW_antitop"      );
	//add2("WWTo2L2Nu"          );
	//add2("WZTo3LNu"           );  
	//add2("ZZTo4L"             ); 
	//add2("WJetsToLNu"         );

	add2("DoubleEG"           );
	add2("DoubleMuon"         );
	add2("MuonEG"             );
	add2("SingleElectron"     );
	add2("SingleMuon"         );
}
///////////////////////////


////////
/// add2
////////

void add2(TString sample){

	////////////////
	/// old branches
	////////////////

	// MC 74
	//TFile *f_old  = TFile::Open("/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox/21Oct_25ns_MC__l2selFix__hadd/latino_" + sample + ".root", "read");

	// MC 76
	//TFile *f_old  = TFile::Open("/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox_76x/22Jan_25ns_mAODv2_MC__l2loose__hadd__bSFL2pTEff__l2tight/latino_" + sample + ".root", "read");

	// 5OCT
	//TFile *f_old  = TFile::Open("/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox/21OctBis_Run2015D_05Oct2015_0553pb__l2sel__hadd/latino_Run2015D_05Oct2015_" + sample + ".root", "read");
	// PromptReco
	TFile *f_old  = TFile::Open("/gpfs/csic_projects/tier3data/LatinosSkims/RunII/cernbox/21OctBis_Run2015D_PromptReco_0716pb__l2sel__hadd/latino_Run2915D_PromptReco_" + sample + ".root", "read");

	TTree *tree_old = (TTree*) f_old -> Get( "latino" );

	float pt1, pt2, mll, MET, pfType1Metphi, njet, channel, dphill, baseW, puW, GEN_weight_SM;
	float dphilmet, dphilmet1, dphilmet2;
	float jetpt1, jeteta1, jetphi1, jetmass1, jetpt2, jeteta2, jetphi2, jetmass2;
	float mtw1, mtw2; 

	vector<float> *std_vector_lepton_pt              ; std_vector_lepton_pt               = 0;
	vector<float> *std_vector_lepton_eta             ; std_vector_lepton_eta              = 0;
	vector<float> *std_vector_lepton_phi             ; std_vector_lepton_phi              = 0;
	vector<float> *std_vector_lepton_flavour         ; std_vector_lepton_flavour          = 0;
        //vector<float> *std_vector_lepton_isMediumMuon    ; std_vector_lepton_isMediumMuon     = 0;
        //vector<float> *std_vector_lepton_BestTrackdxy    ; std_vector_lepton_BestTrackdxy     = 0;
        //vector<float> *std_vector_lepton_BestTrackdz     ; std_vector_lepton_BestTrackdz      = 0;
        //vector<float> *std_vector_lepton_eleIdTight      ; std_vector_lepton_eleIdTight       = 0;
        //vector<float> *std_vector_lepton_chargedHadronIso; std_vector_lepton_chargedHadronIso = 0;
        //vector<float> *std_vector_lepton_photonIso       ; std_vector_lepton_photonIso        = 0;
        //vector<float> *std_vector_lepton_neutralHadronIso; std_vector_lepton_neutralHadronIso = 0;
        //vector<float> *std_vector_lepton_sumPUPt         ; std_vector_lepton_sumPUPt          = 0;

	vector<float> *std_vector_jet_pt      ;  std_vector_jet_pt       = 0;
	vector<float> *std_vector_jet_eta     ;  std_vector_jet_eta      = 0;
	vector<float> *std_vector_jet_phi     ;  std_vector_jet_phi      = 0;
	vector<float> *std_vector_jet_csvv2ivf;  std_vector_jet_csvv2ivf = 0;

	tree_old -> SetBranchAddress( "pt1"          , &pt2           );
	tree_old -> SetBranchAddress( "pt2"          , &pt1           );
	tree_old -> SetBranchAddress( "mll"          , &mll           );
	tree_old -> SetBranchAddress( "pfType1Met"   , &MET           );
	tree_old -> SetBranchAddress( "pfType1Metphi", &pfType1Metphi );
	tree_old -> SetBranchAddress( "njet"         , &njet          );
	tree_old -> SetBranchAddress( "channel"      , &channel       );
	tree_old -> SetBranchAddress( "dphill"       , &dphill        );
	tree_old -> SetBranchAddress( "baseW"        , &baseW         );
	tree_old -> SetBranchAddress( "puW"          , &puW           );
	tree_old -> SetBranchAddress( "GEN_weight_SM", &GEN_weight_SM );

	tree_old -> SetBranchAddress( "dphilmet" , &dphilmet  );
	tree_old -> SetBranchAddress( "dphilmet1", &dphilmet1 );
	tree_old -> SetBranchAddress( "dphilmet2", &dphilmet2 );
	tree_old -> SetBranchAddress( "jetpt1"   , &jetpt1    );
	tree_old -> SetBranchAddress( "jeteta1"  , &jeteta1   );
	tree_old -> SetBranchAddress( "jetphi1"  , &jetphi1   );
	tree_old -> SetBranchAddress( "jetmass1" , &jetmass1  );
	tree_old -> SetBranchAddress( "jetpt2"   , &jetpt2    );
	tree_old -> SetBranchAddress( "jeteta2"  , &jeteta2   );
	tree_old -> SetBranchAddress( "jetphi2"  , &jetphi2   );
	tree_old -> SetBranchAddress( "jetmass2" , &jetmass2  );



	tree_old -> SetBranchAddress( "mtw1"     , &mtw1      );
	tree_old -> SetBranchAddress( "mtw2"     , &mtw2      );

	tree_old -> SetBranchAddress( "std_vector_lepton_pt"              , &std_vector_lepton_pt               );
	tree_old -> SetBranchAddress( "std_vector_lepton_eta"             , &std_vector_lepton_eta              );
	tree_old -> SetBranchAddress( "std_vector_lepton_phi"             , &std_vector_lepton_phi              );
	tree_old -> SetBranchAddress( "std_vector_lepton_flavour"         , &std_vector_lepton_flavour          );
	//tree_old -> SetBranchAddress( "std_vector_lepton_isMediumMuon"    , &std_vector_lepton_isMediumMuon     );
	//tree_old -> SetBranchAddress( "std_vector_lepton_BestTrackdxy"    , &std_vector_lepton_BestTrackdxy     );
	//tree_old -> SetBranchAddress( "std_vector_lepton_BestTrackdz"     , &std_vector_lepton_BestTrackdz      );
	//tree_old -> SetBranchAddress( "std_vector_lepton_eleIdTight"      , &std_vector_lepton_eleIdTight       );
	//tree_old -> SetBranchAddress( "std_vector_lepton_chargedHadronIso", &std_vector_lepton_chargedHadronIso );
	//tree_old -> SetBranchAddress( "std_vector_lepton_photonIso"       , &std_vector_lepton_photonIso        );
	//tree_old -> SetBranchAddress( "std_vector_lepton_neutralHadronIso", &std_vector_lepton_neutralHadronIso );
	//tree_old -> SetBranchAddress( "std_vector_lepton_sumPUPt"         , &std_vector_lepton_sumPUPt          );

	tree_old -> SetBranchAddress( "std_vector_jet_pt"      , &std_vector_jet_pt       );
	tree_old -> SetBranchAddress( "std_vector_jet_eta"     , &std_vector_jet_eta      );
	tree_old -> SetBranchAddress( "std_vector_jet_phi"     , &std_vector_jet_phi      );
	tree_old -> SetBranchAddress( "std_vector_jet_csvv2ivf", &std_vector_jet_csvv2ivf );

	////////////////
	/// new branches
	////////////////

	TFile *f_new =  TFile::Open ( "minitrees/bonsai_" + sample + ".root", "UPDATE" );

	TTree *tree_new = new TTree( "latino", "bonsai tree" );

	Float_t nbjet15, MET_new, METPhi, lep1Pt, lep1Eta, lep1Phi, lep2Pt, lep2Eta, lep2Phi, mll_new, channel_new, njet_new, baseW_new, puW_new, GEN_weight_SM_new;
	Float_t dphilep1lep2, dphilep1jet1, dphilep1jet2, dphilep1MET, dphilep2jet1, dphilep2jet2, dphilep2MET, dphijet1jet2, dphijet1MET, dphijet2MET;
	Float_t dphilep1lep2MET, dphijet1jet2MET; 
	Float_t jetpt1_new, jeteta1_new, jetphi1_new, jetmass1_new, jetpt2_new, jeteta2_new, jetphi2_new, jetmass2_new;
	Float_t mtw1_new, mtw2_new; 

	tree_new -> Branch( "baseW"        , &baseW_new        , "baseW/F"        );
        tree_new -> Branch( "puW"          , &puW_new          , "puW/F"          );
	tree_new -> Branch( "GEN_weight_SM", &GEN_weight_SM_new, "GEN_weight_SM/F");

	tree_new -> Branch( "channel"      , &channel_new      , "channel/F"      );
	tree_new -> Branch( "njet"         , &njet_new         , "njet/F"         );
	tree_new -> Branch( "nbjet15"      , &nbjet15          , "nbjet15/F"      );

	tree_new -> Branch( "mll"          , &mll_new          , "mll/F"          );
	tree_new -> Branch( "mtw1"         , &mtw1_new         , "mtw1/F"         );
	tree_new -> Branch( "mtw2"         , &mtw2_new         , "mtw2/F"         );

	tree_new -> Branch( "MET"          , &MET_new          , "MET/F"          );
	tree_new -> Branch( "METPhi"       , &METPhi           , "METPhi/F"       );	

	tree_new -> Branch( "lep1Pt"       , &lep1Pt           , "lep1Pt/F"       );
	tree_new -> Branch( "lep1Eta"      , &lep1Eta          , "lep1Eta/F"      );
	tree_new -> Branch( "lep1Phi"      , &lep1Phi          , "lep1Phi/F"      );

	tree_new -> Branch( "lep2Pt"       , &lep2Pt           , "lep2Pt/F"       );
	tree_new -> Branch( "lep2Eta"      , &lep2Eta          , "lep2Eta/F"      );
	tree_new -> Branch( "lep2Phi"      , &lep2Phi          , "lep2Phi/F"      );

        tree_new -> Branch( "jetpt1"       , &jetpt1_new       , "jetpt1/F"       );	
        tree_new -> Branch( "jeteta1"      , &jeteta1_new      , "jeteta1/F"      );	
        tree_new -> Branch( "jetphi1"      , &jetphi1_new      , "jetphi1/F"      );	
        tree_new -> Branch( "jetmass1"     , &jetmass1_new     , "jetmass1/F"     );
	
        tree_new -> Branch( "jetpt2"       , &jetpt2_new       , "jetpt2/F"       );	
        tree_new -> Branch( "jeteta2"      , &jeteta2_new      , "jeteta2/F"      );	
        tree_new -> Branch( "jetphi2"      , &jetphi2_new      , "jetphi2/F"      );	
        tree_new -> Branch( "jetmass2"     , &jetmass2_new     , "jetmass2/F"     );

	tree_new -> Branch( "dphilep1lep2", &dphilep1lep2, "dphilep1lep2/F" );
	tree_new -> Branch( "dphilep1jet1", &dphilep1jet1, "dphilep1jet1/F" );
	tree_new -> Branch( "dphilep1jet2", &dphilep1jet2, "dphilep1jet2/F" );
	tree_new -> Branch( "dphilep1MET" , &dphilep1MET , "dphilep1MET/F"  );
	tree_new -> Branch( "dphilep2jet1", &dphilep2jet1, "dphilep2jet1/F" );
	tree_new -> Branch( "dphilep2jet2", &dphilep2jet2, "dphilep2jet2/F" );
	tree_new -> Branch( "dphilep2MET" , &dphilep2MET , "dphilep2MET/F"  );
	tree_new -> Branch( "dphijet1jet2", &dphijet1jet2, "dphijet1jet2/F" );
	tree_new -> Branch( "dphijet1MET" , &dphijet1MET , "dphijet1MET/F"  );
	tree_new -> Branch( "dphijet2MET" , &dphijet2MET , "dphijet2MET/F"  );
	tree_new -> Branch( "dphilep1lep2MET",&dphilep1lep2MET,"dphilep1lep2MET/F" );
	//tree_new -> Branch( "dphijet1jet2MET",&dphijet1jet2MET,"dphijetljet2MET/F" );

	//tree_new -> Branch( "", &_new, "/F" );

	Long64_t nentries = tree_old -> GetEntries();



	cout << "------------------------------------" << endl;
	cout << "  >> accesing to sample:  " << sample << endl;

	for (Long64_t ievt = 0; ievt < nentries; ievt++){
	//for (Long64_t ievt = 0; ievt < 10000; ievt++){		

		if ( ievt % 10000 == 0 ) cout << "    >> processing event " << ievt << " ..." << endl;

		tree_old -> GetEntry(ievt);


		// preselection (e.mail 28i16)
		// ----------------------------------------------------------------------------------------------------------------------------------------
		if ( pt1 < 30 )                                                          continue;       //    lepton1pt > 30 GeV
		if ( pt2 < 10 )                                                          continue;       //    lepton2pt > 10 GeV
		                                                                                         //    exactly 2 leptons
		if ( mll < 20 )                                                          continue;       //    mll > 20 GeV
		if (  ( channel == 0 || channel == 1 )  &&  ( mll > 76 && mll < 106 )  ) continue;       //    | mll - mZ | > 15 GeV for ee and mm channels
		if ( njet <= 0 )                                                         continue;       //    njet > 0
		if ( MET < 50 )                                                          continue;       //    MET > 80 GeV
		// ----------------------------------------------------------------------------------------------------------------------------------------



		// Lepton1 & Lepton2
		// -----------------
		for ( int i = 0; i < 2; i++ ) {

			float pt      = std_vector_lepton_pt      -> at(i);
			float eta     = std_vector_lepton_eta     -> at(i);
			float phi     = std_vector_lepton_phi     -> at(i);
			float flavour = std_vector_lepton_flavour -> at(i);

			Lepton lep;

			lep.index   = i;
			lep.type    = Tight;
    			lep.flavour = flavour;


			float mass = -999;

			if      ( abs(lep.flavour) == ELECTRON_FLAVOUR )  mass = ELECTRON_MASS;
			else if ( abs(lep.flavour) == MUON_FLAVOUR     )  mass = MUON_MASS    ;


			TLorentzVector tlv;

			tlv.SetPtEtaPhiM( pt, eta, phi, mass );

			lep.v = tlv;

			if (i == 0) Lepton1 = lep;
    			if (i == 1) Lepton2 = lep;


		}


		// nbjet15
   		// -------
		nbjet15 = 0;

		//cout << "" << endl;
		//cout << "ievt = " << ievt << endl;
		//cout << "nbjet15 (before) = " << nbjet15 << endl; 

  		int vector_jet_size = std_vector_jet_pt -> size();

	  	for ( int i = 0; i < vector_jet_size; i++ ) {

			float pt  = std_vector_jet_pt  -> at(i);
			float eta = std_vector_jet_eta -> at(i);
			float phi = std_vector_jet_phi -> at(i);

			bool pass = ( pt > 15. && fabs(eta) < 4.7 );

			if ( !pass ) continue;

			TLorentzVector tlv;

			tlv.SetPtEtaPhiM( pt, eta, phi, 0.0 );

			bool is_lepton = false;
			
			if (  tlv.DeltaR( Lepton1.v ) < 0.3 || tlv.DeltaR( Lepton2.v ) < 0.3  ) is_lepton = true;
     			
    			if (is_lepton) continue;

    			if ( std_vector_jet_csvv2ivf -> at(i) > csvv2ivf_looseWP ) nbjet15++;

  		}


		// dphillmet
		// ---------
		TVector2 ET, lep1, lep2, jet1, jet2, lep1lep2, jet1jet2; 

		      ET.SetMagPhi( 1.0, pfType1Metphi                 );
		    lep1.SetMagPhi( 1.0, Lepton1.v.Phi() );
		    lep2.SetMagPhi( 1.0, Lepton2.v.Phi() );	
		    jet1.SetMagPhi( 1.0, jetphi1 );	
		    jet2.SetMagPhi( 1.0, jetphi2 );	
		lep1lep2.SetMagPhi( 1.0, (Lepton1.v + Lepton2.v).Phi() );	
		//jet1jet2.SetMagPhi( 1.0, ().Phi() );	

		MET_new    = MET            ;   
		METPhi     = pfType1Metphi  ;
		lep1Pt     = Lepton1.v.Pt() ;
		lep1Eta    = Lepton1.v.Eta();
		lep1Phi    = Lepton1.v.Phi();
		lep2Pt     = Lepton2.v.Pt() ;
		lep2Eta    = Lepton2.v.Eta();
		lep2Phi    = Lepton2.v.Phi();
		mll_new    = mll            ;     
		channel_new= channel        ;
		njet_new   = njet           ;

		baseW_new  = baseW          ;
		puW_new    = puW            ;
		GEN_weight_SM_new = GEN_weight_SM / abs(GEN_weight_SM); 

		jetpt1_new    = jetpt1      ; 
		jeteta1_new   = jeteta1     ; 
		jetphi1_new   = jetphi1     ; 
		jetmass1_new  = jetmass1    ; 
		jetpt2_new    = jetpt2      ; 
		jeteta2_new   = jeteta2     ; 
		jeteta2_new   = jeteta2     ; 
		jetphi2_new   = jetphi2     ; 
		jetmass2_new  = jetmass2    ; 

		mtw1_new = mtw1; 
		mtw2_new = mtw2; 

		dphilep1lep2 = fabs(  lep1.DeltaPhi(lep2)  );
		dphilep1jet1 = fabs(  lep1.DeltaPhi(jet1)  );				
		dphilep1jet2 = fabs(  lep1.DeltaPhi(jet2)  );				
		dphilep1MET  = fabs(  lep1.DeltaPhi(ET )  );
		dphilep2jet1 = fabs(  lep2.DeltaPhi(jet1)  );				
		dphilep2jet2 = fabs(  lep2.DeltaPhi(jet2)  );				
		dphilep2MET  = fabs(  lep2.DeltaPhi(ET )  );
		dphijet1jet2 = fabs(  jet1.DeltaPhi(jet2)  );				
		dphijet1MET  = fabs(  jet1.DeltaPhi(ET )  );
		dphijet2MET  = fabs(  jet2.DeltaPhi(ET )  );;				
		dphilep1lep2MET  = fabs(  lep1lep2.DeltaPhi(ET )  );
		//dphijet1jet2MET  = fabs(  jet1jet2.DeltaPhi(ET )  );

		tree_new -> Fill();


	}
	
	f_new -> Write();

	f_old  -> Close();

	f_new -> Close(); 
	
	cout << "  >> " <<sample << " bonsai created " << endl;
	cout << "------------------------------------" << endl;

}






