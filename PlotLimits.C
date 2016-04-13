#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TColor.h"  
#include "TLatex.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TInterpreter.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "THStack.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraphPainter.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "cosmetics.h"

void limits2( TString parameter ); 

const int n = 5;   //  number of DM-masses

const int M_star = 100;   // interaction scale

const float exponent = -0.167;   // = -1/6 (dependencia x-section vs interaction scale)

const double     x[n] = { 10, 50, 100, 200, 300 };  // DM-masses

Double_t yellow_low[n]; 
Double_t green_low[n] ; 
Double_t y_exp[n]     ; 
Double_t green_up[n]  ; 
Double_t yellow_up[n] ;
Double_t y_obs[n]     ;

//----------

void PlotLimits() {  

	gInterpreter->ExecuteMacro("GoodStyle.C");


	// ----------------------

	TFile* file = new TFile( "MVAAsymptotic.root", "READ" );

	TTree* tree = (TTree*) file -> Get( "limit" );

	Double_t limit; 

	tree -> SetBranchAddress( "limit", &limit );


	int j = 0;

	for ( long ievt = 0; ievt < tree -> GetEntries(); ievt++ ) { 

		tree -> GetEntry(ievt);

		if (ievt%6==0) yellow_low[j]   = limit; 
		if (ievt%6==1) green_low [j]   = limit;
	 	if (ievt%6==2) y_exp     [j]   = limit;
		if (ievt%6==3) green_up  [j]   = limit;
		if (ievt%6==4) yellow_up [j++] = limit;

	}



	limits2("mu");

	//limits2("M_star");

}



//----------

void limits2(TString parameter){

	TCanvas *c1 = new TCanvas("c", "carioca", 25, 25, 650, 550);

	c1 -> SetTicks(0, 0);

	if ( parameter == "mu"     ) c1 -> SetLogy();
	if ( parameter == "M_star" ) c1 -> SetLogx();

	double a[n], b[n], c[n], d[n], e[n];
	
	for ( int i = 0; i < n; i++ ) {

		a[i] = yellow_low[i];

		b[i] = green_low[i];

		c[i] = y_exp[i];

		d[i] = green_up[i];

		e[i] = yellow_up[i];

	}

	for ( int i = 0; i < n; i++ ) { 


		if ( parameter == "M_star" ) {

			a[i]  = pow( a[i], exponent ) * M_star;

			b[i]  = pow( b[i], exponent ) * M_star;

			c[i]  = pow( c[i], exponent ) * M_star;

			d[i]  = pow( d[i], exponent ) * M_star;

			e[i]  = pow( e[i], exponent ) * M_star;

		}

		a[i] = fabs( a[i] - c[i] );

		b[i] = fabs( b[i] - c[i] );

		d[i] = fabs( d[i] - c[i] );

		e[i] = fabs( e[i] - c[i] );

	}



	TGraph* observed = new TGraph( n, x, y_obs);
	TGraph* expected = new TGraph( n, x, c);
	TGraphAsymmErrors* green  = new TGraphAsymmErrors( n, x, c, 0, 0, b, d );
	TGraphAsymmErrors* yellow = new TGraphAsymmErrors( n, x, c, 0, 0, a, e );

	observed -> SetFillColor(1);	
	observed -> SetLineWidth(3);
	observed -> SetLineStyle(1);

	expected -> SetFillColor(1);	
	expected -> SetLineWidth(3);
	expected -> SetLineStyle(2);

	green    -> SetLineColor(3);
	green    -> SetFillColor(3);
	green    -> SetFillStyle(1001);   // solid
 
	yellow   -> SetLineColor(5); 
	yellow   -> SetFillColor(5); 
	yellow   -> SetFillStyle(1001);   // solid

	yellow   -> Draw("a3");
	green    -> Draw("3");
	expected -> Draw("l same");
	observed -> Draw("c same");		


	// Cosmetics
	//----------------------------------------------------------------------
	//yellow->SetMinimum(0);

	DrawTLatex(61, 0.190, 0.94, 0.055, 11, "CMS" );
	DrawTLatex(52, 0.290, 0.94, 0.030, 11, "Preliminary" );
	DrawTLatex(42, 0.940, 0.94, 0.040, 31, "1.371 fb^{-1} (13 TeV)" );
	if ( parameter == "mu"     ) DrawTLatex(42, 0.700, 0.88, 0.035, 31, "MVA analysis (m_{DM} = 1 GeV) " );
	if ( parameter == "M_star" ) DrawTLatex(42, 0.939, 0.88, 0.035, 31, "cut-and-count analysis, dilepton channel" );
	
	yellow->SetTitle("");

	if ( parameter == "M_star" ) {

		//yellow->GetYaxis()->SetRange(0, 170);
		
		//yellow->SetMaximum(170.0);

		//yellow->SetMinimum(0.0);

	}

	yellow->GetXaxis()->SetTitle("m_{S} [GeV]");
	yellow->GetXaxis()->SetTitleOffset(1.5);

	if ( parameter == "mu"     ) yellow->GetYaxis()->SetTitle("95% CLs upper limit on #sigma/#sigma_{th}");
	if ( parameter == "M_star" ) yellow->GetYaxis()->SetTitle("95% CLs lower limits on M_{*} [GeV]");

	yellow->GetYaxis()->CenterTitle();
	yellow->GetYaxis()->SetTitleOffset(1.35);


	if ( parameter == "mu"     ){
	DrawLegend(0.70, 0.35, (TH1F*)expected, " expected",            "l");
	DrawLegend(0.70, 0.30, (TH1F*)green,    " expected #pm1#sigma", "f");  
	DrawLegend(0.70, 0.25, (TH1F*)yellow,   " expected #pm2#sigma", "f");
	}

	if ( parameter == "M_star"     ){
	DrawLegend(0.22, 0.35, (TH1F*)expected, " expected",            "l");
	DrawLegend(0.22, 0.30, (TH1F*)green,    " expected #pm1#sigma", "f");  
	DrawLegend(0.22, 0.25, (TH1F*)yellow,   " expected #pm2#sigma", "f");
	}

	c1 -> GetFrame() -> DrawClone();
	c1 -> RedrawAxis();

	c1 -> Update();
   	c1 -> Modified();

	c1 -> SaveAs(   "limits_" + parameter + ".pdf"  );
	c1 -> SaveAs(   "limits_" + parameter + ".png"  );

}

