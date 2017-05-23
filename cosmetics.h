	//////////////////////////////
///////////  cosmetics by jonatan  //////////
	//////////////////////////////

const Float_t _yoffset  = 0.045;
const Float_t _yoffset2 = 0.080;

void     DrawTLatex(	Font_t      tfont,
		     	Double_t    x,
		     	Double_t    y,
		     	Double_t    tsize,
		     	Short_t     align,
		     	const char* text,
		     	Bool_t      setndc = true );

void     SetAxis(	TH1F*       h,
		     	TString     xtitle,
		     	TString     units,
		     	Int_t       precision );

TLegend* DrawTLegend(	Float_t     x1,
		     	Float_t     y1,
		     	TH1F*       hist,
		     	TString     label,
		     	TString     option,
		     	Float_t     tsize   = 0.022, //0.03
		     	Float_t     xoffset = 0.20,
		     	Float_t     yoffset = _yoffset );

TLegend* DrawTLegend2(	Float_t     x1,
		     	Float_t     y1,
		     	TGraph*     graph,
		     	TString     label,
		     	TString     option,
		     	Float_t     tsize   = 0.035,
		     	Float_t     xoffset = 0.20,
		     	Float_t     yoffset = _yoffset2 );

void 	 AxisFonts(	TAxis* 	    axis,
	       		TString     coordinate,
	       		TString     title );

void 	 TH2FAxisFonts( TH2F*       h,
		    	TString     coordinate,
		    	TString     title );

void 	 setupGreyScale();

void 	 setupColors();

void 	MoveOverflowBins(  TH1* h,
			   Double_t xmin,
			   Double_t xmax  );

void	DrawLegend(	Float_t x1,
			Float_t y1,
			TH1F*   hist,
			TString label,
			TString option	);

///////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////
/////   DrawLatex   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DrawTLatex( Font_t 	tfont, 
		 Double_t 	x, 
		 Double_t 	y, 
		 Double_t 	tsize, 
		 Short_t 	align, 
		 const char* 	text, 
		 Bool_t	 	setndc ) {
  	
	TLatex* tl = new TLatex(x, y, text);

  	tl->SetNDC      (setndc);
  	tl->SetTextAlign( align);
  	tl->SetTextFont ( tfont);
  	tl->SetTextSize ( tsize);
  	tl->Draw("same");

}/////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////
/////   SetAxis   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void SetAxis(TH1F* h, TString xtitle, TString units, Int_t precision) {

	  TAxis* xaxis = h->GetXaxis();
	  TAxis* yaxis = h->GetYaxis();
	  
	  xaxis->SetLabelFont  (  42);
	  xaxis->SetLabelOffset(0.01);
	  xaxis->SetLabelSize  (0.05);
	  xaxis->SetNdivisions ( 505);
	  xaxis->SetTitleFont  (  42);
	  xaxis->SetTitleOffset( 1.2);
	  xaxis->SetTitleSize  (0.05);

	  yaxis->SetLabelFont  (  42);
	  yaxis->SetLabelOffset(0.01);
	  yaxis->SetLabelSize  (0.05);
	  yaxis->SetNdivisions ( 505);
	  yaxis->SetTitleFont  (  42);
	  yaxis->SetTitleOffset( 1.7);
	  yaxis->SetTitleSize  (0.05);

	  xaxis->SetTitle(xtitle);

	  TString ytitle = Form("events / %s.%df%s", "%", precision, units.Data());

	  yaxis->SetTitle(Form(ytitle.Data(), h->GetBinWidth(0)));

	  if (units.Contains("NULL")) yaxis->SetTitle("events / bin");

}/////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
/////   DrawTLegend   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TLegend* DrawTLegend(	Float_t x1,
		     	Float_t y1,
		     	TH1F*   hist,
		     	TString label,
		     	TString option,
		     	Float_t tsize,
		     	Float_t xoffset,
		     	Float_t yoffset ) {

	TLegend* legend = new TLegend( x1, y1, x1 + xoffset, y1 + yoffset );
  
	legend->SetBorderSize(    0);
	legend->SetFillColor (    0);
	legend->SetTextAlign (   12);
	legend->SetTextFont  (   42);
	legend->SetTextSize  (tsize);

	legend->AddEntry(hist, label.Data(), option.Data());
	legend->Draw();

	return legend;
}///////////////////////////////////////////////////////////////////////////////////////////////////////////// 

////////////////////////////
/////   DrawTLegend2   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TLegend* DrawTLegend2(	Float_t x1,
		     	Float_t y1,
		     	TGraph* graph,
		     	TString label,
		     	TString option,
		     	Float_t tsize,
		     	Float_t xoffset,
		     	Float_t yoffset ) {

	TLegend* legend = new TLegend( x1, y1, x1 + xoffset, y1 + yoffset );
  
	legend->SetBorderSize(    0);
	legend->SetFillColor (    0);
	legend->SetTextAlign (   12);
	legend->SetTextFont  (   42);
	legend->SetTextSize  (tsize);

	legend->AddEntry(graph, label.Data(), option.Data());
	legend->Draw();

	return legend;
}///////////////////////////////////////////////////////////////////////////////////////////////////////////// 

/////////////////////////
/////   AxisFonts   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AxisFonts(TAxis*  axis,
	       TString coordinate,
	       TString title)
{
  axis->SetLabelFont  (   42);
  axis->SetLabelOffset(0.015);
  axis->SetLabelSize  (0.050);
  axis->SetNdivisions (  505);
  axis->SetTitleFont  (   42);
  axis->SetTitleOffset(  1.5);
  axis->SetTitleSize  (0.050);

  if (coordinate == "y") axis->SetTitleOffset(1.6);

  axis->SetTitle(title);
}///////////////////////////////////////////////////////////////////////////////////////////////////////////// 


/////////////////////////////
/////   TH2FAxisFonts   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TH2FAxisFonts(TH2F*   h,
		   TString coordinate,
		   TString title)
{
  TAxis* axis = NULL;

  if (coordinate.Contains("x")) axis = h->GetXaxis();
  if (coordinate.Contains("y")) axis = h->GetYaxis();

  AxisFonts(axis, coordinate, title);
}///////////////////////////////////////////////////////////////////////////////////////////////////////////// 


//////////////////////////////
/////   setupGreyScale   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupGreyScale() 
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red  [NRGBs] = {0.90, 0.65, 0.40, 0.15, 0.00};
  double green[NRGBs] = {0.90, 0.65, 0.40, 0.15, 0.00};
  double blue [NRGBs] = {0.90, 0.65, 0.40, 0.15, 0.00};
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}///////////////////////////////////////////////////////////////////////////////////////////////////////////// 


////////////////////////////
/////   setupColors   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupColors()
{
  const Int_t NRGBs = 5;
//const Int_t NCont = 255;
  const Int_t NCont = 25;
  
  Double_t stops[NRGBs] = {0.00, 0.0625, 0.25, 0.5625, 1.00};
  Double_t red  [NRGBs] = {0.00, 0.00,   0.87, 1.00,   0.51};
  Double_t green[NRGBs] = {0.00, 0.81,   1.00, 0.20,   0.00};
  Double_t blue [NRGBs] = {0.51, 1.00,   0.12, 0.00,   0.00};
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}///////////////////////////////////////////////////////////////////////////////////////////////////////////// 

////////////////////////////////
/////   MoveOverflowBins   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void 	MoveOverflowBins(  TH1* h,
			   Double_t xmin,
			   Double_t xmax  )
{
	UInt_t nbins = h->GetNbinsX();

	TAxis* axis = (TAxis*)h->GetXaxis();

	UInt_t firstBin = (xmin != -999) ? axis->FindBin(xmin) : 1;
	UInt_t lastBin = (xmax != -999) ? axis->FindBin(xmax) : nbins;

        lastBin = (lastBin > nbins) ? nbins : lastBin;

	Double_t firstVal = 0;
	Double_t firstErr = 0;

	Double_t lastVal = 0;
	Double_t lastErr = 0;

	for (UInt_t i=0; i<=nbins+1; i++) {

		if (i <= firstBin) {
			firstVal += h->GetBinContent(i);
			firstErr += (h->GetBinError(i)*h->GetBinError(i));
		}

		if (i >= lastBin) {
			lastVal += h->GetBinContent(i);
			lastErr += (h->GetBinError(i)*h->GetBinError(i));
		}

		if (i < firstBin || i > lastBin) {
			h->SetBinContent(i, 0);
			h->SetBinError (i, 0);
		}

	}

	firstErr = sqrt(firstErr);
	lastErr = sqrt(lastErr);

	h->SetBinContent(firstBin, firstVal);
	h->SetBinError (firstBin, firstErr);

	h->SetBinContent(lastBin, lastVal);
	h->SetBinError (lastBin, lastErr);
}///////////////////////////////////////////////////////////////////////////////////////////////////////////// 


//////////////////////////
/////   DrawLegend   /////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DrawLegend(Float_t x1,
		Float_t y1,
		TH1F*   hist,
		TString label,
		TString option)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + 0.18,
				y1 + 0.05);

  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (0.035);

  legend->AddEntry (hist, label.Data(), option.Data());
  legend->Draw();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
