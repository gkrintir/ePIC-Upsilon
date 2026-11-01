//
// ePIC Style, based on a style file from BaBar, v0.1
//

#include <iostream>
#include "ePICStyle.h"
#include "TROOT.h"

void SetePICStyle ()
{
  static TStyle* epicStyle = 0;
  std::cout << "ePICStyle: Applying nominal settings." << std::endl ;
  if ( epicStyle==0 ) epicStyle = ePICStyle();
  gROOT->SetStyle("ePIC");
  gROOT->ForceStyle();
}

TStyle* ePICStyle() 
{
  TStyle *ePICStyle = new TStyle("ePIC","ePIC style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  ePICStyle->SetFrameBorderMode(icol);
  ePICStyle->SetFrameFillColor(icol);
  ePICStyle->SetCanvasBorderMode(icol);
  ePICStyle->SetCanvasColor(icol);
  ePICStyle->SetPadBorderMode(icol);
  ePICStyle->SetPadColor(icol);
  ePICStyle->SetStatColor(icol);
  //ePICStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  ePICStyle->SetPaperSize(20,26);

  // set margin sizes
  ePICStyle->SetPadTopMargin(0.05);
  ePICStyle->SetPadRightMargin(0.18);
  ePICStyle->SetPadBottomMargin(0.16);
  ePICStyle->SetPadLeftMargin(0.15);

  // set title offsets (for axis label)
  ePICStyle->SetTitleXOffset(1.4);
  ePICStyle->SetTitleYOffset(1.2);
  ePICStyle->SetTitleOffset(1.2,"z"); // Set the offset for Z axis titles expliticly to avoid it being cut off

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  ePICStyle->SetTextFont(font);

  ePICStyle->SetTextSize(tsize);
  ePICStyle->SetLabelFont(font,"x");
  ePICStyle->SetTitleFont(font,"x");
  ePICStyle->SetLabelFont(font,"y");
  ePICStyle->SetTitleFont(font,"y");
  ePICStyle->SetLabelFont(font,"z");
  ePICStyle->SetTitleFont(font,"z");
  
  ePICStyle->SetLabelSize(tsize,"x");
  ePICStyle->SetTitleSize(tsize,"x");
  ePICStyle->SetLabelSize(tsize,"y");
  ePICStyle->SetTitleSize(tsize,"y");
  ePICStyle->SetLabelSize(tsize,"z");
  ePICStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  ePICStyle->SetMarkerStyle(20);
  ePICStyle->SetMarkerSize(1.2);
  ePICStyle->SetHistLineWidth(2.);
  ePICStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //ePICStyle->SetErrorX(0.001);
  // get rid of error bar caps
  ePICStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  ePICStyle->SetOptTitle(0);
  //ePICStyle->SetOptStat(1111);
  ePICStyle->SetOptStat(0);
  //ePICStyle->SetOptFit(1111);
  ePICStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  ePICStyle->SetPadTickX(1);
  ePICStyle->SetPadTickY(1);

  // legend modification
  ePICStyle->SetLegendBorderSize(0);
  ePICStyle->SetLegendFillColor(0);
  ePICStyle->SetLegendFont(font);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
  std::cout << "ePICStyle: ROOT6 mode" << std::endl;
  ePICStyle->SetLegendTextSize(tsize);
  ePICStyle->SetPalette(kBird);
#else
  std::cout << "ePICStyle: ROOT5 mode" << std::endl;
  // color palette - manually define 'kBird' palette only available in ROOT 6
  Int_t alpha = 0;
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
#endif

  ePICStyle->SetNumberContours(80);

  return ePICStyle;

}