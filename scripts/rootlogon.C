// To set this as default, you need a .rootrc file in your home directory,
// containing the following line:
// Rint.Logon: /full/path/to/rootlogon.C

#ifndef ROOTLOGON_C
#define ROOTLOGON_C

#include "TColor.h"
#include "TH1.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"

void rootlogon()
{
  printf("Welcome to the ROOT of all evils \n");

  // Defaults to classic style, but that's OK, we can fix it
  TStyle* novaStyle = new TStyle("novaStyle", "NOvA Style");

  // Centre title
  novaStyle->SetTitleAlign(22);
  novaStyle->SetTitleX(.5);
  novaStyle->SetTitleY(.95);
  novaStyle->SetTitleBorderSize(0);

  // No info box
  //novaStyle->SetOptStat(0);
  //novaStyle->SetOptStat(11);
  //novaStyle->SetStatX(0.98);
  //novaStyle->SetStatY(0.9);
  //novaStyle->SetStatW(0.26);
  //novaStyle->SetStatH(0.2);

  //set the background color to white
  novaStyle->SetFillColor(10);
  novaStyle->SetFrameFillColor(10);
  novaStyle->SetCanvasColor(10);
  novaStyle->SetPadColor(10);
  novaStyle->SetTitleFillColor(0);
  novaStyle->SetStatColor(10);

  // Don't put a colored frame around the plots
  novaStyle->SetFrameBorderMode(0);
  novaStyle->SetCanvasBorderMode(0);
  novaStyle->SetPadBorderMode(0);

  // Set the default line color for a fit function to be red
  novaStyle->SetFuncColor(kRed+1);

  // Marker settings
  //  novaStyle->SetMarkerStyle(kFullCircle);

  // No border on legends
  //novaStyle->SetLegendBorderSize(0);

  // Disabled for violating NOvA style guidelines
  // Scientific notation on axes
  //  TGaxis::SetMaxDigits(3);

  // Axis titles
  novaStyle->SetTitleSize(.055, "xyz");
  novaStyle->SetTitleOffset(.8, "xyz");
  // More space for y-axis to avoid clashing with big numbers
  novaStyle->SetTitleOffset(.9, "y");
  // This applies the same settings to the overall plot title
  novaStyle->SetTitleSize(.055, "");
  novaStyle->SetTitleOffset(.8, "");
  // Axis labels (numbering)
  novaStyle->SetLabelSize(.03, "xyz");
  novaStyle->SetLabelOffset(.005, "xyz");

  // Prevent ROOT from occasionally automatically zero-suppressing
  //novaStyle->SetHistMinimumZero();

  // Thicker lines
  novaStyle->SetHistLineWidth(2);
  novaStyle->SetFrameLineWidth(2);
  novaStyle->SetFuncWidth(2);

  // Set the number of tick marks to show
  novaStyle->SetNdivisions(506, "xyz");

  // Set the tick mark style
  novaStyle->SetPadTickX(1);
  novaStyle->SetPadTickY(1);

  // Fonts
  const int kNovaFont = 42;
  novaStyle->SetStatFont(kNovaFont);
  novaStyle->SetLabelFont(kNovaFont, "xyz");
  novaStyle->SetTitleFont(kNovaFont, "xyz");
  novaStyle->SetTitleFont(kNovaFont, ""); // Apply same setting to plot titles
  novaStyle->SetTextFont(kNovaFont);
  novaStyle->SetLegendFont(kNovaFont);

  // Get moodier colours for colz
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  novaStyle->SetNumberContours(NCont);

  gROOT->SetStyle("novaStyle");

  // Uncomment this line if you want to force all plots loaded from files
  // to use this same style
  gROOT->ForceStyle();
}

// Add a label in top left corner
// Especially useful for "Neutrino Beam" and "Antineutrino Beam" labels
void CornerLabel(std::string Str) {
  TLatex* CornLab = new TLatex(.1, .93, Str.c_str());
  CornLab->SetTextColor(kGray+1);
  CornLab->SetNDC();
  CornLab->SetTextSize (2/30.);
  CornLab->SetTextAlign(11);
  CornLab->Draw();
}

void CenterTitles(TH1* histo)
{
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->CenterTitle();
  histo->GetZaxis()->CenterTitle();  
}

#endif
