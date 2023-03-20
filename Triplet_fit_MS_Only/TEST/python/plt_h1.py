#!/usr/bin/env python

# Standard imports
import sys, os
import argparse

# Scientific imports
import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Typedefs
tH1M = ROOT.ROOT.RDF.TH1DModel
tH2M = ROOT.ROOT.RDF.TH2DModel

def _plt (hist, fout, args):

  if not isinstance(hist, list):
    hist = [hist]

  # Set some histo properties
  for h in hist:
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.8)
    h.SetFillStyle(3244)
    h.SetLineWidth(2)

  # Don't show canvas
  ROOT.gStyle.SetOptStat(0)

  # make canvas and pad
  tc = ROOT.TCanvas("canv", "", 500, 500)
  tp = ROOT.TPad("pad", "", 0, 0, 1, 1)

  # Init
  tp.Draw()
  tp.cd()
  tp.SetTicks(1, 1)

  # Draw histogram
  for i, h in enumerate(hist):

    if i == 0:
      h.Draw("PLC PFC PMC E2 P")
    else:
      h.Draw("PLC PFC PMC E2 P SAME")


  # Save
  tc.SaveAs(fout)




def _run ():

  #===============================
  #
  # Get commandline arguments
  #
  #===============================

  parser = argparse.ArgumentParser(description="Plotting")

  parser.add_argument \
  (
    "--fout",
    type=str,
    default="output.pdf",
    help="Name of output file"
  )
  parser.add_argument \
  (
    "--fin",
    type=str,
    help="Name of input file"
  )
  parser.add_argument \
  (
    "--tree",
    type=str,
    default="tree",
    help="Name of tree"
  )
  parser.add_argument \
  (
    "--config",
    type=str,
    default="Evt_ID,100,0,100",
    help="Name of aranch"
  )
  # Parse arguments
  args = parser.parse_args()

  #===============================
  #
  # Get content of file
  #
  #===============================

  # Go MT
  ROOT.ROOT.EnableImplicitMT()

  # Init RDF
  RDF = ROOT.RDataFrame(args.tree, args.fin)

  # Get the config of the histogram
  branch, nbin, xmin, xmax = args.config.split(",")

  # Container for histograms
  h1 = RDF.Histo1D \
  (
    tH1M(branch, "", int(nbin), float(xmin), float(xmax)),
    branch
  )

  # Trigger main event loop
  h1.Draw()

  # Start plotting
  _plt(h1, args.fout, args)
  

# Main caller
if __name__ == "__main__":

    _run()
