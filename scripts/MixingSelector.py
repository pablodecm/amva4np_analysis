#!/usr/bin/env python 

import ROOT
from ROOT import ExtEvent, DelphesEvent 
from ROOT import TChain, MixingSelector,  vector
from ROOT import TH1

from amva4np.samples.oxford_summer_2016 import delphes

max_events = -100000

TH1.AddDirectory(False)

sub_strs = ['QCD_pp_bbbb_13TeV','pp_hh_bbbb_13TeV']
sub_strs = ['QCD_pp_bbbb_13TeV']
sub_strs = ['pp_hh_bbbb_13TeV_new']
mc_names = delphes.keys()
mc_names=[n for n in mc_names if any(s in n for s in sub_strs)]

o_dir = "../datasets/unmixed_4btag_pt30/" 
p_par = "ofile={}.root;pName={}"

for name in mc_names:
    selector = MixingSelector(ExtEvent(DelphesEvent))(0)
    tchain = TChain("Delphes")
    for f in delphes[name]["files"]:
      tchain.Add(f)
    print "processing {} sample".format(name)
#    print "total number of events",tchain.GetEntries()
    if max_events > 0:
        tchain.Process(selector, p_par.format(o_dir+name, name), max_events)
    else:
        tchain.Process(selector, p_par.format(o_dir+name, name))

