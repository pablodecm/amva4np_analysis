#!/usr/bin/env python 

import ROOT
from ROOT import ThinEvent, ExtEvent
from ROOT import TChain, FrankSelector,  vector
from ROOT import TH1, TFile
from ROOT import TEventList 

max_events = 100000000

TH1.AddDirectory(False)

base_dir = "../datasets/mixed/"
thinEventPath = "../datasets/unmixed/"
name = "QCD_pp_bbbb_13TeV"
mult = 5


tc_hm = TChain("hem_tree")
tc_hm.Add(thinEventPath+name+".root")
tchain = TChain("tree")
tchain.Add(thinEventPath+name+".root")

min_n_CSVM = 4
el = TEventList("el","el")
tchain.Draw(">>el","Sum$(pfjets.getDiscriminatorC(\"BTag\") > 0.5) > {}".format(min_n_CSVM-1))

tc_hm.SetEventList(el)
tchain.SetEventList(el)
selector = FrankSelector(ExtEvent(ThinEvent))(0, tc_hm, mult)
tchain.Process(selector, "ofile="+base_dir+name+".root", max_events)

