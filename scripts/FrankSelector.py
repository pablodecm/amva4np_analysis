#!/usr/bin/env python 

import ROOT
from ROOT import ThinEvent, ExtEvent
from ROOT import TChain, FrankSelector,  vector
from ROOT import TH1, TFile
from ROOT import TEventList 

max_events = 100000000

TH1.AddDirectory(False)

base_dir = "../datasets/mixed/"
thinEventPath = "../datasets/unmixed_4btag/"
names = ["QCD_pp_bbbb_13TeV","pp_hh_bbbb_13TeV"]
fraction = [0.874,5.668e-6]
times_sm = [1., 1000.] 
mult = 5


tc_hm = TChain("hem_tree")
tchain = TChain("tree")
n_evs = [0]
for name in names:
  tc_hm.Add(thinEventPath+name+".root")
  tchain.Add(thinEventPath+name+".root")
  n_evs.append(tchain.GetEntries())

min_n_CSVM = 4
#el = TEventList("el","el")
#tchain.Draw(">>el","Sum$(pfjets.getDiscriminatorC(\"BTag\") > 0.5) > {}".format(min_n_CSVM-1))
els = [TEventList("el"+name,"el"+name) for name in names]
print els
for s, name  in enumerate(names):
  tot_frac = fraction[s]*times_sm[s]
  n_ev = long((n_evs[s+1]-n_evs[s])*tot_frac) 
  print "# {} events to use: {} ".format(name, n_ev)
  print tchain.Draw(">>el{}".format(name),"(Entry$ > {}) && ( Entry$ <= {} )".format(n_evs[s], n_evs[s]+n_ev))

el = TEventList("el","el")
for list_to_add in els:
  el.Add(list_to_add)
print "total number of entries in list" + str(el.GetN()) 
tc_hm.SetEventList(el)
tchain.SetEventList(el)
selector = FrankSelector(ExtEvent(ThinEvent))(0, tc_hm, mult)
print "selector created"
tchain.Process(selector, "ofile="+base_dir+name+"_times_sm_"+str(int(times_sm[1]))+".root", max_events)

