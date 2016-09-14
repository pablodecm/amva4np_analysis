#!/usr/bin/env python 
# coding: utf-8

import argparse

default_nn_vars = ["thrustMayor", "thrustMinor", "sumPz", "invMass"]

parser = argparse.ArgumentParser(description='Create mixed dataset.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--nn_vars', nargs='+', default=default_nn_vars,
                    help="list of vars for metric in NN search")
parser.add_argument('--lumi_factor', type=float, default=1.0,
                    help="percentage of 5fb-1 to use")
parser.add_argument('--hh_times_sm', type=float, default=150.,
                    help="times the expected SM hh contribution")
parser.add_argument('--mult', type=int, default=1,
                    help="number of mixing combinations")
parser.add_argument('--extra_cut', default="",
                    help="TTree::Draw like condition (e.g. && (@pfjets.size() > 5))" )
parser.add_argument('--extra_cut_name', default="",
                    help="label for the extra cut")
args = parser.parse_args()


import ROOT
from ROOT import ThinEvent, ExtEvent
from ROOT import TChain, FrankSelector, vector
from ROOT import TH1, TFile
from ROOT import TEventList 
from check_mixed_data import check_mixed_data

ROOT.gROOT.SetBatch(True)

# input vector for selector
nn_vars = vector("string")()
for nn_var in args.nn_vars: nn_vars.push_back(nn_var)



max_events = 100000000

TH1.AddDirectory(False)

base_dir = "../datasets/mixed/"
thinEventPath = "../datasets/unmixed_4btag_pt30/"
names = ["QCD_pp_bbbb_13TeV","pp_hh_bbbb_13TeV_new"]
fraction = [0.874,5.668e-6]
lumi_frac = 5.0
lumi_factor = args.lumi_factor 
lumi = lumi_frac*lumi_factor
times_sm = [1., args.hh_times_sm] 
mult = args.mult 
extra_cut = args.extra_cut
extra_cut_name = args.extra_cut_name 



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
for s, name  in enumerate(names):
  tot_frac = fraction[s]*times_sm[s]*lumi_factor
  n_ev = long((n_evs[s+1]-n_evs[s])*tot_frac) 
  print "# {} events to use: {} ".format(name, n_ev)
  print tchain.Draw(">>el{}".format(name),"(Entry$ > {}) && ( Entry$ <= {} ) {}".format(n_evs[s], n_evs[s]+n_ev, extra_cut))

el = TEventList("el","el")
for list_to_add in els:
  el.Add(list_to_add)
print "total number of entries in list" + str(el.GetN()) 
tc_hm.SetEventList(el)
tchain.SetEventList(el)
selector = FrankSelector(ExtEvent(ThinEvent))(0, tc_hm, mult, nn_vars)
print "selector created"
ofile = base_dir
ofile+= "mixed_{}_invfb_{}_times_sm_{}_mult_{}_nn_vars{}.root".format(lumi, times_sm[1], mult,
                                                                      "_".join(args.nn_vars),extra_cut_name)
tchain.Process(selector, "ofile="+ofile, max_events)


check_mixed_data(ofile, mult)

