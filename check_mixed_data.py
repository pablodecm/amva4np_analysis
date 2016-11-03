#!/usr/bin/env python
# coding: utf-8

from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D, Canvas, Legend
from ROOT import TLatex
from itertools import combinations
from collections import OrderedDict
import math 
import argparse


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Check mixed data')
  parser.add_argument('file_path')
  parser.add_argument('--mult', type=int, default=1)
  args = parser.parse_args()
  check_mixed_data(args.file_path, args.mult)

def check_mixed_data(file_path, mult = 1):
  f = root_open(file_path)
  title = file_path 

  # memory resident for efficiency
#  f.tree.SetDirectory(0)
#  f.mix_tree.SetDirectory(0)
  
  event_v = ["Sum$(pfjets.E())",
             "sqrt(Sum$(pfjets.E())^2 - (Sum$(pfjets.Px())^2+Sum$(pfjets.Py())^2+Sum$(pfjets.Pz())^2))",
             "dijets[0].M()",
             "dijets[1].M()"]
  event_vname = ["HT",
                 "all_jets_invariant_mass",
                 "leading_dijet_inv_mass",
                 "trailing_dijet_inv_mass"]
  event_vbins = [[50,0.,1000.],
                 [50,0.,1000.],
                 [50,0.,500.],
                 [50,0.,500.]]
  jet_order = ["by_pt","by_pairing"]
  jet_order_name = ["pfjets_by_pt","pfjets"]
  
  bkg_cut = "eventInfo.getPNameC()==\"QCD_pp_bbbb_13TeV\""
  sig_cut = "eventInfo.getPNameC()==\"QCD_pp_bbbb_13TeV\""
  if mult==1:
    nn_cut = ""
  else:
    nn_cut = "((Entry$-0)%{}) == 0".format(mult*mult)
  
  joi = range(4) 
  jet_v = ["{b_name}[{i}].Pt()" ]
  jet_vnames = ["jet_{i}_pt_{order}" ]
  jet_vbins = [[50, 0,500]]
  jet_pair_v = ["abs(ROOT::Math::VectorUtil::Phi_mpi_pi({b_name}[{i0}].Phi()-{b_name}[{i1}].Phi()))",
                "abs({b_name}[{i0}].Eta()-{b_name}[{i1}].Eta())" ]
  
  jet_pair_vnames = ["delta_phi_jet_{i0}_jet_{i1}_{order}",
                     "delta_eta_jet_{i0}_jet_{i1}_{order}"]
  
  jet_pair_vbins = [[30,0.,math.pi],
                    [30,0.,5.]]
  
  draw_dict = OrderedDict() 
  
  for expr, vname, bins in zip(event_v, event_vname, event_vbins):
    draw_dict[vname] = { "expr" : expr,
                         "bins" : bins}
  
  
  
  for order, b_name in zip(jet_order,jet_order_name):
    for i in joi:
      for expr, vname, bins in zip(jet_v, jet_vnames, jet_vbins):
        opt_dict = {"b_name" : b_name, "i" : i, "order":order}
        draw_dict[vname.format(**opt_dict)] = { "expr" : expr.format(**opt_dict),
                                                "bins" : bins }
    for i0, i1 in combinations(joi,2):    
      for expr, vname, bins in zip(jet_pair_v, jet_pair_vnames, jet_pair_vbins):
        opt_dict = {"b_name" : b_name, "i0" : i0, "i1":i1, "order":order}
        draw_dict[vname.format(**opt_dict)] = { "expr" : expr.format(**opt_dict),
                                                "bins" : bins }
  
  out = root_open(file_path.replace(".root","_check.root"),"RECREATE")      
  out.cd()
  h_cols = [2, 1, 4]
  h_alpha = 0.75
  h_names = ["total", "bkg only", "mixing"]
  for k, v in draw_dict.items():
    print "Drawing - {} ".format(k)
    canv = Canvas(2000,1000, name=k+"_canv")
    h_tot = f.tree.draw(v["expr"], hist = Hist(*v["bins"]))
    h_bkg = f.tree.draw(v["expr"], bkg_cut, hist = Hist(*v["bins"]))
    h_mix = f.mix_tree.draw(v["expr"], nn_cut, hist = Hist(*v["bins"]))
    h_mix.Scale(h_bkg.Integral()/h_mix.Integral())
    chi2 = h_bkg.Chi2Test(h_mix, "UUNORM")
    print "chi^2 p={}".format(chi2)
    l = Legend(3)
    # set histogram aesthetics
    for h_i, h in enumerate([h_tot, h_bkg, h_mix]):
      h.SetStats(0)
      h.SetTitle(title)
      h.SetXTitle(k)
      h.SetYTitle("events")
      h.SetLineWidth(2)
      h.SetLineColorAlpha(h_cols[h_i],0.75)
      h.SetMarkerSize(0)
      h.Draw("E1X1same")
      l.AddEntry(h, h_names[h_i], style="LEP")
  
    latex = TLatex(0.25, 0.85,"chi2 prop {:.2e}".format(chi2))
    latex.SetNDC(1)
    latex.Draw("same")
    l.Draw("same")
    canv.Write()
  
  if mult==1:

    f_draw_dict = OrderedDict() 

    f_v = ["tree.dijets[0].M()-mix_tree.dijets[0].M():tree.dijets[1].M()-mix_tree.dijets[1].M()",
           "tree.dijets[0].M():tree.dijets[0].M()-mix_tree.dijets[0].M()",
           "tree.dijets[1].M():tree.dijets[1].M()-mix_tree.dijets[1].M()"]

    f_vxlabel = ["delta_m_true_mixed_dijet_0",
                 "m_dijet_0",
                 "m_dijet_1"]
 
    f_vylabel = ["delta_m_true_mixed_dijet_1",
                 "delta_m_true_mixed_dijet_0",
                 "delta_m_true_mixed_dijet_1"]
    f_vbins = [[20, -100., 100]*2,
               [50, 0., 500., 20, -100., 100],
               [50, 0., 500., 20, -100., 100]]

    cut_exprs = ["",
                 "@pfjets.size() == 4",
                 "@pfjets.size() == 5",
                 "@pfjets.size() > 5"]

    cut_names = ["",
                "cut_4_jets",
                "cut_5_jets",
                "cut__more_than_5_jets"]

    for cut_expr, cut_name in zip(cut_exprs, cut_names):
      for expr, xlabel, ylabel, bins in zip(f_v, f_vxlabel, f_vylabel, f_vbins):
        if f_vylabel == "":
          vname = xlabel
        else:
          vname = xlabel+"_vs_"+ylabel+cut_name
       
        f_draw_dict[vname] = { "expr" : expr,
                               "xlabel" : xlabel,
                               "ylabel" : ylabel,
                               "bins" : bins,
                               "cut_expr" : cut_expr}

      
    f_tree = f.tree
    f_tree.AddFriend(f.mix_tree)

    for k, v in f_draw_dict.items():
      if len(v["bins"]) == 3:
        h = Hist(*v["bins"], name=k)
      elif len(v["bins"]) == 6:  
#        h = Hist2D(*v["bins"], name=k)
        v["expr"] += ">>{}({})".format(k,",".join([str(b) for b in v["bins"]]))
        h = None
      else: 
        h = None
      print "Drawing - {} ".format(k)
      print "expr - {} ".format(v["expr"])
      print "cut - {} ".format(v["cut_expr"])
      h = f_tree.draw(v["expr"], v["cut_expr"], hist = h)
      h.SetTitle(title)
      h.xaxis.SetTitle(v["xlabel"])
      h.yaxis.SetTitle(v["ylabel"])
      h.Write()

  out.Close()  




      


    





