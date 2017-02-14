
#pragma once

#include <algorithm>
#include <math.h>

#include "BaseOperator.h"

#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/MET.h"
#include "mut_framework/mut_dataformats/interface/DiObject.h"

template <class EventClass> class FullWriter : public BaseOperator<EventClass> {

  public:
 
    bool root_;
    std::string dir_;
    // variables to save in branches
    mut::EventInfo * eventInfo = nullptr;
    std::vector<mut::Jet> * pfjets = nullptr;
    std::vector<mut::Candidate> * b_quarks = nullptr;
    std::vector<mut::MET> * mets = nullptr;
    std::vector<mut::Jet> pfjets_by_pt_;
    std::vector<mut::Jet> * pfjets_by_pt = nullptr;

    std::vector<mut::DiObject> * higgs = nullptr;
    std::vector<mut::DiObject> * dihiggs = nullptr;
    std::vector<mut::DiObject> * gen_higgs = nullptr;
    std::vector<mut::DiObject> * gen_dihiggs = nullptr;

    TTree tree_{"tree","Tree using simplified mut::dataformats"};

     FullWriter(bool root = false, std::string dir = "") :
      root_(root),
      dir_(dir),
      pfjets_by_pt(&pfjets_by_pt_) {}
    virtual ~FullWriter() {}

    virtual void init(TDirectory * tdir) {
      if (root_) {
        tdir = tdir->GetFile();
        auto ndir = tdir->mkdir(dir_.c_str());
        if (ndir == 0) {
          tdir = tdir->GetDirectory(dir_.c_str());
        } else {
          tdir = ndir;
        }
      }
      tree_.Branch("eventInfo","mut::EventInfo",
                   &eventInfo, 64000, 1);
      tree_.Branch("pfjets","std::vector<mut::Jet>",
                   &pfjets, 64000, 1);
      tree_.Branch("pfjets_by_pt","std::vector<mut::Jet>",
                   &pfjets_by_pt, 64000, 1);
      tree_.Branch("b_quarks","std::vector<mut::Candidate>",
                   &b_quarks, 64000, 1);
      tree_.Branch("mets","std::vector<mut::MET>",
                   &mets, 64000, 1);

      tree_.Branch("higgs","std::vector<mut::DiObject>",
                   &higgs, 64000, 1);
      tree_.Branch("dihiggs","std::vector<mut::DiObject>",
                   &dihiggs, 64000, 1);
      tree_.Branch("gen_higgs","std::vector<mut::DiObject>",
                   &gen_higgs, 64000, 1);
      tree_.Branch("gen_dihiggs","std::vector<mut::DiObject>",
                   &gen_dihiggs, 64000, 1);

      tree_.SetDirectory(tdir);
      tree_.AutoSave();

   }


    virtual bool process( EventClass & ev ) {

      // fill dijet objects just in case
      ev.higgs_.clear();
      ev.higgs_.emplace_back(ev.jets_.at(0),ev.jets_.at(1));
      ev.higgs_.emplace_back(ev.jets_.at(2),ev.jets_.at(3));
      ev.dihiggs_.clear();
      ev.dihiggs_.emplace_back(ev.higgs_.at(0), ev.higgs_.at(1));

      // copy jet collection to order by pt
      pfjets_by_pt_ = ev.jets_;
      auto pt_comparator = [&](const mut::Jet & a, const mut::Jet & b ){ return a.Pt() > b.Pt(); };
      std::sort(pfjets_by_pt_.begin(), pfjets_by_pt_.end(), pt_comparator);


      // redirect pointers
      eventInfo = dynamic_cast<mut::EventInfo *>(&ev.eventInfo_);
      pfjets = dynamic_cast<std::vector<mut::Jet> *>(&ev.jets_); 
      b_quarks = dynamic_cast<std::vector<mut::Candidate> *>(&ev.b_quarks_); 
      mets = dynamic_cast<std::vector<mut::MET> *>(&ev.mets_); 

      higgs = dynamic_cast<std::vector<mut::DiObject> *>(&ev.higgs_); 
      dihiggs = dynamic_cast<std::vector<mut::DiObject> *>(&ev.dihiggs_); 
      gen_higgs = dynamic_cast<std::vector<mut::DiObject> *>(&ev.gen_higgs_); 
      gen_dihiggs = dynamic_cast<std::vector<mut::DiObject> *>(&ev.gen_dihiggs_); 

      tree_.Fill();

      return true;
    }

    virtual bool output( TFile * tfile) {

      return true;

    }

};
