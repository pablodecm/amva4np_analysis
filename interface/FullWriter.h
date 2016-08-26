
#pragma once

#include <algorithm>
#include <math.h>

#include "BaseOperator.h"

#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/MET.h"

template <class EventClass> class FullWriter : public BaseOperator<EventClass> {

  public:
 
    bool root_;
    std::string dir_;
    // variables to save in branches
    mut::EventInfo * eventInfo = nullptr;
    std::vector<mut::Jet> * pfjets = nullptr;
    std::vector<mut::Candidate> * dijets = nullptr;
    mut::MET * pfmet = nullptr;

    TTree tree_{"tree","Tree using simplified mut::dataformats"};

     FullWriter(bool root = false, std::string dir = "") :
      root_(root),
      dir_(dir) {}
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

      tree_.Branch("dijets","std::vector<mut::Candidate>",
                   &dijets, 64000, 1);

      tree_.SetDirectory(tdir);
      tree_.AutoSave();

   }


    virtual bool process( EventClass & ev ) {

      // set dijets
      ev.dijets_.clear();
      ev.dijets_.emplace_back(ev.jets_.at(0) + ev.jets_.at(1));
      ev.dijets_.emplace_back(ev.jets_.at(2) + ev.jets_.at(3));
      // redirect pointers
      eventInfo = dynamic_cast<mut::EventInfo *>(&ev.eventInfo_);
      pfjets = dynamic_cast<std::vector<mut::Jet> *>(&ev.jets_); 
      dijets = dynamic_cast<std::vector<mut::Candidate> *>(&ev.dijets_); 

      tree_.Fill();

      return true;
    }

    virtual bool output( TFile * tfile) {

      return true;

    }

};
