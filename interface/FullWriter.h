
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
//      tree_.Branch("pfmet","mut::MET",
//                   &pfmet, 64000, 1);

      tree_.SetDirectory(tdir);
      tree_.AutoSave();

   }


    virtual bool process( EventClass & ev ) {


      // to fill tree redirect pointers
      eventInfo = dynamic_cast<mut::EventInfo *>(&ev.eventInfo_);
      pfjets = dynamic_cast<std::vector<mut::Jet> *>(&ev.jets_); 
//      pfmet = dynamic_cast<mut::MET *>(&ev.met_);

      tree_.Fill();

      return true;
    }

    virtual bool output( TFile * tfile) {

      return true;

    }

};
