
#pragma once

#include <algorithm>
#include <math.h>

#include "BaseOperator.h"
#include "BTagJetSelection.h"
#include "DiJetPairSelection.h"

#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/MET.h"


template <class EventClass> class FrankWriter : public BaseOperator<EventClass> {

  public:
 
    // variables to save in branches
    mut::EventInfo * eventInfo = nullptr;
    std::vector<mut::Jet> mix_jets_ = {};
    std::vector<mut::Jet> mix_jets_by_pt_ = {};
    std::vector<mut::Candidate> mix_dijets_ = {};
    // a pointer to avoid undecalred label
    std::vector<mut::Jet> * mix_jets_ptr_ = nullptr;
    std::vector<mut::Jet> * mix_jets_by_pt_ptr_ = nullptr;
    std::vector<mut::Candidate> * mix_dijets_ptr_ = nullptr;

    // hemisphere combinations to save
    std::size_t n_h_mix_;
    std::size_t n_h_skip_;

    // to order jets after mixing
    std::string disc_ = "BTag";
    double d_value_ = 0.5;
    std::size_t n_min_disc_ = 4;

    bool root_;
    std::string dir_;


    TTree tree_{"mix_tree","Tree wth mixed events"};

     FrankWriter(std::size_t n_h_mix = 1, std::size_t n_h_skip = 1,
                bool root = false, std::string dir = "") :
      mix_jets_ptr_(&mix_jets_),  
      mix_jets_by_pt_ptr_(&mix_jets_by_pt_),  
      mix_dijets_ptr_(&mix_dijets_),  
      n_h_mix_(n_h_mix),
      n_h_skip_(n_h_skip), 
      root_(root),
      dir_(dir) {}
    virtual ~FrankWriter() {}

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
                   &mix_jets_ptr_, 64000, 1);
      tree_.Branch("pfjets_by_pt","std::vector<mut::Jet>",
                   &mix_jets_by_pt_ptr_, 64000, 1);
      tree_.Branch("dijets","std::vector<mut::Candidate>",
                   &mix_dijets_ptr_, 64000, 1);


      tree_.SetDirectory(tdir);
      tree_.AutoSave();

   }


    virtual bool process( EventClass & ev ) {


      // most stuff from original event
      eventInfo = dynamic_cast<mut::EventInfo *>(&ev.eventInfo_);

      const auto & bm_hems = ev.best_match_hems_;


      // for each hemisphere i
      for (std::size_t h_i=n_h_skip_; h_i<(n_h_skip_+n_h_mix_); h_i++) {
        // for each hemisphere j
        for (std::size_t h_j=n_h_skip_; h_j<(n_h_skip_+n_h_mix_); h_j++) {

          // rename pName to denote # nn
          ev.eventInfo_.setPName("nn_"+std::to_string(h_i)+"_"+
                                 std::to_string(h_j)+"_"+
                                 ev.eventInfo_.getPName());

          // clear jet collection
          mix_jets_.clear();

          // references for easy access
          const auto jets_i = bm_hems.at(0).at(h_i).jets_;
          const auto jets_j = bm_hems.at(1).at(h_j).jets_;
          mix_jets_.insert(mix_jets_.end(), jets_i.begin(), jets_i.end());
          mix_jets_.insert(mix_jets_.end(), jets_j.begin(), jets_j.end());

          // order by pt for consitency
          auto pt_comparator = [&](const mut::Jet & a, const mut::Jet & b ){ return a.Pt() > b.Pt(); };
          std::sort(mix_jets_.begin(), mix_jets_.end(), pt_comparator);
          // copy collection to keeep order by pt (other would be paired)
          mix_jets_by_pt_ = mix_jets_;

          // order by b-tagging disc
          order_jets_by_disc(mix_jets_, disc_);
          // count n jets which pass discriminator
          std::size_t n_pass_disc = std::count_if(mix_jets_.begin(),
                                                  mix_jets_.end(),
                                                  [&] (const mut::Jet & jet)
    	      {return (jet.getDiscriminator(disc_) > d_value_);});
          // event discarded if not enough tagged jets
          if ( n_pass_disc < n_min_disc_) return false;

          // pair with min mass
          bool pt_order = true; // order by higgs pt
          dijet_pairing_better(mix_jets_, disc_, d_value_, n_min_disc_, n_pass_disc, pt_order);

          // set dijets
          mix_dijets_.clear();
          mix_dijets_.emplace_back(mix_jets_.at(0) + mix_jets_.at(1));
          mix_dijets_.emplace_back(mix_jets_.at(2) + mix_jets_.at(3));
          //  fill tree with combination
          tree_.Fill();
        }
      }



      return true;
    }

    virtual bool output( TFile * tfile) {

      return true;

    }

};
