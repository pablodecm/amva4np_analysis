
#pragma once

#include <algorithm>

#include "Math/GenVector/VectorUtil.h"

#include "BaseOperator.h"
#include "mut_framework/mut_utils/interface/prettyprint.hpp"


template <class EventClass> class GenJetMatcher : public BaseOperator<EventClass> {

  public:

    double max_DeltaR_;

    GenJetMatcher( double max_DeltaR = 0.5 ) :
      max_DeltaR_(max_DeltaR) {}
    virtual ~GenJetMatcher() {}

    virtual bool process( EventClass & ev ) {

      const auto & reco_jets = ev.jets_;
      const auto & gen_quarks = ev.b_quarks_;
      auto & matchs = ev.reco_jet_matchs_;

      matchs.clear();

      // iterate over all reco jets
      for (std::size_t i_r=0; i_r<reco_jets.size(); i_r++) {
        const auto & reco_jet = reco_jets.at(i_r);
        matchs.emplace_back();
        for (std::size_t i_g=0; i_g<gen_quarks.size(); i_g++) { 
          const auto & gen_quark = gen_quarks.at(i_g);
          if ( ROOT::Math::VectorUtil::DeltaR( reco_jet, gen_quark) < max_DeltaR_) {
            matchs.at(i_r).insert(i_g);
          }
        }
      }

      return true;
    }

};

template <class EventClass> class AllGenBMatched : public BaseOperator<EventClass> {

  public:

    AllGenBMatched( ) {}

    virtual ~AllGenBMatched() {}

    virtual bool process( EventClass & ev) {

      const auto & matchs = ev.reco_jet_matchs_;

      std::set<std::size_t> gen_b_is_set{0,1,2,3}; 


      for (const auto & match : matchs) {
        if (match.size() == 1) {
          gen_b_is_set.erase(*match.begin());
        }
      }

      return gen_b_is_set.empty();
    }

    virtual std::string get_name() {
     return std::string{"all_gen_b_matched"};
    }

};

template <class EventClass> class HHJetsMatched : public BaseOperator<EventClass> {

  public:

    std::set<std::set<std::size_t>> higgs_jet_is_; 
    std::vector<std::vector<std::size_t>> dijets_is_; 

    HHJetsMatched( std::set<std::set<std::size_t>> higgs_jet_is ={{0,1},{2,3}},
                   std::vector<std::vector<std::size_t>> dijets_is = {{0,1},{2,3}}) :
      higgs_jet_is_(higgs_jet_is),
      dijets_is_(dijets_is) {}

    virtual ~HHJetsMatched() {}

    virtual bool process( EventClass & ev) {

      const auto & matchs = ev.reco_jet_matchs_;

      // copy of higgs indexes, will be empty if all matched 
      auto higgs_jet_is = higgs_jet_is_;


      for (const auto & pair_is : dijets_is_) {
        auto jet_match_is_0 = matchs.at(pair_is.at(0));
        auto jet_match_is_1 = matchs.at(pair_is.at(1));
        // check if both matched (1 to 1)
        if ( (jet_match_is_0.size() == 1) && (jet_match_is_1.size() == 1) ) {
          std::set<std::size_t> match_set;
          match_set.insert(*jet_match_is_0.begin());
          match_set.insert(*jet_match_is_1.begin());
          higgs_jet_is.erase(match_set);
        } else {
          // not matched because some jets not matched or 1 to many matching
          return false; 
        }
      }
      return higgs_jet_is.empty();
    }

    virtual std::string get_name() {
     return std::string{"both_dijets_matched"};
    }

};

