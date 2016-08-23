
#pragma once

#include <algorithm>
#include <math.h>

#include "Math/GenVector/VectorUtil.h"

#include "BaseOperator.h"
#include "Hemisphere.h"

#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/Candidate.h"

template <class EventClass> class HemisphereProducer : public BaseOperator<EventClass> {

  public:

    HemisphereProducer() {}
    virtual ~HemisphereProducer() {}

    virtual bool process( EventClass & ev ) {

      // get thrust phi
      const auto & t_phi = ev.thrust_phi_;
      auto p_phi = t_phi-M_PI/2.;
     
      // clear old vector and add two hemispheres
      auto & hems = ev.hems_;
      hems.clear();
      hems.emplace_back(ev.eventInfo_, p_phi, false);
      hems.emplace_back(ev.eventInfo_, p_phi, true);

      // loop over all jets and push to hemispheres
      for (std::size_t i = 0; i < ev.jets_.size(); i++ ) {     
        const auto & jet = ev.jets_.at(i);
        auto d_phi = ROOT::Math::VectorUtil::Phi_mpi_pi(jet.Phi() - p_phi);
        // fill and rotate (swap phi if negative)
        if (d_phi < 0) {  
          hems.at(1).jets_.emplace_back(jet);   
          hems.at(1).jets_.back().SetPhi(-d_phi);
        } else {
          hems.at(0).jets_.emplace_back(jet);   
          hems.at(0).jets_.back().SetPhi(d_phi);
        }
      }
      
      // invert eta if sumPz < 0 (symmetry)
      for (auto & hem : hems) {
        if (hem.getSumPz() < 0) {
          hem.sumPz_inv_ = true; // so can be reverted 
          for (auto & j : hem.jets_) {
            j.SetEta(-j.Eta());
          }
        }
      }

      return true;
    }

};
