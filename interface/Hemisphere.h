
#pragma once

#define _USE_MATH_DEFINES

#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/Candidate.h"

#include "Math/GenVector/VectorUtil.h"

class Hemisphere { 
  public:

    mut::EventInfo eventInfo_;
    std::vector<mut::Jet> jets_;

    double p_phi_ = 0.;
    bool sumPz_inv_ = false;
    bool d_phi_inv_ = false;

    Hemisphere() {}
    Hemisphere(const mut::EventInfo & eventInfo, double p_phi, bool d_phi_inv ) : 
      eventInfo_(eventInfo),
      p_phi_(p_phi),
      d_phi_inv_(d_phi_inv) {}

    virtual ~Hemisphere() {}                   


    static double sumPz(const Hemisphere & hem) {
      double sumPz = 0.0;
      for (const auto & j : hem.jets_) {
        sumPz += j.Pz();
      }
      return sumPz;
    }

    static double thrustMayor(const Hemisphere & hem) {
      double thrustMayor = 0.0;
      for (const auto & j : hem.jets_) {
        double d_phi = ROOT::Math::VectorUtil::Phi_mpi_pi(j.Phi()-M_PI/2.);
        thrustMayor += std::abs(j.Pt()* std::cos(d_phi));
      }
      return thrustMayor;
    }

    static double thrustMinor(const Hemisphere & hem) {
      double thrustMinor = 0.0;
      for (const auto & j : hem.jets_) {
        double d_phi = ROOT::Math::VectorUtil::Phi_mpi_pi(j.Phi()-M_PI/2.);
        thrustMinor += std::abs(j.Pt()* std::sin(d_phi));
      }
      return thrustMinor;
    }

    static double invMass(const Hemisphere & hem) {
      mut::PtEtaPhiEVector sum_v; // init null?
      for (const auto & j : hem.jets_) {
        sum_v += j;
      }
      return sum_v.M();
    }

    static int nJets(const Hemisphere & hem) {
      return int(hem.jets_.size());
    }

    static int nTags(const Hemisphere & hem, std::string disc = "CSV",
                     float wp = 0.8) {
      int nTags = 0;
      for (const auto & j : hem.jets_) {
        if (j.getDiscriminator(disc) > wp) nTags++;
      }
      return nTags;
    }


    double getSumPz() {
      return sumPz(*this);
    }

    double getThrustMayor() {
      return thrustMayor(*this);
    }
   
    double getThrustMinor() {
      return thrustMinor(*this);
    }

    double getInvMass() {
      return invMass(*this);
    }

    int getNJets() {
      return nJets(*this);
    }

    int getNTags(std::string disc = "CSV", float wp = 0.8) {
      return nTags(*this,disc,wp);
    }

};

