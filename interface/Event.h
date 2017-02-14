
#pragma once

#include <TTreeReader.h>

#include "JetCollection.h"
#include "EventInfo.h"
#include "Hemisphere.h"
#include "GenParticleCollection.h"
#include "METCollection.h"

#include "mut_framework/mut_dataformats/interface/Reader.h"
#include "mut_framework/mut_dataformats/interface/MET.h"
#include "mut_framework/mut_dataformats/interface/DiObject.h"

class DelphesEvent { 
  public:

    // read from TTree
    EventInfo eventInfo_;
    JetCollection jets_;
    GenParticleCollection b_quarks_;
    METCollection mets_;

    DelphesEvent() {}

    DelphesEvent(TTreeReader & reader) :
     eventInfo_(reader),  
   	 jets_(reader),
     b_quarks_(reader, "Particle", true),
     mets_(reader)
     {}                   

    virtual ~DelphesEvent() {};

    virtual void update() {
      eventInfo_.update();
      jets_.update();
      b_quarks_.update();
      mets_.update();
    }

};

class ThinEvent { 
  public:

    // read from TTree
    mut::Reader<mut::EventInfo> eventInfo_;
    mut::Reader<std::vector<mut::Jet>> jets_;
    mut::Reader<std::vector<mut::Candidate>> b_quarks_;
    mut::Reader<std::vector<mut::MET>> mets_;

    ThinEvent() {}

    ThinEvent(TTreeReader & reader) :
     eventInfo_(reader, "eventInfo" ),  
   	 jets_(reader, "pfjets"),
     b_quarks_(reader, "b_quarks"),
     mets_(reader, "mets")
     {}                   

    virtual ~ThinEvent() {};

    virtual void update() {
      eventInfo_.update();
      jets_.update();
      b_quarks_.update();
      mets_.update();
    }
};

typedef std::vector<mut::Candidate> CandidateCollection;

template <class EventBase> class ExtEvent : public EventBase {
  public:

  std::vector<std::set<std::size_t>> reco_jet_matchs_; 
  // indexes of jet chosen by min mass diff 
  std::vector<std::size_t> free_is_;
  // tranverse thrust phi
  double thrust_phi_ = -1.;
  // hemispheres (rotated and pz positive)
  std::vector<Hemisphere> hems_; 
  // best matching hemispheres ( [first/second][proximity])
  std::vector<std::vector<Hemisphere>> best_match_hems_; 

  std::vector<mut::DiObject> higgs_;
  std::vector<mut::DiObject> dihiggs_;
  std::vector<mut::DiObject> gen_higgs_;
  std::vector<mut::DiObject> gen_dihiggs_;

  // inherit constructors
  using EventBase::EventBase;
  virtual ~ExtEvent() {}

  virtual void update() {
    EventBase::update();
    reco_jet_matchs_.clear();
    free_is_.clear();
    thrust_phi_ = -1.;
    hems_.clear();
    best_match_hems_.clear();

    higgs_.clear();
    dihiggs_.clear();
    gen_higgs_.clear();
    gen_dihiggs_.clear();  

  }

};

