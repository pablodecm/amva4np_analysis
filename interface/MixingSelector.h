
#ifndef MixingSelector_h
#define MixingSelector_h

#include <memory>

#include "BaseSelector.h"
#include "EventCounter.h"
#include "JetSelection.h"
#include "BTagJetSelection.h"
#include "DiJetPairSelection.h"
#include "DiHiggsPlotter.h"
#include "ThrustAxisFinder.h"
#include "HemisphereProducer.h"
#include "FullWriter.h"
#include "HemisphereWriter.h"
#include "GenJetMatcher.h"
#include "DiHiggsPlotter.h"


template <class EventClass> class MixingSelector : public BaseSelector<EventClass> {

  public:


  MixingSelector(TTree * /*tree*/ =0) :
    BaseSelector<EventClass>(0)
    {
      std::size_t n_CSV = 4;
      this->addOperator(new EventCounter<EventClass>());
      this->addOperator(new JetSelection<EventClass>(2.5, 30., 4));
      this->addOperator(new EventCounter<EventClass>());
      this->addOperator(new BTagJetSelection<EventClass>("BTag", 0.5, n_CSV ));
      this->addOperator(new EventCounter<EventClass>());
      this->addOperator(new BetterDiJetPairSelection<EventClass>("BTag",0.5, n_CSV));
      this->addOperator(new EventCounter<EventClass>());
      this->addOperator(new DiHiggsPlotter<EventClass>({}, true));
      this->addOperator(new GenJetMatcher<EventClass>());
      this->addOperator(new FullWriter<EventClass>(true));
      this->addOperator(new ThrustAxisFinder<EventClass>());
      this->addOperator(new HemisphereProducer<EventClass>());
      this->addOperator(new HemisphereWriter<EventClass>(true));
      this->addOperator(new GenJetMatcher<EventClass>());
      this->addOperator(new HHJetsMatched<EventClass>());
      this->addOperator(new FullWriter<EventClass>(true, "matched"));
      this->addOperator(new DiHiggsPlotter<EventClass>({}, true, "matched"));
    }

  virtual ~MixingSelector() {}

};

#endif


