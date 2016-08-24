
#ifndef FrankSelector_h
#define FrankSelector_h

#include <memory>

#include "BaseSelector.h"
#include "EventCounter.h"
#include "FullWriter.h"
#include "ThrustAxisFinder.h"
#include "HemisphereMixer.h"
#include "FrankWriter.h"


template <class EventClass> class FrankSelector : public BaseSelector<EventClass> {

  public:


  FrankSelector(TTree * /*tree*/ =0, TChain * tc_hm = nullptr, std::size_t n_h_mix = 1 ) :
    BaseSelector<EventClass>(0)
    {

      this->addOperator(new FullWriter<EventClass>(true));
      this->addOperator(new ThrustAxisFinder<EventClass>());
      this->addOperator(new HemisphereProducer<EventClass>());
      this->addOperator(new HemisphereMixer<EventClass>(tc_hm));
      this->addOperator(new FrankWriter<EventClass>(n_h_mix, 1, true));

    }

  virtual ~FrankSelector() {}

};

#endif


