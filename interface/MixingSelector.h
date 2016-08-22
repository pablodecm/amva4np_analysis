
#ifndef MixingSelector_h
#define MixingSelector_h

#include <memory>

#include "BaseSelector.h"
#include "EventCounter.h"
#include "JetSelection.h"
#include "BTagJetSelection.h"


template <class EventClass> class MixingSelector : public BaseSelector<EventClass> {

  public:


  MixingSelector(TTree * /*tree*/ =0) :
    BaseSelector<EventClass>(0)
    {
      std::size_t n_CSV = 3;
      this->addOperator(new EventCounter<EventClass>());
      this->addOperator(new JetSelection<EventClass>(2.5, 20., 4));
      this->addOperator(new EventCounter<EventClass>());
      this->addOperator(new BTagJetSelection<EventClass>("BTag", 0.500, n_CSV ));
      this->addOperator(new EventCounter<EventClass>());

    }

  virtual ~MixingSelector() {}

};

#endif


