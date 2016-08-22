
#ifndef MixingSelector_h
#define MixingSelector_h

#include <memory>

#include "BaseSelector.h"
#include "EventCounter.h"


template <class EventClass> class MixingSelector : public BaseSelector<EventClass> {

  public:


  MixingSelector(TTree * /*tree*/ =0) :
    BaseSelector<EventClass>(0)
    {
      this->addOperator(new EventCounter<EventClass>());
    }

  virtual ~MixingSelector() {}

};

#endif


