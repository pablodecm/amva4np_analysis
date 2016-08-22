
#include "amva4np/analysis/interface/BaseSelector.h"
#include "amva4np/analysis/interface/MixingSelector.h"
#include "amva4np/analysis/interface/Event.h"
#include "amva4np/analysis/interface/Hemisphere.h"


namespace {
  struct Delphes_ExtEvent {
    DelphesEvent dummya; 
    ExtEvent<DelphesEvent> dummya1; 
    MixingSelector<ExtEvent<DelphesEvent>> dummya2;
  };
}
