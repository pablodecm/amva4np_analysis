
#pragma once

#include <vector>
#include <memory>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "mut_framework/mut_dataformats/interface/EventInfo.h"


class EventInfo : public mut::EventInfo {

  public:

    
    EventInfo() {} 

    EventInfo(TTreeReader & reader) {}

    ~EventInfo() {}
    
    void update() {
   }
};

