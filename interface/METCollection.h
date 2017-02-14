
#pragma once

#include <vector>
#include <algorithm>
#include <memory>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"

#include "mut_framework/mut_dataformats/interface/MET.h"

  typedef ROOT::Math::PtEtaPhiMVector PtEtaPhiMVector;

  class METCollection : public std::vector<mut::MET> {

    public:

      // had to use smart pointers because TTreeReaderArray is badly implemented
      TTreeReaderArray<float> * met_pts_;
      TTreeReaderArray<float> * met_etas_;
      TTreeReaderArray<float> * met_phis_;

      METCollection() :
        met_pts_(nullptr),
        met_etas_(nullptr),
        met_phis_(nullptr)
      {
      }
      
      METCollection(TTreeReader & reader)
    {
        met_pts_   = new TTreeReaderArray<float>(reader, "MissingET.MET");
        met_etas_  = new TTreeReaderArray<float>(reader, "MissingET.Eta");
        met_phis_  = new TTreeReaderArray<float>(reader, "MissingET.Phi");
    } 

    ~METCollection() {
      delete met_pts_; 
      delete met_etas_; 
      delete met_phis_; 
    }
    
    void update() {

      if (met_pts_ != nullptr) {

        // delete previous elements
        this->clear();

        // use default order
        std::vector<int> order(met_pts_->GetSize());
        std::iota(order.begin(), order.end(), 0); 


        std::vector<int> b_ids;
        // iterate over particles in order 
        for (const auto & i : order) {
          PtEtaPhiMVector met_lv((*met_pts_)[i],
                                 (*met_etas_)[i],
                                 (*met_phis_)[i],
                                 0); 
        // new element using constructor from PtEtaPhiM
      	this->emplace_back(met_lv);
  
      }
    }
  }
    
};

