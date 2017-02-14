
#pragma once

#include <vector>
#include <algorithm>
#include <memory>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"

#include "mut_framework/mut_dataformats/interface/Candidate.h"

  typedef ROOT::Math::PtEtaPhiMVector PtEtaPhiMVector;

  class GenParticleCollection : public std::vector<mut::Candidate> {

    public:

      // had to use smart pointers because TTreeReaderArray is badly implemented
      TTreeReaderArray<float> * gen_particle_pts_;
      TTreeReaderArray<float> * gen_particle_etas_;
      TTreeReaderArray<float> * gen_particle_phis_;
      TTreeReaderArray<float> * gen_particle_masss_;
      TTreeReaderArray<int> * gen_particle_pids_;
      TTreeReaderArray<int> * gen_particle_d1s_;
      TTreeReaderArray<int> * gen_particle_d2s_;

      // flag to use after pythia or tree level objects
      bool after_radiation_ = false;

      GenParticleCollection() :
        gen_particle_pts_(nullptr),
        gen_particle_etas_(nullptr),
        gen_particle_phis_(nullptr),
        gen_particle_masss_(nullptr),
        gen_particle_pids_(nullptr)
      {
      }
      
      GenParticleCollection(TTreeReader & reader, std::string obj_name,
                            bool load = false, bool after_radiation = false) :
        after_radiation_(after_radiation)
    {
      if (load) { 
        gen_particle_pts_   = new TTreeReaderArray<float>(reader, (obj_name+".PT").c_str());
        gen_particle_etas_  = new TTreeReaderArray<float>(reader, (obj_name+".Eta").c_str());
        gen_particle_phis_  = new TTreeReaderArray<float>(reader, (obj_name+".Phi").c_str());
        gen_particle_masss_ = new TTreeReaderArray<float>(reader, (obj_name+".Mass").c_str());
        gen_particle_pids_ = new TTreeReaderArray<int>(reader, (obj_name+".PID").c_str());
        gen_particle_d1s_ = new TTreeReaderArray<int>(reader, (obj_name+".D1").c_str());
        gen_particle_d2s_ = new TTreeReaderArray<int>(reader, (obj_name+".D2").c_str());
      } else {
        gen_particle_pts_   = nullptr;
        gen_particle_etas_  = nullptr;
        gen_particle_phis_  = nullptr;
        gen_particle_masss_ = nullptr;
        gen_particle_pids_ = nullptr;
        gen_particle_d1s_ = nullptr;
        gen_particle_d2s_ = nullptr;
      }
    } 

    ~GenParticleCollection() {
      delete gen_particle_pts_; 
      delete gen_particle_etas_; 
      delete gen_particle_phis_; 
      delete gen_particle_masss_; 
      delete gen_particle_pids_; 
      delete gen_particle_d1s_; 
      delete gen_particle_d2s_; 
    }
    
    void update() {

      if (gen_particle_pts_ != nullptr) {

        // delete previous elements
        this->clear();

        // use default order
        std::vector<int> order(gen_particle_pts_->GetSize());
        std::iota(order.begin(), order.end(), 0); 


        // to keep b bbar b bbar indexes
        std::vector<int> b_ids;
        // iterate over particles in order 
        for (const auto & i : order) {
          
          if ((*gen_particle_pids_)[i] == 25) {
            if ((*gen_particle_d1s_)[i] != (*gen_particle_d2s_)[i])  {
              // always d1 is b and b2 is bbar 
              b_ids.emplace_back((*gen_particle_d1s_)[i]);
              b_ids.emplace_back((*gen_particle_d2s_)[i]);
            }
          }

         if (after_radiation_) {
           // check if element is b or bbar
           std::vector<int>::iterator it = std::find(b_ids.begin(), b_ids.end(), i);
           if (it != b_ids.end())  { // it is a Higgs decay
              if ((*gen_particle_d1s_)[i] != (*gen_particle_d2s_)[i])  {
                *it = i; // is final b/bquark
              } else {
                *it = (*gen_particle_d1s_)[i];
              }
           }
         }
        }

        // iterate over particles in order 
        for (const auto & i : b_ids) {
          PtEtaPhiMVector gen_particle_lv((*gen_particle_pts_)[i],
                                        (*gen_particle_etas_)[i],
                                        (*gen_particle_phis_)[i],
                                        (*gen_particle_masss_)[i]); 
          // new element using constructor from PtEtaPhiM
      	  this->emplace_back(gen_particle_lv);
        }

      }
    }
    
};

