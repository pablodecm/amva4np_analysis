
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

  typedef std::function<int(const Hemisphere &)> FuncI;
  typedef std::vector<std::function<int(const Hemisphere &)>> FuncIVec;
  typedef std::function<double(const Hemisphere &)> FuncD;
  typedef std::vector<std::function<double(const Hemisphere &)>> FuncDVec;
  typedef std::map<std::string,std::function<double(const Hemisphere &)>> FuncDMap;


  FrankSelector(TTree * /*tree*/ =0, TChain * tc_hm = nullptr, std::size_t n_h_mix = 1,
                std::vector<std::string> nn_vars = { "thrustMayor","thrustMinor",
                                                     "sumPz","invMass"}) :
    BaseSelector<EventClass>(0)
    {

        // hemisphere matching variables 
        FuncIVec funcIVec = { FuncI( [] (const Hemisphere & hem) {
                                            int nJets = Hemisphere::nJets(hem);
                                            return ( nJets > 3 ? 4 : nJets);
                                          }),
                                          FuncI( [] (const Hemisphere & hem) {
                                            int nTags = Hemisphere::nTags(hem, "BTag", 0.5);
                                            return ( nTags > 3 ? 4 : nTags);
                                          })};
        FuncDMap funcDMap = {{"thrustMayor", FuncD(&Hemisphere::thrustMayor)},
                             {"thrustMinor",FuncD(&Hemisphere::thrustMinor)},
                             {"sumPz", FuncD(&Hemisphere::sumPz)},
                             {"invMass", FuncD(&Hemisphere::invMass)}};

        FuncDVec funcDVec; 
        for (const auto & nn_var : nn_vars) {
          if (funcDMap.count(nn_var) < 1) {
            std::cout << nn_var << " not present in function map, skipping " << std::endl; 
          } else {
            funcDVec.emplace_back(funcDMap.at(nn_var));
          }
        }
 

        this->addOperator(new FullWriter<EventClass>(true));
        this->addOperator(new ThrustAxisFinder<EventClass>());
        this->addOperator(new HemisphereProducer<EventClass>());
        this->addOperator(new HemisphereMixer<EventClass>(tc_hm, funcIVec,funcDVec));
        this->addOperator(new FrankWriter<EventClass>(n_h_mix, 1, true));

      }

    virtual ~FrankSelector() {}

};

#endif


