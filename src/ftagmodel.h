//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGMODEL_H
#define _FTAGMODEL_H

#include <cstdlib>
#include <string>
#include <vector>
#include <deque>

#include "contextmodel.h"
#include "transitionmodel.h"
#include "ftagtable.h"
#include "ftagstates.h"
#include "transitionestimator.h"
#include "emissionestimator.h"

class FTAGModel
{
public:
   
   FTAGModel();
   FTAGModel(const ContextModel& eModel, const TransitionModel& tModel);
   virtual ~FTAGModel();

   void setSequences(const std::string& seqA, const std::string& seqB,
                     size_t winA = 0,  size_t winB = 0);

   void setEmissionModel(const ContextModel& eModel);
   const ContextModel& getEmissionModel() const;
   
   void setTransitionModel(const TransitionModel& tModel);
   const TransitionModel& getTransitionModel() const;

   void setEmissionEstimator(const EmissionEstimator& ee);
   const EmissionEstimator& getEmissionEstimator() const;

   void setTransitionEstimator(const TransitionEstimator& te);
   const TransitionEstimator& getTransitionEstimator() const;

   double forward();
   double backward();
   double viterbi();
   void bw();

   void viterbiTrace(std::string& outA, std::string& outB, 
                     std::deque<FTAGTableStructs::Trace>& outTrace);

   void computeStateDistribution(std::vector<double>& outDist) const;

   double posterior(const FTAGTableStructs::Trace& trace) const;
   
protected:

   void forwardInit();
   void forwardStepIP(size_t i, size_t p);
   void forwardStepIJKP(size_t i, size_t j, size_t k, size_t p);
   void forwardStepIPQR(size_t i, size_t p, size_t q, size_t r);

   void backwardInit();
   void backwardStepIP(size_t i, size_t p);
   void backwardStepIJKP(size_t i, size_t j, size_t k, size_t p);
   void backwardStepIPQR(size_t i, size_t p, size_t q, size_t r);

   void viterbiInit();
   void viterbiStepIP(size_t i, size_t p);
   void viterbiStepIJKP(size_t i, size_t j, size_t k, size_t p);
   void viterbiStepIPQR(size_t i, size_t p, size_t q, size_t r);

   void bwTransitionInit();
   void bwTransitionStepIP(size_t i, size_t p);
   void bwTransitionStepIJKP(size_t i, size_t j, size_t k, size_t p);
   void bwTransitionStepIPQR(size_t i, size_t p, size_t q, size_t r);

   void bwEmissionInit();
   void bwEmissionStepIP(size_t i, size_t p);
   void bwEmissionStepIJKP(size_t i, size_t j, size_t k, size_t p);
   void bwEmissionStepIPQR(size_t i, size_t p, size_t q, size_t r);

protected:
   
   std::vector<DNA> _seqA;
   std::vector<DNA> _seqB;
   size_t _winA;
   size_t _winB;

   TransitionModel _tm; // transition model
   ContextModel _em; // emission model
   FTAGTable<double, FTAGTableStructs::Sum> _ft; // forward table
   FTAGTable<double, FTAGTableStructs::Sum> _bt; // backward table
   FTAGTable<FTAGTableStructs::VCell, 
             FTAGTableStructs::Max> _vt; // viterbi table
   TransitionEstimator _te; // expectation maximization object for transitions
   EmissionEstimator _ee;

   double _pr;

private:
   FTAGModel(const FTAGModel&);
   FTAGModel& operator=(const FTAGModel&);
};


#endif
