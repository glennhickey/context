//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGTRAIN_H
#define _FTAGTRAIN_H

#include <vector>
#include <string>
#include <ostream>

#include "ftagmodel.h"
#include "transitionestimator.h"
#include "emissionestimator.h"
#include "contextmodel.h"
#include "transitionmodel.h"
#include "agdoptimizer.h"

class FTAGTrain
{
public:
   FTAGTrain();
   ~FTAGTrain();

   void setSequences(const std::vector<std::pair<std::string, std::string> >& 
                     seqs);
   void setSequences(const std::vector<std::pair<std::string, std::string> >& 
                     seqs, size_t winSize);
   void setSequences(const std::vector<std::pair<std::string, std::string> >& 
                     seqs, 
                     const std::vector<std::pair<size_t, size_t> >& wins); 

   void initialize(const AGDOptimizer& opt,
                   const EmissionEstimator& ee, 
                   const TransitionEstimator& te,
                   bool scale = true);

   const TransitionModel& getTransitionModel() const;
   const ContextModel& getEmissionModel() const;
   const FTAGParams& getParams() const;

   typedef std::vector<std::pair<double, FTAGParams> > EMTrace;

   void emLoop(size_t maxIt, double threshold, size_t convRepeats, 
               double ecrThreshold, EMTrace* trace = NULL);

protected:

   std::vector<std::pair<std::string, std::string> > _seqs;
   std::vector<std::pair<size_t, size_t> > _wins;

   FTAGModel _ftag;
   EmissionEstimator _ee;
   TransitionEstimator _te;
   ContextModel _em;
   TransitionModel _tm;
   AGDOptimizer _opt;
   bool _scale;
   size_t _convRepeats;
};

#endif
