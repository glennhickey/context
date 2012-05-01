//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _AGDOPTIMIZER_H
#define _AGDOPTIMIZER_H

#include "contextmodel.h"
#include "transitionmodel.h"

class AGDOptimizer
{
public:

   AGDOptimizer();

   void setModels(const ContextModel& em, const TransitionModel& tm);

   void setLoopParams(size_t nRestarts, size_t maxIt,
                      double threshold, const FTAGParams& offset);

   double optimize();

   const FTAGParams& getParams() const {
      return _params;
   }

   const FTAGParams& getSavedInit() const {
      return _paramsSaved;
   }
   
   void setSavedInit(const FTAGParams& params); 
   void setUseSaved(bool useSaved);
   void setFixed(const FTAGParams& params);
   void setRelative(bool relative);
   void setRandomOrder(bool ro);
   void setStateDistribution(const std::vector<double>& sDist);
   void setEmIt(size_t emit) {
      _emit = emit;
   }
   size_t getEmIt() const {
      return _emit;
   }
   
protected:

   double stepUp(FTAGParams::Param par, double val) const {
      return val + (_relative && val ? val * _offset.get(par) : 
                    _offset.get(par));
   }
   
   double stepDown(FTAGParams::Param par, double val) const {
      return val - (_relative && val ? val * _offset.get(par) : 
                    _offset.get(par));
   }
   
   void initRandom();
   void initLast();
   double mse();
   double runTrial();

   ContextModel _emTarget;
   ContextModel _emCurrent;
   TransitionModel _tmTarget;
   TransitionModel _tmCurrent;
   FTAGParams _params;
   FTAGParams _paramsCurrent;
   FTAGParams _paramsSaved;
   FTAGParams _offset;
   double _mseSaved;
   size_t _nRestarts;
   size_t _maxIt;
   double _threshold;
   bool _useSaved;
   bool _relative;
   size_t _curTrial;
   bool _randomOrder;
   std::vector<FTAGParams::Param> _pvec;
   std::vector<double> _sDist;
   size_t _emit;
};

#endif
