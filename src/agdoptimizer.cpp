//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <cstdlib>
#include <limits>
#include <algorithm>

#include "agdoptimizer.h"
#include "ftagstates.h"

using namespace std;
using namespace FTAGStates;

AGDOptimizer::AGDOptimizer() : _nRestarts(1), _maxIt(250), _threshold(1e-7),
                               _useSaved(false), _relative(false), _curTrial(0),
                               _randomOrder(false), _pvec(FTAGParams::Par_Max)
{
   for (size_t i = 0; i < FTAGParams::Par_Max; ++i)
   {
      _pvec[i] = (FTAGParams::Param)i;
   }
}

void AGDOptimizer::setModels(const ContextModel& em, const TransitionModel& tm)
{
   _emTarget = em;
   _tmTarget = tm;
}

void AGDOptimizer::setLoopParams(size_t nRestarts, size_t maxIt,
                                double threshold, const FTAGParams& offset)
{
   _nRestarts = nRestarts;
   _maxIt = maxIt;
   _threshold = threshold;
   _offset = offset;
}

void AGDOptimizer::setUseSaved(bool useSaved)
{
   _useSaved = useSaved;
}

void AGDOptimizer::setSavedInit(const FTAGParams& params)
{
   // todo: be a little smarter about the whole interface here...
   _params = params;
   _paramsSaved = params;
   _useSaved = true;
}

void AGDOptimizer::setRelative(bool relative)
{
   _relative = relative;
}

void AGDOptimizer::setRandomOrder(bool ro)
{
   _randomOrder = ro;
   if (ro == false)
   {
      for (size_t i = 0; i < FTAGParams::Par_Max; ++i)
      {
         _pvec[i] = (FTAGParams::Param)i;
      }
   }
}

void AGDOptimizer::setStateDistribution(const vector<double>& sDist)
{
   _sDist = sDist;
}

double AGDOptimizer::optimize()
{
   double curError = numeric_limits<double>::max();
   _mseSaved = curError;
   for (_curTrial = 0; _curTrial < _nRestarts; ++_curTrial)
   {
      curError = runTrial();
      if (curError < _mseSaved)
      {
         _mseSaved = curError;
         _paramsSaved = _params;
      }
   }
   _params = _paramsSaved;
   return _mseSaved;
}

void AGDOptimizer::initRandom()
{
   _params.randomizeNonFixed();
}

void AGDOptimizer::initLast()
{
   _params = _paramsSaved;
}

// copy only fixed flags and values into the object
void AGDOptimizer::setFixed(const FTAGParams& params)
{
   _params = params;
   _paramsCurrent = params;
   _paramsSaved = params;
   for (size_t i = 0; i < FTAGParams::Par_Max; ++i)
   {
      if (_params.isFixed((FTAGParams::Param)i) == false)
      {         
         _params.set((FTAGParams::Param)i, 0);
         _paramsCurrent.set((FTAGParams::Param)i, 0);
         _paramsSaved.set((FTAGParams::Param)i, 0);
      }
   }
}

double AGDOptimizer::mse()
{
   _emCurrent.setSingleAndDouble(_params);
   _tmCurrent.setSimplified(_params);
   double es = 0.;
   if (_sDist.size() == S_Max)
   {
      ContextModelTable emCurrentScaled(_emCurrent.getTable());
      ContextModelTable emTargetScaled(_emTarget.getTable());
      emCurrentScaled.scale(_sDist);
      emTargetScaled.scale(_sDist);
      TransitionModel tmCurrentScaled(_tmCurrent);
      TransitionModel tmTargetScaled(_tmTarget);
      tmCurrentScaled.scale(_sDist);
      tmTargetScaled.scale(_sDist);
      es = emCurrentScaled.mse(emTargetScaled);
      es += tmCurrentScaled.getProbMatrix().mse(tmTargetScaled.getProbMatrix());
   }
   else
   {
      es = _emCurrent.getTable().mse(_emTarget.getTable());
      es += _tmCurrent.getProbMatrix().mse(_tmTarget.getProbMatrix());
   }
   return es;
}

double AGDOptimizer::runTrial()
{
   if (_curTrial == 0 && _useSaved == true)
   {
      initLast();
   }
   else
   {
      initRandom();
   }

   double prevError = numeric_limits<double>::max();
   double curError = mse();
   double plusError, minusError;
   FTAGParams parPlus, parMinus;
   
   for (size_t i = 0; i < _maxIt; ++i)
   {
/*      cout << "it(" << _curTrial << "." << i << ") p=" << prevError 
           << " c=" << curError 
           << "   DELTA=" << prevError - curError << "\n"
           << _params << endl;
*/
      _paramsCurrent = _params;
      if (_randomOrder == true)
      {
         random_shuffle(_pvec.begin(), _pvec.end());
      }
      for (size_t j = 0; j < _pvec.size(); ++j)
      {
         FTAGParams::Param par = _pvec[j];
         if (_params.isFixed(par) == false && _params.isDependent(par) == false
             && _params.getPhase(par).first <= _emit 
             && _params.getPhase(par).second >= _emit)
         {
            _paramsCurrent = _params;
            _params.set(par, stepUp(par, _paramsCurrent.get(par)));
            plusError = mse();
            parPlus = _params;
            
            _params = _paramsCurrent;
            _params.set(par, stepDown(par, _paramsCurrent.get(par)));
            minusError = mse();
            parMinus = _params;
            
            if (curError - plusError > _threshold && plusError < minusError)
            {
               _params = parPlus;
               curError = plusError;
            }
            else if (curError - minusError > _threshold)
            {
               _params = parMinus;
               curError = minusError;
            }
            else
            {
               _params = _paramsCurrent;
            }
         }
      }
      if (prevError - curError < _threshold)
      {
         break;
      }
      else
      {
         prevError = curError;
      }
   }
   return curError;
}
