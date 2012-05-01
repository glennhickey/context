//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <algorithm>
#include <ctime>

#include "ftagtrain.h"
#include "ftagstates.h"

using namespace FTAGStates;
using namespace std;


FTAGTrain::FTAGTrain()
{

}

FTAGTrain::~FTAGTrain()
{

}

void FTAGTrain::setSequences(const vector<pair<string, string> >& seqs)
{
   _seqs = seqs;
   _wins.clear();
   for (size_t i = 0; i < _seqs.size(); ++i)
   {
      _wins.push_back(pair<size_t, size_t>(0, 0));
   }
}

void FTAGTrain::setSequences(const vector<pair<string, string> >& seqs,
                             size_t winSize)
{
   _seqs = seqs;
   _wins.clear();
   for (size_t i = 0; i < _seqs.size(); ++i)
   {
      size_t lenA = _seqs[i].first.length();
      size_t lenB = _seqs[i].second.length();
      size_t win = min(winSize, min(lenA, lenB));
      _wins.push_back(pair<size_t, size_t>(win, win));
   }
}

void FTAGTrain::setSequences(const vector<pair<string, string> >& seqs,
                             const vector<pair<size_t, size_t> >& wins)
{
   _seqs = seqs;
   _wins = wins;
}

void FTAGTrain::initialize(const AGDOptimizer& opt,
                           const EmissionEstimator& ee, 
                           const TransitionEstimator& te,
                           bool scale)
{
   _opt = opt;
   _em.setSingleAndDouble(_opt.getParams());
   _tm.setSimplified(_opt.getParams());
   _ee = ee;
   _te = te;
   _scale = scale;
}

const TransitionModel& FTAGTrain::getTransitionModel() const
{
   return _tm;
}

const ContextModel& FTAGTrain::getEmissionModel() const
{
   return _em;
}

const FTAGParams& FTAGTrain::getParams() const
{
   return _opt.getParams();
}

void FTAGTrain::emLoop(size_t maxIt, double threshold, size_t convRepeats,
                       double ecrThreshold, EMTrace* trace)
{
   _ftag.setEmissionModel(_em);
   _ftag.setTransitionModel(_tm);
   assert(_em.getTable().verify() == true);
   assert(_tm.verify() == true);

   size_t curRepeat = 0;
   double prPrev = -1e100;
   for (size_t i = 0; i < maxIt; ++i)
   {      
      _te.init();
      _ee.init();

      _ftag.setTransitionEstimator(_te);
      _ftag.setEmissionEstimator(_ee);
 
      double prf = 0.;
      double prb = 0.;
      vector<double> scaleVec(S_Max, 0.);
      
      // Train with all sequcnes
      for (size_t j = 0; j < _seqs.size(); ++j)
      {
         _ftag.setSequences(_seqs[j].first, _seqs[j].second,
                            _wins[j].first, _wins[j].second);

         prf += log(_ftag.forward()); 
         prb += log(_ftag.backward());
         _ftag.bw();
         if (_scale == true)
         {
            vector<double> sv;
            _ftag.computeStateDistribution(sv);
            for (size_t k = 0; k < S_Max; ++k) scaleVec[k] += sv[k];
         }
      }

      _te = _ftag.getTransitionEstimator();
      _ee = _ftag.getEmissionEstimator();

      _te.finish();
      _ee.finish();

      assert(fabs(prb - prf) < 0.00000001);
      
      // compute delta
      double delta = (prb - prPrev);
      
      if (delta < ecrThreshold)
      {
         ++curRepeat;
      }
      else
      {
         curRepeat = 0;
      }
      
      if (delta < threshold || (curRepeat > 0 && curRepeat == convRepeats))
      {
         break;
      }
      else
      {
         // if not within threshold, invert the model and
         // reiterate
         _tm.setProbMatrix(_te.getMatrix());
         _em.setTable(_ee.getTable());
         assert(_em.getTable().verify() == true);
         assert(_tm.verify() == true);
         _opt.setSavedInit(_opt.getParams());
         if (_scale == true)
         {
            _opt.setStateDistribution(scaleVec);
         }
         _opt.setModels(_em, _tm);

         time_t timer = time(NULL);
         string timeStr = ctime(&timer);
         timeStr.erase(0, 4);
         timeStr.erase(timeStr.length() - 6, 6);
         _opt.setEmIt(i);
         double mse = _opt.optimize();
         cout << "EM iteration " << i << " [" << timeStr
              << "]: param MSE=" << mse 
              << " Pr=" << prb << " params=\n"
              << _opt.getParams() << endl; 
         if (trace)
         {
            trace->push_back(pair<double, FTAGParams>(prb, _opt.getParams()));
         }
         _em.setSingleAndDouble(_opt.getParams());
         _tm.setSimplified(_opt.getParams());
         assert(_em.getTable().verify());
         cout << endl;
         
         _ftag.setEmissionModel(_em);
         _ftag.setTransitionModel(_tm);
         prPrev = prb;
      }
   }

}
