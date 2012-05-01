//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cassert>

#include "ftagmodel.h"

using namespace std;
using namespace FTAGStates;

FTAGModel::FTAGModel() : _pr(0.0)
{

}

FTAGModel::~FTAGModel()
{

}

void FTAGModel::setEmissionModel(const ContextModel& eModel)
{
   _em = eModel;
}

const ContextModel& FTAGModel::getEmissionModel() const
{
   return _em;
}

void FTAGModel::setTransitionModel(const TransitionModel& tModel)
{
   _tm = tModel;
   _ft.setTransitionModel(&_tm);
   _bt.setTransitionModel(&_tm);
   _vt.setTransitionModel(&_tm);
}

const TransitionModel& FTAGModel::getTransitionModel() const
{
   return _tm;
}


void FTAGModel::setEmissionEstimator(const EmissionEstimator& ee)
{
   _ee = ee;
}

const EmissionEstimator& FTAGModel::getEmissionEstimator() const
{
   return _ee;
}

void FTAGModel::setTransitionEstimator(const TransitionEstimator& te)
{
   _te = te;
}

const TransitionEstimator& FTAGModel::getTransitionEstimator() const
{
   return _te;
}

// INDEXING STARTS AT 1, NOT ZERO, TO KEEP IN LINE WITH THE TABLE!!!
// DON'T FORGET THIS!!!
void FTAGModel::setSequences(const string& seqA, const string& seqB,
                             size_t winA, size_t winB)
{
   _seqA.clear();
   _seqB.clear();
   _seqA.push_back(DNA_Max);
   _seqB.push_back(DNA_Max);

   size_t i;
   for (i = 0; i < seqA.length(); ++i)
   {
      if (seqA[i] != '-')
      {
         _seqA.push_back(DNASubModel::charToDNA(seqA[i]));
      }
   }
   for (i = 0; i < seqB.length(); ++i)
   {
      if (seqB[i] != '-')
      {
         _seqB.push_back(DNASubModel::charToDNA(seqB[i]));
      }
   }

   assert(winA <= seqA.size() && winB <= seqB.size());

   _winA = winA > 0 ? winA : seqA.size();
   _winB = winB > 0 ? winB : seqB.size();
}

void FTAGModel::computeStateDistribution(vector<double>& outDist) const
{
   outDist.clear();
   outDist.resize(S_Max, 0.);

   for (size_t i = 0; i < _seqA.size(); ++i)
   {
      for (size_t p = 0; p < _seqB.size(); ++p)
      {
         for (size_t s = 0; s < IPStatesSize; ++s)
         {
            outDist[IPStates[s]] += _ft.get1(i, p, IPStates[s]) * 
               _bt.get1(i, p, IPStates[s]);
         }
         
         for (size_t j = i + 2; j < i + 1 + _winA && j < _seqA.size(); ++j)
         {
            for (size_t k = j; k < i + 1 + _winA && k < _seqA.size(); ++k)
            {
               for (size_t s = 0; s < IJKPStatesSize; ++s)
               {
                  outDist[IJKPStates[s]] += 
                     _ft.get2(i, j, k, p, IJKPStates[s]) * 
                     _bt.get2(i, j, k, p, IJKPStates[s]);
               }
            }
         }

         for (size_t q = p + 2; q < p + 1 + _winB && q < _seqB.size(); ++q)
         {
            for (size_t r = q; r < p + 1 + _winB && r < _seqB.size(); ++r)
            {
               for (size_t s = 0; s < IPQRStatesSize; ++s)
               {
                  outDist[IPQRStates[s]] += 
                     _ft.get3(i, p, q, r, IPQRStates[s]) * 
                     _bt.get3(i, p, q, r, IPQRStates[s]);
               }
            }
         }
      }
   }

   double total = 0;
   for (size_t i = 0; i < outDist.size(); ++i)
   {
      total += outDist[i];
   }
   for (size_t i = 0; i < outDist.size(); ++i)
   {
      outDist[i] /= total;
   }
}

double FTAGModel::posterior(const FTAGTableStructs::Trace& trace) const
{
   if (trace.isIPState() == true)
   {
      return _ft.get1(trace._x1, trace._x2, trace._s) *
         _bt.get1(trace._x1, trace._x2, trace._s);
   }
   else if (trace.isIJKPState() == true)
   {
      return _ft.get2(trace._x1, trace._x2, trace._x3, trace._x4, trace._s) *
         _bt.get2(trace._x1, trace._x2, trace._x3, trace._x4, trace._s);
   }
   else if (trace.isIPQRState() == true)
   {
      return _ft.get3(trace._x1, trace._x2, trace._x3, trace._x4, trace._s) *
         _bt.get3(trace._x1, trace._x2, trace._x3, trace._x4, trace._s);
   }
   throw string("bad trace in posterior");
   return 0;
}

