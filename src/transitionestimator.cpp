//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cstdlib>
#include <cassert>
#include <iostream>

#include "transitionestimator.h"

using namespace std;
using namespace FTAGStates;

TransitionEstimator::TransitionEstimator() :
   _current(S_Max), _bias(S_Max), _numeratorTotal(S_Max), _denomTotal(S_Max),
   _prSeq(0.)
{
   _current.setAll(0.);
   _bias.setAll(0.);
   _numeratorTotal.setAll(0.);
   _denomTotal.setAll(0.);
}

TransitionEstimator::~TransitionEstimator() 
{

}

TransitionEstimator::TransitionEstimator(const TransitionEstimator& te) :
   _current(te._current),
   _bias(te._bias),
   _numeratorTotal(te._numeratorTotal),
   _denomTotal(te._denomTotal),
   _prSeq(te._prSeq)
{

}

TransitionEstimator& TransitionEstimator::operator=(const 
                                                    TransitionEstimator& te)
{
   _current = te._current;
   _bias = te._bias;
   _numeratorTotal = te._numeratorTotal;
   _denomTotal = te._denomTotal;
   _prSeq = te._prSeq;
   return *this;
}

void TransitionEstimator::init()
{
   _current.setAll(0.);
   _numeratorTotal.setAll(0.);
   _denomTotal.setAll(0.);
}

// compute the estimated fractions of reaching each state
void TransitionEstimator::next()
{
   double denomInside[S_Max];
   double denomOutside[S_Max];
   size_t i, j;

   // remember that we have to divide by total probability of the 
   // alignment (which was computed by forward, then stored and passed here)
   _current.multByScalar(1. / _prSeq);
   _current += _bias;
   
   // roll it up to get denominators (incorporating bias as each step)
   for (i = 0; i < S_Max; ++i)
   {
      denomInside[i] = 0.;
      for (j = 0; j < InsideToStatesSize; ++j)
      {
         denomInside[i] += _current.get(i, InsideToStates[j]);
      }
      denomOutside[i] = 0.;
      for (j = 0; j < OutsideToStatesSize; ++j)
      {
         denomOutside[i] += _current.get(i, OutsideToStates[j]);
      }
   }

   // add to numerator
   _numeratorTotal += _current;

   // add to denominator
   for (i = 0; i < S_Max; ++i)
   {
      for (j = 0; j < InsideToStatesSize; ++j)
      {
         _denomTotal.add(i, InsideToStates[j], denomInside[i]);
      }
      for (j = 0; j < OutsideToStatesSize; ++j)
      {
         _denomTotal.add(i, OutsideToStates[j], denomOutside[i]);
      }
   }
   
   // reset 
    _current.setAll(0.);
}

void TransitionEstimator::finish()
{
   _current = _numeratorTotal;
   _current /= _denomTotal;
}

void TransitionEstimator::setBias(double val)
{
   _bias.setAll(0);

   // OUTSIDE

   const State from[] = {S_Start, S_Mx, S_Dx, S_Ix, S_MDx, S_MIx};
   const State to[] = {S_End, S_Mx, S_Dx, S_Ix, S_MDx, S_MIx};

   for (size_t i = 0; i < 6; ++i)
   {
      for (size_t j = 0; j < 6; ++j)
      {
         _bias.set(from[i], to[j], val);
      }
   }

   // INSIDE 

   // Match-Delete Start
   _bias.set(S_MDx, S_Close, val);
   _bias.set(S_MDx, S_MDyM, val);
   
   // Match-Delete Continue
   _bias.set(S_MDyM, S_Close, val);
   _bias.set(S_MDyM, S_MDyM, val);

   // Match-Insert Start
   _bias.set(S_MIx, S_Close, val);
   _bias.set(S_MIx, S_MIzM, val);
   
   // Match-Insert Continue
   _bias.set(S_MIzM, S_Close, val);
   _bias.set(S_MIzM, S_MIzM, val);
}
