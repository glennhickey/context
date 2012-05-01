//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _TRANSITIONESTIMATOR_H
#define _TRANSITIONESTIMATOR_H

#include <cstdlib>

#include "ftagstates.h"
#include "matrix.h"

class TransitionEstimator
{
public:
   TransitionEstimator();
   ~TransitionEstimator();
   TransitionEstimator(const TransitionEstimator&);
   TransitionEstimator& operator=(const TransitionEstimator&);

   void init();
   void next();
   void finish();

   const Matrix<double>& getMatrix() const
   {
      return _current;
   }

   void setBias(double bias);
   void setBias(const Matrix<double>& bias);
   void setSeqProb(double pr)
   {
      _prSeq = pr;
   }

   void update(FTAGStates::State from, FTAGStates::State to, double val)
   {
      if (from == FTAGStates::S_MIzM && to == FTAGStates::S_MDx && val != 0.)
      {
         assert(false);
      }
      _current.add(from, to, val);
   }

protected:

   Matrix<double> _current;
   Matrix<double> _bias;
   Matrix<double> _numeratorTotal;
   Matrix<double> _denomTotal;
   double _prSeq;
};

#endif
