//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include "emissionestimator.h"

using namespace std;

EmissionEstimator::EmissionEstimator()
{
   _current.setAll(0.);
   _bias.setAll(0.);
   _numeratorTotal.setAll(0.);
   _denomTotal.setAll(0.);
   _prSeq = 0.;
}

EmissionEstimator::~EmissionEstimator()
{

}

EmissionEstimator::EmissionEstimator(const EmissionEstimator& ee) :
   _current(ee._current),
   _bias(ee._bias),
   _numeratorTotal(ee._numeratorTotal),
   _denomTotal(ee._denomTotal),
   _prSeq(ee._prSeq)
{

}

EmissionEstimator& EmissionEstimator::operator=(const EmissionEstimator& ee)
{
   _current = ee._current;
   _bias = ee._bias;
   _numeratorTotal = ee._numeratorTotal;
   _denomTotal = ee._denomTotal;
   _prSeq = ee._prSeq;
   return *this;
}


void EmissionEstimator::init()
{
   _current.setAll(0.);
   _numeratorTotal.setAll(0.);
   _denomTotal.setAll(0.);
}

void EmissionEstimator::next()
{
   // total up current and add the proper values to
   // numerator and denominator.  note to self: this method
   // could be made faster without too much trouble
   
   double dD = 0.;
   double dI = 0.;
   double dM = 0.;
   double dMD = 0.;
   double dMI = 0.;

   // remember that we have to divide by total probability of the 
   // alignment (which was computed by forward, then stored and passed here)
   ContextModelTable temp;
   temp.setAll(_prSeq);
   _current /= temp;

   _current += _bias;  // VERIFY THIS IS CORRECT

   // roll it up to get denominators (incorporating bias as each step)
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      dD += _current.getD((DNA)i);
      dI += _current.getI((DNA)i);

      for (size_t j = 0; j < DNA_Max; ++j)
      {
         dM += _current.getM((DNA)i, (DNA)j);
         
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            dMD += _current.getMD((DNA)i, (DNA)j, (DNA)k);
            dMI += _current.getMI((DNA)i, (DNA)j, (DNA)k);
         }
      }
   }

   // add to numerator
   _numeratorTotal += _current;

   // add to denominator
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      _denomTotal.addD((DNA)i, dD);
      _denomTotal.addI((DNA)i, dI);

      for (size_t j = 0; j < DNA_Max; ++j)
      {
         _denomTotal.addM((DNA)i, (DNA)j, dM);
         
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            _denomTotal.addMD((DNA)i, (DNA)j, (DNA)k, dMD);
            _denomTotal.addMI((DNA)i, (DNA)j, (DNA)k, dMI);
         }
      }
   }
   
   // reset 
   _current.setAll(0.);
}

void EmissionEstimator::finish()
{
   _current = _numeratorTotal;
   _current /= _denomTotal;
}
