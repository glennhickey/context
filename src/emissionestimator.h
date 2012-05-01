//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _EMISSIONESTIMATOR_H
#define _EMISSIONESTIMATOR_H

#include "contextmodel.h"

class EmissionEstimator
{
public:

   EmissionEstimator();
   ~EmissionEstimator();
   
   EmissionEstimator(const EmissionEstimator&);
   EmissionEstimator& operator=(const EmissionEstimator& ee);

   void init();
   void next();
   void finish();

   const ContextModelTable& getTable() const
   {
      return _current;
   }

   void setBias(double bias)
   {
      _bias.setAll(bias);
   }
   void setBias(const ContextModelTable& bias)
   {
      _bias = bias;
   }
   void setSeqProb(double pr)
   {
      _prSeq = pr;
   }
   
   void updateD(DNA i, double val) 
   {
      _current.addD(i, val);
   }
   void updateI(DNA p, double val) 
   {
      _current.addI(p, val);
   }
   void updateM(DNA i, DNA p, double val) 
   {
      _current.addM(i, p, val);
   }
   void updateMD(DNA i, DNA k, DNA p, double val)  
   {
      _current.addMD(i, k, p, val);
   }
   void updateMI(DNA i, DNA p, DNA r, double val)  
   {
      _current.addMI(i, p, r, val);
   }
   
protected:
   
   ContextModelTable _current;
   ContextModelTable _bias;
   ContextModelTable _numeratorTotal;
   ContextModelTable _denomTotal;
   double _prSeq;
};

#endif

