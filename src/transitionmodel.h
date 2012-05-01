//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _TRANSITIONMODEL_H
#define _TRANSITIONMODEL_H

#include <ostream>
#include <vector>
#include "ctmarkov.h"
#include "ftagparams.h"

class TransitionModel : public CtMarkov
{
public:
   
   TransitionModel();
   virtual ~TransitionModel();
   TransitionModel(const TransitionModel& tm);
   TransitionModel& operator=(const TransitionModel& tm);
      
   void setSimplified(const FTAGParams& tp);

   void scale(const std::vector<double>& svec);

   // for preliminary debugging purposes only
   void setMatchOnly(double p);

   bool verify(double threshold = 0.00000001) const;
};

std::ostream& operator<<(std::ostream& os, const TransitionModel& tm);

#endif
