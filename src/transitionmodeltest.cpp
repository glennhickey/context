//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <string>
#include <cmath>

#include "transitionmodel.h"
#include "ftagstates.h"
#include "ftagparams.h"

using namespace std;
using namespace FTAGStates;

static double sumToOne(const TransitionModel& tm,
                       State state, bool outside = true)
{
   double rSum = 0.;

   size_t size = outside ? OutsideToStatesSize : InsideToStatesSize;
   const State* states = outside ? OutsideToStates : InsideToStates;

   for (size_t s = 0; s < size; ++s)
   {
      rSum += tm.prTrans(state, states[s]);
   }

   return rSum;
}

int main(int argc, char** argv)
{
   FTAGParams params;
   params.setMu(0.1);
   params.setGamma(0.2);
   params.setRD(0.1);
   params.setRI(0.1);
   params.setT(0.9);
   params.setRMD(0.1);
   params.setRMI(0.1);
   params.setPE(0.09);
   params.setPCD(0.6);
   params.setPCI(0.6);
   params.setKD(0.9);
   params.setKI(0.9);
      
   TransitionModel tm;
   tm.setSimplified(params);

   for (size_t s = 0; s < S_Max; ++s)
   {
      double oSum = sumToOne(tm, (State)s, true);
      double iSum = sumToOne(tm, (State)s, false);

      if (fabs(oSum - 1.) < 0.000000001)
         oSum = 1.;
      if (fabs(iSum - 1.) < 0.000000001)
         iSum = 1.;

      if ((oSum != 0. && oSum != 1.) || (iSum != 0. && iSum != 1.))
         cout << "**";

      cout << s << " oSum=" << oSum << " iSum=" << iSum << endl;
   }
   return 0;
}
