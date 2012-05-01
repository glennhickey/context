//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>

#include "ftagmodel.h"
#include "ftagstates.h"

using namespace FTAGStates;
using namespace FTAGTableStructs;
using namespace std;

double FTAGModel::forward()
{
   forwardInit();

   for (size_t i = 0; i < _seqA.size(); ++i)
   {
      for (size_t p = 0; p < _seqB.size(); ++p)
      {
         forwardStepIP(i, p);
         
         for (size_t j = i + 1; j < i + 1 + _winA && j < _seqA.size(); ++j)
         {
            for (size_t k = j; k < i + 1 + _winA && k < _seqA.size(); ++k)
            {
               forwardStepIJKP(i, j, k, p);
            }
         }

         for (size_t q = p + 1; q < p + 1 + _winB && q < _seqB.size(); ++q)
         {
            for (size_t r = q; r < p + 1 + _winB && r < _seqB.size(); ++r)
            {
               forwardStepIPQR(i, p, q, r);
            }
         }
      }
   }

   double total = 0.;
   for (size_t s = 0; s < IPStatesSize; ++s)
   {
      total += _ft.get1(_seqA.size() - 1, _seqB.size() - 1, IPStates[s]) *
         _tm.prTrans(IPStates[s], S_End);
   }
   
   _pr = total;
   return total;
}

void FTAGModel::forwardInit()
{
   _ft.resize(_seqA.size(), _seqB.size(), _winA, _winB);
   _ft.reset(Sum::seed());

   _ft.set1(0, 0, S_Start, 1);
}

void FTAGModel::forwardStepIP(size_t i, size_t p)
{
   // New "Single" Insert/Match/Delete
   if (i > 0)
   {
      _ft.set1(i, p, S_Dx, _em.prFullD(_seqA[i]) * 
               _ft.all1f(i-1, p, S_Dx, XStates, XStatesSize));
   }
   if (p > 0)
   {
      _ft.set1(i, p, S_Ix, _em.prFullI(_seqB[p]) *
               _ft.all1f(i, p-1, S_Ix, XStates, XStatesSize));
   }
   if (i > 0 && p > 0)
   {
      _ft.set1(i, p, S_Mx, _em.prFullM(_seqA[i], _seqB[p]) *
               _ft.all1f(i-1, p-1, S_Mx, XStates, XStatesSize));         
   }
   
   // Close context and revert back into IP "Single" state
   // (no emission for this type of transition)
   if (i > 1)
   {
      _ft.set1(i, p, S_MDx, 
               _ft.all1fFrom2(i, p, MDxStates, MDxStatesSize));
   }
   if (p > 1)
   {
      _ft.set1(i, p, S_MIx, 
               _ft.all1fFrom3(i, p, MIxStates, MIxStatesSize));
   }
}

void FTAGModel::forwardStepIJKP(size_t i, size_t j, size_t k, size_t p)
{
   assert(i < j && j <= k);

   // Open new context ("outside transition")
   if (j == k && i > 0 && p > 0)
   {
      _ft.set2(i, k, k, p, S_MDx, _em.prFullMD(_seqA[i], _seqA[k], _seqB[p]) *
               _ft.all1f(i-1, p-1, S_MDx, XStates, XStatesSize));
   }

   // Fill in existing context ("inside transition")
   if (i > 0 && k > 0 && k > j && p > 0)
   {
      _ft.set2(i, j, k, p, S_MDyM, _em.prFullMD(_seqA[i], _seqA[k], _seqB[p]) *
               _ft.all2f(i-1, j, k-1, p-1, S_MDyM, IJKPStates, IJKPStatesSize));
   }
}

void FTAGModel::forwardStepIPQR(size_t i, size_t p, size_t q, size_t r)
{
   assert(p < q && q <= r); 

   // Open new context ("outside transition")
   if (q == r && i > 0 && p > 0)
   {
      _ft.set3(i, p, r, r, S_MIx, _em.prFullMI(_seqA[i], _seqB[p], _seqB[r]) *
               _ft.all1f(i-1, p-1, S_MIx, XStates, XStatesSize));
   }

   // Fill in existing context ("inside transition")
   if (i > 0 && p > 0 && r > 0 && r > q)
   {
      _ft.set3(i, p, q, r, S_MIzM, _em.prFullMI(_seqA[i], _seqB[p], _seqB[r]) *
               _ft.all3f(i-1, p-1, q, r-1, S_MIzM, IPQRStates, IPQRStatesSize));
   }
}


