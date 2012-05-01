//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <algorithm>
#include <cmath>

#include "ftagmodel.h"
#include "ftagstates.h"


using namespace FTAGStates;
using namespace FTAGTableStructs;
using namespace std;

double FTAGModel::backward()
{
   backwardInit();
   const int lenA = _seqA.size() - 1;
   const int lenB = _seqB.size() - 1;

   for (int i = lenA; i >= 0; --i)
   {
      for (int p = lenB; p >= 0; --p)
      {
         // don't overwrite init-ed cells
         if (i < lenA || p < lenB)
         {
            backwardStepIP(i, p);
         }

         for (int k = min(i + 1 + _winA, _seqA.size()) - 1; k > i; --k)
         {
            for (int j = k; j > i; --j)
            {
               // don't overwrite init-ed cells
               if (i < lenA || p < lenB || j > i + 1)
               {
                  backwardStepIJKP(i, j, k, p);
               }
            }
         }

         for (int r = min(p + 1 + _winB, _seqB.size()) - 1; r > p; --r)
         {
            for (int q = r; q > p; --q)
            {
               // don't overwrite init-ed cells (maybe not necessary here)
               if (i < lenA || p < lenB || q > p + 1)
               {
                  backwardStepIPQR(i, p, q, r);
               }
            }
         }
      }
   }

   return _bt.get1(0, 0, S_Start);
}

void FTAGModel::backwardInit()
{
   _bt.resize(_seqA.size(), _seqB.size(), _winA, _winB);
   _bt.reset(Sum::seed());

   for (size_t s = 0; s < IPStatesSize; ++s)
   {
      _bt.set1(_seqA.size() - 1, _seqB.size() - 1, IPStates[s], 
               _tm.prTrans(IPStates[s], S_End));
   }
}

void FTAGModel::backwardStepIP(size_t i, size_t p)
{
   size_t s, k, r;
   double pr;
   
   for (s = 0; s < XStatesSize; ++s)
   {
      pr = 0.;
      
      // Back from "single" state
      if (i < _seqA.size() - 1)
      {
         pr += _bt.get1(i+1, p, S_Dx) * _em.prFullD(_seqA[i+1]) * 
            _tm.prTrans(XStates[s], S_Dx);
      }
      if (p < _seqB.size() - 1)
      {
         pr += _bt.get1(i, p+1, S_Ix) * _em.prFullI(_seqB[p+1]) *
            _tm.prTrans(XStates[s], S_Ix);
      }
      if (i < _seqA.size() - 1 && p < _seqB.size() - 1)
      {
         pr += _bt.get1(i+1, p+1, S_Mx) * _em.prFullM(_seqA[i+1], _seqB[p+1]) *
            _tm.prTrans(XStates[s], S_Mx);
      }

      // Back from new context delete
      if (i < _seqA.size() - 1 && p < _seqB.size() - 1)
      {
         for (k = i + 2; k < min(i + 2 + _winA, _seqA.size()); ++k)
         {
            pr += _bt.get2(i+1, k, k, p+1, S_MDx) *
               _em.prFullMD(_seqA[i+1], _seqA[k], _seqB[p+1]) *
               _tm.prTrans(XStates[s], S_MDx);
         }
      }
      // Back from new context insert
      if (i < _seqA.size() - 1 && p < _seqB.size() - 1)
      {
         for (r = p + 2; r < min(p + 2 + _winB, _seqB.size()); ++r)
         {
            pr += _bt.get3(i+1, p+1, r, r, S_MIx) *
               _em.prFullMI(_seqA[i+1], _seqB[p+1], _seqB[r]) *
               _tm.prTrans(XStates[s], S_MIx);
         }
      }
      
      assert(_bt.get1(i, p, XStates[s]) == Sum::seed());
      _bt.set1(i, p, XStates[s], pr);
//      cout << "set bp1 " << i << "," << p << "(s)" << " = " << pr << endl;
   }
}

void FTAGModel::backwardStepIJKP(size_t i, size_t j, size_t k, size_t p)
{
   assert(i < j && j <= k);

   size_t s;
   double pr;

   for (s = 0; s < MDxStatesSize; ++s)
   {
      pr = 0.;

      if (i + 1 == j) 
      {
         // closed filled context MD 
         pr += _bt.get1(k, p, S_MDx) * _tm.prTrans(MDxStates[s], S_Close);
      }

      // continued context
      if (j > i + 1 && k + 1 < min(i + 2 + _winA, _seqA.size()) && 
          p + 1 < _seqB.size())
      {
         pr += _bt.get2(i+1, j, k+1, p+1, S_MDyM) *
            _em.prFullMD(_seqA[i+1], _seqA[k+1], _seqB[p+1]) *
            _tm.prTrans(MDxStates[s], S_MDyM);
      }

      assert(_bt.get2(i, j, k, p, MDxStates[s]) == Sum::seed());
      _bt.set2(i, j, k, p, MDxStates[s], pr);
   }
}

void FTAGModel::backwardStepIPQR(size_t i, size_t p, size_t q, size_t r)
{
   assert(p < q && q <= r); 

   size_t s;
   double pr;

   for (s = 0; s < MIxStatesSize; ++s)
   {
      pr = 0.;
      
      if (p + 1 == q) 
      {
         // closed filled context MI 
         pr += _bt.get1(i, r, S_MIx) * _tm.prTrans(MIxStates[s], S_Close);
      }

      // continued context MI
      if (i + 1 < _seqA.size() && q > p + 1 &&
         r + 1 < min(p + 2 + _winB, _seqB.size()))
      {
         pr += _bt.get3(i+1, p+1, q, r+1, S_MIzM) * 
            _em.prFullMI(_seqA[i+1], _seqB[p+1], _seqB[r+1]) *
            _tm.prTrans(MIxStates[s], S_MIzM);
      }
      
      assert(_bt.get3(i, p, q, r, MIxStates[s]) == Sum::seed());
      _bt.set3(i, p, q, r, MIxStates[s], pr);
   }

}


