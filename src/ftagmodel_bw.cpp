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

void FTAGModel::bw()
{
   const int lenA = _seqA.size() - 1;
   const int lenB = _seqB.size() - 1;

   bwTransitionInit();
   bwEmissionInit();
   
   for (int i = lenA; i >= 0; --i)
   {
      for (int p = lenB; p >= 0; --p)
      {
         // don't recount init-ed cells
         if (i < lenA || p < lenB)
         {
            bwTransitionStepIP(i, p);
            bwEmissionStepIP(i, p);
         }

         for (int k = min(i + 1 + _winA, _seqA.size())-1; k > i; --k)
         {
            for (int j = k; j > i; --j)
            {
               // don't recount init-ed cells
               if (i < lenA || p < lenB || j > i + 1)
               {
                  bwTransitionStepIJKP(i, j, k, p);
                  bwEmissionStepIJKP(i, j, k, p);
               }
            }
         }

         for (int r = min(p + 1 + _winB, _seqB.size()) - 1; r > p; --r)
         {
            for (int q = r; q > p; --q)
            {
               // don't recount init-ed cells 
               if (i < lenA || p < lenB || q > p + 1)
               {
                  bwTransitionStepIPQR(i, p, q, r);
                  bwEmissionStepIPQR(i, p , q, r);
               }
            }
         }
      }
   }

   _te.next();
   _ee.next();
}

void FTAGModel::bwTransitionInit()
{
   _te.setSeqProb(_pr);
   for (size_t s = 0; s < IPStatesSize; ++s)
   {
      _te.update(IPStates[s], S_End,
                 _ft.get1(_seqA.size()-1, _seqB.size()-1, IPStates[s]) *
                 _tm.prTrans(IPStates[s], S_End));
//                      _bt.get1(_seqA.size() - 1, _seqB.size() - 1, S_End));
   }
}

void FTAGModel::bwEmissionInit()
{
   _ee.setSeqProb(_pr);
   // don't think anything needs to go here (transition to end doesn't 
   // emit) but leave for the moment
}

void FTAGModel::bwTransitionStepIP(size_t i, size_t p)
{
   size_t s, k, r;
   
   for (s = 0; s < IPStatesSize; ++s)
   {
      // From "single" state to "single" state
      if (i < _seqA.size() - 1)
      {
         _te.update(IPStates[s], S_Dx, 
                    _ft.get1(i, p, IPStates[s]) *
                    _em.prFullD(_seqA[i+1]) *
                    _tm.prTrans(IPStates[s], S_Dx) *
                    _bt.get1(i+1, p, S_Dx));
      }
            
      if (p < _seqB.size() - 1)
      {
         _te.update(IPStates[s], S_Ix,
                    _ft.get1(i, p, IPStates[s]) *
                    _em.prFullI(_seqB[p+1]) *
                    _tm.prTrans(IPStates[s], S_Ix) *
                    _bt.get1(i, p+1, S_Ix));
      }
      if (i < _seqA.size() - 1 && p < _seqB.size() - 1)
      {
         _te.update(IPStates[s], S_Mx,
                    _ft.get1(i, p, IPStates[s]) *
                    _em.prFullM(_seqA[i+1], _seqB[p+1]) *
                    _tm.prTrans(IPStates[s], S_Mx) *
                    _bt.get1(i+1, p+1, S_Mx));
      }


      // new context delete
      if (i < _seqA.size() - 1 && p < _seqB.size() - 1)
      {
         for (k = i + 2; k < min(i + 2 + _winA, _seqA.size()); ++k)
         {
            _te.update(IPStates[s], S_MDx,
                       _ft.get1(i, p, IPStates[s]) *
                       _em.prFullMD(_seqA[i+1], _seqA[k], _seqB[p+1]) *
                       _tm.prTrans(XStates[s], S_MDx) *
                       _bt.get2(i+1, k, k, p+1, S_MDx));
         }
      }
      
      // new context insert
      if (i < _seqA.size() - 1 && p < _seqB.size() - 1)
      {
         for (r = p + 2; r < min(p + 2 + _winB, _seqB.size()); ++r)
         {
            _te.update(IPStates[s], S_MIx,
                       _ft.get1(i, p, IPStates[s]) *
                       _em.prFullMI(_seqA[i+1], _seqB[p+1], _seqB[r]) *
                       _tm.prTrans(XStates[s], S_MIx) *
                       _bt.get3(i+1, p+1, r, r, S_MIx));
         }
      }
   }
}

void FTAGModel::bwEmissionStepIP(size_t i, size_t p)
{
   _ee.updateD(_seqA[i], _ft.get1(i, p, S_Dx) * _bt.get1(i, p, S_Dx));
   _ee.updateI(_seqB[p], _ft.get1(i, p, S_Ix) * _bt.get1(i, p, S_Ix));
   _ee.updateM(_seqA[i], _seqB[p], 
               _ft.get1(i, p, S_Mx) * _bt.get1(i, p, S_Mx)); 
}

void FTAGModel::bwTransitionStepIJKP(size_t i, size_t j, size_t k, size_t p)
{
   assert(i < j && j <= k);

   size_t s;


   for (s = 0; s < MDxStatesSize; ++s)
   {
      if (i + 1 == j) 
      {
         // closed filled context MD 
         _te.update(MDxStates[s], S_Close,
                    _ft.get2(i, j, k, p, MDxStates[s]) *
                    _tm.prTrans(MDxStates[s], S_Close) *
                    _bt.get1(k, p, S_MDx));
      }
      // continued context
      if (j > i + 1 && k + 1 < min(i + 2 + _winA, _seqA.size()) && 
          p + 1 < _seqB.size())
      {
         _te.update(MDxStates[s], S_MDyM,
                    _ft.get2(i, j, k, p, MDxStates[s]) *
                    _em.prFullMD(_seqA[i+1], _seqA[k+1], _seqB[p+1]) *
                    _tm.prTrans(MDxStates[s], S_MDyM) *
                    _bt.get2(i+1, j, k+1, p+1, S_MDyM));
      }
   }
}

void FTAGModel::bwEmissionStepIJKP(size_t i, size_t j, size_t k, size_t p)
{
   size_t s;
   for (s = 0; s < MDxStatesSize; ++s)
   {
      _ee.updateMD(_seqA[i], _seqA[k], _seqB[p], 
                   _ft.get2(i, j, k, p, MDxStates[s]) * 
                   _bt.get2(i, j, k, p, MDxStates[s]));
   }
}

void FTAGModel::bwTransitionStepIPQR(size_t i, size_t p, size_t q, size_t r)
{
   assert(p < q && q <= r); 

   size_t s;

   for (s = 0; s < MIxStatesSize; ++s)
   {
      if (p + 1 == q) 
      {
         // closed filled context MI
         _te.update(MIxStates[s], S_Close,
                    _ft.get3(i, p, q, r, MIxStates[s]) *
                    _tm.prTrans(MIxStates[s], S_Close) *
                    _bt.get1(i, r, S_MIx));
      }

      // continued context MI
      if (i + 1 < _seqA.size() && q > p + 1 && 
          r + 1 < min(p + 2 + _winB, _seqB.size()))
      {
         _te.update(MIxStates[s], S_MIzM,
                    _ft.get3(i, p, q, r, MIxStates[s]) *
                    _em.prFullMI(_seqA[i+1], _seqB[p+1], _seqB[r+1]) *
                    _tm.prTrans(MIxStates[s], S_MIzM) *
                    _bt.get3(i+1, p+1, q, r+1, S_MIzM));
      }
   }
}

void FTAGModel::bwEmissionStepIPQR(size_t i, size_t p, size_t q, size_t r)
{
   size_t s;
   for (s = 0; s < MIxStatesSize; ++s)
   {
      _ee.updateMI(_seqA[i], _seqB[p], _seqB[r], 
                   _ft.get3(i, p, q, r, MIxStates[s]) * 
                   _bt.get3(i, p, q, r, MIxStates[s]));
   }
}
