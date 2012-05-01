//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <algorithm>
#include <cassert>

#include "ftagmodel.h"
#include "ftagstates.h"
#include "dnasubmodel.h"

using namespace FTAGTableStructs;
using namespace FTAGStates;
using namespace std;

double FTAGModel::viterbi()
{
   viterbiInit();

   for (size_t i = 0; i < _seqA.size(); ++i)
   {
      for (size_t p = 0; p < _seqB.size(); ++p)
      {
         viterbiStepIP(i, p);

         for (size_t j = i + 1; j < i + 1 + _winA && j < _seqA.size(); ++j)
         {
            for (size_t k = j; k < i + 1 + _winA && k < _seqA.size(); ++k)
            {
               viterbiStepIJKP(i, j, k, p);
            }
         }

         for (size_t q = p + 1; q < p + 1 + _winB && q < _seqB.size(); ++q)
         {
            for (size_t r = q; r < p + 1 + _winB && r < _seqB.size(); ++r)
            {
               viterbiStepIPQR(i, p, q, r);
            }
         }
      }
   }

   VCell last = Max::seed();
   for (size_t s = 0; s < IPStatesSize; ++s)
   {
      Max::acc(last, 
               _vt.get1(_seqA.size() - 1, _seqB.size() - 1, IPStates[s]) *
               _tm.prTrans(IPStates[s], S_End),
               Trace(IPStates[s], _seqA.size() - 1, _seqB.size() -1));
   }
   return last._val;
}

void FTAGModel::viterbiTrace(string& outA, string& outB, deque<Trace>& outTrace)
{
   int idx = _seqA.size() + _seqB.size();
   outA.resize(idx);
   outB.resize(idx);
   outTrace.clear();

   VCell last = Max::seed();
   for (size_t s = 0; s < IPStatesSize; ++s)
   {
      Max::acc(last, 
               _vt.get1(_seqA.size() - 1, _seqB.size() - 1, IPStates[s]) *
               _tm.prTrans(IPStates[s], S_End),
               Trace(IPStates[s], _seqA.size() - 1, _seqB.size() - 1));
   }
      
   Trace t = last._trace;
   --idx; 

   while (t._s != S_Start && t._s != S_Max)
   {
      outTrace.push_front(t);

      if (t.isIPState())
      {
         size_t i = t._x1;
         size_t p = t._x2;

         switch(t._s)
         {
         case S_Dx : 
            outA[idx] = DNASubModel::DNAToChar(_seqA[i]);
            outB[idx] = DNASubModel::gap();
            --idx;
            break;
         case S_Ix :
            outA[idx] = DNASubModel::gap();
            outB[idx] = DNASubModel::DNAToChar(_seqB[p]);
            --idx;
            break;
         case S_Mx :
            outA[idx] = DNASubModel::DNAToChar(_seqA[i]);
            outB[idx] = DNASubModel::DNAToChar(_seqB[p]);
            --idx;
            break;
         default :
            break;
         }
         
         t = _vt.get1(i, p, t._s)._trace;
      }
      else if (t.isIJKPState())
      {
         size_t i = t._x1;
         size_t j = t._x2;
         size_t k = t._x3;
         size_t p = t._x4;

         int spanIK = k - i;
         int gapIJ = j - i;
         int lenJK = k - j + 1;
         
         assert(spanIK > 0 && gapIJ > 0 && lenJK > 0);

         outA[idx] = DNASubModel::DNAToChar(_seqA[k]);
         outB[idx] = DNASubModel::gap();
         
         switch(t._s)
         {
         case S_MDx :
         case S_MDyM :
            outA[idx - spanIK] = DNASubModel::DNAToChar(_seqA[i]);
            outB[idx - spanIK] = DNASubModel::DNAToChar(_seqB[p]);
            break;
         default:
            assert(false);
            break;
         }
         
         if (lenJK == 1)
         {
            idx -= spanIK;
         }
         --idx;

         t = _vt.get2(i, j, k, p, t._s)._trace;
      }
      else if (t.isIPQRState())
      {
         size_t i = t._x1;
         size_t p = t._x2;
         size_t q = t._x3;
         size_t r = t._x4;

         int spanPR = r - p;
         int gapPQ = q - p;
         int lenQR = r - q + 1;
         
         assert(spanPR > 0 && gapPQ > 0 && lenQR > 0);

         outA[idx] = DNASubModel::gap();
         outB[idx] = DNASubModel::DNAToChar(_seqB[r]);
         
         switch(t._s)
         {
         case S_MIx :
         case S_MIzM :
            outA[idx - spanPR] = DNASubModel::DNAToChar(_seqA[i]);
            outB[idx - spanPR] = DNASubModel::DNAToChar(_seqB[p]);
            break;
         default:
            assert(false);
            break;
         }
         
         if (lenQR == 1)
         {
            idx -= spanPR;
         }
         --idx;
         t = _vt.get3(i, p, q, r, t._s)._trace;
      }
      else
      {
//         assert(idx == 1);
//         break;
      }
   }

   outA.erase(outA.begin(), outA.begin() + outA.find_first_not_of((char)0));
   outB.erase(outB.begin(), outB.begin() + outB.find_first_not_of((char)0));
}

// should be copy-pasted from forward (until generic version that can
// be reused is made) (except for VCell constructors in init()!)
void FTAGModel::viterbiInit()
{
   _vt.resize(_seqA.size(), _seqB.size(), _winA, _winB);
   _vt.reset(Sum::seed());

   _vt.set1(0, 0, S_Start, VCell(1.));
}

void FTAGModel::viterbiStepIP(size_t i, size_t p)
{
   // New "Single" Insert/Match/Delete
   if (i > 0)
   {
      _vt.set1(i, p, S_Dx, _em.prFullD(_seqA[i]) * 
               _vt.all1f(i-1, p, S_Dx, XStates, XStatesSize));
   }
   if (p > 0)
   {
      _vt.set1(i, p, S_Ix, _em.prFullI(_seqB[p]) *
               _vt.all1f(i, p-1, S_Ix, XStates, XStatesSize));
   }
   if (i > 0 && p > 0)
   {
      _vt.set1(i, p, S_Mx, _em.prFullM(_seqA[i], _seqB[p]) *
               _vt.all1f(i-1, p-1, S_Mx, XStates, XStatesSize));
   }
   
   // Close context and revert back into IP "Single" state
   // (no emission for this type of transition)
   if (i > 1)
   {
      _vt.set1(i, p, S_MDx, 
               _vt.all1fFrom2(i, p, MDxStates, MDxStatesSize));
   }
   if (p > 1)
   {
      _vt.set1(i, p, S_MIx, 
               _vt.all1fFrom3(i, p, MIxStates, MIxStatesSize));
   }
}

void FTAGModel::viterbiStepIJKP(size_t i, size_t j, size_t k, size_t p)
{
   assert(i < j && j <= k);

   // Open new context ("outside transition")
   if (j == k && i > 0 && p > 0)
   {
      _vt.set2(i, k, k, p, S_MDx, _em.prFullMD(_seqA[i], _seqA[k], _seqB[p]) *
               _vt.all1f(i-1, p-1, S_MDx, XStates, XStatesSize));
   }

   // Fill in existing context ("inside transition")
   if (i > 0 && k > 0 && k > j && p > 0)
   {
      _vt.set2(i, j, k, p, S_MDyM, _em.prFullMD(_seqA[i], _seqA[k], _seqB[p]) *
               _vt.all2f(i-1, j, k-1, p-1, S_MDyM, IJKPStates, IJKPStatesSize));
   }
}

void FTAGModel::viterbiStepIPQR(size_t i, size_t p, size_t q, size_t r)
{
   assert(p < q && q <= r); 

   // Open new context ("outside transition")
   if (q == r && i > 0 && p > 0)
   {
      _vt.set3(i, p, r, r, S_MIx, _em.prFullMI(_seqA[i], _seqB[p], _seqB[r]) *
               _vt.all1f(i-1, p-1, S_MIx, XStates, XStatesSize));
   }

   // Fill in existing context ("inside transition")
   if (i > 0 && p > 0 && r > 0 && r > q)
   {
      _vt.set3(i, p, q, r, S_MIzM, _em.prFullMI(_seqA[i], _seqB[p], _seqB[r]) *
               _vt.all3f(i-1, p-1, q, r-1, S_MIzM, IPQRStates, IPQRStatesSize));
   }
}
