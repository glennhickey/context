//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <string>
#include <vector>
#include <list>
#include <cstdlib>
#include <iostream>

#include "ftaggen.h"
#include "ftagtablestructs.h"

using namespace std;
using namespace FTAGStates;
using namespace FTAGTableStructs;

void FTAGGen::setEmissionModel(const ContextModel& eModel)
{
   _em = eModel;
}

void FTAGGen::setTransitionModel(const TransitionModel& tModel)
{
   _tm = tModel;
}


void FTAGGen::genAlignment(string& outA, string& outB, deque<Trace>& outTrace,
                           bool keepGaps, size_t maxSize)
{
   _oq.reset();
   State state = randTransition(S_Start, true);
   DNA val1;
   DNAPair val2;
   DNATriple val3;
   for (size_t size = 0; state != S_End && size < maxSize; ++size)
   {
      switch(state)
      {
      case S_Dx :
         assert(!_oq.isOpen());
         val1 = randEmitSingle(&ContextModel::prFullD);
         _oq.addSingle(state, DNAToChar(val1), gap());
         break;
      case S_Ix :
         assert(!_oq.isOpen());
         val1 = randEmitSingle(&ContextModel::prFullI);
         _oq.addSingle(state, gap(), DNAToChar(val1));
         break;
      case S_Mx :
         assert(!_oq.isOpen());
         val2 = randEmitPair(&ContextModel::prFullM);
         _oq.addSingle(state, DNAToChar(val2.first), DNAToChar(val2.second));
         break;
      case S_DDx : 
         assert(!_oq.isOpen());
         _oq.openContext();
      case S_MDx :
         assert(!_oq.isOpen());
         _oq.openContext();
      case S_MDyM :
         val3 = randEmitTriple(&ContextModel::prFullMD);
         _oq.addDouble(state, DNAToChar(val3.first), DNAToChar(val3.second), 
                       DNAToChar(val3.third), gap());
         break;
      case S_MIx:
         assert(!_oq.isOpen());
         _oq.openContext();
      case S_MIzM:
         val3 = randEmitTriple(&ContextModel::prFullMI);
         _oq.addDouble(state, DNAToChar(val3.first), gap(),
                       DNAToChar(val3.second), DNAToChar(val3.third));
         break;
      default:
         assert(state == S_End);
      }

      State next = randTransition(state, !_oq.isOpen());
      if (next == S_Close)
      {
         assert(_oq.isOpen());
         _oq.closeContext();
         if (findState(state, MDxStates, MDxStatesSize) == true)
         {
            next = S_MDx;
         }
         else if (findState(state, MIxStates, MIxStatesSize) == true)
         {
            next = S_MIx;
         }
         else
         {
            assert(false);
         }
         // effectively at a silent state, so jump one ahead.
         next = randTransition(next, true);
      }
      state = next;
   }
   
   pair<string, string> alignment = _oq.getAlignment(keepGaps);
   outA = alignment.first;
   outB = alignment.second;
   outTrace = _oq.getTrace();
}

State FTAGGen::randTransition(State current, bool outside) const
{
   double rVal = drand48();
   double rSum = 0.;

   size_t size = outside ? OutsideToStatesSize : InsideToStatesSize;
   const State* states = outside ? OutsideToStates : InsideToStates;

   for (size_t s = 0; s < size; ++s)
   {
      rSum += _tm.prTrans(current, states[s]);
      if (rVal < rSum)
      {
         return states[s];
      }
   }
   
   return S_Max;
}

DNA FTAGGen::randEmitSingle(
   double (ContextModel::*prSingle)(DNA) const) const
{
   double rVal = drand48();
   double rSum = 0.;
   
   for (size_t s = 0; s < DNA_Max; ++s)
   {
      rSum += (_em.*prSingle)((DNA)s);
      if (rVal < rSum)
      {
         return (DNA)s;
      }
   }
   throw (string("Emitting DNA_Max"));
   return DNA_Max;
}

FTAGGen::DNAPair FTAGGen::randEmitPair(
   double (ContextModel::*prDouble)(DNA, DNA) const) const
{
   double rVal = drand48();
   double rSum = 0.;
   
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t p = 0; p < DNA_Max; ++p)
      {
         rSum += (_em.*prDouble)((DNA)i, (DNA)p);
         if (rVal < rSum)
         {
            return DNAPair((DNA)i, (DNA)p);
         }
      }
   }
   throw (string("Emitting DNA_Max,DNA_Max"));
   return DNAPair(DNA_Max, DNA_Max);
}

FTAGGen::DNATriple FTAGGen::randEmitTriple(
   double (ContextModel::*prTriple)(DNA, DNA, DNA) const) const
{
   double rVal = drand48();
   double rSum = 0.;
   
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t k = 0; k < DNA_Max; ++k)
      {
         for (size_t p = 0; p < DNA_Max; ++p)
         {
            rSum += (_em.*prTriple)((DNA)i, (DNA)k, (DNA)p);
            if (rVal < rSum)
            {
               return DNATriple((DNA)i, (DNA)k, (DNA)p);
            }
         }
      }
   }
   throw (string("Emitting DNA_Max,DNA_Max,DNA_Max"));
   return DNATriple(DNA_Max, DNA_Max, DNA_Max);
}
