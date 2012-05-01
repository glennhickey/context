//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGGEN_H
#define _FTAGGEN_H

#include "ftagstates.h"
#include "transitionmodel.h"
#include "contextmodel.h"
#include "outputqueue.h"

class FTAGGen
{
public:

   typedef std::pair<DNA, DNA> DNAPair;
   struct DNATriple {
      DNA first;
      DNA second;
      DNA third;
      DNATriple() {}
      DNATriple(DNA f, DNA s, DNA t) : first(f), second(s), third(t) {}
   };

public:
   
   void setEmissionModel(const ContextModel& eModel);
   void setTransitionModel(const TransitionModel& tModel);

   void genAlignment(std::string& outA, std::string& outB, 
                     std::deque<FTAGTableStructs::Trace>& outTrace,
                     bool keepGaps = true, size_t maxSize = 100);

protected:

   // import DNASubModel static functions so they're easier to use
   // todo: put these guys in a namespace in order to use "using"!!!
   static char DNAToChar(DNA dna) {
      return DNASubModel::DNAToChar(dna);
   }
   static char gap() {
      return DNASubModel::gap();
   }

   FTAGStates::State randTransition(FTAGStates::State current,
                                    bool outside) const;

   DNA randEmitSingle(
      double (ContextModel::*prSingle)(DNA) const) const;

   DNAPair randEmitPair(
      double (ContextModel::*prDouble)(DNA, DNA) const) const;

   DNATriple randEmitTriple(
      double (ContextModel::*prTriple)(DNA, DNA, DNA) const) const;

protected:
   
   TransitionModel _tm; // transition model
   ContextModel _em; // emission model
   OutputQueue _oq;
};

#endif
