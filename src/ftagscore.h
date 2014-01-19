//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGSCORE_H
#define _FTAGSCORE_H

#include <cstdlib>
#include <string>
#include <vector>
#include <deque>

#include "ftagmodel.h"

/** This is a class that handles core functionality related to scoring a fixed 
alignment.  IE, the alignment is already known, and we just want to scan it
over to see what it's viterbi score is */

class FTAGScore
{
public:
   FTAGScore(const FTAGModel* model = NULL);
   virtual ~FTAGScore();

   // model must be set before anything gets done
   void setModel(const FTAGModel& model);
   
   // return the (viterbi) log score of a pairwise alignment. performed in 
   // linear time because each indel can be resolved independently and is 
   // either context or not. note that the alignment structure here can 
   // be obtained from an AXTPair object
   double logViterbiScore(const std::pair<std::string, std::string>& alignment);

protected:

   // compute score of a single match column (note no transition taken into
   // account
   void scoreMatch(size_t pos);

   // given the first position of a gap in the indel, scan right computing the
   // score
   void scoreIndel(size_t pos);

   // get the id of a single column.  Either a match, insert, or delete.  Note
   // because we are just looking at one column, we don't consider context
   // events -- just want to know if / where the gap is. 
   FTAGStates::State getColID(size_t pos);

   // convenience for filling out context tables below.  want to be able to
   // quickly get context version of a state if present;
   FTAGStates::State getContext(FTAGStates::State state) const;

   // convenience for copying arrays of size two which get used a fair amount
   template<typename T> void copy2(T a[2], T b[2]) const;
   template<typename T> void init2(T a[2], T b0, T b1) const;
   
protected:

   const FTAGModel* _model;
   const ContextModel* _emModel;
   const TransitionModel* _transModel;

   // length of alignment for convenience
   size_t _N;

   // keep pointer to alignment
   const std::pair<std::string, std::string>* _alignment;

   // the current detected state along the alingment;
   // we always keep track of a non-context (0) and context state (1)
   // in the case of match, we set them both to match
   FTAGStates::State _state[2];

   // if the above state is an indel, we also store its size
   size_t _indelSize;

   // the previous state (before above)
   FTAGStates::State _lastState[2];

   // viterbi dynamic programming table.  stores value for 
   // 0 = non-context, 1 = context.  
   double _table[2];
   // like above but for previous column
   double _lastTable[2];
   
};

inline FTAGStates::State FTAGScore::getContext(FTAGStates::State state) const
{
  if (state == FTAGStates::S_Dx)
  {
    return FTAGStates::S_MDx;
  }
  else if (state == FTAGStates::S_Ix)
  {
    return FTAGStates::S_MIx;
  }
  return state;
}

template<typename T>
inline void FTAGScore::copy2(T a[2], T b[2]) const
{
  a[0] = b[0];
  a[1] = b[1];
}

template<typename T>
inline void FTAGScore::init2(T a[2], T b0, T b1) const
{
  a[0] = b0;
  a[1] = b1;
}

#endif
