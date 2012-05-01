//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGSTATES_H
#define _FTAGSTATES_H

#include <cstdlib>
#include <iostream>

// helper constants to manage states and ranges of states.
// if one enum / array is updated, chances are the rest will have
// to be synchronized by hand!


// **OCT25:  Stripped out all indel-within-context states.  revert to older
//           version (of everything) to get them back.    

namespace FTAGStates
{
   enum State {
      S_Dx = 0, S_Ix, S_MDx, S_MIx, S_MDyM, S_MIzM, S_Mx,
      S_Start, S_Close, S_End, S_Max,
      S_DDx, S_IDx, S_DIx, S_IIx,
   };
   
   // redundant?  maybe delete and replace all usages with IPStates....
   static const State XStates[] = {
      S_Dx, S_Ix, S_MDx, S_MIx, S_Mx, S_Start
   };
   static const size_t XStatesSize = 6;

   static const State YStates[] = {
      S_MDyM
   };
   static const size_t YSTatesSize = 1;


   static const State ZStates[] = {
      S_MIzM
   };
   static const size_t ZStatesSize = 1;
   
   // "type-1" states ie no open context
   static const State IPStates[] = {
      S_Dx, S_Ix, S_MDx, S_MIx, S_Mx, S_Start
   }; 
   static const size_t IPStatesSize = 6;

   // "single" type 1 states states
   static const State IPStrictStates[] = {
      S_Dx, S_Ix, S_Mx
   };
   static const size_t IPStrictStatesSize = 3;
   
   // "type-2" states ie open context deletes
   static const State IJKPStates[] = {
      S_MDx, S_MDyM
   };
   static const size_t IJKPStatesSize = 2;

   // "type-3" states ie open context inserts
   static const State IPQRStates[] = {
      S_MIx, S_MIzM
   };
   static const size_t IPQRStatesSize = 2;
   
   // states that can be "closed" to make a DDx(k,r) state
   static const State DDxStates [] = {
   };
   static const size_t DDxStatesSize = 0;

   // states that can be "closed" to make a IDx(k,r) state
   static const State IDxStates [] = {
   };
   static const size_t IDxStatesSize = 0;

   // states that can be "closed" to make a MDx(k,r) state
   static const State MDxStates [] = {
      S_MDx, S_MDyM
   };
   static const size_t MDxStatesSize = 2;

   // states that can be "closed" to make a DIx(k,r) state
   static const State DIxStates [] = {
   };
   static const size_t DIxStatesSize = 0;

   // states that can be "closed" to make a IIx(k,r) state
   static const State IIxStates [] = {
   };
   static const size_t IIxStatesSize = 0;

   // states that can be "closed" to make a MIx(k,r) state
   static const State MIxStates [] = {
      S_MIx, S_MIzM
   };
   static const size_t MIxStatesSize = 2;

   // states that can be destinations in outside transitions
   static const State OutsideToStates [] = {
      S_Dx, S_Ix, S_MDx, S_MIx, S_Mx, S_Start, S_End
   };
   static const size_t OutsideToStatesSize = 7;

   // states that can be destinations in inside transitions
   static const State InsideToStates [] = {
      S_MDyM, S_MIzM, S_Close
   };
   static const size_t InsideToStatesSize = 3;

   inline bool findState(State s, const State* states, size_t size)
   {
      for (size_t i = 0; i < size; ++i)
      {
         if (states[i] == s)
         {
            return true;
         }
      }
      return false;
   }   
}

inline std::ostream& operator<<(std::ostream& os, FTAGStates::State s)
{
   switch(s) 
   {
   case FTAGStates::S_Dx: os << "Dx"; break;
   case FTAGStates::S_Ix: os << "Ix"; break; 
   case FTAGStates::S_Mx: os << "Mx"; break; 
   case FTAGStates::S_MDx: os << "MDx"; break; 
   case FTAGStates::S_MIx: os << "MIx"; break;
   case FTAGStates::S_MDyM: os << "MDyM"; break;
   case FTAGStates::S_MIzM: os << "MIzM"; break;
   case FTAGStates::S_Start: os << "Start"; break; 
   case FTAGStates::S_Close: os << "Close"; break; 
   case FTAGStates::S_End: os << "End"; break;
   case FTAGStates::S_Max: os <<"Max"; break;
   default:
      os << "wut?";
   }
   return os;
}

#endif
