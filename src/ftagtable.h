//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGTABLE_H
#define _FTAGTABLE_H

#include <cstdlib>
#include <vector>
#include <iostream>

#include "ftagstates.h"
#include "ftagtablestructs.h"
#include "transitionmodel.h"

// T : data type U: object with acc fucntion defined (max or sum)
template<typename T, typename U>
class FTAGTable
{
public:
   typedef FTAGStates::State State;
   
   FTAGTable();
   ~FTAGTable();

   void setTransitionModel(const TransitionModel* tm);
   void resize(size_t lenA, size_t lenB, size_t winA, size_t winB);
   void reset(T val);

   // all indices are global.  ex: jk are not relative to i coming in.
   // so the get/set functions must take care of global / local offsets

   T get1(size_t i, size_t p, State s) const;
   void set1(size_t i, size_t p, State s, T val);

   T get2(size_t i, size_t j, size_t k, size_t p, State s) const;
   void set2(size_t i, size_t j, size_t k, size_t p, State s, T val);

   T get3(size_t i, size_t p, size_t q, size_t r, State s) const;
   void set3(size_t i, size_t p, size_t q, size_t r, State s, T val);
   
   T all1f(size_t i, size_t p, State to,
           const State states[], size_t statesSize) const;

   T all2f(size_t i, size_t j, size_t k, size_t p, State to,
           const State states[], size_t statesSize) const;

   T all3f(size_t i, size_t p, size_t q, size_t r, State to,
           const State states[], size_t statesSize) const;

   T all1fFrom2(size_t i, size_t p, 
                const State states[], size_t statesSize) const;
   
   T all1fFrom3(size_t i, size_t p, 
                const State states[], size_t statesSize) const;

protected:

   void initMaps();

   T* _tab1;
   T* _tab2;
   T* _tab3;
   size_t _lenA;
   size_t _lenB;
   size_t _winA;
   size_t _winB;
   size_t _c1[3];
   size_t _c2[5];
   size_t _c3[5];
   size_t _map1[FTAGStates::S_Max];
   size_t _map2[FTAGStates::S_Max];
   size_t _map3[FTAGStates::S_Max];

   const TransitionModel* _tm;

private:

   FTAGTable(const FTAGTable&);
   FTAGTable& operator=(const FTAGTable&);
};


#include "ftagtable_impl.h"

#endif
