//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cassert>
#include <iostream>
#include <cmath>

template<typename T, typename U>
inline 
FTAGTable<T, U>::FTAGTable() : 
   _tab1(NULL), _tab2(NULL), _tab3(NULL), 
   _lenA(0), _lenB(0),
   _winA(0), _winB(0),
   _tm(NULL)
{
   initMaps();
}

template<typename T, typename U>
inline 
FTAGTable<T, U>::~FTAGTable() 
{
   delete [] _tab1;
   delete [] _tab2;
   delete [] _tab3;
}

template<typename T, typename U>
inline 
void FTAGTable<T, U>::setTransitionModel(const TransitionModel* tm)
{
   _tm = tm;
}

template<typename T, typename U>
inline 
void FTAGTable<T, U>::resize(size_t lenA, size_t lenB, 
                             size_t winA, size_t winB)
{
   if (lenA != _lenA || lenB != _lenB || 
       winA != _winA || winB != _winB)
   {
      _lenA = lenA;
      _lenB = lenB;
      _winA = winA;
      _winB = winB;
      _c1[2] = FTAGStates::IPStatesSize;
      _c1[1] = _lenB * FTAGStates::IPStatesSize;
      _c1[0] = _lenA * _lenB * FTAGStates::IPStatesSize;
      _c2[4] = FTAGStates::IJKPStatesSize;
      _c2[3] = _winA * FTAGStates::IJKPStatesSize;
      _c2[2] = _winA * _winA * FTAGStates::IJKPStatesSize;
      _c2[1] = _lenB * _winA * _winA * FTAGStates::IJKPStatesSize;
      _c2[0] = _lenA * _lenB * _winA * _winA * FTAGStates::IJKPStatesSize;
      _c3[4] = FTAGStates::IPQRStatesSize;
      _c3[3] = _winB * FTAGStates::IPQRStatesSize;
      _c3[2] = _winB * _winB * FTAGStates::IPQRStatesSize;
      _c3[1] = _lenB * _winB * _winB * FTAGStates::IPQRStatesSize;
      _c3[0] = _lenA * _lenB * _winB * _winB * FTAGStates::IPQRStatesSize;

      try{
         // table 1
         delete [] _tab1;
         _tab1 = new T[_c1[0]];
      
         // table 2
         delete [] _tab2;
         _tab2 = new T[_c2[0]];
      
         // table 3
         delete [] _tab3;
         _tab3 = new T[_c3[0]];
      }
      catch (...)
      {
         std::cerr << "ALLOCATION FAIL of \n" 
                   << "T1 : " << (sizeof(T) * _c1[0]) / (1 << 20)+1 << " MB\n"
                   << "T2 : " << (sizeof(T) * _c2[0]) / (1 << 20) << " MB\n"
                   << "T3 : " << (sizeof(T) * _c3[0]) / (1 << 20) << " MB\n" 
                   << std::endl;
         exit(1);
      }
   }
}

template<typename T, typename U>
inline 
void FTAGTable<T, U>::reset(T val) 
{
   size_t i;

   // table 1
   for (i = 0; i < _c1[0]; ++i)
   {
      _tab1[i] = val;
   }
   
   // table 2
   for (i = 0; i < _c2[0]; ++i)
   {
      _tab2[i] = val;
   }

   // table 3
   for (i = 0; i < _c3[0]; ++i)
   {
      _tab3[i] = val;
   }
}

// table 1
template<typename T, typename U>
inline 
T FTAGTable<T, U>::get1(size_t i, size_t p, State s) const
{
   assert(_map1[s] < FTAGStates::IPStatesSize); 
   return _tab1[i * _c1[1] + p * _c1[2] + _map1[s]];
}
                 
template<typename T, typename U>
inline 
void FTAGTable<T, U>::set1(size_t i, size_t p, State s, T val)
{
   assert(_map1[s] < FTAGStates::IPStatesSize); 
   assert(!isnan(val) && !isinf(val));
   _tab1[i * _c1[1] + p * _c1[2] + _map1[s]] = val;
}

// table 2
template<typename T, typename U>
inline 
T FTAGTable<T, U>::get2(size_t i, size_t j, size_t k, size_t p, 
                        State s) const
{
   assert(_map2[s] < FTAGStates::IJKPStatesSize); 
   j -= i + 1;
   k -= i + 1;
   assert(j < _winA && k < _winA);
   return _tab2[i * _c2[1] + p * _c2[2] + j * _c2[3] + k * _c2[4] + _map2[s]];
}
                 
template<typename T, typename U>
inline 
void FTAGTable<T, U>::set2(size_t i, size_t j, size_t k, size_t p, 
                           State s, T val)
{
   assert(_map2[s] < FTAGStates::IJKPStatesSize); 
   assert(!isnan(val) && !isinf(val));
   j -= i + 1;
   k -= i + 1;
   assert(j < _winA && k < _winA);
   _tab2[i * _c2[1] + p * _c2[2] + j * _c2[3] + k * _c2[4] + _map2[s]] = val;
}

// table 3
template<typename T, typename U>
inline 
T FTAGTable<T, U>::get3(size_t i, size_t p, size_t q, size_t r, 
                        State s) const
{
   assert(_map3[s] < FTAGStates::IPQRStatesSize); 
   q -= p + 1;
   r -= p + 1;
   assert(q < _winB && r < _winB);
   return _tab3[i * _c3[1] + p * _c3[2] + q * _c3[3] + r * _c3[4] + _map3[s]];
}
                 
template<typename T, typename U>
inline 
void FTAGTable<T, U>::set3(size_t i, size_t p, size_t q, size_t r, 
                           State s, T val)
{
   assert(_map3[s] < FTAGStates::IPQRStatesSize); 
   assert(!isnan(val) && !isinf(val));
   q -= p + 1;
   r -= p + 1;
   assert(q < _winB && r < _winB);
   _tab3[i * _c3[1] + p * _c3[2] + q * _c3[3] + r * _c3[4] + _map3[s]] = val;
}

template<typename T, typename U>
inline 
T FTAGTable<T, U>::all1f(size_t i, size_t p, State to, const State states[],
                         size_t statesSize) const
{
   T all = U::seed(); 
   for (size_t s = 0; s < statesSize; ++s)
   {
      U::acc(all, get1(i, p, states[s]) * _tm->prTrans(states[s], to), 
             FTAGTableStructs::Trace(states[s], i, p));
   }
   return all;
}

template<typename T, typename U>
inline 
T FTAGTable<T, U>::all2f(size_t i, size_t j, size_t k, size_t p, State to, 
                         const State states[], size_t statesSize) const
{
   T all = U::seed();
   for (size_t s = 0; s < statesSize; ++s)
   {
      U::acc(all, get2(i, j, k, p, states[s]) * _tm->prTrans(states[s], to),
             FTAGTableStructs::Trace(states[s], i, j, k, p));
   }
   return all;
}

template<typename T, typename U>
inline 
T FTAGTable<T, U>::all3f(size_t i, size_t p, size_t q, size_t r, State to, 
                         const State states[], size_t statesSize) const
{
   T all = U::seed();
   for (size_t s = 0; s < statesSize; ++s)
   {
      U::acc(all, get3(i, p, q, r, states[s]) * _tm->prTrans(states[s], to),
             FTAGTableStructs::Trace(states[s], i, p, q, r));
   }
   return all;
}

template<typename T, typename U>
inline 
T FTAGTable<T, U>::all1fFrom2(size_t i, size_t p, const State states[],
                              size_t statesSize) const
{
   T all = U::seed();
   size_t u = i > _winA ? i - _winA : 0;
   for (; u < i; ++u)
   {
      for (size_t s = 0; s < statesSize; ++s)
      {
         U::acc(all, get2(u, u+1, i, p, states[s]) * 
                _tm->prTrans(states[s], FTAGStates::S_Close), 
                FTAGTableStructs::Trace(states[s], u, u+1, i, p));
      }
   }
   return all;
}

template<typename T, typename U>
inline 
T FTAGTable<T, U>::all1fFrom3(size_t i, size_t p, const State states[],
                              size_t statesSize) const
{
   T all = U::seed();
   size_t v = p > _winB ? p - _winB : 0;
   for (; v < p; ++v)
   {
      for (size_t s = 0; s < statesSize; ++s)
      {
         U::acc(all, get3(i, v, v+1, p, states[s]) * 
                _tm->prTrans(states[s], FTAGStates::S_Close), 
                FTAGTableStructs::Trace(states[s], i, v, v+1, p));
      }
   }
   return all;
}

template<typename T, typename U>
inline
void FTAGTable<T, U>::initMaps()
{
   size_t i;
   for (i = 0; i < FTAGStates::S_Max; ++i)
   {
      _map1[i] = FTAGStates::S_Max;
      _map2[i] = FTAGStates::S_Max;
      _map3[i] = FTAGStates::S_Max;
   }
   
   for (i = 0; i < FTAGStates::IPStatesSize; ++i)
   {
      _map1[FTAGStates::IPStates[i]] = i; 
   }
   
   for (i = 0; i < FTAGStates::IJKPStatesSize; ++i)
   {
      _map2[FTAGStates::IJKPStates[i]] = i; 
   }

   for (i = 0; i < FTAGStates::IPQRStatesSize; ++i)
   {
      _map3[FTAGStates::IPQRStates[i]] = i; 
   }
}



