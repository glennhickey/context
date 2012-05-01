//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _FTAGTABLESTRUCTS_H
#define _FTAGTABLESTRUCTS_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cassert>

#include "ftagstates.h"
#include "ftagtablestructs.h"

namespace FTAGTableStructs {

   struct Trace 
   {
      Trace(FTAGStates::State s = FTAGStates::S_Max, 
            unsigned short x1 = 0, unsigned short x2 = 0, 
            unsigned short x3 = 0, unsigned short x4 = 0) : 
      _s(s), _x1(x1), _x2(x2), _x3(x3), _x4(x4) {}
      FTAGStates::State _s;
      unsigned short _x1;
      unsigned short _x2;
      unsigned short _x3;
      unsigned short _x4;
      bool isIPState() const 
      {
         return _x3 == 0 && _x4 == 0 && 
             FTAGStates::findState(_s, FTAGStates::IPStates, 
                                   FTAGStates::IPStatesSize);
      }
      bool isIJKPState() const
      {
         return (_x3 > 0 || _x4 > 0) &&
            FTAGStates::findState(_s, FTAGStates::IJKPStates,
                                  FTAGStates::IJKPStatesSize);
      }
      bool isIPQRState() const
      {
         return (_x3 > 0 || _x4 > 0) &&
            FTAGStates::findState(_s, FTAGStates::IPQRStates,
                                  FTAGStates::IPQRStatesSize);
      }
   };
      
   struct VCell 
   {
      VCell(double v = 0,Trace t = Trace()) 
      : _val(v), _trace(t) {}
      double _val;
      Trace _trace;
   };

   inline VCell operator*(double left, VCell right) 
   {
      return VCell(left * right._val, right._trace);
   }
   inline VCell operator*(VCell left, double right)
   {
      return VCell(left._val * right, left._trace);
   }
      
   struct Sum {
      static void acc(double& current, double other, const Trace& trace) 
      {
         current += other;
      }
      static double seed() 
      {
         return 0.;
      }
   };
   
   struct Max 
   {
      static void acc(VCell& current, VCell other, const Trace& trace) 
      {
         if (log(other._val) > log(current._val) + 1e-10)
         {
            current._val = other._val;
            current._trace = trace;
         }
      }
      static VCell seed()
      {
         return VCell(0., FTAGStates::S_Max);
      }
   };
}

inline bool isnan(FTAGTableStructs::VCell cell) 
{
   return std::isnan(cell._val);
}

inline bool isinf(FTAGTableStructs::VCell cell)
{
   return std::isinf(cell._val);
}

inline std::ostream& operator << (std::ostream& os, 
                                  const FTAGTableStructs::Trace& t)
{
   os << t._s;
   if (t.isIPState() == true)
   {
      os << " <i=" << t._x1 << ", p=" << t._x2 << ">";
   }
   else if (t.isIJKPState() == true)
   {
      os << " <i=" << t._x1 << ", j=" << t._x2 << ", k=" << t._x3 
         << ", p=" << t._x4 << ">";
   }
   else if (t.isIPQRState() == true) 
   {
      os << " <i=" << t._x1 << ", p=" << t._x2 << ", q=" << t._x3 
         << ", r=" << t._x4 << ">";
   }
   else if (t._s != FTAGStates::S_Max || t._s != FTAGStates::S_End ||
            t._s != FTAGStates::S_Close)
   {
      assert(false);
   }
   return os;
}

inline std::ostream& operator << (std::ostream& os, 
                                  const FTAGTableStructs::VCell& v)
{
   os << v._val;
   return os;
}

#endif
