//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <vector>
#include <string>
#include <istream>
#include <ostream>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <limits>
#include <sstream>
#include <cmath>
#include "ftagparams.h"

using namespace std;

FTAGParams::FTAGParams()
{
   for (size_t i = 0; i < Par_Max; ++i)
   {
      _params[i] = 0.;
      _fixed[i] = false;
      _clamps[0][i] = 0.;
      _clamps[1][i] = numeric_limits<double>::max();
      _phases[0][i] = 0;
      _phases[1][i] = numeric_limits<size_t>::max();
   }
   
   _clamps[1][Par_PE] = 1.;
   _clamps[1][Par_PCD] = 1.;
   _clamps[1][Par_PCI] = 1.;
   _clamps[1][Par_KD] = 1.;
   _clamps[1][Par_KI] = 1.;
   _clamps[1][Par_P0] = 1;
   _clamps[1][Par_P1] = 1;
   _clamps[1][Par_P2] = 1;
   _clamps[1][Par_P3] = 1;
   _clamps[1][Par_Q0] = 1;
   _clamps[1][Par_Q1] = 1;
   _clamps[1][Par_Q2] = 1;
   _clamps[1][Par_Q3] = 1;/*
   _clamps[0][Par_P0] = .1;
   _clamps[0][Par_P1] = .1;
   _clamps[0][Par_P2] = .1;
   _clamps[0][Par_P3] = .1;
   _clamps[0][Par_Q0] = .1;
   _clamps[0][Par_Q1] = .1;
   _clamps[0][Par_Q2] = .1;
   _clamps[0][Par_Q3] = .1;*/
   _isSymmetric = true;
   _isUniGap = false;
   _singleF84 = false;
   _doubleF84 = false;
   _isPLinkedGC = false;
   _isPLinkedGA = false;
   _isQLinkedGC = false;
   _isQLinkedGA = false;
   _params[Par_P0] = 0.25;
   _params[Par_P1] = 0.25;
   _params[Par_P2] = 0.25;
   _params[Par_P3] = 0.25;
   _params[Par_Q0] = 0.25;
   _params[Par_Q1] = 0.25;
   _params[Par_Q2] = 0.25;
   _params[Par_Q3] = 0.25;
}

FTAGParams::FTAGParams(const FTAGParams& params)
{
   for (size_t i = 0; i < Par_Max; ++i)
   {
      _params[i] = params._params[i];
      _fixed[i] = params._fixed[i];
      _clamps[0][i] = params._clamps[0][i];
      _clamps[1][i] = params._clamps[1][i];
      _phases[0][i] = params._phases[0][i];
      _phases[1][i] = params._phases[1][i];
   }
   _isSymmetric = params._isSymmetric;
   _isUniGap = params._isUniGap;
   _singleF84 = params._singleF84;
   _doubleF84 = params._doubleF84;
   _isPLinkedGC = params._isPLinkedGC;
   _isPLinkedGA = params._isPLinkedGA;
   _isQLinkedGC = params._isQLinkedGC;
   _isQLinkedGA = params._isQLinkedGA;
}

FTAGParams::~FTAGParams()
{
   
}

FTAGParams& FTAGParams::operator=(const FTAGParams& params)
{
   for (size_t i = 0; i < Par_Max; ++i)
   {
      _params[i] = params._params[i];
      _fixed[i] = params._fixed[i];
      _clamps[0][i] = params._clamps[0][i];
      _clamps[1][i] = params._clamps[1][i];
      _phases[0][i] = params._phases[0][i];
      _phases[1][i] = params._phases[1][i];
   }
   _isSymmetric = params._isSymmetric;
   _isUniGap = params._isUniGap;
   _singleF84 = params._singleF84;
   _doubleF84 = params._doubleF84;
   _isPLinkedGC = params._isPLinkedGC;
   _isPLinkedGA = params._isPLinkedGA;
   _isQLinkedGC = params._isQLinkedGC;
   _isQLinkedGA = params._isQLinkedGA;
   return *this;
}

FTAGParams FTAGParams::operator-(const FTAGParams& params) const
{
   FTAGParams diff(*this);
   assert(params._isSymmetric == _isSymmetric);
   assert(params._singleF84 == _singleF84);
   assert(params._doubleF84 == _doubleF84);
   for (size_t i = 0; i < Par_Max; ++i)
   {
      diff._params[i] = _params[i] - params._params[i];
      assert(_fixed[i] == params._fixed[i]);
      diff._fixed[i] = _fixed[i];
   }
   return diff;
}
   
FTAGParams FTAGParams::operator/(const FTAGParams& params) const
{
   FTAGParams quot(*this);
   assert(params._isSymmetric == _isSymmetric);
   assert(params._singleF84 == _singleF84);
   assert(params._doubleF84 == _doubleF84);
   assert(params._isUniGap == _isUniGap);
   for (size_t i = 0; i < Par_Max; ++i)
   {
      quot._params[i] = _params[i] / params._params[i];
      assert(_fixed[i] == params._fixed[i]);
      quot._fixed[i] = _fixed[i];
   }
   return quot;
}

FTAGParams FTAGParams::operator*(const FTAGParams& params) const
{
   FTAGParams prod(*this);
   assert(params._isSymmetric == _isSymmetric);
   assert(params._singleF84 == _singleF84);
   assert(params._doubleF84 == _doubleF84);
   assert(params._isUniGap == _isUniGap);
   for (size_t i = 0; i < Par_Max; ++i)
   {
      prod._params[i] = _params[i] * params._params[i];
      assert(_fixed[i] == params._fixed[i]);
      prod._fixed[i] = _fixed[i];
   }
   return prod;
}

bool FTAGParams::operator==(const FTAGParams& params) const
{
   // symmetric flag not tested.  should it be? 
   for (size_t i = 0; i < Par_Max; ++i)
   {
      if (_params[i] != params._params[i] ||
          _fixed[i] != params._fixed[i] ||
          _clamps[0][i] != params._clamps[0][i] ||
          _clamps[1][i] != params._clamps[1][i] ||
          _phases[0][i] != params._phases[0][i] ||
          _phases[1][i] != params._phases[1][i])
      {
         return false;
      }
   }
   return true;
}

bool FTAGParams::isDependent(FTAGParams::Param param) const 
{
   return (_isSymmetric && param == Par_KI) ||
      (_isSymmetric && param == Par_RMI) ||
      (_isUniGap && param == Par_RMD) ||
      (_isUniGap && param == Par_RMI) ||
      ((_isPLinkedGC || _isPLinkedGA) && 
       (param == Par_P1 || param == Par_P2 || param == Par_P3)) ||
      ((_isQLinkedGC || _isQLinkedGA) && 
       (param == Par_Q1 || param == Par_Q2 || param == Par_Q3));
}

FTAGParams FTAGParams::abs() const
{
   FTAGParams par = *this;
   for (size_t i = 0; i < Par_Max; ++i)
   {
      par._params[i] = fabs(par._params[i]);
   }
   return par;
}

double FTAGParams::mse(const FTAGParams& params) const
{
   double total = 0;
   double diff;
   for (size_t i = 0; i < Par_Max; ++i)
   {
      diff = _params[i] - params._params[i];
      total += diff * diff;
   }
   return total / (double) Par_Max;
}

void FTAGParams::setA(double a) 
{
   _params[Par_A] = clamp(Par_A, a);
   if (rand() % 5 == 0)
   {
      double delta = a - _params[Par_A];
      _params[Par_B] = clamp(Par_B, _params[Par_B] + delta);
   }
   if (rand() % 25 == 0)
   {
      swap(_params[Par_A], _params[Par_B]);
   }
}

void FTAGParams::setB(double b)
{
   _params[Par_B] = clamp(Par_B, b);
}

void FTAGParams::setC(double c) 
{
   _params[Par_C] = clamp(Par_C, c);
   if (rand() % 5 == 0)
   {
      double delta = c - _params[Par_C];
      _params[Par_D] = clamp(Par_D, _params[Par_D] + delta);
   }
   if (rand() % 25 == 0)
   {
      swap(_params[Par_C], _params[Par_D]);
   }

}

void FTAGParams::setD(double d)
{
   _params[Par_D] = clamp(Par_D, d);
}

void FTAGParams::setP0(double p0) 
{
   if (rand() % 50)
   {
      if (_isPLinkedGC == true)
      {
         _params[Par_P0] = clamp(Par_P0, p0) > 0.5 ? 0.5 : clamp(Par_P0, p0);
         _params[Par_P1] = 0.5 - _params[Par_P0];
         _params[Par_P2] = _params[Par_P1];
         _params[Par_P3] = _params[Par_P0];
      }
      else if (_isPLinkedGA == true)
      {
         _params[Par_P0] = clamp(Par_P0, p0) > 0.5 ? 0.5 : clamp(Par_P0, p0);
         _params[Par_P1] = 0.5 - _params[Par_P0];
         _params[Par_P2] = _params[Par_P0];
         _params[Par_P3] = _params[Par_P1];
      }
      else
      {
         double delta = clamp(Par_P0, p0) - _params[Par_P0];
         double fracDelta;
         _params[Par_P0] += delta;
         size_t sel = rand() % 4;
         
         switch(sel) {
         case 0 : 
            fracDelta = delta / 3. ;
            _params[Par_P1] = clamp(Par_P1, _params[Par_P1] - fracDelta);
            _params[Par_P2] = clamp(Par_P2, _params[Par_P2] - fracDelta);
            _params[Par_P3] = clamp(Par_P3, _params[Par_P3] - fracDelta);
            break;
         case 1: _params[Par_P1] = clamp(Par_P1, _params[Par_P1] - delta); 
            break;
         case 2: _params[Par_P2] = clamp(Par_P2, _params[Par_P2] - delta); 
            break;
         case 3: _params[Par_P3] = clamp(Par_P3, _params[Par_P3] - delta);
            break;
         default : break;
         }
         normalizeP();
      }
   }   
   else if (rand() % 2)
   {
      _params[Par_P0] = .25; _params[Par_P1] = .25;
      _params[Par_P2] = .25; _params[Par_P3] = .25;
   }
   else
   {
      randomSwapP();
   }
}

void FTAGParams::setP1(double p1) 
{
   if (!_isPLinkedGC && ! _isPLinkedGA)
   {
      double delta = clamp(Par_P1, p1) - _params[Par_P1];
      double fracDelta;
      _params[Par_P1] += delta;
      size_t sel = rand() % 4;
      
      switch(sel) {
      case 0 : 
         fracDelta = delta / 3. ;
         _params[Par_P0] = clamp(Par_P0, _params[Par_P0] - fracDelta);
         _params[Par_P2] = clamp(Par_P2, _params[Par_P2] - fracDelta);
         _params[Par_P3] = clamp(Par_P3, _params[Par_P3] - fracDelta);
         break;
      case 1: _params[Par_P0] = clamp(Par_P0, _params[Par_P0] - delta); break;
      case 2: _params[Par_P2] = clamp(Par_P2, _params[Par_P2] - delta); break;
      case 3: _params[Par_P3] = clamp(Par_P3, _params[Par_P3] - delta); break;
      default : break;
      }
      normalizeP();
   }   
}

void FTAGParams::setP2(double p2) 
{
   if (!_isPLinkedGC && ! _isPLinkedGA)
   {
      double delta = clamp(Par_P2, p2) - _params[Par_P2];
      double fracDelta;
      _params[Par_P2] += delta;
      size_t sel = rand() % 4;
      
      switch(sel) {
      case 0 : 
         fracDelta = delta / 3. ;
         _params[Par_P0] = clamp(Par_P0, _params[Par_P0] - fracDelta);
         _params[Par_P1] = clamp(Par_P1, _params[Par_P1] - fracDelta);
         _params[Par_P3] = clamp(Par_P3, _params[Par_P3] - fracDelta);
         break;
      case 1: _params[Par_P0] = clamp(Par_P0, _params[Par_P0] - delta); break;
      case 2: _params[Par_P1] = clamp(Par_P1, _params[Par_P1] - delta); break;
      case 3: _params[Par_P3] = clamp(Par_P2, _params[Par_P3] - delta); break;
      default : break;
      }
      normalizeP();
   }   
}

void FTAGParams::setP3(double p3) 
{
   if (!_isPLinkedGC && ! _isPLinkedGA)
   {
      double delta = clamp(Par_P3, p3) - _params[Par_P3];
      double fracDelta;
      _params[Par_P3] += delta;
      size_t sel = rand() % 4;
      
      switch(sel) {
      case 0 : 
         fracDelta = delta / 3. ;
         _params[Par_P0] = clamp(Par_P0, _params[Par_P0] - fracDelta);
         _params[Par_P1] = clamp(Par_P1, _params[Par_P1] - fracDelta);
         _params[Par_P2] = clamp(Par_P2, _params[Par_P2] - fracDelta);
         break;
      case 1: _params[Par_P0] = clamp(Par_P0, _params[Par_P0] - delta); break;
      case 2: _params[Par_P1] = clamp(Par_P1, _params[Par_P1] - delta); break;
      case 3: _params[Par_P2] = clamp(Par_P2, _params[Par_P2] - delta); break;
      default : break;
      }
      normalizeP();
   }
}

void FTAGParams::setQ0(double q0) 
{
   if (rand() % 50)
   {
      if (_isQLinkedGC == true)
      {
         _params[Par_Q0] = clamp(Par_Q0, q0) > 0.5 ? 0.5 : clamp(Par_Q0, q0);
         _params[Par_Q1] = 0.5 - _params[Par_Q0];
         _params[Par_Q2] = _params[Par_Q1];
         _params[Par_Q3] = _params[Par_Q0];
      }
      else if (_isQLinkedGA == true)
      {
         _params[Par_Q0] = clamp(Par_Q0, q0) > 0.5 ? 0.5 : clamp(Par_Q0, q0);
         _params[Par_Q1] = 0.5 - _params[Par_Q0];
         _params[Par_Q2] = _params[Par_Q0];
         _params[Par_Q3] = _params[Par_Q1];
      }
      else
      {
         double delta = clamp(Par_Q0, q0) - _params[Par_Q0];
         double fracDelta;
         _params[Par_Q0] += delta;
         size_t sel = rand() % 4;
         
         switch(sel) {
         case 0 : 
            fracDelta = delta / 3. ;
            _params[Par_Q1] = clamp(Par_Q1, _params[Par_Q1] - fracDelta);
            _params[Par_Q2] = clamp(Par_Q2, _params[Par_Q2] - fracDelta);
            _params[Par_Q3] = clamp(Par_Q3, _params[Par_Q3] - fracDelta);
            break;
         case 1: _params[Par_Q1] = clamp(Par_Q1, _params[Par_Q1] - delta); 
            break;
         case 2: _params[Par_Q2] = clamp(Par_Q2, _params[Par_Q2] - delta); 
            break;
         case 3: _params[Par_Q3] = clamp(Par_Q3, _params[Par_Q3] - delta); 
            break;
         default : break;
         }
         normalizeQ();
      }
   }   
   else if (rand() % 2)
   {
      _params[Par_Q0] = .25; _params[Par_Q1] = .25;
      _params[Par_Q2] = .25; _params[Par_Q3] = .25;
   }
   else
   {
      randomSwapQ();
   }
}

void FTAGParams::setQ1(double q1) 
{
   if (!_isQLinkedGC && ! _isQLinkedGA)
   {
      double delta = clamp(Par_Q1, q1) - _params[Par_Q1];
      double fracDelta;
      _params[Par_Q1] += delta;
      size_t sel = rand() % 4;
      
      switch(sel) {
      case 0 : 
         fracDelta = delta / 3. ;
         _params[Par_Q0] = clamp(Par_Q0, _params[Par_Q0] - fracDelta);
         _params[Par_Q2] = clamp(Par_Q2, _params[Par_Q2] - fracDelta);
         _params[Par_Q3] = clamp(Par_Q3, _params[Par_Q3] - fracDelta);
         break;
      case 1: _params[Par_Q0] = clamp(Par_Q0, _params[Par_Q0] - delta); break;
      case 2: _params[Par_Q2] = clamp(Par_Q2, _params[Par_Q2] - delta); break;
      case 3: _params[Par_Q3] = clamp(Par_Q3, _params[Par_Q3] - delta); break;
      default : break;
      }
      normalizeQ();
   }   
}

void FTAGParams::setQ2(double q2) 
{
   if (!_isQLinkedGC && ! _isQLinkedGA)
   {
      double delta = clamp(Par_Q2, q2) - _params[Par_Q2];
      double fracDelta;
      _params[Par_Q2] += delta;
      size_t sel = rand() % 4;
      
      switch(sel) {
      case 0 : 
         fracDelta = delta / 3. ;
         _params[Par_Q0] = clamp(Par_Q0, _params[Par_Q0] - fracDelta);
         _params[Par_Q1] = clamp(Par_Q1, _params[Par_Q1] - fracDelta);
         _params[Par_Q3] = clamp(Par_Q3, _params[Par_Q3] - fracDelta);
         break;
      case 1: _params[Par_Q0] = clamp(Par_Q0, _params[Par_Q0] - delta); break;
      case 2: _params[Par_Q1] = clamp(Par_Q1, _params[Par_Q1] - delta); break;
      case 3: _params[Par_Q3] = clamp(Par_Q3, _params[Par_Q3] - delta); break;
      default : break;
      }
      normalizeQ();
   }   
}

void FTAGParams::setQ3(double q3) 
{
   if (!_isQLinkedGC && ! _isQLinkedGA)
   {
      double delta = clamp(Par_Q3, q3) - _params[Par_Q3];
      double fracDelta;
      _params[Par_Q3] += delta;
      size_t sel = rand() % 4;
      
      switch(sel) {
      case 0 : 
         fracDelta = delta / 3. ;
         _params[Par_Q0] = clamp(Par_Q0, _params[Par_Q0] - fracDelta);
         _params[Par_Q1] = clamp(Par_Q1, _params[Par_Q1] - fracDelta);
         _params[Par_Q2] = clamp(Par_Q2, _params[Par_Q2] - fracDelta);
         break;
      case 1: _params[Par_Q0] = clamp(Par_Q0, _params[Par_Q0] - delta); break;
      case 2: _params[Par_Q1] = clamp(Par_Q1, _params[Par_Q1] - delta); break;
      case 3: _params[Par_Q2] = clamp(Par_Q2, _params[Par_Q2] - delta); break;
      default : break;
      }
      normalizeQ();
   }
}

void FTAGParams::normalizeP()
{
   double tot = _params[Par_P0] + _params[Par_P1] + 
      _params[Par_P2] + _params[Par_P3];
   if (tot != 0 && tot != 1.)
   {
      _params[Par_P0] /= tot;
      _params[Par_P1] /= tot;
      _params[Par_P2] /= tot;
      _params[Par_P3] /= tot;
   }
 // hack
   if (_params[Par_P0] < _clamps[0][Par_P0] || 
       _params[Par_P0] > _clamps[1][Par_P0] ||
       _params[Par_P1] < _clamps[0][Par_P1] || 
       _params[Par_P1] > _clamps[1][Par_P1] ||
       _params[Par_P2] < _clamps[0][Par_P2] || 
       _params[Par_P2] > _clamps[1][Par_P2] ||
       _params[Par_P3] < _clamps[0][Par_P3] || 
       _params[Par_P3] > _clamps[1][Par_P3])
   {
      _params[Par_P0] = _params[Par_P1] = 
         _params[Par_P2] = _params[Par_P3] = .25;
   }
}

void FTAGParams::normalizeQ()
{
   double tot = _params[Par_Q0] + _params[Par_Q1] + 
      _params[Par_Q2] + _params[Par_Q3];
   if (tot != 0. && tot != 1.)
   {
      _params[Par_Q0] /= tot;
      _params[Par_Q1] /= tot;
      _params[Par_Q2] /= tot;
      _params[Par_Q3] /= tot;
   }
   // hack
   if (_params[Par_Q0] < _clamps[0][Par_Q0] || 
       _params[Par_Q0] > _clamps[1][Par_Q0] ||
       _params[Par_Q1] < _clamps[0][Par_Q1] || 
       _params[Par_Q1] > _clamps[1][Par_Q1] ||
       _params[Par_Q2] < _clamps[0][Par_Q2] || 
       _params[Par_Q2] > _clamps[1][Par_Q2] ||
       _params[Par_Q3] < _clamps[0][Par_Q3] || 
       _params[Par_Q3] > _clamps[1][Par_Q3])
   {
      _params[Par_Q0] = _params[Par_Q1] = 
         _params[Par_Q2] = _params[Par_Q3] = .25;
   }
}

void FTAGParams::randomSwapP()
{
   if (_isPLinkedGC || _isPLinkedGA)
   {
      swap(_params[Par_P0], _params[Par_P1]);
      swap(_params[Par_P2], _params[Par_P3]);
   }
   else
   {
      size_t i = rand() % 4;
      size_t j = rand() % 4;
      while (j == i) j = rand() % 4;
      swap(_params[Par_P0 + i], _params[Par_P0 + j]);
   }
}

void FTAGParams::randomSwapQ()
{
   if (_isQLinkedGC || _isQLinkedGA)
   {
      swap(_params[Par_Q0], _params[Par_Q1]);
      swap(_params[Par_Q2], _params[Par_Q3]);
   }
   else
   {
      size_t i = rand() % 4;
      size_t j = rand() % 4;
      while (j == i) j = rand() % 4;
      swap(_params[Par_Q0 + i], _params[Par_Q0 + j]);
   }
}

void FTAGParams::setSingleF84(bool sf)
{
   _singleF84 = sf;
   if (sf == true)
   {
      setMu(0);
      setMuFixed(true);
   }
   else
   {
      setPFlatFixed(true);
      setA(0);
      setAFixed(true);
      setB(0);
      setBFixed(true);
   }
}

void FTAGParams::setDoubleF84(bool df)
{
   _doubleF84 = df;
   if (df == true)
   {
      setGamma(0);
      setGammaFixed(true);
   }
   else
   {
      setQFlatFixed(true);
      setC(0);
      setCFixed(true);
      setD(0);
      setDFixed(true);
   }
}

void FTAGParams::set(Param param, double val)
{
   switch (param)
   {
   case Par_T:
      setT(val);
      break;
   case Par_Mu:
      setMu(val);
      break;
   case Par_A:
      setA(val);
      break;
   case Par_B:
      setB(val);
      break;
   case Par_P0:
      setP0(val);
      break;
   case Par_P1:
      setP1(val);
      break;
   case Par_P2:
      setP2(val);
      break;
   case Par_P3:
      setP3(val);
      break;
   case Par_Gamma:
      setGamma(val);
      break;
   case Par_C:
      setC(val);
      break;
   case Par_D:
      setD(val);
      break;
   case Par_Q0:
      setQ0(val);
      break;
   case Par_Q1:
      setQ1(val);
      break;
   case Par_Q2:
      setQ2(val);
      break;
   case Par_Q3:
      setQ3(val);
      break;
   case Par_PE:
      setPE(val);
      break;
   case Par_RD:
      setRD(val);
      break;
   case Par_RI:
      setRI(val);
      break;
   case Par_RMD:
      setRMD(val);
      break;
   case Par_RMI:
      setRMI(val);
      break;
   case Par_PCD:
      setPCD(val);
      break;
   case Par_PCI:
      setPCI(val);
      break;
   case Par_KD:
      setKD(val);
      break;
   case Par_KI:
      setKI(val);
      break;
   case Par_Max:
      assert(false);
      break;
   }
}

void FTAGParams::setAll(double val)
{
   for (size_t i = 0; i < Par_Max; ++i)
   {
      _params[i] = val;
   }
}

void FTAGParams::setFixed(Param param, bool val)
{
   if (!_isSymmetric)
   {
      _fixed[param] = val;
   }
   else
   {
      switch (param)
      {
      case Par_T:
      case Par_Mu:
      case Par_A:
      case Par_B:
      case Par_Gamma:
      case Par_C:
      case Par_D:
      case Par_PE:
         _fixed[param] = val;
         break;
      case Par_RD:
         _fixed[param] = val;
         _fixed[Par_RI] = val;
         break;
      case Par_RMD:
         _fixed[param] = val;
         _fixed[Par_RMI] = val;
         break;
      case Par_PCD:
         _fixed[param] = val;
         _fixed[Par_PCI] = val;
         break;
      case Par_KD:
         _fixed[param] = val;
         _fixed[Par_KI] = val;
         break;
      default:
         break;
      }
   }
   if (_isUniGap)
   {
      _fixed[Par_PCI] = _fixed[Par_KI];
      _fixed[Par_PCD] = _fixed[Par_KD];
   }
}

void FTAGParams::setClamp(Param param, double lower, double upper)
{
   assert(lower < upper);
   _clamps[0][param] = lower;
   _clamps[1][param] = upper;
}

void FTAGParams::randomizeNonFixed()
{
   for (size_t i = 0; i < Par_Max; ++i)
   {
      // redundant since set() should check already
      if (isFixed((Param)i) == false)
      {
         set((Param)i, drand48() * (min(1., _clamps[1][i]) - _clamps[0][i]));
      }
   }
   // seeding with a high PE can be irrecoverable sincel longer 
   // sequences will run out of precision and get pr=0 every time
   // so we prevent with this little hack
   if (!isPEFixed() && getPE() > 0.8)
   {
      setPE(getPE() * 0.8);
   }
   // same goes for indel rates
   // to do: clean this up:
   if (!isRDFixed() && getRD() > 0.1)
   {
      setRD(getRD() * 0.1);
   }
   if (!isRIFixed() && getRI() > 0.1)
   {
      setRI(getRI() * 0.1);
   }
   if (!isKDFixed() && getKD() > 0.5)
   {
      setKD(getKD() * 0.25);
   }
   if (!isKIFixed() && getKI() > 0.5)
   {
      setKI(getKI() * 0.5);
   }
   if (!isRMDFixed() && getRMD() > 0.1)
   {
      setRMD(getRMD() * 0.1);
   }
   if (!isRMIFixed() && getRMI() > 0.1)
   {
      setRMI(getRMI() * 0.1);
   }
   if (!isPCDFixed() && getPCD() > 0.95)
   {
      setPCD(getPCD() * 0.95);
   }
   if (!isPCIFixed() && getPCI() > 0.95)
   {
      setPCI(getPCI() * 0.95);
   }
   if (!isMuFixed() && getMu() > 0.5)
   {
      setMu(getMu() * 0.5);
   }
   if (!isAFixed() && getA() > 0.5)
   {
      setA(getA() * 0.5);
   }
   if (!isBFixed() && getB() > 0.5)
   {
      setB(getB() * 0.5);
   }
   if (!isGammaFixed() && getGamma() > 0.5)
   {
      setGamma(getGamma() * 0.5);
   }
   if (!isCFixed() && getC() > 0.5)
   {
      setC(getC() * 0.5);
   }
   if (!isDFixed() && getD() > 0.5)
   {
      setD(getD() * 0.5);
   }
   // don't randomize priors. 
   if (!isP0Fixed()) _params[Par_P0] = 0.25;
   if (!isP1Fixed()) _params[Par_P1] = 0.25;
   if (!isP2Fixed()) _params[Par_P2] = 0.25;
   if (!isP3Fixed()) _params[Par_P3] = 0.25;
   if (!isQ0Fixed()) _params[Par_Q0] = 0.25;
   if (!isQ1Fixed()) _params[Par_Q1] = 0.25;
   if (!isQ2Fixed()) _params[Par_Q2] = 0.25;
   if (!isQ3Fixed()) _params[Par_Q3] = 0.25;
}

string FTAGParams::asRow() const
{
   stringstream row;
   
   row << getT() << " ";
   row << getMu() << " ";
   row << getA() << " ";
   row << getB() << " ";
   row << getP0() << " ";
   row << getP1() << " ";
   row << getP2() << " ";
   row << getP3() << " ";
   row << getGamma() << " ";
   row << getC() << " ";
   row << getD() << " ";
   row << getQ0() << " ";
   row << getQ1() << " ";
   row << getQ2() << " ";
   row << getQ3() << " ";
   row << getRD() << " ";
   row << getRI() << " ";
   row << getRMD() << " ";
   row << getRMI() << " ";
   row << getPE() << " ";
   row << getPCD() << " ";
   row << getPCI() << " ";
   row << getKD() << " ";
   row << getKI();

   return row.str();
}

ostream& operator<<(ostream& os, const FTAGParams& params)
{
   os << "<";
   if (params.isTFixed()) os << "*";
   os << "t= " << params.getT() << " ";
   if (params.isMuFixed()) os << "*";
   os << "mu= " << params.getMu() << " ";
   if (params.isAFixed()) os << "*";
   os << "a= " << params.getA() << " ";
   if (params.isBFixed()) os << "*";
   os << "b= " << params.getB() << " ";
   if (params.isP0Fixed()) os << "*";
   os << "p0= " << params.getP0() << " ";
   if (params.isP1Fixed()) os << "*";
   os << "p1= " << params.getP1() << " ";
   if (params.isP2Fixed()) os << "*";
   os << "p2= " << params.getP2() << " ";
   if (params.isP3Fixed()) os << "*";
   os << "p3= " << params.getP3() << " ";
   if (params.isGammaFixed()) os << "*";   
   os << "ga= " << params.getGamma() << " ";
   if (params.isCFixed()) os << "*";
   os << "c= " << params.getC() << " ";
   if (params.isDFixed()) os << "*";
   os << "d= " << params.getD() << " ";
   if (params.isQ0Fixed()) os << "*";
   os << "q0= " << params.getQ0() << " ";
   if (params.isQ1Fixed()) os << "*";
   os << "q1= " << params.getQ1() << " ";
   if (params.isQ2Fixed()) os << "*";
   os << "q2= " << params.getQ2() << " ";
   if (params.isQ3Fixed()) os << "*";
   os << "q3= " << params.getQ3() << " ";
   if (params.isRDFixed()) os << "*";
   os << "RD= " << params.getRD() << " ";
   if (params.isRIFixed()) os << "*";
   os << "RI= " << params.getRI() << " ";
   if (params.isRMDFixed()) os << "*";
   os << "RMD= " << params.getRMD() << " ";
   if (params.isRMIFixed()) os << "*";
   os << "RMI= " << params.getRMI() << " ";
   if (params.isPEFixed()) os << "*";
   os << "PE= " << params.getPE() << " ";
   if (params.isPCDFixed()) os << "*";
   os << "PCD= " << params.getPCD() << " ";
   if (params.isPCIFixed()) os << "*";
   os << "PCI= " << params.getPCI() << " ";
   if (params.isKDFixed()) os << "*";
   os << "KD= " << params.getKD() << " ";
   if (params.isKIFixed()) os << "*";
   os << "KI= " << params.getKI() << " >";
   os << endl;
   return os;
}
 
istream& operator>>(istream& is, FTAGParams& params)
{
   string sbuf;
   double dbuf;
   
   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "<*t=" || sbuf == "<t="))
   {
      throw string("Param read error t");
   }
   params.setTFixed(sbuf[1] == '*');
   params.setT(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*mu=" || sbuf == "mu="))
   {
      throw string("Param read error mu");
   }
   params.setMuFixed(sbuf[0] == '*');
   params.setMu(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*a=" || sbuf == "a="))
   {
      throw string("Param read error a");
   }
   params.setAFixed(sbuf[0] == '*');
   params.setA(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*b=" || sbuf == "b="))
   {
      throw string("Param read error b");
   }
   params.setBFixed(sbuf[0] == '*');
   params.setB(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*p0=" || sbuf == "p0="))
   {
      throw string("Param read error p0");
   }
   params.setP0Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_P0, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*p1=" || sbuf == "p1="))
   {
      throw string("Param read error p1");
   }
   params.setP1Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_P1, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*p2=" || sbuf == "p2="))
   {
      throw string("Param read error p2");
   }
   params.setP2Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_P2, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*p3=" || sbuf == "p3="))
   {
      throw string("Param read error p3");
   }
   params.setP3Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_P3, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*ga=" || sbuf == "ga="))
   {
      throw string("Param read error ga: ") + sbuf;
   }
   params.setGammaFixed(sbuf[0] == '*');
   params.setGamma(dbuf);
   
   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*c=" || sbuf == "c="))
   {
      throw string("Param read error c");
   }
   params.setCFixed(sbuf[0] == '*');
   params.setC(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*d=" || sbuf == "d="))
   {
      throw string("Param read error d");
   }
   params.setDFixed(sbuf[0] == '*');
   params.setD(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*q0=" || sbuf == "q0="))
   {
      throw string("Param read error q0");
   }
   params.setQ0Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_Q0, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*q1=" || sbuf == "q1="))
   {
      throw string("Param read error q1");
   }
   params.setQ1Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_Q1, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*q2=" || sbuf == "q2="))
   {
      throw string("Param read error q2");
   }
   params.setQ2Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_Q2, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*q3=" || sbuf == "q3="))
   {
      throw string("Param read error q3");
   }
   params.setQ3Fixed(sbuf[0] == '*');
   params.setDirect(FTAGParams::Par_Q3, dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*RD=" || sbuf == "RD="))
   {
      throw string("Param read error RD");
   }
   params.setRDFixed(sbuf[0] == '*');
   params.setRD(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*RI=" || sbuf == "RI="))
   {
      throw string("Param read error RI");
   }
   params.setRIFixed(sbuf[0] == '*');
   params.setRI(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*RMD=" || sbuf == "RMD="))
   {
      throw string("Param read error RMD");
   }
   params.setRMDFixed(sbuf[0] == '*');
   params.setRMD(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*RMI=" || sbuf == "RMI="))
   {
      throw string("Param read error RMI");
   }
   params.setRMIFixed(sbuf[0] == '*');
   params.setRMI(dbuf);
   
   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*PE=" || sbuf == "PE="))
   {
      throw string("Param read error PE");
   }
   params.setPEFixed(sbuf[0] == '*');
   params.setPE(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*PCD=" || sbuf == "PCD="))
   {
      throw string("Param read error PCD");
   }
   params.setPCDFixed(sbuf[0] == '*');
   params.setPCD(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*PCI=" || sbuf == "PCI="))
   {
      throw string("Param read error PCI");
   }
   params.setPCIFixed(sbuf[0] == '*');
   params.setPCI(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*KD=" || sbuf == "KD="))
   {
      throw string("Param read error KD");
   }
   params.setKDFixed(sbuf[0] == '*');
   params.setKD(dbuf);

   is >> sbuf >> dbuf;
   if (!is || !(sbuf == "*KI=" || sbuf == "KI="))
   {
      throw string("Param read error KI");
   }
   params.setKIFixed(sbuf[0] == '*');
   params.setKI(dbuf);
   
   is >> sbuf;

   return is;
}
