//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cmath>

#include "contextmodeltable.h"
#include "dnasubmodel.h"
#include "ftagstates.h"

using namespace std;
using namespace FTAGStates;

const size_t ContextModelTable::_size =  
   (size_t)IDX1_Max * (size_t)DNA_Max +
   (size_t)IDX2_Max * (size_t)DNA_Max * 
   (size_t)DNA_Max +
   (size_t)IDX3_Max * (size_t)DNA_Max * 
   (size_t)DNA_Max;

// there's room for speed improvements in this implementation but it
// shouldn't be a bottleneck right now..

ContextModelTable::ContextModelTable()
{
   setAll(0.);
}

ContextModelTable::~ContextModelTable()
{

}

ContextModelTable::ContextModelTable(const ContextModelTable& cmt)
{
   *this = cmt;
}

ContextModelTable& ContextModelTable::operator=(const ContextModelTable& cmt)
{
   size_t s;

   for (size_t i = 0; i < DNA_Max; ++i)
   {
      // dear compiler: please unroll me...
      for (s = 0; s < IDX1_Max; ++s)
      {
         _tab1[s][i] = cmt._tab1[s][i];
      }
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (s = 0; s < IDX2_Max; ++s)
         {
            _tab2[s][i][j] = cmt._tab2[s][i][j];
         }
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            for (s = 0; s < IDX3_Max; ++s)
            {
               _tab3[s][i][j][k] = cmt._tab3[s][i][j][k];
            }
         }
      }
   }
   return *this;
}

// NOTE: 0/0 = 0;
ContextModelTable& ContextModelTable::operator/=(const ContextModelTable& cmt)
{
   size_t s;

   for (size_t i = 0; i < DNA_Max; ++i)
   {
      // dear compiler: please unroll me...
      for (s = 0; s < IDX1_Max; ++s)
      {
         assert(cmt._tab1[s][i] != 0. || _tab1[s][i] == 0.);
         if (cmt._tab1[s][i] != 0.)
         {
            _tab1[s][i] /= cmt._tab1[s][i];
         }
      }
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (s = 0; s < IDX2_Max; ++s)
         {
            assert(cmt._tab2[s][i][j] != 0. || _tab2[s][i][j] == 0.);
            if (cmt._tab2[s][i][j] != 0.)
            {
               _tab2[s][i][j] /= cmt._tab2[s][i][j];
            }
         }
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            for (s = 0; s < IDX3_Max; ++s)
            {
               assert(cmt._tab3[s][i][j][k] != 0. || _tab3[s][i][j][k] == 0.);
               if (cmt._tab3[s][i][j][k] != 0.)
               {
                  _tab3[s][i][j][k] /= cmt._tab3[s][i][j][k];
               }
            }
         }
      }
   }
   return *this;
}

ContextModelTable& ContextModelTable::operator+=(const ContextModelTable& cmt)
{
   size_t s;

   for (size_t i = 0; i < DNA_Max; ++i)
   {
      // dear compiler: please unroll me...
      for (s = 0; s < IDX1_Max; ++s)
      {
         _tab1[s][i] += cmt._tab1[s][i];
      }
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (s = 0; s < IDX2_Max; ++s)
         {
            _tab2[s][i][j] += cmt._tab2[s][i][j];
         }
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            for (s = 0; s < IDX3_Max; ++s)
            {
               _tab3[s][i][j][k] += cmt._tab3[s][i][j][k];
            }
         }
      }
   }
   return *this;
}

void ContextModelTable::scale(const vector<double>& svec)
{
   double md = svec[S_MDx] + svec[S_MDyM];
   double mi = svec[S_MIx] + svec[S_MIzM];

   for (size_t i = 0; i < DNA_Max; ++i)
   {
      _tab1[IDX1_D][i] *= svec[S_Dx];
      _tab1[IDX1_I][i] *= svec[S_Ix];

      for (size_t j = 0; j < DNA_Max; ++j)
      {
         _tab2[IDX2_M][i][j] *= svec[S_Mx];
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            _tab3[IDX3_MD][i][j][k] *= md;
            _tab3[IDX3_MI][i][j][k] *= mi;
         }
      }
   }
}

double ContextModelTable::mse(const ContextModelTable& cmt) const
{
   size_t s;
   double diff;
   double total1 = 0;
   double total2 = 0;
   double total3 = 0;
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      // dear compiler: please unroll me...
      for (s = 0; s < IDX1_Max; ++s)
      {
         diff = _tab1[s][i] - cmt._tab1[s][i];
         total1 += diff * diff;
      }
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (s = 0; s < IDX2_Max; ++s)
         {
            diff = _tab2[s][i][j] - cmt._tab2[s][i][j];
            total2 += diff * diff;
         }
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            for (s = 0; s < IDX3_Max; ++s)
            {
               diff = _tab3[s][i][j][k] - cmt._tab3[s][i][j][k];
               total3 += diff * diff;
            }
         }
      }
   }
   return (total1 / 4. + total2 / 8. + total3 / 16.);
}
   
void ContextModelTable::setAll(double val)
{
   size_t s;

   for (size_t i = 0; i < DNA_Max; ++i)
   {
      // dear compiler: please unroll me...
      for (s = 0; s < IDX1_Max; ++s)
      {
         _tab1[s][i] = val;
      }
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (s = 0; s < IDX2_Max; ++s)
         {
            _tab2[s][i][j] = val;
         }
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            for (s = 0; s < IDX3_Max; ++s)
            {
               _tab3[s][i][j][k] = val;
            }
         }
      }
   }
}

bool ContextModelTable::verify(double threshold) const
{
   size_t s;
   double total;

   for (s = 0; s < IDX1_Max; ++s)
   {
      total = 0;
      for (size_t i = 0; i < DNA_Max; ++i)
      {
         total += _tab1[s][i];
      }
      if (fabs(1. - total) > threshold && total != 0.)
      {
         assert(false);
         return false;
      }
   }
   
   for (s = 0; s < IDX2_Max; ++s)
   {
      total = 0;
      for (size_t i = 0; i < DNA_Max; ++i)
      {
         for (size_t j = 0; j < DNA_Max; ++j)
         {
            total += _tab2[s][i][j];
         }
      }
      if (fabs(1. - total) > threshold && total != 0.)
      {
         assert(false);
         return false;
      }
   }

   for (s = 0; s < IDX3_Max; ++s)
   {
      total = 0;
      for (size_t i = 0; i < DNA_Max; ++i)
      {
         for (size_t j = 0; j < DNA_Max; ++j)
         {
            for (size_t k = 0; k < DNA_Max; ++k)
            {
               total += _tab3[s][i][j][k];
            }
         }
      }
      if (fabs(1. - total) > threshold && total != 0.)
      {
         assert(false);
         return false;
      }
   }
   return true;
}

#ifndef NDEBUG
double ContextModelTable::sumTripleSq() const
{
   double total = 0.;
   for (size_t s = 0; s < IDX3_Max; ++s)
   {
      for (size_t i = 0; i < DNA_Max; ++i)
      {
         for (size_t j = 0; j < DNA_Max; ++j)
         {
            for (size_t k = 0; k < DNA_Max; ++k)
            {
               total += _tab3[s][i][j][k] * _tab3[s][i][j][k];
            }
         }
      }
   }
   return total;
}

double ContextModelTable::sumTripleAbs() const
{
   double total = 0.;
   for (size_t s = 0; s < IDX3_Max; ++s)
   {
      for (size_t i = 0; i < DNA_Max; ++i)
      {
         for (size_t j = 0; j < DNA_Max; ++j)
         {
            for (size_t k = 0; k < DNA_Max; ++k)
            {
               total += _tab3[s][i][j][k];
            }
         }
      }
   }
   return total;
}

double ContextModelTable::mse1(const ContextModelTable& cmt) const
{
   size_t s;
   double diff;
   double total = 0;
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      // dear compiler: please unroll me...
      for (s = 0; s < IDX1_Max; ++s)
      {
         diff = _tab1[s][i] - cmt._tab1[s][i];
         total += diff * diff;
      }
   }
   return total;
}

double ContextModelTable::mse2(const ContextModelTable& cmt) const
{
   size_t s;
   double diff;
   double total = 0;
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (s = 0; s < IDX2_Max; ++s)
         {
            diff = _tab2[s][i][j] - cmt._tab2[s][i][j];
            total += diff * diff;
         }
      }
   }
   return total;
}

double ContextModelTable::mse3(const ContextModelTable& cmt) const
{
   size_t s;
   double diff;
   double total = 0;
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            for (s = 0; s < IDX3_Max; ++s)
            {
               diff = _tab3[s][i][j][k] - cmt._tab3[s][i][j][k];
               total += diff * diff;
            }
         }
      }
   }
   return total;
}

#endif
 
ostream& operator<<(ostream& os, const ContextModelTable& tab)
{
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      if (i < 2) 
      {
         os << "Pr[Del i=" << (DNA)i << "] = " << tab.getD((DNA)i) << endl;
         os << "Pr[Ins p=" << (DNA)i << "] = " << tab.getI((DNA)i) << endl;
      }
   }

   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         if (i < 2 && j < 2)
         {
            os << "Pr[M  i=" << (DNA)i << ", p=" << (DNA)j <<" ] = " 
               << tab.getM((DNA)i, (DNA)j) << endl;
/*
            os << "Pr[DD i=" << (DNA)i << ", k=" << (DNA)j <<" ] = " 
               << tab.getDD((DNA)i, (DNA)j) << endl;
            os << "Pr[ID k=" << (DNA)i << ", p=" << (DNA)j <<" ] = " 
               << tab.getID((DNA)i, (DNA)j) << endl;
            os << "Pr[DI i=" << (DNA)i << ", r=" << (DNA)j <<" ] = " 
               << tab.getDI((DNA)i, (DNA)j) << endl;
            os << "Pr[II p=" << (DNA)i << ", r=" << (DNA)j <<" ] = " 
               << tab.getII((DNA)i, (DNA)j) << endl;
*/
         }
      }
   }
   
   for (size_t i = 0; i < DNA_Max; ++i)
   {
      for (size_t j = 0; j < DNA_Max; ++j)
      {
         for (size_t k = 0; k < DNA_Max; ++k)
         {
            if (i < 2 && j < 2 && k < 2)
            {
            os << "Pr[MD i=" << (DNA)i << ", k=" << (DNA)j <<", p=" << (DNA)k
               << "] = " << tab.getMD((DNA)i, (DNA)j, (DNA)k) << endl;
            os << "Pr[MI i=" << (DNA)i << ", p=" << (DNA)j <<", r=" << (DNA)k
               << "] = " << tab.getMI((DNA)i, (DNA)j, (DNA)k) << endl;
            }
         }
      }
   }
   return os;
}
