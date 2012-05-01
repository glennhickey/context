//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _CONTEXTMODELTABLE_H
#define _CONTEXTMODELTABLE_H

#include <ostream>
#include <cmath>
#include "dnasubmodel.h"

class ContextModelTable
{
public:
      
   enum IDX1D { IDX1_D = 0, IDX1_I, IDX1_Max };
   enum IDX2D { IDX2_M = 0, IDX2_Max };
   enum IDX3D { IDX3_MD = 0, IDX3_MI, IDX3_Max};

   
   ContextModelTable();
   ContextModelTable(const ContextModelTable& cmt);
   ContextModelTable& operator=(const ContextModelTable& cmt);
   ~ContextModelTable();

   ContextModelTable& operator/=(const ContextModelTable& cmt);
   ContextModelTable& operator+=(const ContextModelTable& cmt);
   void scale(const std::vector<double>& svec);
   double mse(const ContextModelTable& cmt) const;
   
   void setAll(double val);

   double getD(DNA i) const
   {
      return _tab1[IDX1_D][i];
   }
   double getI(DNA p) const
   {
      return _tab1[IDX1_I][p];
   }
   double getM(DNA i, DNA p) const
   {
      return _tab2[IDX2_M][i][p];
   }
   double getMD(DNA i, DNA k, DNA p) const 
   {
      return _tab3[IDX3_MD][i][k][p];
   }
   double getMI(DNA i, DNA p, DNA r) const 
   {
      return _tab3[IDX3_MI][i][p][r];
   }

   void setD(size_t i, double val) 
   {
      assert(!std::isnan(val) && !std::isinf(val));
      _tab1[IDX1_D][i] = val;
   }
   void setI(size_t p, double val) 
   {
      assert(!std::isnan(val) && !std::isinf(val));
      _tab1[IDX1_I][p] = val;
   }
   void setM(size_t i, size_t p, double val) 
   {
      assert(!std::isnan(val) && !std::isinf(val));
      _tab2[IDX2_M][i][p] = val;
   }
   void setMD(DNA i, DNA k, DNA p, double val)  
   {
      assert(!std::isnan(val) && !std::isinf(val));
      _tab3[IDX3_MD][i][k][p] = val;
   }
   void setMI(DNA i, DNA p, DNA r, double val)  
   {
      assert(!std::isnan(val) && !std::isinf(val));
      _tab3[IDX3_MI][i][p][r] = val;
   }

   void addD(size_t i, double val) 
   {
      _tab1[IDX1_D][i] += val;
   }
   void addI(size_t p, double val) 
   {
      _tab1[IDX1_I][p] += val;
   }
   void addM(size_t i, size_t p, double val) 
   {
      _tab2[IDX2_M][i][p] += val;
   }
   void addMD(DNA i, DNA k, DNA p, double val)  
   {
      _tab3[IDX3_MD][i][k][p] += val;
   }
   void addMI(DNA i, DNA p, DNA r, double val)  
   {
      _tab3[IDX3_MI][i][p][r] += val;
   }

   // debug
   bool verify(double threshold = 0.0001) const;
   double sumTripleSq() const;
   double sumTripleAbs() const;
   double mse1(const ContextModelTable&) const;
   double mse2(const ContextModelTable&) const;
   double mse3(const ContextModelTable&) const;

protected:

   double _tab1[IDX1_Max][DNA_Max];
   double _tab2[IDX2_Max][DNA_Max][DNA_Max];
   double _tab3[IDX3_Max][DNA_Max][DNA_Max][DNA_Max];

   static const size_t _size;
};

std::ostream& operator<<(std::ostream& os, const ContextModelTable& tab);

#endif
