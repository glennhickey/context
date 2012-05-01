//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cassert>
#include "dnasubmodel.h"

using namespace std;

DNASubModel::DNASubModel() : CtMarkov(DNA_Max)
{

}

DNASubModel::~DNASubModel()
{

}

DNASubModel::DNASubModel(const DNASubModel& model) : CtMarkov(model._size)
{
   *this = model;
}

DNASubModel& DNASubModel::operator=(const DNASubModel& model)
{
   this->CtMarkov::operator=(model);
   return *this;
   // nothing class specific right now
}


DNA DNASubModel::charToDNA(char c)
{
   if (c == 'A' || c == 'a') return DNA_A;
   else if (c == 'C' || c == 'c') return DNA_C;
   else if (c == 'G' || c == 'g') return DNA_G;
   else if (c == 'T' || c == 't') return DNA_T;
   else 
   {
      assert(false);
      return DNA_Max;
   }
}

char DNASubModel::DNAToChar(DNA dna)
{
   if (dna == DNA_A) return 'A';
   else if (dna == DNA_C) return 'C';
   else if (dna == DNA_G) return 'G';
   else if (dna == DNA_T) return 'T';
   else
   {
      assert(false);
      return 'N';
   }
}

void DNASubModel::setJukesCantor(double mu, double t)
{
   _t = t;
   /*
   double tran = mu / (double)_size;
   double stay = -tran * (double)(_size - 1);
   double val;
   for (size_t i = 0; i < _size; ++i)
   {
      for (size_t j = 0; j < _size; ++j)
      {
         val = i != j ? tran : stay;
         _Q0.set(i, j, val);
      }
      _pi[i] = 1. / (double)_size;
   }
   computeProbMatrix();
   */

   double ds = (double)_size;
   double tran = (1. / ds) - (1. / ds) * exp(-t * mu);
   double stay = (1. / ds) + ((ds - 1.) / ds) * exp(-t * mu);
   double val;
   for (size_t i = 0; i < _size; ++i)
   {
      for (size_t j = 0; j < _size; ++j)
      {
         val = i != j ? tran : stay;
         _Qt.set(i, j, val);
      }
      _pi[i] = 1. / ds;
   }
}

void DNASubModel::setF84(double p0, double p1, double p2, double p3, double a, 
                         double b, double t)
{
   _pi[0] = p0;
   _pi[1] = p1;
   _pi[2] = p2;
   _pi[3] = p3;
   _t = t;
/*
   double pr = p0 + p2;
   double py = p1 + p3;

   double b = 1. / (2. * pr * py * (1 + r));
   double a = (pr * py * r - p0 * p2 - p1 * p3) / 
      (2. * (1. + r) * (py * p0 * p2  + pr * p1 * p3));
*/ 
   double wk[4][4];
   for (size_t i = 0; i < 4; ++i)
   {
      for (size_t j = 0; j < 4; ++j)
      {
         wk[i][j] = (int)abs((int)i - (int)j) % 2 == 0 ? 1. : 0.;
      }
   }

   for (size_t i = 0; i < 4; ++i)
   {
      for (size_t j = 0; j < 4; ++j)
      {
         
         double dij = i == j ? 1. : 0.;

         double val = exp(-(a + b) * t) * dij +
            exp(-b * t) * (1 - exp(-a * t)) * (_pi[j] * wk[i][j] / 
                                         (wk[j][0] * _pi[0] +
                                          wk[j][1] * _pi[1] +
                                          wk[j][2] * _pi[2] +
                                          wk[j][3] * _pi[3])) +
            (1. - exp(-b * t)) * _pi[j];
         _Qt.set(i, j, val);
      }
   }
}
