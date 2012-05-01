//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _CTMARKOV_H
#define _CTMARKOV_H

#include <cstdlib>
#include <vector>

#include "matrix.h"


class CtMarkov
{
public:

   CtMarkov(size_t size);
   virtual ~CtMarkov();

   CtMarkov(const CtMarkov& model);
   CtMarkov& operator=(const CtMarkov& model);

   void setAllZero();

   void setRateMatrix(const Matrix<double>& rm);
   void setProbMatrix(const Matrix<double>& pm);
   void setStationaryDist(const std::vector<double>& pi);
   void setT(double t);

   void computeProbMatrix();
   void computeStationaryDist();

   double rate(size_t from, size_t to) const { return _Q0.get(from, to); }
   double prTrans(size_t from, size_t to) const { return _Qt.get(from, to); }
   double pr(size_t state) const { return _pi[state]; }
   double prFull(size_t from, size_t to) const { 
      return prTrans(from, to) * pr(from);
   }
   size_t size() const { return _size; }
   double getT() const { return _t; }

   const Matrix<double>& getRateMatrix() const { return _Q0; }
   const Matrix<double>& getProbMatrix() const { return _Qt; }
   const std::vector<double> getStationaryDist() const { return _pi; }

   bool verifyRateMatrix() const;
   bool verifyProbMatrix() const;
   bool verifyStationaryDist() const; 

protected:

   void allocMatrix();
   void deleteMatrix();

protected:

   size_t _size; // Alphabet Size
   std::vector<double> _pi; // equilibrium probs
   Matrix<double> _Q0; // Rate Matrix
   Matrix<double> _Qt; // Sub Matrix at t
   double _t; // time

   CtMarkov();
};

#endif
