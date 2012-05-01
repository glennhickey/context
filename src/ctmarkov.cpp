//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

//#include <gsl/gsl_linalg.h>

#include "ctmarkov.h"

using namespace std;

CtMarkov::CtMarkov(size_t size) : _size(size), _pi(size), _Q0(size), _Qt(size)
{
   setAllZero();
}

CtMarkov::~CtMarkov()
{

}

CtMarkov::CtMarkov(const CtMarkov& model) : _size(model._size), 
                                            _pi(model._pi), 
                                            _Q0(model._Q0),
                                            _Qt(model._Qt),
                                            _t(model._t)
{

}

CtMarkov& CtMarkov::operator=(const CtMarkov& model)
{
   _size = model._size; 
   _pi = model._pi; 
   _Q0 = model._Q0;
   _Qt = model._Qt;
   _t = model._t;

   return *this;
}

void CtMarkov::setAllZero()
{
   // reset the probability matrix
   for (size_t from = 0; from < _size; ++from)
   {
      _pi[from] = 0.; // reset stationary (though it's not used)
   }
   _Q0.setAll(0);
   _Qt.setAll(0);
}

void CtMarkov::setRateMatrix(const Matrix<double>& rm)
{
   _Q0 = rm;
}

void CtMarkov::setProbMatrix(const Matrix<double>& pm)
{
   _Qt = pm;
}

void CtMarkov::setStationaryDist(const vector<double>& pi)
{
   _pi = pi;
}

void CtMarkov::setT(double t)
{
   _t = t;
}

// Compute the substution probability matrix given the instantaneous
// rate matrix (Q0) and t, using exponentiation.  
/*
void CtMarkov::computeProbMatrix()
{
   // multiply by t
   Matrix<double> temp = _Q0;
   temp.multByScalar(_t);

   // gsl wants double* for interface
   size_t numBytes = _size * _size * sizeof(double);
   double* resBuf = (double*)malloc(numBytes);
   double* buf = (double*)malloc(numBytes);
   memcpy(buf, temp.getPointer(), numBytes); 

   gsl_matrix_view view = gsl_matrix_view_array(buf, _size, _size);
   gsl_matrix_view resView = gsl_matrix_view_array(resBuf, _size, _size);
   
   gsl_linalg_exponential_ss(&view.matrix, &resView.matrix, .01);

   for (size_t i = 0; i < _size; ++i) 
   {
      for (size_t j = 0; j < _size; ++j)
      {
         _Qt.set(i, j, gsl_matrix_get(&resView.matrix, i, j));
      }
   }
   
   free(buf);
   free(resBuf);
}*/

void CtMarkov::computeStationaryDist()
{
   // find pi such that pi x Q0 = pi
   // pi is therefore normalized eigenvector for 1, or something
}

bool CtMarkov::verifyRateMatrix() const
{
   for (size_t i = 0; i < _size; ++i)
   {
      double sum = 0.;
      for (size_t j = 0; j < _size; ++j)
      {
         sum += _Q0.get(i, j);
      }
      if (sum != 0.0) // precision ?
         return false;
   }
   return true;
}

bool CtMarkov::verifyProbMatrix() const
{
   for (size_t i = 0; i < _size; ++i)
   {
      double sum = 0.;
      for (size_t j = 0; j < _size; ++j)
      {
         sum += _Qt.get(i, j);
      }
      if (sum != 1.0) // precision ?
         return false;
   }
   return true;
}
