//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _MATRIX_H
#define _MATRIX_H

#include <ostream>

template<class T>
class Matrix
{
public:
   Matrix();
   virtual ~Matrix();
   Matrix(const Matrix& mat);
   Matrix(size_t width, size_t height);
   Matrix(size_t width);

   void resize(size_t width, size_t height);
   size_t getWidth() const;
   size_t getHeight() const;
   
   Matrix& operator=(const Matrix<T>& mat);

   Matrix& operator+=(const Matrix<T>& mat);
   Matrix& operator/=(const Matrix<T>& mat);
   double mse(const Matrix<T>& mat) const;

   T sumAllElems() const;
   T get(size_t row, size_t col) const;
   void set(size_t row, size_t col, T val);
   void add(size_t row, size_t col, T val);
   void multByScalar(T val);
   void setAll(T val);

   const T* getPointer() const;
   
protected:
   
   size_t _h;
   size_t _w;
   T* _mat;
};

template<class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat);

#include "matrix_impl.h"

#endif
