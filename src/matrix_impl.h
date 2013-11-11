//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cmath>
#include <cassert>
#include <cstring>

template<class T> inline
Matrix<T>::Matrix() : _h(0), _w(0), _mat(NULL)
{

}

template<class T> inline
Matrix<T>::~Matrix()
{
   delete [] _mat;
}

template<class T> inline
Matrix<T>::Matrix(const Matrix<T>& mat) : _h(mat._h), _w(mat._w)
{
   size_t size = _h * _w;
   _mat = new T[size];
   for (size_t i = 0; i < size; ++i)
   {
      _mat[i] = mat._mat[i];
   }
}

template<class T> inline
Matrix<T>::Matrix(size_t width) : _h(width), _w(width)
{
   _mat = new T[_h * _w];
}

template<class T> inline
Matrix<T>::Matrix(size_t width, size_t height) : _h(height), _w(width)
{
   _mat = new T[_h * _w];
}

template<class T> inline
void Matrix<T>::resize(size_t width, size_t height)
{
   if (width != _w || height != _h)
   {
      delete [] _mat;
      _w = width;
      _h = height;
      _mat = new T[_h * _w];
      std::memset(_mat, 0, sizeof(T) * _h * _w);
   }
}

template<class T> inline
size_t Matrix<T>::getWidth() const
{
   return _w;
}

template<class T> inline
size_t Matrix<T>::getHeight() const
{
   return _h;
}

template<class T> inline
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat)
{
   resize(mat._w, mat._h);
   size_t size = _h * _w;
   for (size_t i = 0; i < size; ++i)
   {
      _mat[i] = mat._mat[i];
   }
   return *this;
}

template<class T> inline
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& mat)
{
   assert(_w == mat._w && _h == mat._h);
   size_t size = _h * _w;
   for (size_t i = 0; i < size; ++i)
   {
      _mat[i] += mat._mat[i];
   }
   return *this;
}

// term by term divide.  note: 0 /= 0 = 0
template<class T> inline
Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& mat)
{
   assert(_w == mat._w && _h == mat._h);
   size_t size = _h * _w;
   for (size_t i = 0; i < size; ++i)
   {
      assert(!_mat[i] || mat._mat[i]);
      if (_mat[i])
      {
         _mat[i] /= mat._mat[i];
         assert(!std::isnan(_mat[i]));
      }
   }
   return *this;
}

template<class T> inline
double Matrix<T>::mse(const Matrix<T>& mat) const
{
   assert(_w == mat._w && _h == mat._h);
   size_t size = _h * _w;
   double diff;
   double total = 0;
   for (size_t i = 0; i < size; ++i)
   {
      diff = _mat[i] - mat._mat[i];
      total += diff * diff;
   }
   return total;
}

template<class T> inline
T Matrix<T>::sumAllElems() const
{
   size_t size = _h * _w;
   T total = (T)0;
   for (size_t i = 0; i < size; ++i)
   {
      total += _mat[i];
   }
   return total;
}


template<class T> inline
T Matrix<T>::get(size_t row, size_t col) const
{
   return _mat[row * _w + col];
}

template<class T> inline
void Matrix<T>::set(size_t row, size_t col, T val)
{
   assert(!std::isnan(val));
   _mat[row * _w + col] = val;
}

template<class T> inline
void Matrix<T>::add(size_t row, size_t col, T val)
{
   _mat[row * _w + col] += val;
}

template<class T> inline
void Matrix<T>::multByScalar(T val)
{
   assert(!std::isnan(val));
   size_t size = _h * _w;
   for (size_t i = 0; i < size; ++i)
   {
      _mat[i] *= val;
   }
}

template<class T> inline
void Matrix<T>::setAll(T val)
{
   assert(!std::isnan(val));
   size_t size = _h * _w;
   for (size_t i = 0; i < size; ++i)
   {
      _mat[i] = val;
   }
}

template<class T> inline
const T* Matrix<T>::getPointer() const
{
   return _mat;
}

template<class T> inline
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
{
   for (size_t row = 0; row < mat.getHeight(); ++row)
   {
      os << "| ";
      for (size_t col = 0; col < mat.getWidth(); ++col)
      {
         os << mat.get(row, col) << " ";
      }
      os << "|\n";
   }
   return os;
}
