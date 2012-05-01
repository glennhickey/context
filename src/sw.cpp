//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <cassert>
#include <limits>
#include <iostream>

#include "sw.h"
#include "pairalignment.h"


using namespace std;

const double SW::InsStartCost = 4;
const double SW::InsContCost = 1;
const double SW::DelStartCost = 4;
const double SW::DelContCost = 1;
const double SW::MisMatchCost = 1;
const double SW::MatchCost = 0;

SW::SW() : _matrix(NULL) 
{

}

SW::~SW()
{
}

void SW::align(const string& seqA, const string& seqB, PairAlignment& alignment)
{
   init(seqA, seqB);
   fillTable();
   getAlignment(alignment);
}

void SW::init(const string& seqA, const string& seqB)
{
   _seqA = seqA;
   _seqB = seqB;
   _matrix.resize(_seqA.length() );
   for (size_t i = 0; i < _matrix.size(); ++i)
   {
      _matrix[i].resize(_seqB.size());
      for (size_t j = 0; j < _matrix[i].size(); ++j)
      {
         _matrix[i][j]._score = 0.;
      }
   }
}

void SW::fillTable()
{
   size_t i, j;
   for (i = 1; i < _matrix.size(); ++i)
   {
      for (j = 1; j < _matrix[i].size(); ++j)
      {
         fillCell(i, j);
      }
   }
}

void SW::fillCell(size_t i, size_t j)
{
   size_t k;
   Cell& cell = _matrix[i][j];
   cell._score = numeric_limits<double>::max();
   double score; 
   
   //match
   if (i == 0 && j == 0)
   {
      cell._score = scoreMatch(i, j);
      cell._trace = 0;
   }
   if (i > 0 && j > 0)
   {
      score = _matrix[i-1][j-1]._score + scoreMatch(i, j);
      if (score < cell._score)
      {
         cell._score = score;
         cell._trace = 0;
      }
   }
   //del
   for (k = 0; k < i; ++k)
   {
      score = _matrix[i-k][j]._score + scoreDelete(i, j, i-k); 
      if (score < cell._score)
      {
         cell._score = score;
         cell._trace = i-k;
      }
   }
   //ins
   for (k = 0; k < j; ++k)
   {
      score = _matrix[i][j-k]._score + scoreInsert(i, j, j-k);
      if (score < cell._score)
      {
         cell._score = score;
         cell._trace = -(j-k);
      }
   }
   if (cell._trace)
   cout << "(" << i << "," << j << ") <- " << cell._trace << endl;
}

// score would be its own class is a full implementaiton
double SW::scoreDelete(size_t i, size_t j, size_t length)
{
   return DelStartCost + DelContCost * length;
}

double SW::scoreInsert(size_t i, size_t j, size_t length)
{
   return InsStartCost + InsContCost * length;
}

double SW::scoreMatch(size_t i, size_t j)
{
   return _seqA[i] == _seqB[j] ? MatchCost : MisMatchCost;
}

void SW::getAlignment(PairAlignment& alignment)
{
   alignment._seqA.clear();
   alignment._seqB.clear();
   int i = _seqA.length() - 1;
   int j = _seqB.length() - 1;
   
   while (i >= 0 && j >= 0)
   {
      int trace = _matrix[i][j]._trace;
      int k;
      //del
      for (k = 0; k < trace; ++k)
      {
         alignment._seqA.push_front(_seqA[i--]);
         alignment._seqB.push_front('-');
      }
      //ins
      for (k = 0; k > trace; --k)
      {
         alignment._seqA.push_front('-');
         alignment._seqB.push_front(_seqB[j--]);
      }
      //match
      if (trace == 0)
      {
         alignment._seqA.push_front(_seqA[i--]);
         alignment._seqB.push_front(_seqB[j--]);
      }
   }
   if (i == 0 && j == 0)
   {
      alignment._seqA.push_front(_seqA[i--]);
      alignment._seqB.push_front(_seqB[j--]);

   }

   while (i >= 0)
   {
      alignment._seqA.push_front(_seqA[i--]);
      alignment._seqB.push_front('-');
   }
   while (j >= 0)
   {
      alignment._seqA.push_front('-');
      alignment._seqB.push_front(_seqB[j--]);
   }
}

