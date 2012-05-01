//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _SW_H
#define _SW_H

#include <string>
#include <vector>
struct PairAlignment;

class SW
{
public:
   SW();
   ~SW();
   void align(const std::string& seqA, const std::string& seqB, 
              PairAlignment& alignment);

   static const double InsStartCost;
   static const double InsContCost;
   static const double DelStartCost;
   static const double DelContCost;
   static const double MisMatchCost;
   static const double MatchCost;
   
protected:

   void init(const std::string& seqA, const std::string& seqB);
   void fillTable();
   void fillCell(size_t i, size_t j);

   double scoreDelete(size_t i, size_t j, size_t length);
   double scoreInsert(size_t i, size_t j, size_t length);
   double scoreMatch(size_t i, size_t j);

   void getAlignment(PairAlignment& alignment);
   
   struct Cell 
   {
      int _trace; // neg: insert, pos: delete, 0: match
      double _score;
   };

   std::vector<std::vector<Cell> > _matrix;
   std::string _seqA;
   std::string _seqB;
};

#endif
