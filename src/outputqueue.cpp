//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>

#include "outputqueue.h"
#include "dnasubmodel.h"

using namespace std;
using namespace FTAGTableStructs;
using namespace FTAGStates;


OutputQueue::OutputQueue() : _isOpen(false)
{
   
}

OutputQueue::~OutputQueue() 
{
   reset();
}

void OutputQueue::reset()
{
   for (Alignment::iterator i = _align.begin(); i != _align.end(); ++i)
   {
      delete *i;
   }
   _align.clear();
   _eventSeq.clear();
	_isOpen = false;
}

void OutputQueue::openContext()
{
   assert(_isOpen == false);
   _insideInsPoint = _align.end();
}

void OutputQueue::closeContext()
{
   assert(_isOpen == true);
   _isOpen = false;
}

void OutputQueue::addSingle(State state, char a, char b)
{
   // only add single to outside
   assert(isOpen() == false);

   Column* col = new Column(a, b);
   Event event(state, col, NULL);

   _align.push_back(col);

   if (_insideInsPoint == _align.end())
   {
      --_insideInsPoint;
   }

   _eventSeq.push_back(event);
}

void OutputQueue::addDouble(State state, char a1, char a2, char b1, char b2)
{
   Column* col1 = new Column(a1, b1);
   Column* col2 = new Column(a2, b2);
   Event event(state, col1, col2);

   _align.insert(_insideInsPoint, col1);
   _align.push_back(col2); 

   if (_insideInsPoint == _align.end())
   {
      --_insideInsPoint;
   }

   _eventSeq.push_back(event);
}

pair<string, string> OutputQueue::getAlignment(bool keepGaps) const
{
   pair<string, string> stringPair;
   Alignment::const_iterator i;
   for (i = _align.begin(); i != _align.end(); ++i)
   {
      if (keepGaps || (*i)->_a != DNASubModel::gap())
      {
         stringPair.first += (*i)->_a;
      }
      if (keepGaps || (*i)->_b != DNASubModel::gap())
      {
         stringPair.second += (*i)->_b;
      }
   }
   assert(!keepGaps || stringPair.first.length() == stringPair.second.length());

   return stringPair;
}

deque<Trace> OutputQueue::getTrace() const
{
   deque<Trace> traceBack;
   for (size_t i = 0; i < _eventSeq.size(); ++i)
   {
      Trace trace;
      bool onFirst = true;
      trace._s = _eventSeq[i]._state;
      size_t pos = 0;
      Alignment::const_iterator j;
      for (j = _align.begin(); j != _align.end(); ++j)
      {
         if (onFirst && *j == _eventSeq[i]._col1)
         {
            trace._x1 = pos;
            onFirst = false;
         }
         else if (!onFirst && *j == _eventSeq[i]._col2)
         {
            trace._x2 = pos;
            trace._x3 = pos; // hack <-- to fix.
            break;
         }
         ++pos;
      }
      traceBack.push_back(trace);
   }
   return traceBack;
}
