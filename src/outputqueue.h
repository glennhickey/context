//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#ifndef _OUTPUTQUEUE
#define _OUTPUTQUEUE

#include <list>
#include <deque>
#include <vector>
#include <string>

#include "ftagstates.h"
#include "ftagtablestructs.h"

class OutputQueue
{
public:
   struct Column 
   {
      Column(char a, char b) : _a(a), _b(b) {}
      char _a;
      char _b;
   };

   typedef std::list<Column*> Alignment;

   struct Event
   {
      Event(FTAGStates::State state, Column* col1, Column* col2) :
      _state(state), _col1(col1), _col2(col2) {}
      FTAGStates::State _state;
      Column* _col1;
      Column* _col2;
   };

   typedef std::vector<Event> EventSequence;

public:

   OutputQueue();
   ~OutputQueue();

   void reset();
   void openContext();
   void closeContext();
   bool isOpen() const { return  _isOpen; }

   void addSingle(FTAGStates::State state, char a, char b);
   void addDouble(FTAGStates::State state, char a1, char a2, char b1, char b2);

   std::pair<std::string, std::string> getAlignment(bool keepGaps = true) const;
   std::deque<FTAGTableStructs::Trace> getTrace() const;
   
protected:

   Alignment _align;
   EventSequence _eventSeq;
   Alignment::iterator _insideInsPoint;
	bool _isOpen;
};

#endif
