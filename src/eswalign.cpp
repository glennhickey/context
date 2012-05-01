//Copyright (C) 2010 by Glenn Hickey
//
//Released under the MIT license, see LICENSE.txt

#include <iostream>
#include <gsl/gsl_linalg.h>
#include "sw.h"
#include "pairalignment.h"

using namespace std;

int main(int argc, char** argv)
{
   if (argc == 3)
   {
      SW sw;
      PairAlignment pa;
      size_t i;
      sw.align(argv[1], argv[2], pa);
      cout << pa << endl;
   }
      
   return 0;
}
