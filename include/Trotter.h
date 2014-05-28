#ifndef TROTTER_H
#define TROTTER_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 13-05-2014\n\n
 * This class Trotter contains the two-site gates for the imaginary time evolution in the trotter decomposition
 */
class Trotter {

   public:

      static void set(double tau);

      //!actual operators: left
      static DArray<3> LO;

      //!actual operators: right
      static DArray<3> RO;

      //!timestep
      static double tau;

   private:

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
