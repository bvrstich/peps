#ifndef PEPS_H
#define PEPS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class PEPS is a class written for the construction of projected entangled pair states on a rectangular lattice
 */
template<typename T,class Q>
class PEPS : public vector< QSTArray<T,5,Q> > {

   public:

      //constructor
      PEPS(int,int,int,int);

      //copy constructor
      PEPS(const PEPS &);

      //destructor
      virtual ~PEPS();

      int gLx() const;

      int gLy() const;

      int gd() const;

      int gD() const;

   private:

      //!x dimension
      int Lx;

      //!y dimension
      int Ly;

      //!physical dimension
      int d;

      //!cutoff virtual dimension
      int D;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
