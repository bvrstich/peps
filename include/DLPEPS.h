#ifndef DLPEPS_H
#define DLPEPS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class DLPEPS is a helper class, written for the contraction of two PEPS. DL stands for double layer
 */
template<typename T>
class DLPEPS : public vector< TArray<T,4> > {

   public:

      //constructor
      DLPEPS(const PEPS<T> &,const PEPS<T> &);

      //copy constructor
      DLPEPS(const DLPEPS &);

      //destructor
      virtual ~DLPEPS();

      int gD() const;

      const TArray<T,4> &operator()(int,int) const;

      TArray<T,4> &operator()(int,int);

   private:

      //!cutoff virtual dimension
      int D;


};

/**
 * output stream operator overloaded for DLPEPS<T> 
 */
template<typename T>
ostream &operator<<(ostream &output,const DLPEPS<T> &dlpeps_p){

   for(int r = 0;r < PEPS<T>::lat.gLy();++r)
      for(int c = 0;c < PEPS<T>::lat.gLx();++c){

         output << std::endl;
         output << "Tensor on site (" << r << "," << c << ")\t" << std::endl;
         output << std::endl;
         output << dlpeps_p(r,c) << std::endl;

      }

   return output;

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
