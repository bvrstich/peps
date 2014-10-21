/**
 * @mainpage 
 * This is an implementation of an imaginary time evolution algorithm for the optimization of PEPS.
 * @author Brecht Verstichel
 * @date 25-03-2014
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);//dimension of the lattice: LxL
   int d = atoi(argv[2]);//physical dimension
   int D = atoi(argv[3]);//virtual dimension
   int D_aux = atoi(argv[4]);//auxiliary dimension for the contraction

   bool update = true;

   double tau = 0.001;

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,tau);

   double f = 0.74;

   PEPS<double> peps;
   peps.initialize_jastrow(f);
   peps.normalize();

   global::env.calc('A',peps);
   global::env.test();

   char filename[200];

   if(update)
      sprintf(filename,"output/%dx%d/D=%d/full/D_aux=%d.txt",L,L,D,D_aux);
   else 
      sprintf(filename,"output/%dx%d/D=%d/simple/D_aux=%d.txt",L,L,D,D_aux);

   ofstream out(filename);
   out.precision(16);

   for(int i = 0;i < 10000;++i){

      propagate::step(update,peps);

      cout << i << endl;

      if(i % 100 == 0){

         global::env.calc('A',peps);
         out << i << "\t" << peps.energy() << endl;

         char peps_dir[200];

         if(update)
            sprintf(peps_dir,"output/%dx%d/D=%d/full/peps",L,L,D);
         else
            sprintf(peps_dir,"output/%dx%d/D=%d/simple/peps",L,L,D);

         peps.save(peps_dir);

      }

   }

   return 0;

}
