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
   int D_aux = atoi(argv[4]);//virtual dimension

   bool update = true;
   double tau = 0.01;
   int n_steps = 10;

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,tau);

   double f = 0.74;

   PEPS<double> peps;
   peps.initialize_jastrow(f);
   peps.normalize();

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d/full/D_aux=%d.txt",L,L,D,D_aux);

   ofstream out(filename);
   out.precision(16);

   global::env.calc('A',peps);

   double ener = peps.energy();

   out << 0 << "\t" << ener << endl;

   for(int i = 1;i < 1000;++i){

      propagate::step(true,peps,n_steps);

      if(i % 10 == 0){

         peps.normalize();
         global::env.calc('A',peps);

         double tmp = peps.energy();

         if(tmp > ener)
            break;
         else
            ener = tmp;

         out << i << "\t" << tmp << endl;

         char dir_name[200];
         sprintf(dir_name,"output/%dx%d/D=%d/full/D_aux=%d.txt",L,L,D,D_aux);

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i = 1000;i < 2000;++i){

      propagate::step(true,peps,n_steps);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);

         double tmp = peps.energy();

         if(tmp > ener)
            break;
         else
            ener = tmp;

         out << i << "\t" << tmp << endl;

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i = 2000;i < 3000;++i){

      propagate::step(true,peps,n_steps);

      if(i % 10 == 0){

         peps.normalize();
         global::env.calc('A',peps);

         double tmp = peps.energy();

         if(tmp > ener)
            break;
         else
            ener = tmp;

         out << i << "\t" << peps.energy() << endl;

      }

   }

   return 0;

}
