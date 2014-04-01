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

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);//dimension of the lattice: LxL
   int d = atoi(argv[2]);//physical dimension
   int D = atoi(argv[3]);//virtual dimension
/*
   //initialize the dimensions
   PEPS<double>::lat.set(L,L,d);

   PEPS<double> peps1(D);
   PEPS<double> peps2(D);

   MPS<double> mps(peps1,peps2);
   MPO<double> mpo(1,peps1,peps2);

   mps.gemv('L',mpo);
*/
   
   TArray<double,3> A(10,2,10);
   A.generate(PEPS<double>::rgen);

   ofstream out("orig.out");
   out.precision(15);
   out << A << endl;

   TArray<double,2> R;

   Geqrf(A,R);

   TArray<double,3> tmp(10,2,10);

   enum {j,k,l,m};

   Contract(1.0,A,shape(j,k,l),R,shape(l,m),0.0,tmp,shape(j,k,m));

   cout << tmp << endl;

}
