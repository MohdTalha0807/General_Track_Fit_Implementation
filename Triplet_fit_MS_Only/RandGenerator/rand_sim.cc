#include "Libraries.h"
#include <bits/stdc++.h>
#define PI 3.1415926
#define B_avg 2.0   // in Tesla


using namespace std;


int main()
{
   gROOT->SetStyle("Plain");
   
   double x,y,z,phi;
   int nevent = 10000;
   
   TRandom3  *rand1 = new TRandom3(0);          // 0 is there for random seeds
   TRandom3  *rand2 = new TRandom3(0);   
   TRandom3  *rand3 = new TRandom3(0);   
  // TRandom3  *rand4 = new TRandom3(0);   
  // TRandom3  *rand5 = new TRandom3(0);   
   
   ofstream fout;
   fout.open ("triplet_sim.txt");
   //fout<<"x"<<"\t"<<"\t"<<"y"<<"\t"<<"\t"<<"z"<<"\t"<<"\t"<<"phi"<<endl;       
   
   
   for(int i=1; i <= nevent; i++)
  { 
   phi = rand1->Rndm() * 2.0 * M_PI;      // all generated in the units of mm
   x = -2.0 + rand1-> Rndm() * 4.0;
   y = -2.0 + rand2-> Rndm() * 4.0;
   z = -50.0 + rand3-> Rndm() * 100.0;
   
   fout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<phi<<endl; 
  }

  fout.close();  
   
  return 0;
}


 /*  phi = rand1->Rndm() * 2.0 * M_PI;      // all generated in the units of mm
   x = -2.0 + rand1-> Rndm() * 4.0;
   y = -2.0 + rand2-> Rndm() * 4.0;
   z = -50 + rand3-> Rndm() * 100.0;
  */ 
