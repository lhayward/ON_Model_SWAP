/**********************************************************************************************
************************************ OPE COEFFICIENTS CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   IsingSpins.cpp
*
* Ising spins will be stored as boolean variables (0 or 1) but returned as integers 
* (-1 or +1, respectively).
**********************************************************************************************/

#include <iostream>
#include "IsingSpins.h"

/**************** IsingSpins(uint alpha, uint Ltau, uint Nspat) (constructor) *****************
* Input: alpha (number of replicas)
*        Ltau  (number of imaginary time slices)
*        Nspat (number of spins in each spatial slice in each replica)
* This constructor initializes the array of spin degrees of freedom 
**********************************************************************************************/
IsingSpins::IsingSpins(uint alpha, uint Ltau, uint Nspat)
{
  alpha_ = alpha;
  Ltau_  = Ltau;
  Nspat_ = Nspat;
  
  spins_ = new bool**[alpha_];
  for( uint a=0; a<alpha_; a++ )
  { 
    spins_[a] = new bool*[Ltau_];
    for( uint t=0; t<Ltau_; t++ )
    {
      spins_[a][t] = new bool[Nspat_];
      for( uint i=0; i<Nspat_; i++ )
      { spins_[a][t][i] = 0; }
    } //t
  } //a
}

/********************************* ~IsingSpins() (destructor) ********************************/
IsingSpins::~IsingSpins()
{ 
  //delete the spins_ array:
  for(uint a=0; a<alpha_; a++)
  { 
    for( uint t=0; t<Ltau_; t++ )
    { 
      if( spins_[a][t] != NULL )
      { delete[] spins_[a][t]; }
      spins_[a][t] = NULL;
    } //t
    
    if( spins_[a] != NULL )
    { delete[] spins_[a]; }
    spins_[a] = NULL; 
  } //a
  
  if( spins_ != NULL )
  { delete[] spins_; }
  spins_ = NULL;
}

/****************************** flipSpin(uint a, uint t, uint i) *****************************/
void IsingSpins::flipSpin(uint a, uint t, uint i)
{
  spins_[a][t][i] = !spins_[a][t][i];
}

/****************************** getSpin(uint a, uint t, uint i) ******************************/
int IsingSpins::getSpin(uint a, uint t, uint i)
{
  return (2*spins_[a][t][i] - 1);
}

/****************************************** print() ******************************************/
void IsingSpins::print()
{
  for( uint a=0; a<alpha_; a++ )
  {
    std::cout << "Replica #" << (a+1) << ": ";
    for( uint t=0; t<Ltau_; t++ )
    {
      for( uint i=0; i<Nspat_; i++ )
      { 
        //print an extra space if spin i is in the +1 state:
        if( spins_[a][t][i] )
        {  std::cout << " "; }
    
        std::cout << (2*spins_[a][t][i] - 1) << " "; 
      } //i
      std::cout << std::endl;
    } //t
    std::cout << std::endl;
  }//a
} //print method

/******************************** randomize(MTRand* randomGen) *******************************/
void IsingSpins::randomize(MTRand &randomGen)
{
  for( uint a=0; a<alpha_; a++ )
  {
    for( uint t=0; t<Ltau_; t++ )
    {
      for( uint i=0; i<Nspat_; i++ )
      { spins_[a][t][i] = randomGen.randInt(1); } //end of loop over i
    } //t
  } //a
} //randomize method
