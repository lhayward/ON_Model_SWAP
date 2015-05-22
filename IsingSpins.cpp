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

/************************ IsingSpins(uint alpha, uint N) (constructor) ************************
* Input: alpha (number of replicas)
*        N     (number of spins)
* This constructor initializes the array of spin degrees of freedom 
**********************************************************************************************/
IsingSpins::IsingSpins(uint alpha, uint N)
{
  alpha_ = alpha;
  N_     = N;
  
  spins_ = new bool*[alpha_];
  for( uint a=0; a<alpha_; a++ )
  { 
    spins_[a] = new bool[N_];
    for( uint i=0; i<N_; i++ )
    { spins_[a][i] = 0; }
  }
}

/********************************* ~IsingSpins() (destructor) ********************************/
IsingSpins::~IsingSpins()
{   
  //delete the spins_ array:
  for(uint a=0; a<alpha_; a++)
  { 
    if( spins_[a] != NULL )
    { delete[] spins_[a]; }
    spins_[a] = NULL; 
  }
  if( spins_ != NULL )
  { delete[] spins_; }
  spins_ = NULL;
}

/************************************** flipSpin(uint i) *************************************/
void IsingSpins::flipSpin(uint a, uint i)
{
  spins_[a][i] = !spins_[a][i];
}

/*************************************** getSpin(int i) **************************************/
int IsingSpins::getSpin(uint a, uint i)
{
  return (2*spins_[a][i] - 1);
}

/****************************************** print() ******************************************/
void IsingSpins::print()
{
  for( uint a=0; a<alpha_; a++ )
  {
    std::cout << "Replica #" << (a+1) << ": ";
    for( uint i=0; i<N_; i++ )
    { 
      //print an extra space if spin i is in the +1 state:
      if( spins_[a][i] )
      {  std::cout << " "; }
    
      std::cout << (2*spins_[a][i] - 1) << " "; 
    } //i
    std::cout << std::endl;
  }//a
} //print method

/******************************** randomize(MTRand* randomGen) *******************************/
void IsingSpins::randomize(MTRand &randomGen)
{
  for( uint a=0; a<alpha_; a++ )
  {
    for( uint i=0; i<N_; i++ )
    { spins_[a][i] = randomGen.randInt(1); } //end of loop over i
  } //a
} //randomize method
