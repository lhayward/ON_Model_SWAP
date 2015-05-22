/**********************************************************************************************
************************************ OPE COEFFICIENTS CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   VectorSpins.cpp
**********************************************************************************************/

#include <iostream>
#include "Vector_NDim.h"
#include "VectorSpins.h"

/**************** VectorSpins(uint alpha, uint N, uint spinDim) (constructor) *****************
* Input: alpha (number of replicas) 
*        N (number of spins)
*        spinDim (dimension of each spin, i.e. N of the O(N) model)
* This constructor initializes the array of vector spin degrees of freedom 
**********************************************************************************************/
VectorSpins::VectorSpins(uint alpha, uint N, uint spinDim)
{
  alpha_   = alpha;
  N_       = N;
  spinDim_ = spinDim;
  
  spins_ = new Vector_NDim**[alpha_];
  for( uint a=0; a<alpha_; a++ )
  { 
    spins_[a] = new Vector_NDim*[N];
    for( uint i=0; i<N_; i++ )
    { spins_[a][i] = NULL; }
  }
}

/******************************** ~VectorSpins() (destructor) ********************************/
VectorSpins::~VectorSpins()
{ 
  //delete the spins_ array:
  for(uint a=0; a<alpha_; a++)
  { 
    for( uint i=0; i<N_; i++ )
    { 
      if( spins_[a][i] != NULL )
      { delete spins_[a][i]; }
      spins_[a][i] = NULL;
    } //i
    if( spins_[a] != NULL )
    { delete[] spins_[a]; }
    spins_[a] = NULL; 
  } //i
  if( spins_ != NULL )
  { delete[] spins_; }
  spins_ = NULL;
}

/********************************** getSpin(uint a, uint i) **********************************/
Vector_NDim* VectorSpins::getSpin(uint a, uint i)
{
  /*Vector_NDim* result = NULL;
  
  if( (spins_ != NULL) && a<alpha_ && i<N_ )
  { result = spins_[a][i]; }
  else
  { 
    std::cout << "ERROR in VectorSpins::getSpin(uint a, uint i): NULL spins_ array or index "
              << "out of bounds" << std::endl; 
  }
  
  return result;*/
  return spins_[a][i];
}

/****************************************** print() ******************************************/
void VectorSpins::print()
{
  for( uint a=0; a<alpha_; a++ )
  {
    std::cout << "Replica #" << (a+1) << ":\n";
    for( uint i=0; i<N_; i++ )
    {
      std::cout << "  Spin " << i << ": ";
      spins_[a][i]->print();
    } //i
    std::cout << std::endl;
  } //a
} //print method

/******************************** randomize(MTRand* randomGen) *******************************/
void VectorSpins::randomize(MTRand &randomGen)
{
  for( uint a=0; a<alpha_; a++ )
  {
    for( uint i=0; i<N_; i++ )
    {
      //if there is already a Vector_NDim object in this location, then delete it first to
      //avoid memory leaks:
      if( spins_[a][i] != NULL )
      { delete spins_[a][i]; }
    
      spins_[a][i] = new Vector_NDim(spinDim_,randomGen); 
    } //i
  } //a
} //randomize method

/*********************** setSpin(uint a, uint i, Vector_NDim* newSpin) ***********************/
void VectorSpins::setSpin(uint a, uint i, Vector_NDim* newSpin)
{
  /*if( (newSpin != NULL) && a<alpha_ && i<N_ )
  { spins_[a][i] = newSpin; }
  else
  { 
    std::cout << "ERROR in VectorSpins::setSpin(uint a, uint i, Vector_NDim* newSpin): NULL"
              << " Vector_NDim object passed or index out of bounds" << std::endl; 
  }*/
  spins_[a][i] = newSpin;
}
