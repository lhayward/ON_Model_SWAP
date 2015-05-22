/**********************************************************************************************
************************************ OPE COEFFICIENTS CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   VectorSpins.h
**********************************************************************************************/

#ifndef VECTOR_SPINS_H
#define VECTOR_SPINS_H

#include "MersenneTwister.h"
#include "Vector_NDim.h"

class VectorSpins 
{
  public:
    typedef unsigned int uint;
  
  private:
    uint           alpha_;   //number of replicas
    uint           N_;       //number of spins
    uint           spinDim_; //dimensionality of each vector
    Vector_NDim*** spins_;   //array of the vector spin degrees of freedom
    
  public:
    VectorSpins(uint alpha, uint N, uint spinDim);
    virtual ~VectorSpins();
    
    Vector_NDim* getSpin(uint a, uint i);
    void         print();
    void         randomize(MTRand &randomGen);
    void         setSpin(uint a, uint i, Vector_NDim* newSpin);
};

#endif  // VECTOR_SPINS_H
