/**********************************************************************************************
************************************ OPE COEFFICIENTS CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   IsingSpins.h
**********************************************************************************************/

#ifndef ISING_SPINS_H
#define ISING_SPINS_H

#include "MersenneTwister.h"

class IsingSpins 
{
  public:
    typedef unsigned int uint;
  
  private:
    uint   alpha_;   //number of replicas
    uint   N_;       //number of spins
    bool** spins_;   //array of the spin degrees of freedom
    
  public:
    IsingSpins(uint alpha, uint N);
    virtual ~IsingSpins();
    
    void flipSpin(uint a, uint i);
    int  getSpin(uint a, uint i);
    void print();
    void randomize(MTRand &randomGen);
};

#endif  // ISING_SPINS_H
