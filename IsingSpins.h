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
    uint    alpha_;   //number of replicas
    uint    Ltau_;    //number of imaginary time slices
    uint    Nspat_;   //number of spins in each spatial slice in each replica
    bool*** spins_;   //array of the spin degrees of freedom
    
  public:
    IsingSpins(uint alpha, uint Ltau, uint Nspat);
    virtual ~IsingSpins();
    
    void flipSpin(uint a, uint t, uint i);
    int  getSpin (uint a, uint t, uint i);
    void print();
    void randomize(MTRand &randomGen);
};

#endif  // ISING_SPINS_H
