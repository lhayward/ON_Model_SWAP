/**********************************************************************************************
************************************ OPE COEFFICIENTS CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   ON_Model.h (Abstract Class)
**********************************************************************************************/

#ifndef ON_MODEL_H
#define ON_MODEL_H

#include <stdlib.h>
#include <string>
#include "Hyperrectangle.h"
#include "Measure.h"
#include "MersenneTwister.h"

class ON_Model 
{ 
  public:
    typedef unsigned int  uint;
    
  protected:
    uint            spinDim_; //dimension of the spins on the lattice (i.e. N of the O(N) model)
    double          J_;       //coupling 
    double          h_;       //field
    double          T_;       //current temperature
    Measure         measures; //observables to record
    bool            warmupDone; //has the system finished warming up?
    std::ofstream   fout;
    Hyperrectangle* hrect_; //the hyperrectangular lattice on which the d.o.f. live
    bool*           inCluster_;  //boolean array for wolffUpdate indicating whether or not each
                                 //spin is in the cluster (redundant information to cluster 
                                 //vec, but stored for efficiency purposes)
    bool***         inClusterNew_;  //new version of inCluster_ array: stores all three 
                                    //coordinates of spins in cluster
    
    //parameters of the lattice:
    uint    N_;       //total number of spins
    uint    alpha_;   //number of replicas (Renyi index)
    uint    Ltau_;    //lattice's linear length in the imaginary time direction
    uint*   L_;       //lattice's linear length in each spatial dimension
    uint    Dspat_;   //spatial dimension
    uint    Nspat_;   //number of spins living on the spatial hyperrectangle
    
    //acceptance rates for the two kinds of updates:
    int numAccept_local_;
    int numAccept_clust_;
    
    //variables related to storing the cluster sizes to generate a histogram:
    static const bool writeClusts_ = false;  //should we write data for the cluster histograms?
    uint*             clustSizes_;
    uint*             clustSizes_accepted_;
    uint*             clustSizes_rejected_;
    
  public:
    ON_Model(std::ifstream* fin, std::string outFileName, Hyperrectangle* lattice);
    virtual ~ON_Model();
    
    //methods implemented in ON_Model class:
    void writeBin           (int binNum, int numMeas, int sweepsPerMeas);
    void writeClustHistoData(std::string fileName);
    void zeroMeasurements   ();
    
    //methods that can be overwritten by child classes:
    virtual void markWarmupDone  ();
    virtual void printParams     ();
    virtual void changeT         (double newT);
    
    //pure virtual methods (to be implemented by all child classes):
    virtual void localUpdate        (MTRand &randomGen) = 0;
    virtual void makeMeasurement    () = 0;
    virtual void printSpins         () = 0;
    virtual void randomizeLattice   (MTRand &randomGen) = 0;
    virtual void sweep              (MTRand &randomGen, bool pr) = 0;
    
  protected:
    void clearCluster        (std::vector<uint> &cluster);
    uint getNeighRep_minusTau(uint a, uint t, uint i);
    uint getNeighRep_plusTau (uint a, uint t, uint i);
    void randomizeCoords     (MTRand &randomGen, uint& a, uint& t, uint& i);
    uint uintPower           (uint base, uint exp);
};

#endif  // ON_MODEL_H
