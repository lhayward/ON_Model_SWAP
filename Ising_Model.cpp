/**********************************************************************************************
************************************ OPE COEFFICIENTS CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   Ising_Model.cpp
**********************************************************************************************/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>
#include "FileReading.h"
#include "Ising_Model.h"

//typdef needed because uint is a return types:
typedef Ising_Model::uint uint;

/**************** Ising_Model(std::ifstream* fin, std::string outFileName, ... ****************
*****************             Hyperrectangle* lattice, MTRand* randomGen)      ****************
*****************                        (constructor)                         ***************/
Ising_Model::Ising_Model(std::ifstream* fin, std::string outFileName, Hyperrectangle* lattice, 
                         MTRand &randomGen)
  : ON_Model(fin, outFileName, lattice)
{
  std::cout.precision(15);
  
  spinDim_ = 1;
  spins_ = new IsingSpins(alpha_, Ltau_, Nspat_);
  randomizeLattice(randomGen);
}

/******************************** ~Ising_Model() (destructor) ********************************/
Ising_Model::~Ising_Model()
{
  //delete the IsingSpins object:
  if( spins_ != NULL )
  { delete spins_; }
  spins_ = NULL;
} //destructor

/**************************************** flipCluster *****************************************
* This function flips all spins in the passed cluster.
**********************************************************************************************/
void Ising_Model::flipCluster(std::vector<uint> &cluster)
{
  uint clustSize = (uint)cluster.size();
  
  for( uint i=0; i<clustSize; i++ )
  { spins_->flipSpin(0, 0, cluster[i]); }  //!!!Change zeros
} //flipCluster

/**************************************** getEnergy() ****************************************/
double Ising_Model::getEnergy()
{
  double  energyJ       = 0;
  double  energyh       = 0;
  int     currSpin;
  int     neighbour;
  uint    nextTau;
  
  //loops over replicas:
  for( uint a=0; a<alpha_; a++ )
  {
    //loop over imaginary time:
    for( uint t=0; t<Ltau_; t++ )
    {
      //loop over spatial indices:
      for( uint i=0; i<Nspat_; i++ )
      { 
        currSpin    = spins_->getSpin(a,t,i);
    
        //nearest neighbour term in spatial direction:
        for( uint d=0; d<Dspat_; d++ )
        {
          neighbour = spins_->getSpin(a,t,hrect_->getNeighbour(i,d)); //nearest neighbour along
                                                                      //d direction
          energyJ   += currSpin*neighbour;
        } //d
        
        nextTau = (t+1)%Ltau_;
        neighbour = spins_->getSpin( getNeighRep_plusTau(a,t,i),nextTau,i );
        energyJ   += currSpin*neighbour;
    
        //field term:
        energyh += currSpin;
      } //i
    } //t
  } //a
  
  return ((-1.0*J_*energyJ) + (-1.0*h_*energyh)) ;
} //getEnergy()

/*************************************** getSigmaTot() ***************************************/
int Ising_Model::getSigmaTot()
{
  int sigmaTot = 0;
  
  for( uint i=0; i<Nspat_; i++ )
  { sigmaTot += spins_->getSpin(0, 0, i); } //!!!Change zeros
  
  return sigmaTot;
}

/******************************* localUpdate(MTRand* randomGen) ******************************/
void Ising_Model::localUpdate(MTRand &randomGen)
{ 
  //randomly selected spin location:
  uint index_a;   //replica index
  uint index_t;   //tau index
  uint index_i;   //spatial index
  
  uint neighTau;  //tau index for neighbour
  uint neighRep;  //replica index for neighbour
  double deltaE;
  int    spin_old;    //previous state of spin (at randomly selected lattice site)
  int    nnSum = 0;
  
  //randomly select a spin on the lattice:
  randomizeCoords(randomGen, index_a, index_t, index_i);
  spin_old = spins_->getSpin(index_a, index_t, index_i);
  
  //Calculate the nearest neighbour sum:
  //Loop over spatial neighbours:
  for( uint d=0; d<Dspat_; d++ )
  { 
    nnSum += ( spins_->getSpin( index_a, index_t, hrect_->getNeighbour(index_i, d       ) ) 
             + spins_->getSpin( index_a, index_t, hrect_->getNeighbour(index_i, d+Dspat_) ) ); 
  }
  //Add in the contribution from +tau neighbour:
  neighRep = getNeighRep_plusTau (index_a, index_t, index_i);
  neighTau = (index_t+1)%Ltau_;
  nnSum += spins_->getSpin( neighRep, neighTau, index_i );
  
  //Add in the contribution from -tau neighbour:
  neighRep = getNeighRep_minusTau(index_a, index_t, index_i);
  neighTau = (index_t-1+Ltau_)%Ltau_;
  nnSum += spins_->getSpin( neighRep, neighTau, index_i );
  
  //calculate the energy change for the proposed move:
  deltaE = 2.0*J_*spin_old*nnSum + 2.0*h_*spin_old;
  
  //if the move is accepted:
  if( deltaE<=0 || randomGen.randDblExc() < exp(-deltaE/T_) )
  { 
    spins_->flipSpin( index_a, index_t, index_i );
    numAccept_local_++;
  }
}

/************************************* makeMeasurement() *************************************/
void Ising_Model::makeMeasurement()
{
  double energyPerSpin = getEnergy()/(1.0*N_);
  
  measures.accumulate( "E",     energyPerSpin ) ;
  measures.accumulate( "ESq",   pow(energyPerSpin,2) );
}

/**************************************** printSpins() ***************************************/
void Ising_Model::printSpins()
{ spins_->print(); }

/**************************** randomizeLattice(MTRand* randomGen) ****************************/
void Ising_Model::randomizeLattice(MTRand &randomGen)
{ 
  spins_->randomize(randomGen);
}

/********************************** sweep(MTRand* randomGen) *********************************/
void Ising_Model::sweep(MTRand &randomGen, bool pr)
{ 
  uint N1 = N_/2;     //number of local updates before Wolff step
  uint N2 = N_ - N1;  //number of local updates after Wolff step
  
  for( uint i=0; i<N1; i++ )
  { localUpdate(randomGen); }
  
  //wolffUpdate(randomGen, pr);
  
  for( uint i=0; i<N2; i++ )
  { localUpdate(randomGen); }
}

/******************************* wolffUpdate(MTRand* randomGen) ******************************/
void Ising_Model::wolffUpdate(MTRand &randomGen, bool pr)
{
  //randomly selected spin location:
  uint index_a;   //replica index
  uint index_t;   //tau index
  uint index_i;   //spatial index
  
  uint               latticeSite;
  uint               neighSite;
  uint               clustSize;
  int                clusterState;
  double             PAdd = 1.0 - exp(-2.0*J_/T_);
  double             deltaE_onsite;
  std::vector<uint>  buffer;  //indices of spins to try to add to cluster (will loop until 
                              //buffer is empty)
  std::vector<uint>  cluster; //vector storing the indices of spins in the cluster
  
  buffer.reserve(N_);
  cluster.reserve(N_);
  
  //randomly select a spin on the lattice:
  randomizeCoords(randomGen, index_a, index_t, index_i);
  clusterState = spins_->getSpin(index_a, index_t, index_i);
  inClusterNew_[index_a][index_t][index_i] = 1;
  cluster.push_back(latticeSite); //!
  buffer.push_back(latticeSite);  //!
  
  while( !buffer.empty() )
  {
    latticeSite = buffer.back();
    buffer.pop_back();
    
    //loop over spatial neighbours:
    for( uint d=0; d<(2*Dspat_); d++ )
    {
      neighSite = hrect_->getNeighbour( index_i, d );
      
      if ( !( inClusterNew_[index_a][index_t][neighSite] ) 
           && ( spins_->getSpin(index_a, index_t, neighSite) == clusterState ) 
           && (randomGen.randDblExc() < PAdd) )
      { 
        inClusterNew_[index_a][index_t][neighSite] = 1;
        cluster.push_back( neighSite ); //!
        buffer.push_back( neighSite );  //!
      }
    }  //for loop over spatial neighbours
    
    //!add loop over tau neighbours:
  } //while loop for buffer
  
  clustSize = cluster.size();
  deltaE_onsite = 2*h_*clusterState*clustSize;
  
  if( writeClusts_ )
  { clustSizes_[clustSize-1]++; }
  
  //If the onsite energy diff. is negative, accept the cluster move. If the onsite energy 
  //diff. is positive, accept the cluster move with probability exp^(-beta*onsiteEnery_diff).  
  if( deltaE_onsite <= 0 || randomGen.randDblExc() < exp(-1.0*deltaE_onsite/T_) )
  {
    flipCluster(cluster);
    numAccept_clust_++;
    
    if( writeClusts_ )
    { clustSizes_accepted_[clustSize-1]++; }
    
    if(pr)
    {
      if( deltaE_onsite <=0 )
      { std::cout << "  onsite <= 0" << std::endl; }
      else
      { std::cout << "  PAccept = " << exp(-1.0*deltaE_onsite/T_) << std::endl; }
      std::cout << "  size = " << clustSize << std::endl << std::endl;
    }
  } //if
  else if( writeClusts_ )
  { clustSizes_rejected_[clustSize-1]++; }
  
  clearCluster(cluster);
} //wolffUpdate()
