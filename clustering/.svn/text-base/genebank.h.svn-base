#ifndef GENEBANK_H
#define GENEBANK_H

#include "genotype.h"

#include <map>
#include <string>

class Genebank {
private:
  map<int, Genotype*> m_genotypeMap;
  Genebank( const Genebank & );
  const Genebank & operator=( const Genebank & );
public:
  Genebank();
  virtual ~Genebank();

  /**
   * Resets the genebank to an empty state.
   **/
  void clear();

  void addGenotype( Genotype* g );
  void removeGenotype( Genotype* g );

  /**
   * Creates a new genotype with the given characteristics. The genotype
   * does not need to be added.
   **/
  Genotype* createGenotype( int id, int parentId, int treeDepth, int
			      birth, int death, string genes );
  
  /**
   * Sets up the correct parent pointers.
   **/
  void setupParentPointers();

  /**
   * Reliably finds the last coalescent genotype if the given genotype is from
   * the current generation.
   **/
  bool checkCoalescence( Genotype *g );

  /**
   * Calculates the Hamming distance between two genotypes.
   **/
  int calcTreeDist( const Genotype *g1, const Genotype *g2 ) const;
  int calcHammingDist( const Genotype *g1, const Genotype *g2 ) const;
  int calcMRCADist( const Genotype *g1, const Genotype *g2 ) const;

  Genotype* getGenotype( int id ) const {
    return (*m_genotypeMap.find( id )).second; }

  
  /**
   * Print out the state of the genebank.
   **/
  void print() const;
};

#endif
