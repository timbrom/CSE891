// genebank.cpp

#include "genebank.h"

#include "genotype.h"

#include <assert.h>
#include <iostream>
#include <stdlib.h>

#include <map>
#include <string>


Genebank::Genebank()
{
  clear();
}

Genebank::~Genebank()
{
  clear();
}

void Genebank::clear()
{
  map<int, Genotype*>::iterator it = m_genotypeMap.begin();
  for ( ; it!=m_genotypeMap.end(); it++ )
    delete (*it).second;
  m_genotypeMap.clear();
}

void Genebank::addGenotype( Genotype* g )
{
  g->incrementCount();
}

void Genebank::removeGenotype( Genotype *g )
{
  if ( g==0 ) return;

  if ( g->decrementCount() ) {
    // this genotyp has lost all its references. We have to remove it.
    m_genotypeMap.erase( g->getId() );

    // recursively remove all reference counts from parents
    removeGenotype( g->getParent() );
    delete g;
  }
}

Genotype* Genebank::createGenotype( int id, int parentId, int
				    treeDepth, int birth, int
				    death, string genes )
{
  Genotype *g = new Genotype( 0, parentId, treeDepth, id, birth,
			      death, genes );
  m_genotypeMap[id] = g;

  return g;
}

void Genebank::setupParentPointers()
{
  //  cout << "void Genebank::setupParentPointers() called\n";
  map<int, Genotype*>::iterator it = m_genotypeMap.begin();
  map<int, Genotype*>::iterator parent;

  for ( ; it != m_genotypeMap.end(); it++ ){
    int parentId = (*it).second->getParentId();
    if ( parentId >= 0 ){
      parent = m_genotypeMap.find( parentId );
      (*it).second->setParent( (*parent).second );
      (*parent).second->incrementCount();
    }
  }
}


bool Genebank::checkCoalescence( Genotype *g )
{
  if ( g == 0 || g->isCoalescent() ){
    //    cout << g->getId() << " is already coalescent." << endl;
    return true;
  }
  //else
  //  cout << g->getId() << " is not coalescent. Checking parent..." << endl;


  Genotype *parent = g->getParent();

  assert( g != 0 );

  if ( checkCoalescence( parent ) ){
    // if parent is coalescent and has only a count of 1, then we are
    // coalescent as well

    //  cout << "parent " << parent->getId() << " is coalescent, checking count..." << endl;
    if ( parent == 0 || parent->getCount() == 1 ){
      // cout << "count is 1. " << g->getId() << " is coalescent." << endl;
      g->setCoalescent();
      return true;
    }
    else{
      //cout << "count is > 1. " << g->getId() << " is not coalescent." << endl;
     return false;
    }
  }
  else{
    //cout << "parent " << parent->getId() << " is not coalescent. We are done." << endl;
    return false;
  }
}


int Genebank::calcTreeDist( const Genotype *g1, const Genotype *g2 ) const
{
  const Genotype *g, *gmut;

  assert( g1!=0 );
  assert( g2!=0 );

  // first we tag all parents of g1 till the last coalescent ancestor
  g = g1;
  while( !g->isCoalescent() ){
    //    cout << "Tagging " << g->getId() << endl;
    g->setTagged();
    g = g->getParent();
    assert( g!= 0 );
  }
  //  cout << "Tagging " << g->getId() << endl;
  g->setTagged();

  // now we search from g2 backwards till we find a tagged ancestor
  g = g2;
  while( !g->isTagged() ){
    //    cout << "Adding dist " << g->getParentDist() << " from " << g->getId() << endl;
    g = g->getParent();
    assert( g!= 0 );
  }
  gmut = g;
  
  //  cout << "Common ancestor of " << g1->getId() << " and " << g2->getId()
  //       << " is " << gmut->getId() << endl;

  // now we untag everything again
  g = g1;
  while( !g->isCoalescent() ){
    //    cout << "Untag " << g->getId() << endl;
    g->setTagged( false );
    g = g->getParent();
    assert( g!= 0 );
  }
  g->setTagged( false );
  //  cout << "Untag " << g->getId() << endl;
  return g1->getTreeDepth() + g2->getTreeDepth() - 2*gmut->getTreeDepth();
}

int Genebank::calcHammingDist( const Genotype *g1, const Genotype *g2 ) const
{
  //  const Genotype *g, *gmut;

  assert( g1!=0 );
  assert( g2!=0 );

  if ( g1->getLength() != g2->getLength() ) {
    cout << "Error: length mismatch while calculating Hamming distance.";
    return -1;
  }

  int dist = 0;
  for ( int i = 0; i < g1->getLength(); i++ ) {
    if ( g1->getGene( i ) != g2->getGene( i ) ) {
      dist++;
    }
  }

  return dist;
}


int Genebank::calcMRCADist( const Genotype *g1, const Genotype *g2 ) const
{
  const Genotype *g, *gmut, *gmrca;
  int birth1 = 0;
  int birth2 = 0;
  gmrca = g1;

  assert( g1!=0 );
  assert( g2!=0 );

  // first we tag all parents of g1 till the last coalescent ancestor
  g = g1;
  while( !g->isCoalescent() ){
    //    cout << "Tagging " << g->getId() << endl;
    g->setTagged();
    g = g->getParent();
    assert( g!= 0 );
  }
  //  cout << "Tagging " << g->getId() << endl;
  g->setTagged();

  // now we search from g2 backwards till we find a tagged ancestor
  g = g2;
  birth2 = g->getBirth();
  while( !g->isTagged() ){
    //    cout << "Adding dist " << g->getParentDist() << " from " <<
    //    g->getId() << endl;
    birth2 = g->getBirth();
    g = g->getParent();
    assert( g!= 0 );
  }
  gmut = g;
  gmrca = g;
  
  //  cout << "Common ancestor of " << g1->getId() << " and " << g2->getId()
  //       << " is " << gmut->getId() << endl;

  // now we find the offspring of the MRCA
  g = g1;
  birth1 = g->getBirth();
  while( g != gmrca ){
    //    cout << "Untag " << g->getId() << endl;
    birth1 = g->getBirth();
    g = g->getParent();
  }

  // now we untag everything again
  g = g1;
  while( !g->isCoalescent() ){
    //    cout << "Untag " << g->getId() << endl;
    g->setTagged( false );
    g = g->getParent();
    assert( g!= 0 );
  }
  g->setTagged( false );
  //  cout << "Untag " << g->getId() << endl;
  // return (g1->getBirth() + g2->getBirth() - 2*(birth1 <? birth2)); THB - This is the original code, but I have no idea what the <? operator is supposed to do.
  return (g1->getBirth() + g2->getBirth() - birth1 - birth2); // THB - I think this is the correct replacement for the above code

  //  return g1->getTreeDepth() + g2->getTreeDepth() - 2*gmut->getTreeDepth();
}

void Genebank::print() const
{
  cout << "--- Genebank ---\nId parentId treeDepth Count" << endl;

  map<int, Genotype*>::const_iterator it = m_genotypeMap.begin();
  for ( ; it!=m_genotypeMap.end(); it++ ){
    int parentId;
    Genotype *g = (*it).second;
    if ( g->getParent() == 0 )
      parentId = -1;
    else
      parentId = g->getParent()->getId();
    cout << g->getId() << " " << parentId << " " << g->getTreeDepth() << " " << g->getCount() << " ";
    if ( g->isCoalescent() )
      cout << "coalescent";
    cout << endl;
  }
  cout << endl;
}



