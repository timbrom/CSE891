#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <string>

using namespace std;

class Genotype {
private:
  Genotype *m_parent;
  const int m_parentId;
  const int m_treeDepth;
  const int m_id;
  const int m_birth;
  const int m_death;
  int m_count;
  bool m_coalescent;
  mutable bool m_tagged;
  string m_genes;
  
  Genotype();
  Genotype( const Genotype & g );
  const Genotype & operator=( const Genotype &g );
public:
  Genotype( Genotype* parent, int parentId, int treeDepth, int id,
	    int birth, int death, string genes ) :
    m_parent( parent ), m_parentId( parentId ), m_treeDepth( treeDepth ), 
    m_id( id ), m_birth( birth ), m_death( death ), m_count( 0 ),
    m_coalescent( false ), m_tagged( false ), m_genes ( genes ) {}

  void incrementCount() {
    m_count += 1; }

  bool decrementCount() {
    m_count -= 1;
    return m_count == 0;
  }

  Genotype* getParent() const {
    return m_parent; }

  void setParent( Genotype* parent ) {
    m_parent = parent; }
  
  void setCoalescent( bool c = true ) {
    m_coalescent = c; }

  void setTagged( bool c = true ) const {
    m_tagged = c; }

  int getParentId() const {
    return m_parentId; }
  
  int getTreeDepth() const {
    return m_treeDepth; }

  int getId() const {
    return m_id; }
  
  int getBirth() const {
    return m_birth; }
  
  int getDeath() const {
    return m_death; }
  
  char getGene( int pos ) const {
    return m_genes[pos]; }

  int getLength() const {
    return static_cast<int>( m_genes.size() ); }

  int getCount() const {
    return m_count; }

  bool isCoalescent() const {
    return m_coalescent; }

  bool isTagged() const {
    return m_tagged; }
};

#endif
