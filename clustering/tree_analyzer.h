#ifndef TREE_ANALYZER_H
#define TREE_ANALYZER_H

#include <vector>
#include "genebank.h"

class TreeAnalyzer {
private:
  Genebank m_genebank;
  
  vector<int> m_finalPop;
  vector<vector<int> > m_distanceMatrix;
  vector<int> m_picked; // new
  int m_finalPopSize;
  int m_maxTreeDepth;
  double m_aveDistance;
  int m_maxDistance;
  
public:
  TreeAnalyzer();
  ~TreeAnalyzer();
  
  //  void loadData();
  void loadData( const char *gzDetailFile, const char *gzHistoricFile);
  void calculateDistanceMatrix();
  void calculateHammingDistanceMatrix();
  void calculateMRCADistanceMatrix();
  void doClusteringAnalysis( const char *clusterDataFile, int cutoff = 1 );
  void sortData( const char *sortedDataFile );
};

#endif
