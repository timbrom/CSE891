#include "tree_analyzer.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>

int main( int argc, char **argv )
{
  if ( argc < 4 ) {
    cout << "wrong number of arguments --\n";
    cout << "format: treeCS <gzipped detail file> <gzipped historic file> "; 
    cout << "<clustering data output file> ";
    cout << "[<cluster cutoff value>]\n";
    return 0;
  }
  TreeAnalyzer t;


  ifstream test( argv[3] );
  if ( test.good() ){
    cout << "cluster data file " << argv[3] << " exists already. ";
    cout << "Skipping directory..." << endl;

    return 0;
  }
  
  test.close();
  
  for (int i = 0; i < argc; i++)
    cout << argv[i] << " ";
  cout << endl;

  t.loadData( argv[1], argv[2] );
  t.calculateDistanceMatrix();
  // t.calculateMRCADistanceMatrix();
  // t.calculateHammingDistanceMatrix();
  if ( argc == 4 )
    t.doClusteringAnalysis( argv[3] );
  else
    t.doClusteringAnalysis( argv[3], atoi( argv[4] ) );

  t.sortData( argv[3] );
  cout << "done!\n\n";

  return 0;
}

