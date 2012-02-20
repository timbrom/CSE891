#include "tree_analyzer.h"

#include "genebank.h"
#include "genotype.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string> 

TreeAnalyzer::TreeAnalyzer() 
  : m_maxTreeDepth( 0 ), m_aveDistance( 0 ), m_maxDistance( 0 )
{
}


TreeAnalyzer::~TreeAnalyzer()
{
}

void TreeAnalyzer::loadData( const char *gzDetail, const char *gzHistoric )
{
  cout << "loading data" << endl;

  const string gzExt = ".gz";
  string gzDetailFile;
  string detailFile;
  string gzHistoricFile;
  string historicFile;
  
  const string gunzipCommand = "gunzip ";
  const string gzipCommand = "gzip ";
  string systemCommand;

  gzDetailFile = gzDetail;
  gzHistoricFile = gzHistoric;
  detailFile = gzDetailFile;
  detailFile.resize(gzDetailFile.size() - gzExt.size());
  historicFile = gzHistoricFile;
  historicFile.resize(gzHistoricFile.size() - gzExt.size());

  systemCommand = gunzipCommand + gzDetailFile;
  system( systemCommand.c_str() );
  
  systemCommand = gunzipCommand + gzHistoricFile;
  system( systemCommand.c_str() );
  
  ifstream final( detailFile.c_str() );
  if ( final.fail() ){
    cerr << "error opening " << detailFile << ". Exiting" << endl;
    exit( -1 );
  }
  
  int id, parentId, treeDepth, birth, death;
  double dummy;
  string genome; // new
  m_finalPopSize = 0;

  // read in the data we are interested in
  while ( final >> id >> parentId >> dummy >> dummy >> dummy >> dummy
	  >> dummy >> dummy >> dummy >> birth >> death >> treeDepth >>
	  genome ) {     // note: death = -1 for alive organisms

    if ( treeDepth > m_maxTreeDepth )
      m_maxTreeDepth = treeDepth;
    
    m_genebank.createGenotype( id, parentId, treeDepth, birth, death,
			       genome );
    m_finalPop.push_back( id );
    m_finalPopSize += 1;
  }
  final.close();

  cout << "final pop files opened ok" << endl;

  ifstream historic( historicFile.c_str() );
  if ( historic.fail() ) {
    cerr << "error opening " << historicFile << ". Exiting" << endl;
    exit( -1 );
  }
  
  while ( historic >> id >> parentId >> dummy >> dummy >> dummy >>
	  dummy >> dummy >> dummy >> dummy >> birth >> death >>
	  treeDepth >> genome ) {

    m_genebank.createGenotype( id, parentId, treeDepth, birth, death,
			       genome );
  }
  historic.close();
    
  cout << "historic files opened ok" << endl;

  m_genebank.setupParentPointers();
  //  Genotype *g = m_genebank.getGenotype( m_finalPop[0] );
  m_genebank.checkCoalescence( m_genebank.getGenotype( m_finalPop[0] ) );
  //  m_genebank.print();
  
  systemCommand = gzipCommand + detailFile;
  system( systemCommand.c_str() );
  
  systemCommand = gzipCommand + historicFile;
  system( systemCommand.c_str() );  
}


void TreeAnalyzer::calculateDistanceMatrix()
{
  cout << "calculating distance matrix" << endl;
  
  // reserve enough space
  m_distanceMatrix.resize( m_finalPopSize );
  for ( int i=0; i<m_finalPopSize; i++ )
    m_distanceMatrix[i].resize( m_finalPopSize );
  
  for ( int i=0; i<m_finalPopSize; i++ ){
    //    cout <<  "[" << i+1 << "/" << m_finalPopSize << "] ";
    //    cout.flush();
    for ( int j=i; j<m_finalPopSize; j++ ){
	int ID1 = m_finalPop[i];
	int ID2 = m_finalPop[j];
	int dist = m_genebank.calcTreeDist( m_genebank.getGenotype( ID1 ),
					    m_genebank.getGenotype( ID2 ) );
	m_distanceMatrix[i][j] = dist;
	m_distanceMatrix[j][i] = dist;
	if ( i==j )
	  m_aveDistance += dist;
	else
	  m_aveDistance += dist*2;
	
	if ( dist > m_maxDistance )
	  m_maxDistance = dist;
    }
  }

  m_aveDistance /= static_cast<double>(m_finalPopSize*m_finalPopSize);
}


void TreeAnalyzer::calculateHammingDistanceMatrix()
{
  cout << "calculating Hamming distance matrix" << endl;
  
  // reserve enough space
  m_distanceMatrix.resize( m_finalPopSize );
  for ( int i=0; i<m_finalPopSize; i++ )
    m_distanceMatrix[i].resize( m_finalPopSize );
  
  for ( int i=0; i<m_finalPopSize; i++ ){
    for ( int j=i; j<m_finalPopSize; j++ ){
	int ID1 = m_finalPop[i];
	int ID2 = m_finalPop[j];
	int dist = m_genebank.calcHammingDist( m_genebank.getGenotype( ID1 ),
					       m_genebank.getGenotype( ID2 ) );
	m_distanceMatrix[i][j] = dist;
	m_distanceMatrix[j][i] = dist;
	if ( i==j ) {
	  assert ( dist == 0 );
	  m_aveDistance += dist;
	}
	else
	  m_aveDistance += dist*2;
	
	if ( dist > m_maxDistance )
	  m_maxDistance = dist;
    }
  }

  m_aveDistance /= static_cast<double>(m_finalPopSize*m_finalPopSize);
}

void TreeAnalyzer::calculateMRCADistanceMatrix()
{
  cout << "calculating depth to MRCA distance matrix" << endl;
  
  // reserve enough space
  m_distanceMatrix.resize( m_finalPopSize );
  for ( int i=0; i<m_finalPopSize; i++ )
    m_distanceMatrix[i].resize( m_finalPopSize );
  
  for ( int i=0; i<m_finalPopSize; i++ ){
    for ( int j=i; j<m_finalPopSize; j++ ){
	int ID1 = m_finalPop[i];
	int ID2 = m_finalPop[j];
	int dist = m_genebank.calcMRCADist( m_genebank.getGenotype( ID1 ),
					    m_genebank.getGenotype( ID2 ) );
	assert(dist >= 0);
	m_distanceMatrix[i][j] = dist;
	m_distanceMatrix[j][i] = dist;
	if ( i==j ) {
	  assert ( dist == 0 );
	  m_aveDistance += dist;
	}
	else
	  m_aveDistance += dist*2;
	
	if ( dist > m_maxDistance )
	  m_maxDistance = dist;
    }
  }

  m_aveDistance /= static_cast<double>(m_finalPopSize*m_finalPopSize);
}

void TreeAnalyzer::doClusteringAnalysis( const char *clusterDataFile, int cutoff )
{
  cout << "doing clustering analysis" << endl;
  
  vector<int> nearestDistToPicked;
  vector<int> nonPickedGenotypes;
  
  nearestDistToPicked.reserve( m_finalPopSize );
  nonPickedGenotypes.reserve( m_finalPopSize );

  // initialize values
  for ( int i = 0; i < m_finalPopSize; i++ ){
    nearestDistToPicked[i] = 2*m_maxDistance; // was 2*m_maxTreeDepth (why?)
    nonPickedGenotypes[i] = i;
  }	
  
  int maxNonPicked = m_finalPopSize - 1;

  // output file

  ofstream out( clusterDataFile );
  if ( out.fail() ) {
    cerr << "cannot open cluster data file " << clusterDataFile << endl;
    exit( -1 );
  }
  out << "#Maximum distance between organisms: " << m_maxDistance << endl;
  out << "#Average distance between organisms: " << m_aveDistance << endl;
  out << "#Cutoff pick value: " << cutoff << endl;
  out <<  "#<organism ID> <pick value>\n";

  // now find best to pick
  int pick;
  int pickValue;
  int pickI;
  int pickOrg = 0;
  for ( pick = 0; pick < m_finalPopSize; pick++ ){ // for each possible pick
    pickValue = 0;
    pickI = 0;
    //    cout << "[" << pick+1 << "/" << m_finalPopSize << "] ";
    for ( int i = 0; i <= maxNonPicked; i++ ){ // for each organism

      // calculate drop in distance ( = value )
      int value = 0;
      for ( int j = 0; j <= maxNonPicked; j++ ){
	if ( nearestDistToPicked[nonPickedGenotypes[j]]
	     > m_distanceMatrix[nonPickedGenotypes[i]][nonPickedGenotypes[j]] ){
	  value += nearestDistToPicked[nonPickedGenotypes[j]]
		   - m_distanceMatrix[nonPickedGenotypes[i]][nonPickedGenotypes[j]];
	}
      }
      if ( value > pickValue ){
	pickValue = value;
	pickI = i;
	pickOrg = nonPickedGenotypes[i];
      }
    }	

    // now actually adjust distances
    for ( int j = 0; j <= maxNonPicked; j++ ){
      if ( nearestDistToPicked[nonPickedGenotypes[j]]
	   > m_distanceMatrix[pickOrg][nonPickedGenotypes[j]] ){
	nearestDistToPicked[nonPickedGenotypes[j]]
	  = m_distanceMatrix[pickOrg][nonPickedGenotypes[j]];
      }
    }

    // and remove picked
    nonPickedGenotypes[pickI] = nonPickedGenotypes[maxNonPicked];
    maxNonPicked -= 1;

    // cout << "Picked " << m_finalPop[pickOrg] << " at " <<  pickValue << "\n";
    out << m_finalPop[pickOrg] << " " <<  pickValue << "\n";

    m_picked.push_back( pickI ); // new

    // 2002/11/22: adjusted so that m_picked holds at least one value
    // before: m_picked holds only those _above_ cutoff for clarity.
    
    if ( m_picked.size() > 1 && pickValue < cutoff ) {
      m_picked.pop_back();
      break;
    }

    // we don't have to pick all organisms, just a couple
    if ( pick>=99 )
      break;
  }

  out.close();
}

void TreeAnalyzer::sortData( const char *sortedDataFile )
{
  cout << "sorting data\n";

  ofstream out( sortedDataFile, ios::app ); // reopen data file and append data
  if ( out.fail() ) {
    cerr << "cannot reopen cluster data file " << sortedDataFile << endl;
    exit( -1 );
  }

  // correction - 2002/10/10
  // there is always at least one cluster, even if its score is less
  // than the cutoff ( so m_picked.size() is 0 )
  if (m_picked.size() < 1) {
    cout << "Error! #Number of Clusters: " << m_picked.size() << "\n";
    out << "Error! #Number of Clusters: " << m_picked.size() << "\n";
  }
  else
    out << "#Number of Clusters: " << m_picked.size() << "\n";
    
  int minDist;
  int closest;
  int cluster[m_finalPopSize][2];
  if ( m_picked.size() > 1 ) { // if more than one cluster, then sort
    for ( int i = 0; i <= m_finalPopSize; i++ ) {
      cluster[i][0] = i;
      cluster[i][1] = -1;
    }
  }
  out << "#<genotype> <closest picked genotype> <distance>\n";
  for ( int i = 0; i < m_finalPopSize; i++ ) {
    //    minDist = 2 * m_maxTreeDepth;
    minDist = 2 * m_maxDistance;
    closest = -1;
    vector<int>::iterator pickedIter;
    for ( pickedIter = m_picked.begin(); pickedIter != m_picked.end();
	  pickedIter++ ) {
      if ( m_distanceMatrix[i][*pickedIter] < minDist ) {
	minDist = m_distanceMatrix[i][*pickedIter];
	closest = *pickedIter;
      }
    }
    assert ( closest != -1 );
    cluster[i][1] = closest;
    out << m_finalPop[i] << " " << m_finalPop[closest] << " " << minDist << endl;
  }

  out.close();
}










