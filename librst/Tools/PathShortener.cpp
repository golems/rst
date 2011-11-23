/**
 * @file Smoothie.cpp
 * @brief Implementation of Smoothie class for Shortening-Smoothing
 */
#include "PathShortener.h"

PathShortener::PathShortener()
{}

PathShortener::~PathShortener()
{}

/**
 * @function Initialize
 */
void PathShortener::Initialize(  World &world, 
                            int robotId,
                            std::vector<int> linksId )
{
   mWorld = &world;
   mRobotId = robotId;
   mLinksId = linksId;
}


/**
 * @function BruteForce
 * @brief Random smoother preserving points separated by stepSize (avoids the ugly cornered profiles)
 */
void PathShortener::BruteForce( std::list< Eigen::VectorXd > &path, double stepSize )
{
   srand(time(NULL));

   mStepSize = stepSize;

   int node1Index; int node2Index; int nodeAuxIndex;
   int nodeMinIndex; int nodeMaxIndex;
   int numPoints;
   int numChecks;
   list<Eigen::VectorXd>::iterator node1Iter;
   list<Eigen::VectorXd>::iterator node2Iter;

   numPoints = path.size();
   numChecks = numPoints * 100;

   // Number of checks
   for( int count = 0; count < numChecks; count++ )
   {
      if( path.size() < 3 ) {    mBruteForcePath = path;  return; } //-- No way we can reduce something leaving out the extremes

      nodeMinIndex = 0;
      nodeMaxIndex = path.size() - 1;

      do
      { 
        node1Index = (int) RAND12( nodeMinIndex, nodeMaxIndex );
        node2Index = (int) RAND12( nodeMinIndex, nodeMaxIndex ); 
      } while( node2Index == node1Index || abs(node1Index - node2Index) < 2 );

	  if(node1Index == nodeMaxIndex) {
		  cout << "Hurra\n";
	  }

      if( node2Index < node1Index ) 
      {  nodeAuxIndex = node1Index;
         node1Index = node2Index;
         node2Index = nodeAuxIndex; }
      
      //-- Check
      node1Iter = path.begin(); node2Iter = path.begin();
      advance( node1Iter, node1Index );
      advance( node2Iter, node2Index );
      Eigen::VectorXd node1 = *node1Iter;
      Eigen::VectorXd node2 = *node2Iter;

      bool result = CheckPathSegment( *node1Iter, *node2Iter );
      //printf("[%d] Start path size: %d \n", count, path.size() );
      if( result == true )
      { 
        int oldMidNodes = node2Index - node1Index - 1;
        int newMidNodes = (int)( (*node2Iter - *node1Iter).norm() / mStepSize - 1 );
        int n = newMidNodes + 1;
        //printf("-- Old nodes: %d new mid nodes: %d \n", oldMidNodes, newMidNodes );

        if( newMidNodes <= oldMidNodes )
        {
           int times = node2Index- node1Index - 1;
           for( int j = 0; j < times; j++ )
           {  list<Eigen::VectorXd>::iterator temp = path.begin(); 
              advance( temp, node1Index + 1 );
              path.erase( temp );  
           }
           
      
           for( int j = 1; j <= newMidNodes; j++ )
           {  list<Eigen::VectorXd>::iterator temp = path.begin(); 
              advance( temp, node1Index + j );
              Eigen::VectorXd conf = (double)(n - j)/(double)n * node1 + (double)(j)/(double)n * node2;
              path.insert( temp, conf );
           }
           
        }
        else
        {  if( newMidNodes == oldMidNodes )
           { printf("Okay old and new are equal, I can live with it \n"); }
           else
           { printf("--(!) What the heck? This is weird -- Check the BruteForce Shortener! \n");} 
        }
      //printf("[%d] End path size: %d \n", count, path.size() );
      } //-- end if result true

   }   

   mBruteForcePath = path;   

}


/**
 * @function CheckPathSegment
 */
bool PathShortener::CheckPathSegment( Eigen::VectorXd config1, Eigen::VectorXd config2 ) const {
	int n = (int)((config2 - config1).norm() / mStepSize);
	for(int i = 0; i < n; i++) {
		Eigen::VectorXd conf = (double)(n - i)/(double)n * config1 + (double)(i)/(double)n * config2;
		mWorld->robots[mRobotId]->setConf( mLinksId, conf, true);
		if( mWorld->checkCollisions()) {
			return false;
		}
	}
	return true;
}


