#include "PathShortener.h"
#include "Robot.h"
#include "RRT.h"

using namespace std;
using namespace Eigen;

extern int numCollisionChecks;

#define RAND12(N1,N2) N1 + ((N2-N1) * ((double)rand() / ((double)RAND_MAX + 1))) // random # between N&M

PathShortener::PathShortener() {}

PathShortener::PathShortener(World* world, int robotId, const std::vector<int> &linkIds, double stepSize) :
   world(world),
   robotId(robotId),
   linkIds(linkIds),
   stepSize(stepSize)
{}

PathShortener::~PathShortener()
{}

void PathShortener::shortenPath(list<VectorXd> &path)
{
	printf("--> Start Brute Force Shortener \n"); 
	srand(time(NULL));

	const int numShortcuts = path.size() * 20;
	
	// Number of checks
	for( int count = 0; count < numShortcuts; count++ ) {
		if( path.size() < 3 ) { //-- No way we can reduce something leaving out the extremes
			return;
		}
		
		int node1Index;
		int node2Index;
		do {
			node1Index = (int) RAND12(0, path.size());
			node2Index = (int) RAND12(0, path.size());
		} while(node2Index <= node1Index + 1);
		
		list<VectorXd>::iterator node1Iter = path.begin();
		advance(node1Iter, node1Index);
		list<VectorXd>::iterator node2Iter = node1Iter;
		advance(node2Iter, node2Index - node1Index);

		list<VectorXd> intermediatePoints;
		if(localPlanner(intermediatePoints, node1Iter, node2Iter)) {
			list<VectorXd>::iterator node1NextIter = node1Iter;
			node1NextIter++;
			path.erase(node1NextIter, node2Iter);
			path.splice(node2Iter, intermediatePoints);
        }
	}
	printf("End Brute Force Shortener \n");
}

bool PathShortener::localPlanner(list<VectorXd> &intermediatePoints, list<VectorXd>::const_iterator it1, list<VectorXd>::const_iterator it2) {
	return segmentCollisionFree(intermediatePoints, *it1, *it2);
}

// true iff collision-free
// does not check endpoints
// interemdiatePoints are only touched if collision-free
bool PathShortener::segmentCollisionFree(list<VectorXd> &intermediatePoints, const VectorXd &config1, const VectorXd &config2) {
	const double length = (config1 - config2).norm();
	if(length <= stepSize) {
		return true;
	}

	const int n = (int)(length / stepSize) + 1; // number of intermediate segments
	int n1 = n / 2;
	int n2 = n / 2;
	if(n % 2 == 1) {
		n2 += 1;
	}

	VectorXd midpoint = (double)n2 / (double)n * config1 + (double)n1 / (double)n * config2;
	list<VectorXd> intermediatePoints1, intermediatePoints2;
	numCollisionChecks++;
	world->robots[robotId]->setConf(linkIds, midpoint);
	if(!world->checkCollisions() && segmentCollisionFree(intermediatePoints1, config1, midpoint)
			&& segmentCollisionFree(intermediatePoints2, midpoint, config2))
	{
		intermediatePoints.clear();
		intermediatePoints.splice(intermediatePoints.end(), intermediatePoints1);
		intermediatePoints.push_back(midpoint);
		intermediatePoints.splice(intermediatePoints.end(), intermediatePoints2);
		return true;
	}
	else {
		return false;
	}
}
