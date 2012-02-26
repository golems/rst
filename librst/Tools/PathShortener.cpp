#include "PathShortener.h"
#include "Robot.h"
#include "RRT.h"

using namespace std;
using namespace Eigen;

#define RAND12(N1,N2) N1 + ((N2-N1) * ((double)rand() / ((double)RAND_MAX + 1))) // random # between N&M

PathShortener::PathShortener() {}

PathShortener::PathShortener(World* world, int robotId, std::vector<int> linkIds) :
   world(world),
   robotId(robotId),
   linkIds(linkIds)
{}

PathShortener::~PathShortener()
{}

void PathShortener::shortenPath(list<VectorXd> &path, double stepSize)
{
	printf("--> Start Brute Force Shortener \n");  
	srand(time(NULL));

	this->stepSize = stepSize;

	const int numShortcuts = path.size() * 100;
	
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
	return !checkSegment(intermediatePoints, *it1, *it2);
}

// true iff collision
// interemdiatePoints are only touched if collision-free
bool PathShortener::checkSegment(list<VectorXd> &intermediatePoints, const VectorXd &config1, const VectorXd &config2) {
	list<VectorXd> tempIntermediatePoints;
	const int n = (int)((config1 - config2).norm() / stepSize) + 1; // number of intermediate segments
	for(int i = 1; i < n; i++) {
		const VectorXd config = (double)(n - i) / (double)n * config1 + (double)i / (double)n * config2;
		world->robots[robotId]->setConf(linkIds, config);
		if(world->checkCollisions()) {
			return true;
		}
		tempIntermediatePoints.push_back(config);
	}
	intermediatePoints.clear();
	intermediatePoints.splice(intermediatePoints.begin(), tempIntermediatePoints);
	return false;
}