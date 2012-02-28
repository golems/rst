#include <Eigen/Core>
#include <vector>
#include <list>
#include "Link.h"
#include "World.h"
#include "RRT.h"
#include <iostream>
#include "Robot.h"

#ifndef _RST_PATH_PLANNER_
#define _RST_PATH_PLANNER_

template <class R = RRT>
class PathPlanner {
public:
	PathPlanner(); // You need to call one of the other constructors before the object is usable.
	PathPlanner(World& world, bool copyWorld = true, double stepsize = 0.1 );
	~PathPlanner();
	double stepsize;
	World* world;
	bool planPath(int robotId, const std::vector<int> &linkIds, const Eigen::VectorXd &start, const Eigen::VectorXd &goal, std::list<Eigen::VectorXd> &path, bool bidirectional = true, bool connect = true, unsigned int maxNodes = 0) const;
	bool planPath(int robotId, const std::vector<int> &linkIds, const std::vector<Eigen::VectorXd> &start, const std::vector<Eigen::VectorXd> &goal, std::list<Eigen::VectorXd> &path, bool bidirectional = true, bool connect = true, unsigned int maxNodes = 0) const;
	bool checkPathSegment(int robotId, const std::vector<int> &linkIds, const Eigen::VectorXd &config1, const Eigen::VectorXd &config2) const;
	static inline double randomInRange(double min, double max);
private:
	bool copyWorld;
	bool planSingleTreeRrt(int robot, const std::vector<int> &links, const std::vector<Eigen::VectorXd> &start, const Eigen::VectorXd &goal, std::list<Eigen::VectorXd> &path, bool connect, unsigned int maxNodes) const;
	bool planBidirectionalRrt(int robot, const std::vector<int> &links, const std::vector<Eigen::VectorXd> &start, const std::vector<Eigen::VectorXd> &goal, std::list<Eigen::VectorXd> &path, bool connect, unsigned int maxNodes) const;
};


template <class R>
PathPlanner<R>::PathPlanner() : copyWorld(false), world(NULL) {}

template <class R>
PathPlanner<R>::PathPlanner(World& world, bool copyWorld, double stepSize) {
	this->copyWorld = copyWorld;
	if(copyWorld) {
		this->world = new World(world);
	}
	else {
		this->world = &world;
	}
	this->stepsize = stepSize;
}

template <class R>
PathPlanner<R>::~PathPlanner() {
	if(this->copyWorld) {
		delete this->world;
	}
}

// true iff collision-free
// end points are not checked
template <class R>
bool PathPlanner<R>::checkPathSegment(int robotId, const vector<int> &linkIds, const Eigen::VectorXd &config1, const Eigen::VectorXd &config2) const {
	int n = (int)((config2 - config1).norm() / stepsize) + 1;
	for(int i = 1; i < n; i++) {
		Eigen::VectorXd conf = (double)(n - i)/(double)n * config1 + (double)i/(double)n * config2;
		world->robots[robotId]->setConf(linkIds, conf, true);
		if(world->checkCollisions()) {
			return false;
		}
	}
	return true;
}

template <class R>
inline double PathPlanner<R>::randomInRange(double min, double max) {
	return min + ((max-min) * ((double)rand() / ((double)RAND_MAX + 1)));
}

template <class R>
bool PathPlanner<R>::planPath(int robotId, const std::vector<int> &links, const Eigen::VectorXd &start, const Eigen::VectorXd &goal, std::list<Eigen::VectorXd> &path, bool bidirectional, bool connect, unsigned int maxNodes) const {
	world->robots[robotId]->setConf(links, start);
	if(world->checkCollisions())
		return false;
	world->robots[robotId]->setConf(links, goal);
	if(world->checkCollisions())
		return false;

	std::vector<Eigen::VectorXd> startVector, goalVector;
	startVector.push_back(start);
	goalVector.push_back(goal);
	return planPath(robotId, links, startVector, goalVector, path, bidirectional, connect, maxNodes);
}

template <class R>
bool PathPlanner<R>::planPath(int robotId, const std::vector<int> &links, const std::vector<Eigen::VectorXd> &start, const std::vector<Eigen::VectorXd> &goal, std::list<Eigen::VectorXd> &path, bool bidirectional, bool connect, unsigned int maxNodes) const {
	bool result;
	if(bidirectional)
		result = planBidirectionalRrt(robotId, links, start, goal, path, connect, maxNodes);
	else
		result = planSingleTreeRrt(robotId, links, start, goal.front(), path, connect, maxNodes);
	return result;
}

template <class R>
bool PathPlanner<R>::planSingleTreeRrt(int robot, const std::vector<int> &links, const std::vector<Eigen::VectorXd> &start, const Eigen::VectorXd &goal, std::list<Eigen::VectorXd> &path, bool connect, unsigned int maxNodes) const {

	R rrt(world, robot, links, start, stepsize);
	typename R::StepResult result = R::STEP_PROGRESS;
	double smallestGap = DBL_MAX;
	while (result != RRT::STEP_REACHED) {
		if(connect) {
			rrt.connect();
			if(rrt.connect(goal)) {
				result = RRT::STEP_REACHED;
			}
		}
		else {
			rrt.tryStep();
			result = rrt.tryStep(goal);
		}
		if(maxNodes > 0 && rrt.getSize() > maxNodes)
			return false;
		double gap = rrt.getGap(goal);
		if(gap < smallestGap) {
			smallestGap = gap;
			cout << "Gap: " << smallestGap << "    Tree size: " << rrt.configVector.size() << endl;
		}
	}

	rrt.tracePath(rrt.activeNode, path);
	return true;
}

template <class R>
bool PathPlanner<R>::planBidirectionalRrt(int robot, const std::vector<int> &links, const std::vector<Eigen::VectorXd> &start, const std::vector<Eigen::VectorXd> &goal, std::list<Eigen::VectorXd> &path, bool connect, unsigned int maxNodes) const {
	
	R start_rrt(world, robot, links, start, stepsize);
	R goal_rrt(world, robot, links, goal, stepsize);
	R* rrt1 = &start_rrt;
	R* rrt2 = &goal_rrt;
	
	double smallestGap = DBL_MAX;
	bool connected = false;
	while(!connected) {
		R* temp = rrt1;
		rrt1 = rrt2;
		rrt2 = temp;

		if(connect) {
			//rrt1->connect();
			rrt1->tryStep();
			connected = rrt2->connect(rrt1->configVector[rrt1->activeNode]);
		}
		else {
			if(rrt1->tryStep() != RRT::STEP_COLLISION) {
				connected = (RRT::STEP_REACHED == rrt2->tryStep(rrt1->configVector[rrt1->activeNode]));
			}
		}

		if(maxNodes > 0 && rrt1->getSize() + rrt2->getSize() > maxNodes) {
			cout << "Max number of nodes reached." << endl;
			return false;
		}

		double gap = rrt2->getGap(rrt1->configVector[rrt1->activeNode]);
		if(gap < smallestGap) {
			smallestGap = gap;
			cout << "Gap: " << smallestGap << "  Sizes: " << start_rrt.configVector.size() << "/" << goal_rrt.configVector.size() << endl;
			cout << "  Coll.: " << start_rrt.numCollisions + goal_rrt.numCollisions
				<< "  Progress: " << start_rrt.numNoProgress + goal_rrt.numNoProgress
				<< "  Step: " << start_rrt.numStepTooLarge + goal_rrt.numStepTooLarge
				<< "  Task error: " << start_rrt.numErrorIncrease + goal_rrt.numErrorIncrease
				<< endl;
		}
	}
	
	start_rrt.tracePath(start_rrt.activeNode, path);
	goal_rrt.tracePath(goal_rrt.activeNode, path, true);
	
	return true;
}



#endif /** _RST_PATH_PLANNER_ */
