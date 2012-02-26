#pragma once

#include <list>
#include <vector>
#include <Eigen/Core>
#include "World.h"


class PathShortener
{
public:
	PathShortener();
	PathShortener(World* world, int robotId, std::vector<int> linksId);
	~PathShortener();
	virtual void shortenPath(std::list<Eigen::VectorXd> &rawPath, double stepSize = 0.1);
	bool checkSegment(std::list<Eigen::VectorXd> &waypoints, const Eigen::VectorXd &config1, const Eigen::VectorXd &config2);
protected:
	World* world;
	int robotId;
	std::vector<int> linkIds;
	double stepSize;
	virtual bool localPlanner(std::list<Eigen::VectorXd> &waypoints, std::list<Eigen::VectorXd>::const_iterator it1, std::list<Eigen::VectorXd>::const_iterator it2);
};