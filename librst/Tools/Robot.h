#ifndef ROBOT_H
#define ROBOT_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <queue>
#include <list>
#include <string>

using namespace std;

class World;

#define config vec<double>

#include "Link.h"

class Robot {
public:
	Eigen::Vector3d COM;

	Robot();
	Robot(Robot &);
	~Robot();

	vector<Link*> links;
	vector<Link*> activeLinks;
	Link* baseLink;
	World* world;

	int idNum;
	double mass;

	string pathname;
	string name;

	int findLink(string name);
	int Load(string, World*);
	void setConf(const vector<double> &conf);
	void getConf(vector<double> &conf);

	void updateCOM();
	void drawCOM();
	void Draw();
	void Info();
};

#endif
