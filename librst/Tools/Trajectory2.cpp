#include "Trajectory2.h"
#include <limits>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;

const double Trajectory2::timeStep = 0.001;
const double Trajectory2::eps = 0.000001;

Trajectory2::Trajectory2(const Path &path, const VectorXd &maxVelocity, const VectorXd &maxAcceleration) :
	path(path),
	maxVelocity(maxVelocity),
	maxAcceleration(maxAcceleration),
	n(maxVelocity.size()),
	valid(true),
	cachedTime(numeric_limits<double>::max())
{
	// debug
	{
	ofstream file("maxVelocity.txt");
	for(double s = 0.0; s < path.getLength(); s += 0.0001) {
		file << s << "  " << getMaxPathVelocity(s) << endl;
	}
	file.close();
	}

	//{
	//ofstream dfile("maxVelocityDebug.txt");
	//for(double s = 0.0; s < path.getLength(); s += 0.001) {
	//	dfile << s << " ";
	//	VectorXd maxPathVelocity = getMaxPathVelocities(s);
	//	for(int i = 0; i < maxPathVelocity.size(); i++) {
	//		dfile << "  " << maxPathVelocity[i];
	//	}
	//	dfile << endl;
	//}
	//dfile.close();
	//}


	list<TrajectoryStep> startTrajectory;
	startTrajectory.push_back(TrajectoryStep(0.0, 0.0));
	double afterAcceleration = getMinMaxPathAcceleration(0.0, 0.0, true);
	while(!integrateForward(startTrajectory, afterAcceleration)) {
		double beforeAcceleration;
		list<TrajectoryStep> trajectory = getNextSwitchingPoint(startTrajectory.back().pathPos, beforeAcceleration, afterAcceleration);
		integrateBackward(trajectory, startTrajectory, beforeAcceleration);
	}

	list<TrajectoryStep> endTrajectory;
	endTrajectory.push_front(TrajectoryStep(path.getLength(), 0.0));
	double beforeAcceleration = getMinMaxPathAcceleration(path.getLength(), 0.0, false);
	integrateBackward(endTrajectory, startTrajectory, beforeAcceleration);
	
	this->trajectory = startTrajectory;

	// calculate timing
	list<TrajectoryStep>::iterator previous = trajectory.begin();
	list<TrajectoryStep>::iterator it = previous;
	it->time = 0.0;
	it++;
	while(it != trajectory.end()) {
		it->time = previous->time + (it->pathPos - previous->pathPos) / ((it->pathVel + previous->pathVel) / 2.0);
		previous = it;
		it++;
	}

	// debug
	ofstream file("trajectory.txt");
	for(list<TrajectoryStep>::iterator it = trajectory.begin(); it != trajectory.end(); it++) {
		file << it->pathPos << "  " << it->pathVel << endl;
	}
	file.close();
}

Trajectory2::~Trajectory2(void) {
}	

list<Trajectory2::TrajectoryStep> Trajectory2::getNextSwitchingPoint(double pathPos, double &beforeAcceleration, double &afterAcceleration) {
	const double eps = 0.000001;

	double switchingPathPos = pathPos - eps;

	double switchingPathVel;
	double afterPathVel;
	double afterPathPos;
	double beforePathVel;
	double beforePathPos;
	do {
		bool discontinuity;
		switchingPathPos = path.getNextSwitchingPoint(switchingPathPos + eps, discontinuity);
		
		if(discontinuity) {
			switchingPathVel = min(getMaxPathVelocity(switchingPathPos - eps), getMaxPathVelocity(switchingPathPos + eps));
			beforeAcceleration = getMinMaxPathAcceleration(switchingPathPos - eps, switchingPathVel, false);
			afterAcceleration = getMinMaxPathAcceleration(switchingPathPos + eps, switchingPathVel, true);
		}
		else {
			switchingPathVel = getMaxPathVelocity(switchingPathPos);
			beforeAcceleration = 0.0;
			afterAcceleration = 0.0;
		}

		beforePathVel = switchingPathVel - timeStep * beforeAcceleration;
		beforePathPos = switchingPathPos - timeStep * 0.5 * (switchingPathVel + beforePathVel);
		
		afterPathVel = switchingPathVel + timeStep * afterAcceleration;
		afterPathPos = switchingPathPos + timeStep * 0.5 * (switchingPathVel + afterPathVel);
	} while(beforePathVel > getMaxPathVelocity(beforePathPos) || afterPathVel > getMaxPathVelocity(afterPathPos));
	
	list<Trajectory2::TrajectoryStep> trajectory;
	trajectory.push_back(TrajectoryStep(switchingPathPos, switchingPathVel));
	return trajectory;
}

bool Trajectory2::integrateForward(list<TrajectoryStep> &trajectory, double acceleration) {
	
	double pathPos = trajectory.back().pathPos;
	double pathVel = trajectory.back().pathVel;
	
	while(true)
	{
		double oldPathVel = pathVel;
		pathVel += timeStep * acceleration;
		pathPos += timeStep * 0.5 * (oldPathVel + pathVel);
		trajectory.push_back(TrajectoryStep(pathPos, pathVel));
		acceleration = getMinMaxPathAcceleration(pathPos, pathVel, true);

		if(pathPos > path.getLength()) {
			return true;
		}
		else if(pathVel < 0.0) {
			ofstream file("trajectory.txt");
			for(list<TrajectoryStep>::iterator it = trajectory.begin(); it != trajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			file.close();
			valid = false;
			cout << "error " << path.getLength() << endl;
			return true;
		}
		else if(pathVel >= getMaxPathVelocity(pathPos)) {
			// find more accurate intersection with max-velocity curve using bisection
			TrajectoryStep overshoot = trajectory.back();
			trajectory.pop_back();
			double slope = getSlope(trajectory.back(), overshoot);
			double before = trajectory.back().pathPos;
			double after = overshoot.pathPos;
			while(after - before > 0.00001) {
				const double midpoint = 0.5 * (before + after);
				const double midpointPathVel = trajectory.back().pathVel + slope * (midpoint - trajectory.back().pathPos);
				if(midpointPathVel > getMaxPathVelocity(midpoint))
					after = midpoint;
				else
					before = midpoint;
			}
			trajectory.push_back(TrajectoryStep(before, trajectory.back().pathVel + slope * (before - trajectory.back().pathPos)));
			return false;
		}
	}
}


void Trajectory2::integrateBackward(list<TrajectoryStep> &trajectory, list<TrajectoryStep> &startTrajectory, double acceleration) {
	list<TrajectoryStep>::reverse_iterator before = startTrajectory.rbegin();
	double pathPos = trajectory.front().pathPos;
	double pathVel = trajectory.front().pathVel;

	while(true)
	{
		double oldPathVel = pathVel;
		pathVel -= timeStep * acceleration;
		pathPos -= timeStep * 0.5 * (oldPathVel + pathVel);
		trajectory.push_front(TrajectoryStep(pathPos, pathVel));
		acceleration = getMinMaxPathAcceleration(pathPos, pathVel, false);

		if(pathVel < 0.0 || pathPos < 0.0) {
			ofstream file("trajectory.txt");
			for(list<TrajectoryStep>::iterator it = trajectory.begin(); it != trajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			for(list<TrajectoryStep>::iterator it = startTrajectory.begin(); it != startTrajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			file.close();
			valid = false;
			cout << "error " << pathPos << " " << pathVel << endl;
			return;
		}

		while(before != startTrajectory.rend() && before->pathPos > pathPos) {
			before++;
		}

		if(before != startTrajectory.rbegin() && pathVel >= before->pathVel + getSlope(before.base()) * (pathPos - before->pathPos)) {
			TrajectoryStep overshoot = trajectory.front();
			trajectory.pop_front();
			list<TrajectoryStep>::iterator after = before.base();
			TrajectoryStep intersection = getIntersection(startTrajectory, after, overshoot, trajectory.front());

			if(after != startTrajectory.end()) {
				startTrajectory.erase(after, startTrajectory.end());
				startTrajectory.push_back(intersection);
			}
			startTrajectory.splice(startTrajectory.end(), trajectory);

			return;
		}
		else if(pathVel >= getMaxPathVelocity(pathPos)) {
			ofstream file("trajectory.txt");
			for(list<TrajectoryStep>::iterator it = trajectory.begin(); it != trajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			for(list<TrajectoryStep>::iterator it = startTrajectory.begin(); it != startTrajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			file.close();
			cout << "error" << endl;
			valid = false;
			return;
		}
	}
}

inline double Trajectory2::getSlope(const TrajectoryStep &point1, const TrajectoryStep &point2) {
	return (point2.pathVel - point1.pathVel) / (point2.pathPos - point1.pathPos);
}

inline double Trajectory2::getSlope(list<TrajectoryStep>::const_iterator lineEnd) {
	list<TrajectoryStep>::const_iterator lineStart = lineEnd;
	lineStart--;
	return getSlope(*lineStart, *lineEnd);
}

Trajectory2::TrajectoryStep Trajectory2::getIntersection(const list<TrajectoryStep> &trajectory, list<TrajectoryStep>::iterator &it, const TrajectoryStep &linePoint1, const TrajectoryStep &linePoint2) {
	
	const double lineSlope = getSlope(linePoint1, linePoint2);
	it--;

	double factor = 1.0;
	if(it->pathVel > linePoint1.pathVel + lineSlope * (it->pathPos - linePoint1.pathPos))
		factor = -1.0;
	it++;
	
	while(it != trajectory.end() && factor * it->pathVel < factor * (linePoint1.pathVel + lineSlope * (it->pathPos - linePoint1.pathPos))) {
		it++;
	}

	if(it == trajectory.end()) {
		return TrajectoryStep(0.0, 0.0);
	}
	else {
		const double trajectorySlope = getSlope(it);
		const double intersectionPathPos = (it->pathVel - linePoint1.pathVel + lineSlope * linePoint1.pathPos - trajectorySlope * it->pathPos)
			/ (lineSlope - trajectorySlope);
		const double intersectionPathVel = linePoint1.pathVel + lineSlope * (intersectionPathPos - linePoint1.pathPos);
		return TrajectoryStep(intersectionPathPos, intersectionPathVel);
	}
}


double Trajectory2::getMinMaxPathAcceleration(double pathPos, double pathVel, bool max) {
	VectorXd configDeriv = path.getConfigDeriv(pathPos);
	VectorXd configDeriv2 = path.getConfigDeriv2(pathPos);
	double factor = max ? 1.0 : -1.0;
	double maxPathAcceleration = numeric_limits<double>::max();
	for(unsigned int i = 0; i < n; i++) {
		if(configDeriv[i] != 0.0) {
			maxPathAcceleration = min(maxPathAcceleration,
				maxAcceleration[i]/abs(configDeriv[i]) - factor * configDeriv2[i] * pathVel*pathVel / configDeriv[i]);
		}
	}
	return factor * maxPathAcceleration;
}

double Trajectory2::getMaxPathVelocity(double pathPos) {
	double maxPathVelocity = numeric_limits<double>::max();
	const VectorXd configDeriv = path.getConfigDeriv(pathPos);
	const VectorXd configDeriv2 = path.getConfigDeriv2(pathPos);
	for(unsigned int i = 0; i < n; i++) {
		if(configDeriv[i] != 0.0) {
			for(unsigned int j = i + 1; j < n; j++) {
				if(configDeriv[j] != 0.0) {
					double A_ij = configDeriv2[i] / configDeriv[i] - configDeriv2[j] / configDeriv[j];
					if(A_ij != 0.0) {
						maxPathVelocity = min(maxPathVelocity,
							sqrt((maxAcceleration[i] / abs(configDeriv[i]) + maxAcceleration[j] / abs(configDeriv[j]))
							/ abs(A_ij)));
					}
				}
			}
		}
		else if(configDeriv2[i] != 0.0) {
			maxPathVelocity = min(maxPathVelocity, sqrt(maxAcceleration[i] / abs(configDeriv2[i])));
		}
		maxPathVelocity = min(maxPathVelocity, maxVelocity[i] / abs(configDeriv[i]));
	}
	return maxPathVelocity;
}


// debug
VectorXd Trajectory2::getMaxPathVelocities(double pathPos) {
	VectorXd maxPathVelocity(n * (n-1) / 2 + n);
	const VectorXd configDeriv = path.getConfigDeriv(pathPos);
	const VectorXd configDeriv2 = path.getConfigDeriv2(pathPos);
	int k = 0;
	for(unsigned int i = 0; i < n; i++) {
		//if(configDeriv[i] != 0.0) {
			for(unsigned int j = i + 1; j < n; j++) {
				maxPathVelocity[k] = 100.0;
				if(configDeriv[i] != 0.0 || configDeriv[j] != 0.0) {
					double A_ij = configDeriv2[i] / configDeriv[i] - configDeriv2[j] / configDeriv[j];
					if(A_ij != 0.0) {
						maxPathVelocity[k] = sqrt((maxAcceleration[i] / abs(configDeriv[i]) + maxAcceleration[j] / abs(configDeriv[j])) / abs(A_ij));
					}
				}
				k++;
			}
		//}
		//else if(configDeriv2[i] != 0.0) {
		//	maxPathVelocity = min(maxPathVelocity, sqrt(maxAcceleration[i] / abs(configDeriv2[i])));
		//}

		//maxPathVelocity = min(maxPathVelocity, 5.0);
	}
	for(unsigned int i = 0; i < n; i++) {
		maxPathVelocity[k] = 100.0;
		if(configDeriv[i] != 0.0) {
			maxPathVelocity[k] = maxVelocity[i] / abs(configDeriv[i]);
		}
		k++;
	}

	return maxPathVelocity;
}



bool Trajectory2::isValid() const {
	return valid;
}

double Trajectory2::getDuration() const {
	return trajectory.back().time;
}

list<Trajectory2::TrajectoryStep>::const_iterator Trajectory2::getTrajectorySegment(double time) const {
	if(time >= trajectory.back().time) {
		list<TrajectoryStep>::const_iterator last = trajectory.end();
		last--;
		return last;
	}
	else {
		if(time < cachedTime) {
			cachedTrajectorySegment = trajectory.begin();
		}
		while(time >= cachedTrajectorySegment->time) {
			cachedTrajectorySegment++;
		}
		cachedTime = time;
		return cachedTrajectorySegment;
	}
}

VectorXd Trajectory2::getPosition(double time) const {
	list<TrajectoryStep>::const_iterator it = getTrajectorySegment(time);
	list<TrajectoryStep>::const_iterator previous = it;
	previous--;
	
	const double pathPos = previous->pathPos + (time - previous->time) * (previous->pathVel + it->pathVel) / 2.0;
	return path.getConfig(pathPos);
}