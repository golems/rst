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
	list<TrajectoryStep> endTrajectory;

	startTrajectory.push_back(TrajectoryStep(0.0, 0.0));
	endTrajectory.push_front(TrajectoryStep(path.getLength(), 0.0));


	integrateBackward(endTrajectory, startTrajectory, true);
	while(!integrateForward(startTrajectory, endTrajectory)) {
		list<TrajectoryStep> trajectory;
		trajectory.push_front(getNextTangentPoint(startTrajectory.back().pathPos));
		integrateBackward(trajectory, startTrajectory, false);
	}
	this->trajectory = endTrajectory;

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
	

	ofstream file("trajectory.txt");
	for(list<TrajectoryStep>::iterator it = endTrajectory.begin(); it != endTrajectory.end(); it++) {
		file << it->pathPos << "  " << it->pathVel << endl;
	}
	file.close();
}

Trajectory2::~Trajectory2(void) {
}

static inline double squared(double d) {
	return d * d;
}


Trajectory2::TrajectoryStep Trajectory2::getNextTangentPoint(double before) {
	//double switchingPathPos = path.getNextSwitchingPoint(before);
	//return TrajectoryStep(switchingPathPos, getMaxPathVelocity(switchingPathPos - 0.00000001));


	const double stepSize = 0.001;
	const double accuracy = 0.00001;

	before += accuracy;
	double after = before;

	{
		double maxPathVel = getMaxPathVelocity(after);
		double nextMaxPathVel = getMaxPathVelocity(after + stepSize);
		while(getMaxPathVelocityDeriv(after, maxPathVel) - (nextMaxPathVel - maxPathVel) / stepSize > 0.0) {
			before = after;
			after += stepSize;
			maxPathVel = nextMaxPathVel;
			nextMaxPathVel = getMaxPathVelocity(after + stepSize);
		}
	}

	double pathPos[5];
	pathPos[0] = before;
	pathPos[4] = after + stepSize;
	double maxPathVel[5];
	while(pathPos[4] - pathPos[0] > accuracy) {
		for(int i = 1; i < 4; i++) {
			pathPos[i] = pathPos[0] + (double)i / 4 * (pathPos[4] - pathPos[0]);
		}
		
		for(int i = 0; i < 5; i++) {
			maxPathVel[i] = getMaxPathVelocity(pathPos[i]);
		}

		int j = 0;
		while(j < 3 && getMaxPathVelocityDeriv(pathPos[j], maxPathVel[j]) - (maxPathVel[j+1] - maxPathVel[j]) / (pathPos[j+1] - pathPos[j]) > 0.0) {
			j++;
		}
		if(j > 0) {
			pathPos[0] = pathPos[j-1];
		}
		pathPos[4] = pathPos[j+1];
	}
	
	return TrajectoryStep(pathPos[0], getMaxPathVelocity(pathPos[0]));
}


bool Trajectory2::integrateForward(list<TrajectoryStep> &trajectory, list<TrajectoryStep> &endTrajectory) {
	
	bool checkFeasibility = false;
	list<TrajectoryStep>::iterator after = endTrajectory.begin();
	list<TrajectoryStep>::iterator after2nd = endTrajectory.begin();
	double pathPos = trajectory.back().pathPos;
	double pathVel = trajectory.back().pathVel;
	
	while(true)
	{
		if(pathVel < getMaxPathVelocity(pathPos))
			checkFeasibility = true;
		double oldPathVel = pathVel;

		if(trajectory.size() == 1)
			pathVel += timeStep * getMaxPathAcceleration(pathPos + eps, pathVel);
		else
			pathVel += timeStep * getMaxPathAcceleration(pathPos, pathVel);
		pathPos += timeStep * 0.5 * (oldPathVel + pathVel);

		//Vector2d state = rungeKutta4(Vector2d(pathPos, pathVel), timeStep, true, checkFeasibility);
		//pathPos = state[0];
		//pathVel = state[1];

		trajectory.push_back(TrajectoryStep(pathPos, pathVel));
		after2nd = after;
		while(after->pathPos < pathPos) {
			after++;
		}


		if(after != endTrajectory.begin() && pathVel >= after->pathVel + getSlope(after) * (pathPos - after->pathPos)) {
			TrajectoryStep overshoot = trajectory.back();
			trajectory.pop_back();
			TrajectoryStep intersection = getIntersection(endTrajectory, after2nd, trajectory.back(), overshoot);
			endTrajectory.erase(endTrajectory.begin(), after2nd);
			endTrajectory.push_front(intersection);
			endTrajectory.splice(endTrajectory.begin(), trajectory);
			return true;
		}
		else if(pathVel < 0.0 || pathPos > path.getLength()) {
			ofstream file("trajectory.txt");
			for(list<TrajectoryStep>::iterator it = trajectory.begin(); it != trajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			for(list<TrajectoryStep>::iterator it = endTrajectory.begin(); it != endTrajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			file.close();
			valid = false;
			cout << "error " << path.getLength() << endl;
			return true;
		}
		else if(checkFeasibility && pathVel >= getMaxPathVelocity(pathPos)) {
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
			trajectory.push_back(TrajectoryStep(after, trajectory.back().pathVel + slope * (after - trajectory.back().pathPos)));
			return false;
		}
	}
}


void Trajectory2::integrateBackward(list<TrajectoryStep> &trajectory, list<TrajectoryStep> &startTrajectory, bool checkFeasibility) {
	list<TrajectoryStep>::reverse_iterator before = startTrajectory.rbegin();
	double pathPos = trajectory.front().pathPos;
	double pathVel = trajectory.front().pathVel;


	while(true)
	{
		double oldPathVel = pathVel;

		if(trajectory.size() == 1)
			pathVel -= timeStep * getMinPathAcceleration(pathPos - eps, pathVel);
		else
			pathVel -= timeStep * getMinPathAcceleration(pathPos, pathVel);
		pathPos -= timeStep * 0.5 * (oldPathVel + pathVel);
		
		//Vector2d state = rungeKutta4(Vector2d(pathPos, pathVel), timeStep, false, checkFeasibility);
		//pathPos = state[0];
		//pathVel = state[1];

		trajectory.push_front(TrajectoryStep(pathPos, pathVel));
		
		if(pathVel < 0.0 || pathPos < 0.0) {
			if(!checkFeasibility) {
				valid = false;
				cout << "error " << pathPos << " " << pathVel << endl;
			}
			return;
		}

		while(before != startTrajectory.rend() && before->pathPos > pathPos) {
			before++;
		}

		// debug
		if(before == startTrajectory.rend()) {
			ofstream file("trajectory.txt");
			for(list<TrajectoryStep>::iterator it = startTrajectory.begin(); it != startTrajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			for(list<TrajectoryStep>::iterator it = trajectory.begin(); it != trajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			file.close();
		}
		

		if(before != startTrajectory.rbegin() && pathVel >= before->pathVel + getSlope(before.base()) * (pathPos - before->pathPos)) {
			TrajectoryStep overshoot = trajectory.front();
			trajectory.pop_front();
			list<TrajectoryStep>::iterator after = before.base();
			TrajectoryStep intersection = getIntersection(startTrajectory, after, overshoot, trajectory.front());

			//double startEnd = after->pathPos;
			//double trajectoryStart = trajectory.front().pathPos;
			//if(abs(startEnd - trajectoryStart) > 0.1) {
			//	cout << "error " << startEnd << " " << trajectoryStart << " " << startTrajectory.back().pathPos << endl;
			//	valid = false;
			//}

			if(after != startTrajectory.end()) {
				startTrajectory.erase(after, startTrajectory.end());
				startTrajectory.push_back(intersection);
			}
			startTrajectory.splice(startTrajectory.end(), trajectory);

			if(!valid) {
				ofstream file("trajectory.txt");
				for(list<TrajectoryStep>::iterator it = startTrajectory.begin(); it != startTrajectory.end(); it++) {
					file << it->pathPos << "  " << it->pathVel << endl;
				}
				file.close();
			}
			return;
		}
		else if(checkFeasibility && pathVel >= getMaxPathVelocity(pathPos)) {
			return;
		}

	}
}


inline Vector2d Trajectory2::getMinMaxStateDerivs(const Vector2d& state, bool max) {
	return Vector2d(state[1], getMinMaxPathAcceleration(state[0], state[1], max));
}

inline Vector2d Trajectory2::rungeKutta4(const Vector2d &y, double step, bool max, bool checkFeasibility) {
	
	const Vector2d k1 = step * getMinMaxStateDerivs(y, max);
	if(checkFeasibility && k1[1] > getMaxPathVelocity(k1[0]))
		return k1;
	const Vector2d k2 = step * getMinMaxStateDerivs(y + 0.5 * k1, max);
	if(checkFeasibility && k2[1] > getMaxPathVelocity(k2[0]))
		return k2;
	const Vector2d k3 = step * getMinMaxStateDerivs(y + 0.5 * k2, max);
	if(checkFeasibility && k3[1] > getMaxPathVelocity(k3[0]))
		return k3;
	const Vector2d k4 = step * getMinMaxStateDerivs(y + k3, max);
	if(checkFeasibility && k4[1] > getMaxPathVelocity(k4[0]))
		return k4;
	return y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
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
		//if(it == trajectoryDebug.end()) {
		//	valid = false;
		//	cout << "error " << endl;
		//	cout << trajectoryParam->pathPos << " " << trajectoryParam->pathVel << endl;
		//	cout << trajectoryDebug.back().pathPos << " " << trajectoryDebug.back().pathVel << endl;
		//	cout << linePoint1.pathPos << " " << linePoint1.pathVel << endl;
		//	cout << linePoint2.pathPos << " " << linePoint2.pathVel << endl;
		//	trajectory = trajectoryDebug.end();
		//	return TrajectoryStep(0.0, 0.0);
		//}
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
				maxAcceleration[i]/abs(configDeriv[i]) - factor * configDeriv2[i] * squared(pathVel) / configDeriv[i]);
		}
	}
	return factor * maxPathAcceleration;
}

double Trajectory2::getMinPathAcceleration(double pathPos, double pathVel) {
	//return (-maxAcceleration.cwiseQuotient(path.getConfigDeriv(pathPos).cwiseAbs())
	//	- (path.getConfigDeriv2(pathPos) * squared(pathVel)).cwiseQuotient(path.getConfigDeriv(pathPos))).maxCoeff();

	VectorXd configDeriv = path.getConfigDeriv(pathPos);
	VectorXd configDeriv2 = path.getConfigDeriv2(pathPos);
	double minPathAcceleration = -numeric_limits<double>::max();
	for(unsigned int i = 0; i < n; i++) {
		if(configDeriv[i] != 0.0) {
			minPathAcceleration = max(minPathAcceleration,
				-maxAcceleration[i]/abs(configDeriv[i]) - configDeriv2[i] * squared(pathVel) / configDeriv[i]);
		}
	}
	return minPathAcceleration;
}

double Trajectory2::getMaxPathAcceleration(double pathPos, double pathVel) {
	//return (maxAcceleration.cwiseQuotient(path.getConfigDeriv(pathPos).cwiseAbs())
	//	- (path.getConfigDeriv2(pathPos) * squared(pathVel)).cwiseQuotient(path.getConfigDeriv(pathPos))).minCoeff();

	VectorXd configDeriv = path.getConfigDeriv(pathPos);
	VectorXd configDeriv2 = path.getConfigDeriv2(pathPos);
	double maxPathAcceleration = numeric_limits<double>::max();
	for(unsigned int i = 0; i < n; i++) {
		if(configDeriv[i] != 0.0) {
			maxPathAcceleration = min(maxPathAcceleration,
				maxAcceleration[i]/abs(configDeriv[i]) - configDeriv2[i] * squared(pathVel) / configDeriv[i]);
		}
	}
	return maxPathAcceleration;
}


double Trajectory2::getMinPathVelocityDeriv(double pathPos, double pathVel) {
	return getMinPathAcceleration(pathPos, pathVel) / pathVel;
}

double Trajectory2::getMaxPathVelocityDeriv(double pathPos, double pathVel) {
	return getMaxPathAcceleration(pathPos, pathVel) / pathVel;
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

		//maxPathVelocity = min(maxPathVelocity, 5.0);
	}
	return maxPathVelocity;
}

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


double Trajectory2::getMaxPathVelocity_Deriv(double pathPos) {
	return (getMaxPathVelocity(pathPos + 0.5*eps) - getMaxPathVelocity(pathPos - 0.5*eps)) / eps;
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