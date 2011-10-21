#include "Trajectory.h"

using namespace std;
using namespace Eigen;

// Partially based on
// John J. Craig. Introduction to Robotics - Mechanics and Control, 3rd Edition. Chapter 7.3

Trajectory::Trajectory(const list<VectorXd> &_path, const VectorXd &maxVelocity, const VectorXd &maxAcceleration, bool slowDown, bool removeWayPoints) :
	path(_path.begin(), _path.end()),
	velocities(path.size() - 1),
	accelerations(path.size()),
	durations(path.size() - 1),
	blendDurations(path.size()),
	duration(0.0)
{
	// calculate time between waypoints and initial velocities of linear segments
	for(int i = 0; i < path.size() - 1; i++) {
		durations[i] = 0.0;
		for(int j = 0; j < path[i].size(); j++) {
			durations[i] = max(durations[i], abs(path[i+1][j] - path[i][j]) / maxVelocity[j]);
		}
		velocities[i] = (path[i+1] - path[i]) / durations[i];
	}

	if(slowDown) {
		
		vector<double> slowDownFactors(path.size() - 1, 1.0);

		for(int i = 0; i < path.size(); i++) {
			// calculate initial blend duration
			VectorXd previousVelocity = (i == 0) ? VectorXd::Zero(path[i].size()) : velocities[i-1];
			VectorXd nextVelocity = (i == path.size() - 1) ? VectorXd::Zero(path[i].size()) : velocities[i];
			blendDurations[i] = 0.0;
			for(int j = 0; j < path[i].size(); j++) {
				blendDurations[i] = max(blendDurations[i], abs(nextVelocity[j] - previousVelocity[j]) / maxAcceleration[j]);
			}

			// calculate maximum allowable blend duration for initial linear velocities
			double maxDuration = numeric_limits<double>::max();
			if(i > 0) {
				maxDuration = min(maxDuration, durations[i-1]);
			}
			if(i < path.size() - 1) {
				maxDuration = min(maxDuration, durations[i]);
			}
			
			// calculate slow down factors for neighboring linear velocities such that
			// the blend phase replaces at most half of the neighboring linear segments
			if(blendDurations[i] > maxDuration) {
				double slowDownFactor = sqrt(maxDuration / blendDurations[i]);
				if(i > 0) {
					slowDownFactors[i-1] = min(slowDownFactors[i-1], slowDownFactor);
				}
				if(i < path.size() - 1) {
					slowDownFactors[i] = slowDownFactor;
				}
			}
		}

		// apply slow down factors to linear segments
		for(int i = 0; i < path.size() - 1; i++) {
			velocities[i] *= slowDownFactors[i];
			durations[i] /= slowDownFactors[i];
		}

	}



	// calculate final blend durations
	valid = true;
	for(int i = 0; i < path.size(); i++) {
		VectorXd previousVelocity = (i == 0) ? VectorXd::Zero(path[i].size()) : velocities[i-1];
		VectorXd nextVelocity = (i == path.size() - 1) ? VectorXd::Zero(path[i].size()) : velocities[i];
		blendDurations[i] = 0.0;
		for(int j = 0; j < path[i].size(); j++) {
			blendDurations[i] = max(blendDurations[i], abs(nextVelocity[j] - previousVelocity[j]) / maxAcceleration[j]);
		}
		if((i > 0 && blendDurations[i] > durations[i-1] + 0.000001)
			|| (i < path.size() - 1 && blendDurations[i] > durations[i] + 0.000001))
		{
			valid = false;
		}

		accelerations[i] = (nextVelocity - previousVelocity) / blendDurations[i];
	}


	// calculate total time of trajectory
	for(int i = 0; i < path.size() - 1; i++) {
		duration += durations[i];
	}
	duration += 0.5 * blendDurations.front();
	duration += 0.5 * blendDurations.back();
}

bool Trajectory::isValid() {
	return valid;
}

VectorXd Trajectory::getPosition(double time) const {
	if(time > duration) {
		return path.back();
	}
	double t = time;
	if(t <= 0.5 * blendDurations[0]) {
		return path[0] + 0.5 * t * t * accelerations[0];
	}
	else {
		t -= 0.5 * blendDurations[0];
	}
	int i = 0;
	while(i < path.size() - 1 && t > durations[i]) {
		t -= durations[i];
		i++;
	}
	if(i == path.size() - 1) {
		t = 0.5 * blendDurations.back() - t;
		return path.back() + 0.5 * t * t * accelerations.back();
	}

	double switchingTime1 = 0.5 * blendDurations[i];
	double switchingTime2 = durations[i] - 0.5 * blendDurations[i+1];

	if(t < switchingTime1) {
		t = switchingTime1 - t;
		return path[i] + switchingTime1 * velocities[i] - t * velocities[i] + 0.5 * t * t * accelerations[i];
	}
	else if(t > switchingTime2) {
		t -= switchingTime2;
		return path[i] + switchingTime2 * velocities[i] + t * velocities[i] + 0.5 * t * t * accelerations[i+1];
	}
	else {
		return path[i] + t * velocities[i];
	}
}


VectorXd Trajectory::getVelocity(double time) const {
	if(time > duration) {
		return VectorXd::Zero(path.back().size());
	}
	double t = time;
	if(t <= 0.5 * blendDurations[0]) {
		return t * accelerations[0];
	}
	else {
		t -= 0.5 * blendDurations[0];
	}
	int i = 0;
	while(i < path.size() - 1 && t > durations[i]) {
		t -= durations[i];
		i++;
	}
	if(i == path.size() - 1) {
		t = 0.5 * blendDurations.back() - t;
		return - t * accelerations.back();
	}

	double switchingTime1 = 0.5 * blendDurations[i];
	double switchingTime2 = durations[i] - 0.5 * blendDurations[i+1];

	if(t < switchingTime1) {
		t = switchingTime1 - t;
		return velocities[i] - t * accelerations[i];
	}
	else if(t > switchingTime2) {
		t -= switchingTime2;
		return velocities[i] + t * accelerations[i+1];
	}
	else {
		return velocities[i];
	}
}


double Trajectory::getDuration() const {
	return duration;
}