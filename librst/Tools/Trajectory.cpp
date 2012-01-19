/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author: Tobias Kunz <tobias@gatech.edu>
 * Date: 10/2011
 *
 * Humanoid Robotics Lab      Georgia Institute of Technology
 * Director: Mike Stilman     http://www.golems.org
 *
 * Algorithm details and publications:
 * http://www.golems.org/node/1570
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include "Trajectory.h"
#include <limits>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;

const double Trajectory::timeStep = 0.001;
const double Trajectory::eps = 0.000001;

Trajectory::Trajectory(const Path &path, const VectorXd &maxVelocity, const VectorXd &maxAcceleration) :
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
		double maxVelocity = getAccelerationMaxPathVelocity(s);
		if(maxVelocity == numeric_limits<double>::infinity())
			maxVelocity = 10.0;
		file << s << "  " << maxVelocity << "  " << getVelocityMaxPathVelocity(s) << endl;
	}
	file.close();
	}

	list<TrajectoryStep> startTrajectory;
	startTrajectory.push_back(TrajectoryStep(0.0, 0.0));
	double afterAcceleration = getMinMaxPathAcceleration(0.0, 0.0, true);
	while(!integrateForward(startTrajectory, afterAcceleration)) {
		double beforeAcceleration;
		TrajectoryStep switchingPoint;
		if(getNextSwitchingPoint(startTrajectory.back().pathPos, switchingPoint, beforeAcceleration, afterAcceleration)) {
			break;
		}
		list<TrajectoryStep> trajectory;
		trajectory.push_back(switchingPoint);
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

Trajectory::~Trajectory(void) {
}


// returns true if end of path is reached.
bool Trajectory::getNextSwitchingPoint(double pathPos, TrajectoryStep &nextSwitchingPoint, double &beforeAcceleration, double &afterAcceleration) {
	TrajectoryStep accelerationSwitchingPoint(pathPos, 0.0);
	double accelerationBeforeAcceleration, accelerationAfterAcceleration;
	bool accelerationReachedEnd;
	do {
		accelerationReachedEnd = getNextAccelerationSwitchingPoint(accelerationSwitchingPoint.pathPos, accelerationSwitchingPoint, accelerationBeforeAcceleration, accelerationAfterAcceleration);
		double test = getVelocityMaxPathVelocity(accelerationSwitchingPoint.pathPos);
	} while(!accelerationReachedEnd && accelerationSwitchingPoint.pathVel > getVelocityMaxPathVelocity(accelerationSwitchingPoint.pathPos));
	
	TrajectoryStep velocitySwitchingPoint(pathPos, 0.0);
	double velocityBeforeAcceleration, velocityAfterAcceleration;
	bool velocityReachedEnd;
	do {
		velocityReachedEnd = getNextVelocitySwitchingPoint(velocitySwitchingPoint.pathPos, velocitySwitchingPoint, velocityBeforeAcceleration, velocityAfterAcceleration);
	} while(!velocityReachedEnd && velocitySwitchingPoint.pathPos <= accelerationSwitchingPoint.pathPos
		&& (velocitySwitchingPoint.pathVel > getAccelerationMaxPathVelocity(velocitySwitchingPoint.pathPos - eps)
		|| velocitySwitchingPoint.pathVel > getAccelerationMaxPathVelocity(velocitySwitchingPoint.pathPos + eps)));

	if(accelerationReachedEnd && velocityReachedEnd) {
		return true;
	}
	else if(!accelerationReachedEnd && (velocityReachedEnd || accelerationSwitchingPoint.pathPos <= velocitySwitchingPoint.pathPos)) {
		nextSwitchingPoint = accelerationSwitchingPoint;
		beforeAcceleration = accelerationBeforeAcceleration;
		afterAcceleration = accelerationAfterAcceleration;
		return false;
	}
	else {
		nextSwitchingPoint = velocitySwitchingPoint;
		beforeAcceleration = velocityBeforeAcceleration;
		afterAcceleration = velocityAfterAcceleration;
		return false;
	}
}

bool Trajectory::getNextAccelerationSwitchingPoint(double pathPos, TrajectoryStep &nextSwitchingPoint, double &beforeAcceleration, double &afterAcceleration) {
	double switchingPathPos = pathPos;
	double switchingPathVel;
	while(true) {
		bool discontinuity;
		switchingPathPos = path.getNextSwitchingPoint(switchingPathPos, discontinuity);

		if(switchingPathPos > path.getLength() - eps) {
			return true;
		}
		
		if(discontinuity) {
			const double beforePathVel = getAccelerationMaxPathVelocity(switchingPathPos - eps);
			const double afterPathVel = getAccelerationMaxPathVelocity(switchingPathPos + eps);
			switchingPathVel = min(beforePathVel, afterPathVel);
			beforeAcceleration = getMinMaxPathAcceleration(switchingPathPos - eps, switchingPathVel, false);
			afterAcceleration = getMinMaxPathAcceleration(switchingPathPos + eps, switchingPathVel, true);
			
			if((beforePathVel > afterPathVel
				|| getMinMaxPhaseSlope(switchingPathPos - eps, switchingPathVel, false) > getAccelerationMaxPathVelocityDeriv(switchingPathPos - 2.0*eps))
				&& (beforePathVel < afterPathVel
				|| getMinMaxPhaseSlope(switchingPathPos + eps, switchingPathVel, true) < getAccelerationMaxPathVelocityDeriv(switchingPathPos + 2.0*eps)))
			{
				break;
			}
		}
		else {
			switchingPathVel = getAccelerationMaxPathVelocity(switchingPathPos);
			beforeAcceleration = 0.0;
			afterAcceleration = 0.0;

			if(getAccelerationMaxPathVelocityDeriv(switchingPathPos - eps) < 0.0 && getAccelerationMaxPathVelocityDeriv(switchingPathPos + eps) > 0.0) {
				break;
			}
		}
	}
	
	nextSwitchingPoint = TrajectoryStep(switchingPathPos, switchingPathVel);
	return false;
}

bool Trajectory::getNextVelocitySwitchingPoint(double pathPos, TrajectoryStep &nextSwitchingPoint, double &beforeAcceleration, double &afterAcceleration) {
	const double stepSize = 0.001;
	const double accuracy = 0.00001;

	bool start = false;
	pathPos -= stepSize;
	do {
		pathPos += stepSize;
		if(getMinMaxPhaseSlope(pathPos, getVelocityMaxPathVelocity(pathPos), false) >= getVelocityMaxPathVelocityDeriv(pathPos)) {
			start = true;
		}
	} while((!start || getMinMaxPhaseSlope(pathPos, getVelocityMaxPathVelocity(pathPos), false) > getVelocityMaxPathVelocityDeriv(pathPos))
		&& pathPos < path.getLength());

	if(pathPos >= path.getLength()) {
		return true; // end of trajectory reached
	}

	double beforePathPos = pathPos - stepSize;
	double afterPathPos = pathPos;
	while(afterPathPos - beforePathPos > accuracy) {
		pathPos = (beforePathPos + afterPathPos) / 2.0;
		if(getMinMaxPhaseSlope(pathPos, getVelocityMaxPathVelocity(pathPos), false) > getVelocityMaxPathVelocityDeriv(pathPos)) {
			beforePathPos = pathPos;
		}
		else {
			afterPathPos = pathPos;
		}
	}

	beforeAcceleration = getMinMaxPathAcceleration(beforePathPos, getVelocityMaxPathVelocity(beforePathPos), false);
	afterAcceleration = getMinMaxPathAcceleration(afterPathPos, getVelocityMaxPathVelocity(afterPathPos), true);
	nextSwitchingPoint = TrajectoryStep(afterPathPos, getVelocityMaxPathVelocity(afterPathPos));
	return false;
}

bool Trajectory::integrateForward(list<TrajectoryStep> &trajectory, double acceleration) {
	
	double pathPos = trajectory.back().pathPos;
	double pathVel = trajectory.back().pathVel;
	
	list<pair<double, bool> > switchingPoints = path.getSwitchingPoints();
	list<pair<double, bool> >::iterator nextDiscontinuity = switchingPoints.begin();

	while(true)
	{
		if(pathPos > 1.304) {
			int test = 48;
		}

		while(nextDiscontinuity != switchingPoints.end() && (nextDiscontinuity->first <= pathPos || !nextDiscontinuity->second)) {
			nextDiscontinuity++;
		}

		double oldPathPos = pathPos;
		double oldPathVel = pathVel;
		
		pathPos += timeStep * pathVel;
		pathVel += timeStep * acceleration;

		if(nextDiscontinuity != switchingPoints.end() && pathPos > nextDiscontinuity->first + eps) {
			pathVel = oldPathVel + (nextDiscontinuity->first + eps - oldPathPos) * (pathVel - oldPathVel) / (pathPos - oldPathPos);
			pathPos = nextDiscontinuity->first + eps;
		}

		//pathVel += timeStep * acceleration;
		//pathPos += timeStep * 0.5 * (oldPathVel + pathVel);

		if(pathPos > path.getLength()) {
			return true;
		}
		else if(pathVel < 0.0) {
			valid = false;
			cout << "error" << endl;
			return true;
		}

		double test1 = getMinMaxPhaseSlope(oldPathPos, getVelocityMaxPathVelocity(oldPathPos), false);
		double test2 = getVelocityMaxPathVelocityDeriv(oldPathPos);

		if(pathVel > getVelocityMaxPathVelocity(pathPos)
			&& getMinMaxPhaseSlope(oldPathPos, getVelocityMaxPathVelocity(oldPathPos), false) <= getVelocityMaxPathVelocityDeriv(oldPathPos))
		{
			pathVel = getVelocityMaxPathVelocity(pathPos);
		}

		trajectory.push_back(TrajectoryStep(pathPos, pathVel));
		acceleration = getMinMaxPathAcceleration(pathPos, pathVel, true);

		if(pathVel > getAccelerationMaxPathVelocity(pathPos) || pathVel > getVelocityMaxPathVelocity(pathPos)) {
			// find more accurate intersection with max-velocity curve using bisection
			TrajectoryStep overshoot = trajectory.back();
			trajectory.pop_back();
			double slope = getSlope(trajectory.back(), overshoot);
			double before = trajectory.back().pathPos;
			double after = overshoot.pathPos;
			while(after - before > 0.00001) {
				const double midpoint = 0.5 * (before + after);
				double midpointPathVel = trajectory.back().pathVel + slope * (midpoint - trajectory.back().pathPos);

				if(midpointPathVel > getVelocityMaxPathVelocity(midpoint)
					&& getMinMaxPhaseSlope(before, getVelocityMaxPathVelocity(before), false) <= getVelocityMaxPathVelocityDeriv(before))
				{
					midpointPathVel = getVelocityMaxPathVelocity(midpoint);
				}

				if(midpointPathVel > getAccelerationMaxPathVelocity(midpoint) || midpointPathVel > getVelocityMaxPathVelocity(midpoint))
					after = midpoint;
				else
					before = midpoint;
			}
			trajectory.push_back(TrajectoryStep(before, trajectory.back().pathVel + slope * (before - trajectory.back().pathPos)));

			if(getAccelerationMaxPathVelocity(after) < getVelocityMaxPathVelocity(after)) {
				if(after > nextDiscontinuity->first) {
					return false;
				}
				else if(getMinMaxPhaseSlope(trajectory.back().pathPos, trajectory.back().pathVel, true) > getAccelerationMaxPathVelocityDeriv(trajectory.back().pathPos)) {
					return false;
				}
			}
			else {
				if(getMinMaxPhaseSlope(trajectory.back().pathPos, trajectory.back().pathVel, false) > getVelocityMaxPathVelocityDeriv(trajectory.back().pathPos)) {
					return false;
				}
			}
			
		}
	}
}


void Trajectory::integrateBackward(list<TrajectoryStep> &trajectory, list<TrajectoryStep> &startTrajectory, double acceleration) {
	list<TrajectoryStep>::reverse_iterator before = startTrajectory.rbegin();
	double pathPos = trajectory.front().pathPos;
	double pathVel = trajectory.front().pathVel;

	while(true)
	{
		pathPos -= timeStep * pathVel;
		pathVel -= timeStep * acceleration;

		//double oldPathVel = pathVel;
		//pathVel -= timeStep * acceleration;
		//pathPos -= timeStep * 0.5 * (oldPathVel + pathVel);

		trajectory.push_front(TrajectoryStep(pathPos, pathVel));
		acceleration = getMinMaxPathAcceleration(pathPos, pathVel, false);

		if(pathVel < 0.0 || pathPos < 0.0) {
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
		else if(pathVel > getAccelerationMaxPathVelocity(pathPos) || pathVel > getVelocityMaxPathVelocity(pathPos)) {
			ofstream file("trajectory.txt");
			for(list<TrajectoryStep>::iterator it = startTrajectory.begin(); it != startTrajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			for(list<TrajectoryStep>::iterator it = trajectory.begin(); it != trajectory.end(); it++) {
				file << it->pathPos << "  " << it->pathVel << endl;
			}
			file.close();
			cout << "error" << endl;
			valid = false;
			return;
		}
	}
}

inline double Trajectory::getSlope(const TrajectoryStep &point1, const TrajectoryStep &point2) {
	return (point2.pathVel - point1.pathVel) / (point2.pathPos - point1.pathPos);
}

inline double Trajectory::getSlope(list<TrajectoryStep>::const_iterator lineEnd) {
	list<TrajectoryStep>::const_iterator lineStart = lineEnd;
	lineStart--;
	return getSlope(*lineStart, *lineEnd);
}

Trajectory::TrajectoryStep Trajectory::getIntersection(const list<TrajectoryStep> &trajectory, list<TrajectoryStep>::iterator &it, const TrajectoryStep &linePoint1, const TrajectoryStep &linePoint2) {
	
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


double Trajectory::getMinMaxPathAcceleration(double pathPos, double pathVel, bool max) {
	VectorXd configDeriv = path.getTangent(pathPos);
	VectorXd configDeriv2 = path.getCurvature(pathPos);
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

double Trajectory::getMinMaxPhaseSlope(double pathPos, double pathVel, bool max) {
	return getMinMaxPathAcceleration(pathPos, pathVel, max) / pathVel;
}

double Trajectory::getAccelerationMaxPathVelocity(double pathPos) {
	double maxPathVelocity = numeric_limits<double>::infinity();
	const VectorXd configDeriv = path.getTangent(pathPos);
	const VectorXd configDeriv2 = path.getCurvature(pathPos);
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
	}
	return maxPathVelocity;
}


double Trajectory::getVelocityMaxPathVelocity(double pathPos) {
	const VectorXd tangent = path.getTangent(pathPos);
	double maxPathVelocity = numeric_limits<double>::max();
	for(unsigned int i = 0; i < n; i++) {
		maxPathVelocity = min(maxPathVelocity, maxVelocity[i] / abs(tangent[i]));
	}
	return maxPathVelocity;
}

double Trajectory::getAccelerationMaxPathVelocityDeriv(double pathPos) {
	return (getAccelerationMaxPathVelocity(pathPos + eps) - getAccelerationMaxPathVelocity(pathPos - eps)) / (2.0 * eps);
}

double Trajectory::getVelocityMaxPathVelocityDeriv(double pathPos) {
	const VectorXd tangent = path.getTangent(pathPos);
	double maxPathVelocity = numeric_limits<double>::max();
	unsigned int activeConstraint;
	for(unsigned int i = 0; i < n; i++) {
		const double thisMaxPathVelocity = maxVelocity[i] / abs(tangent[i]);
		if(thisMaxPathVelocity < maxPathVelocity) {
			maxPathVelocity = thisMaxPathVelocity;
			activeConstraint = i;
		}
	}
	return - (maxVelocity[activeConstraint] * path.getCurvature(pathPos)[activeConstraint])
		/ (tangent[activeConstraint] * abs(tangent[activeConstraint]));
}

bool Trajectory::isValid() const {
	return valid;
}

double Trajectory::getDuration() const {
	return trajectory.back().time;
}

list<Trajectory::TrajectoryStep>::const_iterator Trajectory::getTrajectorySegment(double time) const {
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

VectorXd Trajectory::getPosition(double time) const {
	list<TrajectoryStep>::const_iterator it = getTrajectorySegment(time);
	list<TrajectoryStep>::const_iterator previous = it;
	previous--;
	
	const double pathPos = previous->pathPos + (time - previous->time) * (previous->pathVel + it->pathVel) / 2.0;
	return path.getConfig(pathPos);
}