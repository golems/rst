#pragma once

#include <Eigen/Core>
#include "Tools/Path.h"

class Trajectory2
{
public:
	Trajectory2(const Path &path, const Eigen::VectorXd &maxVelocity, const Eigen::VectorXd &maxAcceleration);
	~Trajectory2(void);

	bool isValid() const;
	double getDuration() const;
	Eigen::VectorXd getPosition(double time) const;

private:
	struct TrajectoryStep {
		TrajectoryStep() {}
		TrajectoryStep(double pathPos, double pathVel) :
			pathPos(pathPos),
			pathVel(pathVel)
		{}
		double pathPos;
		double pathVel;
		double time;
	};

	std::list<TrajectoryStep> getNextSwitchingPoint(double pathPos, double &beforeAcceleration, double &afterAcceleration);
	bool integrateForward(std::list<TrajectoryStep> &trajectory, double acceleration);
	void integrateBackward(std::list<TrajectoryStep> &trajectory, std::list<TrajectoryStep> &startTrajectory, double acceleration);
	double getMinMaxPathAcceleration(double pathPosition, double pathVelocity, bool max);
	double getMaxPathVelocity(double pathPos);
	Eigen::VectorXd getMaxPathVelocities(double pathPos);
	
	TrajectoryStep getIntersection(const std::list<TrajectoryStep> &trajectory, std::list<TrajectoryStep>::iterator &it, const TrajectoryStep &linePoint1, const TrajectoryStep &linePoint2);
	inline double getSlope(const TrajectoryStep &point1, const TrajectoryStep &point2);
	inline double getSlope(std::list<TrajectoryStep>::const_iterator lineEnd);
	
	std::list<TrajectoryStep>::const_iterator getTrajectorySegment(double time) const;
	
	Path path;
	Eigen::VectorXd maxVelocity;
	Eigen::VectorXd maxAcceleration;
	unsigned int n;
	bool valid;
	std::list<TrajectoryStep> trajectory;

	static const double eps;
	static const double timeStep;

	mutable double cachedTime;
	mutable std::list<TrajectoryStep>::const_iterator cachedTrajectorySegment;
};