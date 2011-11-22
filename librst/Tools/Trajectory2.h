#pragma once

#include <Eigen/Core>
#include "Tools/Path.h"



class Trajectory2
{
public:
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

	Trajectory2(const Path &path, const Eigen::VectorXd &maxVelocity, const Eigen::VectorXd &maxAcceleration);
	~Trajectory2(void);
	double getMinPathAcceleration(double pathPosition, double pathVelocity);
	double getMaxPathAcceleration(double pathPosition, double pathVelocity);
	double getMinMaxPathAcceleration(double pathPosition, double pathVelocity, bool max);
	double getMaxPathVelocity(double pathPos);
	double getMinPathVelocityDeriv(double pathPos, double pathVel);
	double getMaxPathVelocityDeriv(double pathPos, double pathVel);
	Eigen::VectorXd getMaxPathVelocities(double pathPos);
	double getMaxPathVelocity_Deriv(double pathPos);
	TrajectoryStep getNextTangentPoint(double pathPos);
	bool integrateForward(std::list<TrajectoryStep> &trajectory, std::list<TrajectoryStep> &endTrajectory);
	void integrateBackward(std::list<TrajectoryStep> &trajectory, std::list<TrajectoryStep> &startTrajectory, bool checkFeasibility);
	TrajectoryStep getIntersection(const std::list<TrajectoryStep> &trajectory, std::list<TrajectoryStep>::iterator &it, const TrajectoryStep &linePoint1, const TrajectoryStep &linePoint2);
	bool isValid() const;
	double getDuration() const;
	Eigen::VectorXd getPosition(double time) const;
private:
	inline double getSlope(const TrajectoryStep &point1, const TrajectoryStep &point2);
	inline double getSlope(std::list<TrajectoryStep>::const_iterator lineEnd);
	inline Eigen::Vector2d getMinMaxStateDerivs(const Eigen::Vector2d& state, bool max);
	inline Eigen::Vector2d rungeKutta4(const Eigen::Vector2d &y, double step, bool max, bool checkFeasibility);
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