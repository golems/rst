#pragma once

#include <list>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>

class PathSegment
{
public:
	PathSegment(double length = 0.0) :
		length(length)
	{
	}
	
	virtual ~PathSegment() {}

	double getLength() const {
		return length;
	}
	virtual Eigen::VectorXd getConfig(double s) const = 0;
	virtual Eigen::VectorXd getConfigDeriv(double s) const = 0;
	virtual Eigen::VectorXd getConfigDeriv2(double s) const = 0;
	virtual PathSegment* clone() const = 0;
protected:
	double length;
};

class LinearPathSegment : public PathSegment
{
public:
	LinearPathSegment(const Eigen::VectorXd &start, const Eigen::VectorXd &end) :
		start(start),
		end(end),
		PathSegment((end-start).norm())
	{
	}

	Eigen::VectorXd getConfig(double s) const {
		s /= length;
		s = std::max(0.0, std::min(1.0, s));
		return (1.0 - s) * start + s * end;
	}

	Eigen::VectorXd getConfigDeriv(double s) const {
		return (end - start) / length;
	}

	Eigen::VectorXd getConfigDeriv2(double s) const {
		return Eigen::VectorXd::Zero(start.size());
	}

	LinearPathSegment* clone() const {
		return new LinearPathSegment(*this);
	}

private:
	Eigen::VectorXd start;
	Eigen::VectorXd end;
};

class CircularPathSegment : public PathSegment
{
public:
	CircularPathSegment(const Eigen::VectorXd &start, const Eigen::VectorXd &intersection, const Eigen::VectorXd &end, double maxDeviation) {
		
		Eigen::VectorXd startDirection = (intersection - start).normalized();
		Eigen::VectorXd endDirection = (end - intersection).normalized();

		double distance = std::min((start - intersection).norm(), (end - intersection).norm());
		double angle = acos(startDirection.dot(endDirection));
		distance = std::min(distance, maxDeviation * sin(0.5 * angle) / (1.0 - cos(0.5 * angle)));  // enforce max deviation

		radius = distance / tan(0.5 * angle);
		length = angle * radius;

		center = intersection + (endDirection - startDirection).normalized() * radius / cos(0.5 * angle);
		x = (intersection - distance * startDirection - center).normalized();
		y = end - center;
		y = (y - x * y.dot(x)).normalized();
	}

	Eigen::VectorXd getConfig(double s) const {
		const double angle = s / radius;
		return center + radius * (x * cos(angle) + y * sin(angle));
	}

	Eigen::VectorXd getConfigDeriv(double s) const {
		const double angle = s / radius;
		return - x * sin(angle) + y * cos(angle);
	}

	Eigen::VectorXd getConfigDeriv2(double s) const {
		const double angle = s / radius;
		return - 1.0 / radius * (x * cos(angle) + y * sin(angle));
	}

	CircularPathSegment* clone() const {
		return new CircularPathSegment(*this);
	}

private:
	double radius;
	Eigen::VectorXd center;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
};


class Path
{
public:
	Path(const std::list<Eigen::VectorXd> &path, double maxDeviation = 0.0);
	Path(const Path &path);
	~Path();
	double getLength() const;
	Eigen::VectorXd getConfig(double s) const;
	Eigen::VectorXd getConfigDeriv(double s) const;
	Eigen::VectorXd getConfigDeriv2(double s) const;
private:
	PathSegment* getPathSegment(double &s) const;
	
	double length;
	std::list<PathSegment*> pathSegments;
};