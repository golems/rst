#pragma once

#include <list>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>

#include <iostream>

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
	virtual double getNextSwitchingPoint(double s, bool &discontinuity) const = 0;
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

	Eigen::VectorXd getConfigDeriv(double /* s */) const {
		return (end - start) / length;
	}

	Eigen::VectorXd getConfigDeriv2(double /* s */) const {
		return Eigen::VectorXd::Zero(start.size());
	}

	double getNextSwitchingPoint(double /* s */, bool &discontinuity) const {
		discontinuity = true;
		return length;
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
		if((intersection - start).norm() < 0.000001 || (end - intersection).norm() < 0.000001) {
			length = 0.0;
			radius = 1.0;
			center = intersection;
			x = Eigen::VectorXd::Zero(start.size());
			y = Eigen::VectorXd::Zero(start.size());
			return;
		}

		const Eigen::VectorXd startDirection = (intersection - start).normalized();
		const Eigen::VectorXd endDirection = (end - intersection).normalized();

		if((startDirection - endDirection).norm() < 0.000001) {
			length = 0.0;
			radius = 1.0;
			center = intersection;
			x = Eigen::VectorXd::Zero(start.size());
			y = Eigen::VectorXd::Zero(start.size());
			return;
		}

		double distance = std::min((start - intersection).norm(), (end - intersection).norm());
		const double angle = acos(startDirection.dot(endDirection));

		distance = std::min(distance, maxDeviation * sin(0.5 * angle) / (1.0 - cos(0.5 * angle)));  // enforce max deviation

		radius = distance / tan(0.5 * angle);
		length = angle * radius;

		center = intersection + (endDirection - startDirection).normalized() * radius / cos(0.5 * angle);
		x = (intersection - distance * startDirection - center).normalized();
		y = intersection - center;
		y = (y - x * y.dot(x)).normalized();

		// calculate switching points
		const double dim = start.size();
		for(unsigned int i = 0; i < dim; i++) {
			double switchingAngle = atan2(y[i], x[i]);
			if(switchingAngle < 0.0) {
				switchingAngle += M_PI;
			}
			const double switchingPoint = switchingAngle * radius;
			if(switchingPoint < length) {
				switchingPoints.push_back(switchingPoint);
			}
		}
		switchingPoints.sort();

		//debug
		double dotStart = startDirection.dot((intersection - getConfig(0.0)).normalized());
		double dotEnd = endDirection.dot((getConfig(length) - intersection).normalized());
		if(abs(dotStart - 1.0) > 0.0001 || abs(dotEnd - 1.0) > 0.0001) {
			std::cout << "Error\n";
		}
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

	double getNextSwitchingPoint(double s, bool &discontinuity) const {
		discontinuity = false;
		for(std::list<double>::const_iterator it = switchingPoints.begin(); it != switchingPoints.end(); it++) {
			if(*it > s) {
				return *it;
			}
		}
		discontinuity = true;
		return length;
	}

	CircularPathSegment* clone() const {
		return new CircularPathSegment(*this);
	}

private:
	double radius;
	Eigen::VectorXd center;
	Eigen::VectorXd x;
	Eigen::VectorXd y;
	std::list<double> switchingPoints;
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
	double getNextSwitchingPoint(double s, bool &discontinuity) const;
private:
	PathSegment* getPathSegment(double &s) const;
	
	double length;
public:
	std::list<PathSegment*> pathSegments;
};