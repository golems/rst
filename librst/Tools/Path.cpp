#include "Path.h"
#include <fstream>


using namespace std;
using namespace Eigen;

Path::Path(const list<VectorXd> &path, double maxDeviation) :
	length(0.0)
{
	if(path.size() < 2)
		return;
	list<VectorXd>::const_iterator config1 = path.begin();
	list<VectorXd>::const_iterator config2 = config1;
	config2++;
	list<VectorXd>::const_iterator config3;
	VectorXd startConfig = *config1;
	while(config2 != path.end()) {
		config3 = config2;
		config3++;
		if(maxDeviation > 0.0 && config3 != path.end()) {
			CircularPathSegment* blendSegment = new CircularPathSegment(0.5 * (*config1 + *config2), *config2, 0.5 * (*config2 + *config3), maxDeviation);
			VectorXd endConfig = blendSegment->getConfig(0.0);
			pathSegments.push_back(new LinearPathSegment(startConfig, endConfig));
			pathSegments.push_back(blendSegment);
			startConfig = blendSegment->getConfig(blendSegment->getLength());

			//debug
			if(((endConfig - *config1).norm() > 0.000001 && (*config2 - endConfig).norm() > 0.000001
				&& abs((endConfig - *config1).normalized().dot((*config2 - endConfig).normalized()) - 1.0) > 0.000001)
				|| ((startConfig - *config2).norm() > 0.000001 && (*config3 - startConfig).norm() > 0.000001
				&& abs((startConfig - *config2).normalized().dot((*config3 - startConfig).normalized()) - 1.0) > 0.000001)) {
					cout << "error" << endl;
			}
		}
		else {
			pathSegments.push_back(new LinearPathSegment(startConfig, *config2));
			startConfig = *config2;
		}
		config1 = config2;
		config2++;
	}
	for(list<PathSegment*>::const_iterator it = pathSegments.begin(); it != pathSegments.end(); it++) {
		length += (*it)->getLength();
	}
}

Path::Path(const Path &path) :
	length(path.length)
{
	for(list<PathSegment*>::const_iterator it = path.pathSegments.begin(); it != path.pathSegments.end(); it++) {
		pathSegments.push_back((*it)->clone());
	}
}

Path::~Path() {
	for(list<PathSegment*>::iterator it = pathSegments.begin(); it != pathSegments.end(); it++) {
		delete *it;
	}
}

double Path::getLength() const {
	return length;
}

PathSegment* Path::getPathSegment(double &s) const {
	list<PathSegment*>::const_iterator it = pathSegments.begin();
	while(it != pathSegments.end() && s >= (*it)->getLength()) {
		s -= (*it)->getLength();
		it++;
	}
	if(it == pathSegments.end()) {
		it--;
	}
	return *it;
}

VectorXd Path::getConfig(double s) const {
	const PathSegment* pathSegment = getPathSegment(s);
	return pathSegment->getConfig(s);
}

VectorXd Path::getConfigDeriv(double s) const {
	const PathSegment* pathSegment = getPathSegment(s);
	return pathSegment->getConfigDeriv(s);
}

VectorXd Path::getConfigDeriv2(double s) const {
	const PathSegment* pathSegment = getPathSegment(s);
	return pathSegment->getConfigDeriv2(s);
}

double Path::getNextSwitchingPoint(double s, bool &discontinuity) const {
	double localPathPos = s;
	const PathSegment* pathSegment = getPathSegment(localPathPos);
	return pathSegment->getNextSwitchingPoint(localPathPos, discontinuity) + (s - localPathPos);
}