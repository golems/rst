/*
 * Copyright (c) 2010, Georgia Tech Research Corporation
 * 
 * Humanoid Robotics Lab      Georgia Institute of Technology
 * Director: Mike Stilman     http://www.golems.org
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *     * Redistributions of source code must retain the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     * Neither the name of the Georgia Tech Research Corporation nor
 *       the names of its contributors may be used to endorse or
 *       promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GEORGIA TECH RESEARCH CORPORATION ''AS
 * IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GEORGIA
 * TECH RESEARCH CORPORATION BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef LINK_H
#define LINK_H

#include <Eigen/Geometry>
#include <Tools/Object.h>
#include <Tools/Constants.h>
#include <vector>
class Robot;

class Link: public Object {
public:

	Link()
		: attachedObject(NULL) {
	}

	Eigen::Transform<double, 3, Eigen::Affine> pose;
	Eigen::Transform<double, 3, Eigen::Affine> jTrans;
	Eigen::Vector3d jAxis;

	Robot *robot;
	Link *parent;
	std::vector<Link*> children;
	Object* attachedObject;
	Eigen::Transform<double, 3, Eigen::Affine> attachedObjectPose;

	enum JointType {
		REVOL, PRISM, FIXED, FREE
	};

	int index;
	JointType jType;
	double jVal;
	double jMin, jMax;
        //-- Added on Fri, Sept 9th
        double accMax;
        double velMax;
        //--

	void updateRelPose();
	void updateAbsPose(bool updateCollisionModels = true);

	void updateParentPose();
	void updateParentPoseRecursive(bool fromJoints = false, bool collisions = false);

	void recursiveSetAncestry(Robot*, Link*);
	void updateRecursiveCOM();

	bool attachObject(Object* object);
	void releaseObject();
};

#endif
