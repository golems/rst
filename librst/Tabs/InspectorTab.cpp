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

#include <iostream>
#include <wx/wx.h>

#include "InspectorTab.h"

#include <GUI/Viewer.h>
#include <GUI/GUI.h>
#include <GUI/RSTSlider.h>
#include <GUI/RSTFrame.h>


using namespace std;
using namespace Eigen;

//Give each slider a number so we recognize them
enum sliderNames{
	J_SLIDER = 1000,
	X_SLIDER = 1001,
	Y_SLIDER = 1002,
	Z_SLIDER = 1003,
	ROLL_SLIDER = 1004,
	PITCH_SLIDER = 1005,
	YAW_SLIDER = 1006
};

//Add a handler for slider changes
BEGIN_EVENT_TABLE(InspectorTab, wxPanel)
EVT_COMMAND (wxID_ANY, wxEVT_RST_SLIDER_CHANGE, InspectorTab::OnSlider)
END_EVENT_TABLE()

// Class constructor for the tab: Each tab will be a subclass of wxPanel
IMPLEMENT_DYNAMIC_CLASS(InspectorTab, RSTTab)
InspectorTab::InspectorTab(wxWindow *parent, const wxWindowID id,
                       const wxPoint& pos, const wxSize& size,
                       long style)
					   : RSTTab(parent, id, pos, size, style)
 {
		sizerFull = new wxBoxSizer(wxHORIZONTAL);

		wxStaticBox* jointBox = new wxStaticBox(this,-1,wxT("Item Information"));
		wxStaticBox* posBox = new wxStaticBox(this,-1,wxT("Position"));
		wxStaticBox* rotBox = new wxStaticBox(this,-1,wxT("Orientation"));

		wxStaticBoxSizer* jointBoxS = new wxStaticBoxSizer(jointBox, wxVERTICAL);
		wxStaticBoxSizer* posBoxS = new wxStaticBoxSizer(posBox, wxVERTICAL);
		wxStaticBoxSizer* rotBoxS = new wxStaticBoxSizer(rotBox, wxVERTICAL);

		sizerFull->Add(jointBoxS, 1, wxEXPAND | wxALL, 6);
		sizerFull->Add(posBoxS, 1, wxEXPAND | wxALL, 6);
		sizerFull->Add(rotBoxS, 1, wxEXPAND | wxALL, 6);
		SetSizer(sizerFull);

		// Add the Joint Information and Joint Slider
		itemName = new wxStaticText(this,-1,wxT("Item: (none)"),wxDefaultPosition,wxDefaultSize,wxALIGN_CENTRE);
		parentName = new wxStaticText(this,-1,wxT("Parent: (none)"),wxDefaultPosition,wxDefaultSize,wxALIGN_CENTRE);

		jSlider = new RSTSlider("Joint",-180,180,2000,0,1000,2000,this,J_SLIDER);
		jointBoxS->Add(parentName,1,wxEXPAND | wxALL, 6);
		jointBoxS->Add(itemName,1,wxEXPAND | wxALL, 6);
		jointBoxS->Add(jSlider,1,wxEXPAND | wxALL, 6);


		xSlider = new RSTSlider("X",-10,10,500,0,100,500,this,X_SLIDER);
		ySlider = new RSTSlider("Y",-10,10,500,0,100,500,this,Y_SLIDER);
		zSlider = new RSTSlider("Z",-10,10,500,0,100,500,this,Z_SLIDER);
		posBoxS->Add(xSlider,1,wxEXPAND | wxALL, 6);
		posBoxS->Add(ySlider,1,wxEXPAND | wxALL, 6);
		posBoxS->Add(zSlider,1,wxEXPAND | wxALL, 6);


		rollSlider = new RSTSlider("R",-180,180,500,0,100,500,this,ROLL_SLIDER);
		pitchSlider = new RSTSlider("P",-180,180,500,0,100,500,this,PITCH_SLIDER);
		yawSlider = new RSTSlider("Y",-180,180,500,0,100,500,this,YAW_SLIDER);
		rotBoxS->Add(rollSlider,1,wxEXPAND | wxALL, 6);
		rotBoxS->Add(pitchSlider,1,wxEXPAND | wxALL, 6);
		rotBoxS->Add(yawSlider,1,wxEXPAND | wxALL, 6);
 }

//Handle slider changes
void InspectorTab::OnSlider(wxCommandEvent &evt){
	Robot *r;
	Object* o;
	Link* l;
	Matrix3d rot;
	Vector3d tempTrans;
	int slnum = evt.GetId();
	double pos = *(double*)evt.GetClientData();
	char numBuf[1000];
    numBuf[0] = '\0';
    //sprintf(numBuf,"");

	if(selectedTreeNode==NULL){ return; }

	int selected = selectedTreeNode->dType;
	if(selected == Return_Type_Object){
		o = (Object*)(selectedTreeNode->data);
	}
	else if(selected == Return_Type_Link){
		l = (Link*)(selectedTreeNode->data);
		o = (Object*)(selectedTreeNode->data);
	}
	else if(selected == Return_Type_Robot){
		r = (Robot*)(selectedTreeNode->data);
		o = (Object*)r->baseLink;
	}

	//TODO ask mike if the logic for fromJoints is legit wrt reverseLinkORder
	if(slnum == J_SLIDER && selected == Return_Type_Link){
		if(l->jType == Link::REVOL)
			l->jVal = DEG2RAD(pos);
		else
			l->jVal = pos;
		if(reverseLinkOrder){
			l->updateRelPose();
			l->updateParentPoseRecursive(true, check_for_collisions);
		}else{
			l->updateRelPose();
			l->updateAbsPose();
			//l->updateRecursive(true, check_for_collisions);
			//world->updateRobot(l->robot); // was commented in 2.1 for some reason
		}
		sprintf(numBuf,"Joint Change: %7.4f", pos);
	}else{
		switch(slnum){
			case X_SLIDER:
				o->absPose(0,3) = pos;
				sprintf(numBuf,"X Change: %7.4f", pos);
				break;
			case Y_SLIDER:
				o->absPose(1,3) = pos;
				sprintf(numBuf,"Y Change: %7.4f", pos);
				break;
			case Z_SLIDER:
				o->absPose(2,3) = pos;
				sprintf(numBuf,"Z Change: %7.4f", pos);
				break;
			case ROLL_SLIDER:
			case PITCH_SLIDER:
			case YAW_SLIDER:
				rot = AngleAxisd(DEG2RAD(yawSlider->pos),Vector3d::UnitZ()) * AngleAxisd(DEG2RAD(pitchSlider->pos), Vector3d::UnitY()) * AngleAxisd(DEG2RAD(rollSlider->pos), Vector3d::UnitX());
				tempTrans = o->absPose.translation();
				o->absPose = rot;
				o->absPose.translation() = tempTrans;
				sprintf(numBuf,"Angle Change: %7.4f", pos);
				break;
			default:
				return;
		}
		if(selected == Return_Type_Link){
			l = (Link*)o;
			l->updateParentPoseRecursive(false, check_for_collisions);
		}
		if(selected == Return_Type_Robot){
			l = (Link*)o;
			//l->updateRecursive(false, check_for_collisions);
		}

		if(selected == Return_Type_Object && check_for_collisions){
			world->updateCollisionModel(o);
		}
	}

	world->updateAllCollisionModels(); // added because the collision updates above weren't catching everything
//	static int cdCount = 0;
//	if (cdCount % 5 == 0) { //  && check_for_collisions
	world->clearCollisions();
	if(check_for_collisions){
		world->detectCollisions(); // should be calling this periodically in a thread somewhere, rather than this dumb hack
	}
//	}
//	cdCount++;

	if(frame!=NULL)	frame->SetStatusText(wxString(numBuf,wxConvUTF8));
	viewer->UpdateCamera();
}


// This function is called when an object is selected in the Tree View
// the sliders are shown / hidden and set according to the type of object
// that is selected and its properties
void InspectorTab::RSTStateChange() {
	if(selectedTreeNode==NULL){
		itemName->SetLabel(wxString("Item: None",wxConvUTF8));
		parentName->Hide();
		jSlider->Hide();
		return;
	}

	Object* o;
	Robot* r;
	Link* l;
	string statusBuf;
	string buf,buf2;
	int selected = selectedTreeNode->dType;

	// Get the absPose data from each type of element

	if(selected == Return_Type_Object){
		o = (Object*)(selectedTreeNode->data);
	}
	else if(selected == Return_Type_Link){
		l = (Link*)(selectedTreeNode->data);
		o = (Object*)l;
	}
	else if(selected == Return_Type_Robot){
		r = (Robot*)(selectedTreeNode->data);
		o = (Object*)r->baseLink;
	}

	// Everything can be treated as an object now
	double roll=atan2(o->absPose(2,1), o->absPose(2,2));
	double pitch=-asin(o->absPose(2,0));
	double yaw=atan2(o->absPose(1,0), o->absPose(0,0));

	xSlider->setValue(o->absPose(0,3),false);
	ySlider->setValue(o->absPose(1,3),false);
	zSlider->setValue(o->absPose(2,3),false);

	rollSlider->setValue(RAD2DEG(roll),false);
	pitchSlider->setValue(RAD2DEG(pitch),false);
	yawSlider->setValue(RAD2DEG(yaw),false);


	if(selected == Return_Type_Object){
		statusBuf = " Selected Object: " + o->name;
		buf = "Object: " + o->name;
		itemName->SetLabel(wxString(buf.c_str(),wxConvUTF8));
		parentName->Show();
		parentName->SetLabel(wxString("",wxConvUTF8));
		jSlider->Hide();
	}

	if(selected == Return_Type_Robot){
		statusBuf = " Selected Robot: " + r->name;
		buf = "Robot: " + r->name;
		itemName->SetLabel(wxString(buf.c_str(),wxConvUTF8));
		parentName->Show();
		parentName->SetLabel(wxString("",wxConvUTF8));
		jSlider->Hide();
	}

	if(selected == Return_Type_Link){
		//cout << world->robots[0]->links[4] << " " << world->robots[0]->links[4]->name << " " << world->robots[0]->links[4]->parent << endl;
		//cout << l << " " << l->name << " " << l->parent << endl;

		statusBuf = " Selected Link: " + l->name + " of Robot: " + l->robot->name;
		buf = "Link: " + l->name;
		itemName->SetLabel(wxString(buf.c_str(),wxConvUTF8));
		if(l->parent != NULL){
			buf2 = "Parent Link: " + l->parent->name + "   Robot: " + l->robot->name;
			jSlider->setRange(RAD2DEG(l->jMin),RAD2DEG(l->jMax));
			if(l->jType == Link::REVOL){
				jSlider->setValue(RAD2DEG(l->jVal));
			}else{
				jSlider->setValue(l->jVal);
			}
			jSlider->Show();
		}else{
			buf2 = " (Root Link) ";
			jSlider->Hide();
		}
		parentName->SetLabel(wxString(buf2.c_str(),wxConvUTF8));
		parentName->Show();
	}

	//frame->SetStatusText(wxString(statusBuf.c_str(),wxConvUTF8));
	//sizerFull->Layout();
}
