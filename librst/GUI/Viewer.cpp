//---------------------------------------------------------------------
//  Copyright (c) 2008 Mike Stilman
//  All Rights Reserved.
//
//  Permission to duplicate or use this software in whole or in part
//  is only granted by consultation with the author.
//
//    Mike Stilman              mstilman@cc.gatech.edu
//
//	  Robotics and Intelligent Machines
//    Georgia Tech
//--------------------------------------------------------------------

#include "Viewer.h"
#include "../Tools/World.h"
#include "wx/sizer.h"
#include "GUI.h"
#include <wx/glcanvas.h>

using namespace Eigen;

void Viewer::shown(wxShowEvent& WXUNUSED(evt)){
    int w, h;
    GetClientSize(&w, &h);
	#ifndef __WXMOTIF__
    if (GetContext())
	#endif
    {
        //SetCurrent();
        glViewport(0, 0, (GLint) w, (GLint) h);
    }
	UpdateCamera();
}

void Viewer::UpdateCamera(void){

	SetCurrent();
	camT = camRotT;
	camT.translate(Vector3d(camRadius, 0, 0));
//	for(int i = 0; i < 4; i++) {
//		for(int j = 0; j < 4; j++) {
//			printf("%lf\t", camT(i, j));
//		}
//		printf("\n");
//	}
//	printf("CamRadius  %lf\n", camRadius);
	DrawGLScene();
}

void Viewer::ResetCamera(void) {
	ResetGL();
}

int Viewer::DrawGLScene()
{
	redrawFlag = true;
	glLoadIdentity();
	glPolygonMode (GL_FRONT, GL_FILL);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glHint(GL_FOG_HINT,GL_NICEST);
//	glHint(GL_BLEND_HINT,GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);

	Vector3d zvec(0,0,1);
	upV = camT*zvec;

	gluLookAt(camT(0,3),camT(1,3),camT(2,3),
			  targT(0,3),targT(1,3),targT(2,3),
			  upV[0],upV[1],upV[2]);

//	printf("%lf %lf %lf -> %lf %lf %lf\n", camT(0,3),camT(1,3),camT(2,3),
//			  targT(0,3),targT(1,3),targT(2,3));

	float position[]= {camT(0,3),camT(1,3),camT(2,3), 1.0};
//	float position2[]= {camT.pos.x+10, camT.pos.y, camT.pos.z+10, 1.0};
	glEnable(GL_TEXTURE_2D);
	glColor4f(1.0f,1.0f,1.0f,0.0f);

	glPushMatrix();


	glEnable(GL_LIGHTING);
	glDisable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
	glPushMatrix();
	//float direction[] = {0.0, 0.5, -10.0, 0.0};
	float no_mat[] = {0.0, 0.0, 0.0, 1.0};
	float diffuse[] = {.7, .7, .7, .7};
	float specular[] = {0.5, 0.5, 0.5, 1.0f};
	//float yellowAmbientDiffuse[] = {1.0f, 1.0f, 1.0f, 1.0f};
	glLightfv(GL_LIGHT1, GL_POSITION, position);
	glLightfv(GL_LIGHT1, GL_AMBIENT, no_mat);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);

	GLfloat HeadlightAmb[4] = {0.0f, 0.0f, 0.0f, 1.0f};
	GLfloat HeadlightDif[4] = {1.0f, 1.0f, 1.0f, 0.0f};
	GLfloat HeadlightSpc[4] = {0.0f, 0.0f, 0.0f, 0.0f};
	glLightfv(GL_LIGHT2, GL_AMBIENT, HeadlightAmb);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, HeadlightDif);
	glLightfv(GL_LIGHT2, GL_SPECULAR, HeadlightSpc);
	glLightfv(GL_LIGHT2, GL_POSITION, position);


	glPopMatrix();

	glTranslated(worldT(0, 3),worldT(1, 3),worldT(2, 3));

	float ambRefl[] = {0.2f, 0.2f, 0.2f, 1.0f}; // default
	glMaterialfv(GL_FRONT, GL_AMBIENT, ambRefl);
	float diffRefl[] = {0.8f, 0.8f, 0.8f, 1.0f}; // default
	glMaterialfv(GL_FRONT, GL_DIFFUSE, diffRefl);
	float specRefl[] = {1.0f, 1.0f, 1.0f, 1.0f}; // default
	glMaterialfv(GL_FRONT, GL_SPECULAR, specRefl);

	glPushMatrix();
	if (gridActive){  addGrid(); }
	glPopMatrix();

	//USUALLY BAD
	//glDisable(GL_LIGHTING);
	glDisable(GL_FOG);

	glPushMatrix();
	if(world!=NULL) world->Draw();
	glPopMatrix();

	glPopMatrix();

	glFlush();
	SwapBuffers();
	return TRUE;
}


void Viewer::resized(wxSizeEvent& evt){
	if(!IsShown() || !GetParent()->IsShown()) return;
	wxGLCanvas::OnSize(evt);
    int w, h;
    GetClientSize(&w, &h);
	#ifndef __WXMOTIF__
    if (GetContext())
	#endif
    {
        SetCurrent();
        glViewport(0, 0, (GLint) w, (GLint) h);
    }
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPolygonMode (GL_FRONT, GL_FILL);
	gluPerspective(45.0f,(GLdouble)w/(GLdouble)h,0.1f,100.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	DrawGLScene();
}

void Viewer::mouseWheelMoved(wxMouseEvent& evt){
	SetFocus();
	int wheelRot = evt.GetWheelRotation();

	double camR = camRadius - (double)wheelRot / 500.f;
	if(camR > .05){
		camRadius = camR;
		UpdateCamera();
	}
}

void Viewer::OnCaptureLost(wxMouseCaptureLostEvent& WXUNUSED(evt)){
	mouseCaptured = false;
}

void Viewer::OnMouse(wxMouseEvent& evt){
	evt.GetPosition(&x,&y);
	if(evt.ButtonUp() && mouseCaptured){
		ReleaseMouse();
		mouseCaptured = false;
	}

	if(evt.ButtonDown() && !mouseCaptured){
		prevCamT = camT;
		prevWorldT = worldT;
		xInit = x;
		yInit = y;
		SetFocus();
		CaptureMouse();
		mouseCaptured = true;
		return;
	}
	double dx = (double)(x-xInit);
	double dy = (double)(y-yInit);

	if (evt.LeftIsDown() && evt.RightIsDown()){
		xInit = x;
		yInit = y;
		double camR = camRadius + dy *CAMERASPEED* 2.0f;
		if(camR > .05)
			camRadius = camR;
		UpdateCamera();
	}else if(evt.LeftIsDown()){
		double theta,phi;
		theta = -(dx/90.f);
		phi = -(dy/90.f);
//		Transform<double, 3, Affine> drot,dtrans,dT;
//		Mat33 tm,pm;
//		tm.zRot(theta);
//		pm.yRot(phi);
//		drot.rot = prevCamT.rot;
//		drot.rot = drot.rot*tm;
//		drot.rot = drot.rot*pm;
		camRotT.rotate(AngleAxisd(phi, Vector3d::UnitY())); //.rotate(AngleAxisd(phi, Vector3d::UnitY()));
		//camRotT = camRotT*prevCamT;
		existsUpdate = true;
		UpdateCamera();
	}else if(evt.MiddleIsDown() || evt.RightIsDown()){
//		Vector3d dispV = camRotT*Vector3d(0,(double)dx * CAMERASPEED , -(double)dy * CAMERASPEED);
//		worldT.translate(dispV);
		existsUpdate = false;
		UpdateCamera();
	}
}

void Viewer::render(wxPaintEvent& WXUNUSED(evt)){
    if(!IsShown()) return;
    wxGLCanvas::SetCurrent();
    wxPaintDC(this); // only to be used in paint events. use wxClientDC to paint outside the paint event
	DrawGLScene();
}

void Viewer::InitGL(){
	existsUpdate=false;
	doCollisions=false;
	Move=false;
	pflag=false;
	threadCounter=0;
	gridActive=true;
	redrawFlag = true;
	gridActive = true;
	mouseLDown = false;
	mouseRDown = false;
	mouseMDown = false;
	mouseCaptured = false;

	xInit=0; yInit=0; x=0; y=0;

	redrawCount=0;


    GetClientSize(&w, &h);
	#ifndef __WXMOTIF__
    if (GetContext())
	#endif
    {
		Show();
		SetCurrent();
        glViewport(0, 0, (GLint) w, (GLint) h);
    }

	glShadeModel(GL_SMOOTH);
	//glShadeModel(GL_FLAT);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearDepth(100.0f);
	glEnable(GL_DEPTH_TEST);
	//glDisable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	//glEnable(GL_FOG);
	float FogCol[3]={0.0f,0.0f,0.0f};
	glFogfv(GL_FOG_COLOR,FogCol);
	glFogf(GL_FOG_START, 10.f);
	glFogf(GL_FOG_END, 15.f);
	glFogi(GL_FOG_MODE, GL_LINEAR);


	camT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	worldT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	targT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	camRotT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	prevCamT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	prevWorldT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	camRadius = defaultCamRadius;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPolygonMode (GL_FRONT, GL_FILL);
	gluPerspective(45.0f,(GLdouble)w/(GLdouble)h,0.1f,15.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	UpdateCamera();
}

void Viewer::ResetGL(){
    GetClientSize(&w, &h);
	#ifndef __WXMOTIF__
    if (GetContext())
	#endif
    {
		SetCurrent();
        glViewport(0, 0, (GLint) w, (GLint) h);
    }

	camT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	worldT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	targT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	camRotT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	prevCamT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	prevWorldT = AngleAxis<double> (0, Vector3d(0, 0, 0));
	camRadius = defaultCamRadius;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPolygonMode (GL_FRONT, GL_FILL);
	gluPerspective(45.0f,(GLdouble)w/(GLdouble)h,0.1f,15.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	UpdateCamera();
}

void Viewer::setClearColor(double r,double g, double b, double a){
	glClearColor(r,g,b,a);
	float FogCol[3]={r,g,b};
	glFogfv(GL_FOG_COLOR,FogCol);
}

void Viewer::addGrid(){
	glEnable(GL_FOG);
	double sizeX=100.0f;
	double sizeY=100.0f;
	double inc=.5f;

	double startX=-sizeX/2.f;
	double endX  =sizeX/2.f;
	double startY=-sizeY/2.f;
	double endY  =sizeY/2.f;

	glLineWidth(0.9f);
	glEnable(GL_COLOR_MATERIAL);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(.2f,.2f,0.0f);
	for(double i=startX; i<=endX; i+=inc){ glVertex3d( i, startY, 0.0f); glVertex3d( i, endY, 0.0f);	}
	for(double i=startY; i<=endY; i+=inc){ glVertex3d( startX, i, 0.0f); glVertex3d( endX, i, 0.0f);	}
	glEnd();
	glColor3f(1.0f,1.0f,1.0f);
	glEnable(GL_LIGHTING);
	glDisable(GL_FOG);
}

BEGIN_EVENT_TABLE(Viewer, wxGLCanvas)
	EVT_SHOW(Viewer::shown)
    EVT_SIZE(Viewer::resized)
    EVT_PAINT(Viewer::render)
	EVT_MOUSE_CAPTURE_LOST(Viewer::OnCaptureLost)
	EVT_MOUSEWHEEL(Viewer::mouseWheelMoved)
    EVT_MOTION(Viewer::OnMouse)
	EVT_RIGHT_DOWN(Viewer::OnMouse)
	EVT_MIDDLE_DOWN(Viewer::OnMouse)
	EVT_LEFT_DOWN(Viewer::OnMouse)
	EVT_LEFT_UP(Viewer::OnMouse)
	EVT_RIGHT_UP(Viewer::OnMouse)
	EVT_MIDDLE_UP(Viewer::OnMouse)
END_EVENT_TABLE()
