/*
 * @file example.cpp
 * @Brief examples simulations of homogenization with VegaFEM
 * @author Fei Zhu, Jing Zhao
 * 
 * Copyright (C) 2015 Fei Zhu.
 *
 * This Source Code Form is subject to the terms of the GNU General Public License v2.0. 
 * If a copy of the GPL was not distributed with this file, you can obtain one at:
 * http://www.gnu.org/licenses/gpl-2.0.html
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <ctime>
#include <math.h>
using namespace std;
using std::fstream;
using std::stringstream;
using std::vector;
using std::string;

#ifdef WIN32
#include <windows.h>
#endif

#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

//VegaFEM includes
#include "VegaHeaders.h"

// graphics 
const int string_length=4096;//sizeof(GLUI_String);
char windowTitleBase[string_length] = "Real-time simulation";
int windowID;
int windowWidth = 800;
int windowHeight = 600;
double zNear=0.01;
double zFar=10.0;
double cameraRadius;
double focusPositionX, focusPositionY, focusPositionZ;
double cameraLongitude, cameraLattitude;
SphericalCamera * camera = NULL;
int g_iLeftMouseButton=0, g_iMiddleMouseButton=0, g_iRightMouseButton=0;
int g_vMousePos[2] = {0,0};
int shiftPressed=0;
int altPressed=0;
int ctrlPressed=0;
int renderAxes=0;
int renderWireframe=0;
int renderVertices=0;
int renderVelocity = 0;
int renderDeformableObject=1;
int renderMaterialColors = 0;
int useRealTimeNormals = 0;
int renderGroundPlane = 0;
int renderFixedVertices = 1;
int renderSprings = 0;
int lockScene=0;
int pauseSimulation=1;
Lighting * lighting = NULL;
SceneObjectDeformable * volumetricSurfaceMesh = NULL;
SceneObjectDeformable * renderingMesh = NULL;
SceneObjectDeformable ** objectRenderSurfaceMesh = NULL;
SceneObjectDeformable **extraObjectsInScene=NULL;
SceneObject * extraSceneGeometry = NULL;
int renderVolumetricSurface=1;
char groundPlaneString[128];
double groundPlaneHeight;
double groundPlaneLightHeight = 10.0;
double groundPlaneSize = 15.0;
GLuint displayListGround;
Planes *planesInScene=NULL;
int planeNumber=0;
//configFile
string configFilename;
char renderingMeshFilename[string_length];
char volumetricSurfaceMeshFilename[string_length];
char ** objectRenderSurfaceMeshFilename;
char ** objectRenderSurfaceMeshInterpolationFilename;
char extraObjectsFileNameBase[string_length];
int extraObjectsNum=0;

//interpolation to embedded object render surface mesh
int objectRenderSurfaceMeshFileNum=1;
int * objectRenderSurfaceMeshInterpolationElementVerticesNum;
int ** objectRenderSurfaceMeshInterpolationVertices=NULL;
double ** objectRenderSurfaceMeshInterpolationWeights=NULL;



char deformableObjectMethod[string_length];
char fixedVerticesFilename[string_length];
char massMatrixFilename[string_length];
char invertibleMaterialString[string_length] = "__none";
char initialPositionFilename[string_length];
char initialVelocityFilename[string_length];
char forceLoadsFilename[string_length];
char elasticTensorFilename[string_length];
//outputFile
char outputFilename[string_length];
char outputFileNameBase[string_length];
int timestepPerOutputFile=20;
int outputFileIndex=1;//index of current output file
bool saveMeshToFile=false;

char planesFilename[string_length]="__none";
int corotationalLinearFEM_warp = 1;
const int max_corotationalLinearFEM_warp = 2;
char implicitSolverMethod[string_length];
char solverMethod[string_length];
char extraSceneGeometryFilename[string_length];
char lightingConfigFilename[string_length];
//set dampingMassCoef=dampingStiffnessCoef=0.0 if the material is anistropic material
float dampingMassCoef=0.0;
float dampingStiffnessCoef=0.0;
float dampingLaplacianCoef = 0.0;
float deformableObjectCompliance = 1.0;
float baseFrequency = 1.0;
int maxIterations;
double epsilon;
char backgroundColorString[string_length] = "220 220 220";
int numInternalForceThreads;
int numSolverThreads;

// simulation
int syncTimestepWithGraphics=1;
float timeStep = 1.0;
int timeStepCounter=0;
int totalSteps=0;
double frame_rate=1.0;
int totalFrame = 0;
int simulation_frame_count=0;
double total_simulation_time=0.0;
double total_save_file_time=0.0;

float newmarkBeta = 0.25;
float newmarkGamma = 0.5;
int stepsPerFrame = 1;
double principalStretchThreshold;
int enableTextures = 0;
int staticSolver = 0;
int pulledVertex = -1;
int forceNeighborhoodSize = 5;
int dragStartX, dragStartY;

int test_case=0;
int numFixedVertices;
int * fixedVertices;
int numForceLoads = 0;
double * forceLoads = NULL;
IntegratorBase * integratorBase = NULL;
ImplicitNewmarkSparse * implicitNewmarkSparse = NULL;
IntegratorBaseSparse * integratorBaseSparse = NULL;
int centralDifferencesTangentialDampingUpdateMode = 1;
int positiveDefinite = 0;

ForceModel * forceModel = NULL;
StVKInternalForces * stVKInternalForces = NULL;
StVKStiffnessMatrix * stVKStiffnessMatrix = NULL;
StVKForceModel * stVKForceModel = NULL;
AnisotropicInternalForces * anisotropicInternalForces = NULL;
AnisotropicStiffnessMatrix * anisotropicStiffnessMatrix = NULL;
AnisotropicForceModel * anisotropicForceModel = NULL;
CorotationalAnisotropicFEMForceModel * corotationalAnisotropicFEMForceModel = NULL;
CorotationalLinearFEMForceModel * corotationalLinearFEMForceModel = NULL;
bool addGravity=false;
double g=0;
unsigned int output_mode=1;
VolumetricMesh * volumetricMesh = NULL;
TetMesh * tetMesh = NULL;
Graph * meshGraph = NULL;
enum deformableObjectType { STVK, COROTLINFEM, LINFEM, INVERTIBLEFEM, ANISOTROPIC, COROTANISTROPIC,UNSPECIFIED} deformableObject = UNSPECIFIED;
enum invertibleMaterialType { INV_STVK, INV_NEOHOOKEAN, INV_MOONEYRIVLIN, INV_NONE } invertibleMaterial = INV_NONE;
enum solverType { IMPLICITNEWMARK, IMPLICITBACKWARDEULER, EULER, SYMPLECTICEULER, CENTRALDIFFERENCES, UNKNOWN } solver = UNKNOWN;

SparseMatrix * massMatrix = NULL;
SparseMatrix * LaplacianDampingMatrix = NULL;
int simulation_vertice_num;
double **uRenderSurface=NULL;//displacement of the object render surface mesh
double * u = NULL;
double * f_ext = NULL;
double * f_col=NULL;//collision force
double * f_extBase = NULL;
double * uInitial = NULL;
double * velInitial = NULL;
double * gravity_force = NULL;
int * fixedDOFs;
// glui
GLUI * glui;
GLUI_Spinner * timeStep_spinner;
GLUI_StaticText * systemSolveStaticText;
GLUI_StaticText * forceAssemblyStaticText;

//test_case_1:bar twist
double * initial_pos;
//function declaration
void addGravitySwitch(bool addGravity);
void displayFunction(void);
void saveCurrentObjectSurface(int code);
void changeSimulationMode(int code);
void initFunction(int test_case_);
void simulationFunction(int test_case_);
void testStiffnessMatrix(void);
//font is, for example, GLUT_BITMAP_9_BY_15
void print_bitmap_string(float x, float y, float z, void * font, char * s)
{
	glRasterPos3f(x,y,z);
	if (s && strlen(s)) 
	{
		while (*s) 
		{
			glutBitmapCharacter(font, *s);
			s++;
		}
	}
}

//save current object render surface  to an obj file
void saveCurrentObjectSurface(int code)
{
	PerformanceCounter save_file_counter;
	double each_frame_save_file=0.0;
	save_file_counter.StartCounter();
	stringstream adaptor;
	string outputFileName[objectRenderSurfaceMeshFileNum];
	string outputFileIndexStr;
	string outputObjectFileName,object_idx_str;
	adaptor<<outputFileIndex++;
	adaptor>>outputFileIndexStr;
	
	for(int i=0;i<objectRenderSurfaceMeshFileNum;++i)
	{
		adaptor.str("");
		adaptor.clear();
		adaptor<<i;
		adaptor>>object_idx_str;
		outputObjectFileName=object_idx_str;
		outputFileName[i]=outputFileNameBase+object_idx_str+"_"+outputFileIndexStr+".obj";
		ObjMesh *mesh=objectRenderSurfaceMesh[i]->GetMesh();
		mesh->save(outputFileName[i],0);
		cout<<outputFileName[i]<<" saved.\n";
	}
	save_file_counter.StopCounter();
	each_frame_save_file=save_file_counter.GetElapsedTime();
	total_save_file_time+=each_frame_save_file;
	std::cout<<"Current frame save file time is "<<each_frame_save_file<<" s.\n";
	if(timeStepCounter==(totalSteps-1))
	{
		std::cout<<"Total save files time is "<<total_save_file_time<<" s; Average: "<<total_save_file_time/totalFrame<<" s/frame.\n";
	}
}

//***********************************************graphics loop function***************************************
// graphics loop function.
void displayFunction(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	// setup model transformations
	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity();

	camera->Look();
	// set OpenGL lighting 
	renderingMesh->SetLighting(lighting);

	glEnable(GL_LIGHTING);

	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

	glStencilFunc(GL_ALWAYS, 0, ~(0u));

	if(extraObjectsNum>0)//render the extra objects in scene
	{
		for(int i=0;i<extraObjectsNum;++i)
			extraObjectsInScene[i]->Render();
	}

	glStencilFunc(GL_ALWAYS,1,~(0u));
	//render the exterior surface of volumetric mesh
	//only when the embedded render mesh is rendered shall we render the volumetric surface
	if(renderVolumetricSurface)
	{
		//render transparent volumetric surface
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0,1.0);
		glDrawBuffer(GL_NONE);
		volumetricSurfaceMesh->Render();
		glDisable(GL_POLYGON_OFFSET_FILL);
		glDrawBuffer(GL_BACK);
		glEnable(GL_LIGHTING);
		glColor3f(1.0,0.0,0.0);
		volumetricSurfaceMesh->Render();
		if(renderVertices)//render vertices of the volumetric surface
		{
			glDisable(GL_LIGHTING);
			glColor3f(0.5,0,0);
			glPointSize(8.0);
			volumetricSurfaceMesh->RenderVertices();
			glEnable(GL_LIGHTING);
		}
		if(renderVelocity)//render velocity of vertices
		{
			glDisable(GL_LIGHTING);
			for(unsigned int i=0;i<simulation_vertice_num;++i)
			{
				Vec3d vert_pos;
				Vec3d vert_new_pos;
				for(unsigned int j=0;j<3;++j)
				{
					if(timeStepCounter>=(simulation_frame_count/frame_rate)/timeStep)
					{
						vert_pos[j]=(*volumetricMesh->getVertex(i))[j]+integratorBase->Getq()[3*i+j];
						vert_new_pos[j]=vert_pos[j]+integratorBase->Getqvel()[3*i+j];
					}
					else
					{
						vert_pos[j]=(*volumetricMesh->getVertex(i))[j]+u[3*i+j];
						vert_new_pos[j]=vert_pos[j]+velInitial[3*i+j];
					}
				}
				glColor3f(1.0,0.3,0);
				glLineWidth(1);
				glBegin(GL_LINES);
				glVertex3f(vert_pos[0],vert_pos[1],vert_pos[2]);
				glVertex3f(vert_new_pos[0],vert_new_pos[1],vert_new_pos[2]);
				glEnd();
			}
			glEnable(GL_LIGHTING);
		}
		if(renderWireframe)//render wireframe of the volumetric surface
		{
			glDisable(GL_LIGHTING);
			glColor3f(0.1,0.1,0.1);
			renderingMesh->RenderEdges();
			glEnable(GL_LIGHTING);
		}
		glDisable(GL_BLEND);
	}

	glStencilFunc(GL_ALWAYS,0,~(0u));
	glDisable(GL_TEXTURE_2D);
	//render planes
	if(planeNumber>0)
	{
		planesInScene->render();
	}

	glDisable(GL_LIGHTING);
	glStencilFunc(GL_ALWAYS,0,~(0u));
	glColor3f(1.0,0.0,0.0);
	glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
	//render axis
	if(renderAxes)
	{
		glLineWidth(1.0);
		drawAxes(1.0);
	}
	// render the currently pulled vertex
	if (pulledVertex >= 0)
	{
		glColor3f(0,1,0);
		double pulledVertexPos[3];
		volumetricSurfaceMesh->GetSingleVertexPositionFromBuffer(pulledVertex,
			&pulledVertexPos[0], &pulledVertexPos[1], &pulledVertexPos[2]);
		glEnable(GL_POLYGON_OFFSET_POINT);
		glPolygonOffset(-1.0,-1.0);
		glPointSize(8.0);
		glBegin(GL_POINTS);
		glVertex3f(pulledVertexPos[0], pulledVertexPos[1], pulledVertexPos[2]);
		glEnd();
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	// render model fixed vertices
	if (renderFixedVertices)
	{
		for(int i=0; i<numFixedVertices; i++)
		{
			glColor3f(1,0,0);
			double fixedVertexPos[3];
			volumetricSurfaceMesh->GetSingleVertexRestPosition(fixedVertices[i],
				&fixedVertexPos[0], &fixedVertexPos[1], &fixedVertexPos[2]);

			glEnable(GL_POLYGON_OFFSET_POINT);
			glPolygonOffset(-1.0,-1.0);
			glPointSize(8.0);
			glBegin(GL_POINTS);
			glVertex3f(fixedVertexPos[0], fixedVertexPos[1], fixedVertexPos[2]);
			glEnd();
			glDisable(GL_POLYGON_OFFSET_FILL);
		}
	}
	glStencilFunc(GL_ALWAYS,0,~(0u));
	if(renderMaterialColors)
	{
		int materialNum=volumetricMesh->getNumMaterials();
		int regionNum=volumetricMesh->getNumRegions();
		int setNum=volumetricMesh->getNumSets();
		for(int i = 1; i < setNum-1; ++i)
		{
			VolumetricMesh::Set * elementSet=volumetricMesh->getSet(i);
			std::set<int> elements;
			elementSet->getElements(elements);
			glColor3f(fabs(setNum-i*2.0)/setNum,(setNum-i*1.0)/setNum,i/setNum);
			for(set<int>::iterator iter=elements.begin();iter!=elements.end();++iter)
			{
				for(int dim=0;dim<volumetricMesh->getNumElementVertices();++dim)
				{
					int vert=volumetricMesh->getVertexIndex(*iter-1,dim);
					double materialVertexPos[3];
					volumetricSurfaceMesh->GetSingleVertexPositionFromBuffer(vert,
						&materialVertexPos[0], &materialVertexPos[1], &materialVertexPos[2]);
					glEnable(GL_POLYGON_OFFSET_POINT);
					glPolygonOffset(-1.0,-1.0);
					glPointSize(8.0);
					glBegin(GL_POINTS);
					glVertex3f(materialVertexPos[0], materialVertexPos[1], materialVertexPos[2]);
					glEnd();
					glDisable(GL_POLYGON_OFFSET_FILL);
				}	
			}
		}
	}
	glutSwapBuffers();
}
/***********************************************************************************************/
//output files loop function
void outputFilesLoop(void)
{
	// apply any scripted force loads
	PerformanceCounter each_frame_performance_counter;
	double each_frame_time=0.0;
	each_frame_performance_counter.StartCounter();
	while(timeStepCounter<totalSteps)
	{
		memcpy(f_ext, f_extBase, sizeof(double) * 3 * simulation_vertice_num);
		if (timeStepCounter < numForceLoads)
		{
			printf("  External forces read from the binary input file.\n");
			for(int i=0; i<3*simulation_vertice_num; i++)
				f_ext[i] += forceLoads[ELT(3*simulation_vertice_num, i, timeStepCounter)];
		}
		//apply the collision management with planes in scene
		if(planeNumber>0)
		{
			//planesInScene->resolveContact(volumetricSurfaceMesh->GetMesh(),f_col);

			velInitial=integratorBase->Getqvel();
			planesInScene->resolveContact(volumetricSurfaceMesh->GetMesh(),velInitial);
			for(unsigned int i=0;i<simulation_vertice_num;++i)
			{				
				u[3*i]=volumetricSurfaceMesh->GetMesh()->getPosition(i)[0]-(*volumetricMesh->getVertex(i))[0];
				u[3*i+1]=volumetricSurfaceMesh->GetMesh()->getPosition(i)[1]-(*volumetricMesh->getVertex(i))[1];
				u[3*i+2]=volumetricSurfaceMesh->GetMesh()->getPosition(i)[2]-(*volumetricMesh->getVertex(i))[2];
			}
			integratorBase->SetState(u,velInitial);
		}
		if(deformableObject == COROTLINFEM)
		{
			
			for(int i=0;i<3*simulation_vertice_num;++i)
			{
				f_ext[i]+=gravity_force[i];
			}
		}
		integratorBaseSparse->SetExternalForces(f_ext);
		simulationFunction(test_case);
		int code = integratorBase->DoTimestep();
		cout<<"timeStep "<<timeStepCounter<<" begins \n";

		memcpy(u, integratorBase->Getq(), sizeof(double) * 3 * simulation_vertice_num);	
		volumetricSurfaceMesh->SetVertexDeformations(u);	

		for(int i=0;i<objectRenderSurfaceMeshFileNum;++i)
		{
			VolumetricMesh::interpolate(u,uRenderSurface[i],objectRenderSurfaceMesh[i]->Getn(),objectRenderSurfaceMeshInterpolationElementVerticesNum[i],
				objectRenderSurfaceMeshInterpolationVertices[i],objectRenderSurfaceMeshInterpolationWeights[i]);
			objectRenderSurfaceMesh[i]->SetVertexDeformations(uRenderSurface[i]);	
		}		
		volumetricSurfaceMesh->BuildNormals();
		renderingMesh->BuildNormals();
		
		//save object surface mesh to files			
		if(((timeStepCounter+1)%(int)(1.0/(frame_rate*timeStep))==0)&&(timeStepCounter>0))
		{
			cout<<"frame "<<(timeStepCounter+1)%(int)(1.0/(frame_rate*timeStep))<<" begins \n";
			each_frame_performance_counter.StopCounter();
			each_frame_time=each_frame_performance_counter.GetElapsedTime();
			std::cout<<"Current frame simulation time is "<<each_frame_time<<" s.\n";
			total_simulation_time+=each_frame_time;
			saveCurrentObjectSurface(0);
			if(timeStepCounter==(totalSteps-1))
			{
				std::cout<<"Total simulation time is "<<total_simulation_time<<" s; Average: "<<total_simulation_time/totalFrame<<" s/frame.\n";
			}	
			each_frame_performance_counter.StartCounter();
		}
		timeStepCounter++;
	}
	each_frame_performance_counter.StopCounter();
}
// called periodically by GLUT:
void idleFunction(void)
{
	glutSetWindow(windowID);
	// reset external forces (usually to zero)
	memcpy(f_ext, f_extBase, sizeof(double) * 3 * simulation_vertice_num);
	if ((!lockScene) && (!pauseSimulation))
	{
		// determine force in case user is pulling on a vertex
		if (g_iLeftMouseButton) 
		{
			if (pulledVertex != -1)
			{
				double forceX = (g_vMousePos[0] - dragStartX);
				double forceY = -(g_vMousePos[1] - dragStartY);
				double externalForce[3];
				camera->CameraVector2WorldVector_OrientationOnly3D(forceX, forceY, 0, externalForce);
				for(int j=0; j<3; j++)
				{
					externalForce[j] *= deformableObjectCompliance;
					cout<<externalForce[j]<<",deformableObjectCompliance"<<deformableObjectCompliance<<"----";
				}
				cout<<pulledVertex<<" fx: "<<forceX<<" fy: "<<forceY<<" | "<<externalForce[0]<<" "<<externalForce[1]<<" "<<externalForce[2]<<"\n";
				// register force on the pulled vertex
				f_ext[3*pulledVertex+0] += externalForce[0];
				f_ext[3*pulledVertex+1] += externalForce[1];
				f_ext[3*pulledVertex+2] += externalForce[2];
				// distribute force over the neighboring vertices
				set<int> affectedVertices;
				set<int> lastLayerVertices;
				affectedVertices.insert(pulledVertex);
				lastLayerVertices.insert(pulledVertex);
				for(int j=1; j<forceNeighborhoodSize; j++)
				{
					// linear kernel
					double forceMagnitude = 1.0 * (forceNeighborhoodSize - j) / forceNeighborhoodSize;
					set<int> newAffectedVertices;
					for(set<int> :: iterator iter = lastLayerVertices.begin(); iter != lastLayerVertices.end(); iter++)
					{
						// traverse all neighbors and check if they were already previously inserted
						int vtx = *iter;
						int deg = meshGraph->GetNumNeighbors(vtx);
						for(int k=0; k<deg; k++)
						{
							int vtxNeighbor = meshGraph->GetNeighbor(vtx, k);
							if (affectedVertices.find(vtxNeighbor) == affectedVertices.end())
							{
								// discovered new vertex
								newAffectedVertices.insert(vtxNeighbor);
							}
						}
					}
					lastLayerVertices.clear();
					for(set<int> :: iterator iter = newAffectedVertices.begin(); iter != newAffectedVertices.end(); iter++)
					{
						// apply force
						f_ext[3* *iter + 0] += forceMagnitude * externalForce[0];
						f_ext[3* *iter + 1] += forceMagnitude * externalForce[1];
						f_ext[3* *iter + 2] += forceMagnitude * externalForce[2];
						// generate new layers
						lastLayerVertices.insert(*iter);
						affectedVertices.insert(*iter);
					}
				}
			}
		}
		// apply any scripted force loads
		if (timeStepCounter < numForceLoads)
		{
			printf("  External forces read from the binary input file.\n");
			for(int i=0; i<3*simulation_vertice_num; i++)
				f_ext[i] += forceLoads[ELT(3*simulation_vertice_num, i, timeStepCounter)];
		}
		//apply the penalty collision forces with planes in scene
		if(planeNumber>0)
		{
			/*planesInScene->resolveContact(volumetricSurfaceMesh->GetMesh(),f_col);*/
			//for(int i=0;i<3*simulation_vertice_num;++i)
			//{
			////	f_ext[i]+=f_col[i];
			//	
			//}	
			velInitial=integratorBase->Getqvel();
			planesInScene->resolveContact(volumetricSurfaceMesh->GetMesh(),velInitial);
			for(unsigned int i=0;i<simulation_vertice_num;++i)
			{				
				u[3*i]=volumetricSurfaceMesh->GetMesh()->getPosition(i)[0]-(*volumetricMesh->getVertex(i))[0];
				u[3*i+1]=volumetricSurfaceMesh->GetMesh()->getPosition(i)[1]-(*volumetricMesh->getVertex(i))[1];
				u[3*i+2]=volumetricSurfaceMesh->GetMesh()->getPosition(i)[2]-(*volumetricMesh->getVertex(i))[2];
			}
			integratorBase->SetState(u,velInitial);
		}

		if(deformableObject == COROTLINFEM)
		{
			for(int i=0;i<3*simulation_vertice_num;++i)
			{
				f_ext[i]+=gravity_force[i];
			}
		}
		// set forces to the integrator
		integratorBaseSparse->SetExternalForces(f_ext);
		if(timeStepCounter < totalSteps)
		{
			//simulationFunction(test_case);
			int code = integratorBase->DoTimestep();
			printf("."); 
			timeStepCounter++;
		}
	}
	memcpy(u, integratorBase->Getq(), sizeof(double) * 3 * simulation_vertice_num);
	volumetricSurfaceMesh->SetVertexDeformations(u);	

	for(int i=0;i<objectRenderSurfaceMeshFileNum;++i)
	{
		VolumetricMesh::interpolate(u,uRenderSurface[i],objectRenderSurfaceMesh[i]->Getn(),objectRenderSurfaceMeshInterpolationElementVerticesNum[i],
			objectRenderSurfaceMeshInterpolationVertices[i],objectRenderSurfaceMeshInterpolationWeights[i]);
		objectRenderSurfaceMesh[i]->SetVertexDeformations(uRenderSurface[i]);	
	}
	//save object surface mesh to files
	if((!lockScene)&&(!pauseSimulation)&&saveMeshToFile&&(timeStepCounter%timestepPerOutputFile==0))
		saveCurrentObjectSurface(0);
	if (useRealTimeNormals)
	{
		// recompute normals
		volumetricSurfaceMesh->BuildNormals();
		renderingMesh->BuildNormals();
	}
	glutPostRedisplay();
}

// reacts to pressed keys
void keyboardFunction(unsigned char key, int x, int y)
{
	double cameraX,cameraY,cameraZ;

	switch (key)
	{
	case 27:
		exit(0);
	case 32://space button, pause the simulation
		pauseSimulation=!pauseSimulation;
		break;
	case 'a'://render axis
		renderAxes=!renderAxes;
		break;
	case 'b'://render fixed vertices
		renderFixedVertices=!renderFixedVertices;
		break;
	case 'f'://switch on/off output file
		saveMeshToFile=!saveMeshToFile;
		break;
	case 'g'://switch on/off gravity
		addGravity=!addGravity;
		addGravitySwitch(addGravity);
		break;
		// in this mode, can move the camera while the object's deformations are frozen
	case 'l':
		lockScene = !lockScene;
		if (lockScene)
		{
			camera->PushPosition();
		}
		else
		{
			camera->PopPosition();
		}
		break;
	case 'r'://reset camera
		camera->Reset();
		break;
	case 's'://render the exterior surface of the volumetric mesh
		renderVolumetricSurface=!renderVolumetricSurface;
		break;
	case 'v':
		renderVertices = !renderVertices;
		break;
	case 'V':
		renderVelocity = !renderVelocity;
		break;
	case 'm':
		renderMaterialColors = !renderMaterialColors;
		break;
	case 'n'://use real time normals when rendering
		useRealTimeNormals=!useRealTimeNormals;
		break;
	case 'w':
		renderWireframe = !renderWireframe;
		break;
	}
}

// reacts to pressed "special" keys
void specialFunction(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_LEFT:
		break;
	case GLUT_KEY_RIGHT:
		break;
	case GLUT_KEY_DOWN:
		break;
	case GLUT_KEY_UP:
		break;
	case GLUT_KEY_PAGE_UP:
		break;
	case GLUT_KEY_PAGE_DOWN:
		break;
	case GLUT_KEY_HOME:
		break;
	case GLUT_KEY_END:
		break;
	case GLUT_KEY_INSERT:
		break;
	default:
		break;
	}
}

void reshape(int x, int y)
{
	glViewport(0,0,x,y);
	windowWidth = x;
	windowHeight = y;
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	// compute the window aspect ratio 
	gluPerspective(45.0f, 1.0 * windowWidth / windowHeight, zNear, zFar);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

// reacts to mouse motion
void mouseMotionFunction(int x, int y)
{
	int mouseDeltaX = x-g_vMousePos[0];
	int mouseDeltaY = y-g_vMousePos[1];
	g_vMousePos[0] = x;
	g_vMousePos[1] = y;
	if (g_iLeftMouseButton) // left mouse button 
	{
	}
	if (g_iRightMouseButton) // right mouse button handles camera rotations
	{
		const double factor = 0.2;
		camera->MoveRight(factor * mouseDeltaX);
		camera->MoveUp(factor * mouseDeltaY);
	}
	if ((g_iMiddleMouseButton) || (g_iLeftMouseButton && altPressed)) // handle zoom in/out
	{
		const double factor = 0.1;
		camera->ZoomIn(cameraRadius * factor * mouseDeltaY);
	}
}

// reacts to pressed mouse buttons
void mouseButtonActivityFunction(int button, int state, int x, int y)
{
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		g_iLeftMouseButton = (state==GLUT_DOWN);
		shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
		altPressed = (glutGetModifiers() == GLUT_ACTIVE_ALT);
		ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
		if ((g_iLeftMouseButton) && (!shiftPressed) && (!ctrlPressed))//user pulled vertex, apply force
		{
			// apply force to vertex
			GLdouble model[16];
			glGetDoublev (GL_MODELVIEW_MATRIX, model);
			GLdouble proj[16];
			glGetDoublev (GL_PROJECTION_MATRIX, proj);
			GLint view[4];
			glGetIntegerv (GL_VIEWPORT, view);
			int winX = x;
			int winY = view[3]-1-y;
			float zValue;
			glReadPixels(winX,winY,1,1, GL_DEPTH_COMPONENT, GL_FLOAT, &zValue); 
			GLubyte stencilValue;
			glReadPixels(winX, winY, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_BYTE, &stencilValue);
			GLdouble worldX, worldY, worldZ;
			gluUnProject(winX, winY, zValue, model, proj, view, &worldX, &worldY, &worldZ);
			if (stencilValue == 1)
			{
				dragStartX = x;
				dragStartY = y;
				Vec3d pos(worldX, worldY, worldZ);
				pulledVertex = volumetricSurfaceMesh->GetClosestVertex(pos);
				printf("Clicked on vertex: %d (0-indexed)\n", pulledVertex);
			}
			else
			{
				printf("Clicked on empty stencil: %d.\n", stencilValue);
			}
		}
		if (!g_iLeftMouseButton)
		{
			pulledVertex = -1;
		}
		break;
	case GLUT_MIDDLE_BUTTON:
		g_iMiddleMouseButton = (state==GLUT_DOWN);
		shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
		ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
		break;
	case GLUT_RIGHT_BUTTON:
		g_iRightMouseButton = (state==GLUT_DOWN);
		shiftPressed = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
		ctrlPressed = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
		break;
	}
	g_vMousePos[0] = x;
	g_vMousePos[1] = y;
}

//************************************************************************************************************
//switch on/off the gravity
void addGravitySwitch(bool addGravity)
{
	if(deformableObject==STVK||deformableObject==LINFEM)
		stVKInternalForces->SetGravity(addGravity);
	/*if(deformableObject==INVERTIBLEFEM)
		isotropicHyperelasticFEM->SetGravity(addGravity);*/
	if(deformableObject == ANISOTROPIC)
		anisotropicInternalForces->SetGravity(addGravity);
	if(addGravity)
		cout<<"Gravity switched on.\n";
	else
		cout<<"Gravity switched off.\n";
}

//call back function when changing which mesh to render
void updateRenderingMesh(int code)
{
	renderingMesh=volumetricSurfaceMesh/*objectRenderSurfaceMesh0*/;
	if(enableTextures)
		renderingMesh->SetUpTextures(SceneObject::MODULATE,SceneObject::NOMIPMAP);  
	renderingMesh->ResetDeformationToRest();
	renderingMesh->BuildNeighboringStructure();
	renderingMesh->BuildNormals();
	renderingMesh->resetRenderMode();
}
// program initialization
void initSimulation()
{
	 totalSteps=(int)((1.0/timeStep)/frame_rate)*totalFrame;
	// init lighting
	try
	{
		lighting = new Lighting(lightingConfigFilename);
	}
	catch(int exceptionCode)
	{
		printf("Error (%d) reading lighting information from %s .\n",exceptionCode, lightingConfigFilename);
		exit(1);
	}
	// init camera
	delete(camera);
	double virtualToPhysicalPositionFactor = 1.0;
	initCamera(cameraRadius, cameraLongitude, cameraLattitude,
		focusPositionX, focusPositionY, focusPositionZ,
		1.0 / virtualToPhysicalPositionFactor,
		&zNear, &zFar, &camera);
	//init planes in the scene
	if(planeNumber>0)
	{
		if(strcmp(planesFilename,"__none")!=0)
			planesInScene=new Planes(planesFilename,planeNumber);
		else
		{
			cout<<"Error: configuration file for planes in scene not specified.\n";
			exit(1);
		}
	}
	volumetricMesh = NULL;
	// set deformable material type
	if (strcmp(volumetricSurfaceMeshFilename, "__none") != 0)
	{
		if (strcmp(deformableObjectMethod, "StVK") == 0)
			deformableObject = STVK;
		if (strcmp(deformableObjectMethod, "Anisotropic") == 0)
			deformableObject = ANISOTROPIC;
		if (strcmp(deformableObjectMethod, "CLAnistropic") ==0)
			deformableObject = COROTANISTROPIC;
		if (strcmp(deformableObjectMethod, "CLFEM") == 0)
			deformableObject = COROTLINFEM;
		if (strcmp(deformableObjectMethod, "LinearFEM") == 0)
			deformableObject = LINFEM;
		if (strcmp(deformableObjectMethod, "InvertibleFEM") == 0)
			deformableObject = INVERTIBLEFEM;
	}
	if (deformableObject == UNSPECIFIED)
	{
		printf("Error: no deformable model specified.\n");
		exit(1);
	}
	// load mesh
	if ((deformableObject == STVK) || (deformableObject == ANISOTROPIC) || (deformableObject == COROTLINFEM) || (deformableObject == COROTANISTROPIC) || (deformableObject == INVERTIBLEFEM))
	{
		printf("Loading volumetric mesh from file %s...\n", volumetricSurfaceMeshFilename);
		volumetricMesh = VolumetricMeshLoader::load(volumetricSurfaceMeshFilename);
		if (volumetricMesh == NULL)
		{
			printf("Error: unable to load the volumetric mesh from %s.\n", volumetricSurfaceMeshFilename);
			exit(1);
		}
		simulation_vertice_num = volumetricMesh->getNumVertices();
		printf("Num vertices: %d. Num elements: %d\n",simulation_vertice_num, volumetricMesh->getNumElements());
		meshGraph = GenerateMeshGraph::Generate(volumetricMesh);
		// load mass matrix
		if (strcmp(massMatrixFilename, "__none") == 0)
		{
			printf("Error: mass matrix for the StVK deformable model not specified (%s).\n", massMatrixFilename);
			exit(1);
		}
		printf("Loading the mass matrix from file %s...\n", massMatrixFilename);
		// get the mass matrix
		SparseMatrixOutline * massMatrixOutline;
		try
		{
			massMatrixOutline = new SparseMatrixOutline(massMatrixFilename, 1); // We get the massMatrix from LargeModelDeformationFactory, which is already 3n*3n,so here the parameter,we get 1.
		}
		catch(int exceptionCode)
		{
			printf("Error loading mass matrix %s.\n", massMatrixFilename);
			exit(1);
		}
		massMatrix = new SparseMatrix(massMatrixOutline);
		delete(massMatrixOutline);
		if (deformableObject == STVK || deformableObject == LINFEM)  //LINFEM constructed from stVKInternalForces
		{
			unsigned int loadingFlag = 0; // 0 = use low-memory version, 1 = use high-memory version
			StVKElementABCD * precomputedIntegrals = StVKElementABCDLoader::load(volumetricMesh, loadingFlag);
			if (precomputedIntegrals == NULL)
			{
				printf("Error: unable to load the StVK integrals.\n");
				exit(1);
			}
			printf("Generating internal forces and stiffness matrix models...\n"); fflush(NULL);
			if (numInternalForceThreads == 0)
			{
				stVKInternalForces = new StVKInternalForces(volumetricMesh, precomputedIntegrals, addGravity, g);
				stVKStiffnessMatrix = new StVKStiffnessMatrix(stVKInternalForces);
			}
			else
			{
				stVKInternalForces = new StVKInternalForcesMT(volumetricMesh, precomputedIntegrals, addGravity, g, numInternalForceThreads);
				stVKStiffnessMatrix = new StVKStiffnessMatrixMT(stVKInternalForces, numInternalForceThreads);
			}
		}
		else if(deformableObject == ANISOTROPIC)
		{
			if(numInternalForceThreads == 0)
			{
				anisotropicInternalForces = new AnisotropicInternalForces(volumetricMesh, addGravity, g);
				anisotropicStiffnessMatrix = new AnisotropicStiffnessMatrix(anisotropicInternalForces);
			}	
			//load elastic tensor
			if(strcmp(elasticTensorFilename, "__none")!=0)
			{
				anisotropicInternalForces->loadElasticTensorOnCoaseMesh(elasticTensorFilename);
				anisotropicStiffnessMatrix->loadElasticTensorOnCoaseMesh(elasticTensorFilename);
			}
			else
			{
				cout<<"Error: need to load elastic tensor file name.\n";
				exit(1);
			}
		}
	}
	else
	{
		cout<<"Error: Unsupported material type.\n";
		exit(1);
	}
	//convert different material elements of volumetric mesh to volumetric_surface_mesh
	int scaleRows = 1;
	meshGraph->GetLaplacian(&LaplacianDampingMatrix, scaleRows);
	LaplacianDampingMatrix->ScalarMultiply(dampingLaplacianCoef);
	// initialize the rendering mesh for the volumetric mesh
	if (strcmp(renderingMeshFilename, "__none") == 0)
	{
		printf("Error: rendering mesh was not specified.\n");
		exit(1);
	}
	volumetricSurfaceMesh = new SceneObjectDeformable(renderingMeshFilename);
	if (enableTextures)
		volumetricSurfaceMesh->SetUpTextures(SceneObject::MODULATE, SceneObject::NOMIPMAP);
	volumetricSurfaceMesh->ResetDeformationToRest();
	volumetricSurfaceMesh->BuildNeighboringStructure();
	volumetricSurfaceMesh->BuildNormals(); 
	volumetricSurfaceMesh->SetMaterialAlpha(0.5);
	if(objectRenderSurfaceMeshFileNum>0)
	{
		objectRenderSurfaceMesh = new SceneObjectDeformable*[objectRenderSurfaceMeshFileNum];
		objectRenderSurfaceMeshInterpolationElementVerticesNum=new int[objectRenderSurfaceMeshFileNum];
		objectRenderSurfaceMeshInterpolationVertices=new int*[objectRenderSurfaceMeshFileNum];
		objectRenderSurfaceMeshInterpolationWeights=new double*[objectRenderSurfaceMeshFileNum];
		uRenderSurface=new double*[objectRenderSurfaceMeshFileNum];
		for(int i=0;i<objectRenderSurfaceMeshFileNum;++i)
		{
			objectRenderSurfaceMesh[i]=new SceneObjectDeformable(objectRenderSurfaceMeshFilename[i]);
			uRenderSurface[i]=(double*)calloc(3*objectRenderSurfaceMesh[i]->Getn(),sizeof(double));//allocate space for the displacements of surface vertices
			objectRenderSurfaceMeshFilename[i]=new char[string_length];
		}	
	}
	//default rendering mesh is the object render surface mesh
	delete(renderingMesh);
	renderingMesh=volumetricSurfaceMesh/*objectRenderSurfaceMesh0*/;
	if(enableTextures)
		renderingMesh->SetUpTextures(SceneObject::MODULATE,SceneObject::NOMIPMAP);  
	renderingMesh->ResetDeformationToRest();
	renderingMesh->BuildNeighboringStructure();
	renderingMesh->BuildNormals();

	//load interpolation structure for object render surface
	for(int i=0;i<objectRenderSurfaceMeshFileNum;++i)
	{
		if(strcmp(objectRenderSurfaceMeshInterpolationFilename[i],"__none")==0)
		{
			cout<<"Error: no object render surface mesh interpolation filename specified.\n";
			exit(1);
		}
		objectRenderSurfaceMeshInterpolationElementVerticesNum[i]=VolumetricMesh::getNumInterpolationElementVertices(objectRenderSurfaceMeshInterpolationFilename[i]);
		if(objectRenderSurfaceMeshInterpolationElementVerticesNum[i]<0)
		{
			cout<<"Error: unable to open file "<<objectRenderSurfaceMeshInterpolationFilename[i]<<".\n";
			exit(1);
		}
		cout<<"Num interpolation element vertices: "<<objectRenderSurfaceMeshInterpolationElementVerticesNum[i]<<".\n";
		VolumetricMesh::loadInterpolationWeights(objectRenderSurfaceMeshInterpolationFilename[i],objectRenderSurfaceMesh[i]->Getn(),objectRenderSurfaceMeshInterpolationElementVerticesNum[i],
			&objectRenderSurfaceMeshInterpolationVertices[i],&objectRenderSurfaceMeshInterpolationWeights[i]);
	}

	// read the fixed vertices
	// 1-indexed notation
	if (strcmp(fixedVerticesFilename, "__none") == 0)
	{
		numFixedVertices = 0;
		fixedVertices = NULL;
	}
	else
	{
		if (LoadList::load(fixedVerticesFilename, &numFixedVertices,&fixedVertices) != 0)
		{
			printf("Error reading fixed vertices.\n");
			exit(1);
		}
		LoadList::sort(numFixedVertices, fixedVertices);
	}
	printf("Loaded %d fixed vertices. They are:\n",numFixedVertices);
	LoadList::print(numFixedVertices,fixedVertices);
	// create 0-indexed fixed DOFs
	int numFixedDOFs = 3 * numFixedVertices;
	/*int **/ fixedDOFs = (int*) malloc (sizeof(int) * numFixedDOFs);
	for(int i=0; i<numFixedVertices; i++)
	{
		fixedDOFs[3*i+0] = 3*fixedVertices[i]-3;
		fixedDOFs[3*i+1] = 3*fixedVertices[i]-2;
		fixedDOFs[3*i+2] = 3*fixedVertices[i]-1;
	}
	for(int i=0; i<numFixedVertices; i++)
		fixedVertices[i]--;
	printf("Boundary vertices processed.\n");
	// make room for deformation and force vectors
	u = (double*) calloc (3*simulation_vertice_num, sizeof(double));
	f_ext = (double*) calloc (3*simulation_vertice_num, sizeof(double));
	f_extBase = (double*) calloc (3*simulation_vertice_num, sizeof(double));
	f_col=(double*)calloc(3*simulation_vertice_num,sizeof(double));
	gravity_force=(double*)calloc(3*simulation_vertice_num,sizeof(double));
	if(test_case==1)
	{
		initial_pos=(double*) calloc (3*simulation_vertice_num, sizeof(double));
	}
	// load initial condition
	if (strcmp(initialPositionFilename, "__none") != 0)
	{
		int m1, n1;
		ReadMatrixFromDisk_(initialPositionFilename, &m1, &n1, &uInitial);
		if ((m1 != 3*simulation_vertice_num) || (n1 != 1))
		{
			printf("Error: initial position matrix size mismatch.\n");
			exit(1);
		}
	}
	else
		uInitial = (double*) calloc (3*simulation_vertice_num, sizeof(double));
	
	// load initial velocity
	if (strcmp(initialVelocityFilename, "__none") != 0)
	{
		int m1, n1;
		ReadMatrixFromDisk_(initialVelocityFilename, &m1, &n1, &velInitial);
		if ((m1 != 3*simulation_vertice_num) || (n1 != 1))
		{
			printf("Error: initial position matrix size mismatch.\n");
			exit(1);
		}
	}
	else
		velInitial = (double*) calloc (3*simulation_vertice_num, sizeof(double));
	for(unsigned int i=0;i<3*simulation_vertice_num;++i)
		velInitial[i]=0.0;
	// load force loads
	if (strcmp(forceLoadsFilename, "__none") != 0)
	{
		int m1;
		ReadMatrixFromDisk_(forceLoadsFilename, &m1, &numForceLoads, &forceLoads);
		if (m1 != 3*simulation_vertice_num)
		{
			printf("Mismatch in the dimension of the force load matrix.\n");
			exit(1);
		}
	}
	// create force models, to be used by the integrator
	printf("Creating force models...\n");
	if (deformableObject == STVK)
	{
		stVKForceModel = new StVKForceModel(stVKInternalForces, stVKStiffnessMatrix);
		forceModel = stVKForceModel;
		stVKForceModel->GetInternalForce(uInitial, u);
	}
	if (deformableObject == ANISOTROPIC)
	{
		anisotropicForceModel = new AnisotropicForceModel(anisotropicInternalForces, anisotropicStiffnessMatrix);
		forceModel = anisotropicForceModel;
		anisotropicForceModel->GetInternalForce(uInitial,u);

	}
	if (deformableObject == COROTANISTROPIC)
	{
		TetMesh *tetMesh = dynamic_cast<TetMesh*>(volumetricMesh);
		if(tetMesh == NULL)
		{
			cout<<"Error: the input mesh is not a tet mesh(Corotational Anisotropic deformable model).\n";
			exit(1);
		}
		CorotationalAnisotropicFEM * corotationalAnisotropicFEM;
		corotationalAnisotropicFEM = new CorotationalAnisotropicFEM(tetMesh, addGravity, g);
		corotationalAnisotropicFEMForceModel = new CorotationalAnisotropicFEMForceModel(corotationalAnisotropicFEM);
		//load elastic tensor
		if(strcmp(elasticTensorFilename, "__none")!=0)
		{
			corotationalAnisotropicFEM->loadElasticTensorOnCoaseMesh(elasticTensorFilename);
		}
		else
		{
			cout<<"Error: need to load elastic tensor file name.\n";
			exit(1);
		}
		forceModel = corotationalAnisotropicFEMForceModel;
		
	}
	if (deformableObject == COROTLINFEM)
	{
		TetMesh * tetMesh = dynamic_cast<TetMesh*>(volumetricMesh);
		if (tetMesh == NULL)
		{
			printf("Error: the input mesh is not a tet mesh (CLFEM deformable model).\n");
			exit(1);
		}
		CorotationalLinearFEM * corotationalLinearFEM;

		if (numInternalForceThreads == 0)
			corotationalLinearFEM = new CorotationalLinearFEM(tetMesh);
		else
			corotationalLinearFEM = new CorotationalLinearFEMMT(tetMesh, numInternalForceThreads);

		corotationalLinearFEMForceModel = new CorotationalLinearFEMForceModel(corotationalLinearFEM, corotationalLinearFEM_warp);
		forceModel = corotationalLinearFEMForceModel;
		volumetricMesh->computeGravity(gravity_force,g,addGravity);
		
	}
	if (deformableObject == LINFEM)
	{
		LinearFEMForceModel * linearFEMForceModel = new LinearFEMForceModel(stVKInternalForces);
		forceModel = linearFEMForceModel;
	}
	// initialize the integrator
	printf("Initializing the integrator, n = %d...\n", simulation_vertice_num);
	printf("Solver type: %s\n", solverMethod);
	integratorBaseSparse = NULL;
	if (solver == IMPLICITNEWMARK)
	{
		implicitNewmarkSparse = new ImplicitNewmarkSparse(3*simulation_vertice_num, timeStep, massMatrix, forceModel, positiveDefinite, numFixedDOFs, fixedDOFs,
			dampingMassCoef, dampingStiffnessCoef, maxIterations, epsilon, newmarkBeta, newmarkGamma, numSolverThreads);
		integratorBaseSparse = implicitNewmarkSparse;
	}
	else if (solver == IMPLICITBACKWARDEULER)
	{
		implicitNewmarkSparse = new ImplicitBackwardEulerSparse(3*simulation_vertice_num, timeStep, massMatrix, forceModel, positiveDefinite, numFixedDOFs, fixedDOFs,
			dampingMassCoef, dampingStiffnessCoef, maxIterations, epsilon, numSolverThreads);
		integratorBaseSparse = implicitNewmarkSparse;
	}
	else if (solver == EULER)
	{
		int symplectic = 0;
		integratorBaseSparse = new EulerSparse(3*simulation_vertice_num, timeStep, massMatrix, forceModel, symplectic, numFixedDOFs, fixedDOFs, dampingMassCoef);
	}
	else if (solver == SYMPLECTICEULER)
	{
		int symplectic = 1;
		integratorBaseSparse = new EulerSparse(3*simulation_vertice_num, timeStep, massMatrix, forceModel, symplectic, numFixedDOFs, fixedDOFs, dampingMassCoef);
	}
	else if (solver == CENTRALDIFFERENCES)
	{
		integratorBaseSparse = new CentralDifferencesSparse(3*simulation_vertice_num, timeStep, massMatrix, forceModel, numFixedDOFs, fixedDOFs, dampingMassCoef, dampingStiffnessCoef, centralDifferencesTangentialDampingUpdateMode, numSolverThreads);
	}
	integratorBase = integratorBaseSparse;
	if (integratorBase == NULL)
	{
		printf("Error: failed to initialize numerical integrator.\n");
		exit(1);
	}

	/***************************************************************/
	// set integration parameters
	integratorBaseSparse->SetDampingMatrix(LaplacianDampingMatrix);
	integratorBase->ResetToRest();
	integratorBase->SetState(uInitial, velInitial);
	integratorBase->SetTimestep(timeStep);
	if (implicitNewmarkSparse != NULL)
	{
		implicitNewmarkSparse->UseStaticSolver(staticSolver);
		if (velInitial != NULL)
			implicitNewmarkSparse->SetState(implicitNewmarkSparse->Getq(), velInitial);
	}
	// load any external geometry file (e.g. some static scene for decoration; usually there will be none)
	if (strcmp(extraSceneGeometryFilename,"__none") != 0)
	{
		extraSceneGeometry = new SceneObject(extraSceneGeometryFilename);
		extraSceneGeometry->BuildNormals(85.0);
	}
	else
		extraSceneGeometry = NULL;	

	if(output_mode==1)
	{
		// set background color
		int colorR, colorG, colorB;
		sscanf(backgroundColorString, "%d %d %d", &colorR, &colorG, &colorB);
		glClearColor(1.0 * colorR / 255, 1.0 * colorG / 255, 1.0 * colorB / 255, 0.0);
		glui->sync_live();
	}	
}

// set up the configuration file
void initConfigurations()
{
	printf("Parsing configuration file %s...\n", configFilename.c_str());
	ConfigFile configFile;
	// specify the entries of the config file
	//graphics configuration
	configFile.addOptionOptional("windowWidth", &windowWidth, windowWidth);
	configFile.addOptionOptional("windowHeight", &windowHeight, windowHeight);
	configFile.addOptionOptional("cameraRadius", &cameraRadius, 17.5);
	configFile.addOptionOptional("focusPositionX", &focusPositionX, 0.0);
	configFile.addOptionOptional("focusPositionY", &focusPositionY, 0.0);
	configFile.addOptionOptional("focusPositionZ", &focusPositionZ, 0.0);
	configFile.addOptionOptional("cameraLongitude", &cameraLongitude, -60.0);
	configFile.addOptionOptional("cameraLattitude", &cameraLattitude, 20.0);
	configFile.addOption("lightingConfigFilename", lightingConfigFilename);

	configFile.addOptionOptional("deformableObjectMethod", deformableObjectMethod, "StVK");
	configFile.addOptionOptional("initialPositionFilename", initialPositionFilename, "__none");
	configFile.addOptionOptional("initialVelocityFilename", initialVelocityFilename, "__none");
	configFile.addOptionOptional("outputFilename", outputFilename, "__none");
	configFile.addOptionOptional("outputFileNameBase",outputFileNameBase,"output");
	configFile.addOptionOptional("timestepPerOutputFile",&timestepPerOutputFile,timestepPerOutputFile);
	configFile.addOptionOptional("saveMeshToFile",&saveMeshToFile,saveMeshToFile); 
	configFile.addOptionOptional("volumetricMeshFilename", volumetricSurfaceMeshFilename, "__none");
	configFile.addOptionOptional("renderingMeshFilename", renderingMeshFilename, "__none");
	std::cout<<"-----------------------------------------\n";
	std::cout<<"renderingMeshFilename:"<<renderingMeshFilename<<"\n";
	std::cout<<"-----------------------------------------\n";
	configFile.addOptionOptional("objectRenderSurfaceMeshFileNum", &objectRenderSurfaceMeshFileNum, objectRenderSurfaceMeshFileNum);
	objectRenderSurfaceMeshFilename=(char**)malloc(sizeof(char)*objectRenderSurfaceMeshFileNum*string_length);
	objectRenderSurfaceMeshInterpolationFilename=(char**)malloc(sizeof(char)*objectRenderSurfaceMeshFileNum*string_length);
	for(int i=0;i<objectRenderSurfaceMeshFileNum;++i)
	{
		objectRenderSurfaceMeshFilename[i]=(char*)malloc(sizeof(char)*string_length);
		objectRenderSurfaceMeshInterpolationFilename[i]=(char*)malloc(sizeof(char)*string_length);
	}
	for(int i=0;i<objectRenderSurfaceMeshFileNum;++i)
	{
		char option_name[128];
		char index_ch='0'+i;
		sprintf(option_name,"objectRenderSurfaceMeshFilename%c",index_ch);
		configFile.addOptionOptional(option_name,objectRenderSurfaceMeshFilename[i],"__none");
		sprintf(option_name,"objectRenderSurfaceMeshInterpolationFilename%c",index_ch);
		configFile.addOptionOptional(option_name,objectRenderSurfaceMeshInterpolationFilename[i],"__none");
	}
	configFile.addOptionOptional("fixedVerticesFilename", fixedVerticesFilename, "__none");
	configFile.addOptionOptional("massMatrixFilename", massMatrixFilename, "__none");
	configFile.addOptionOptional("forceLoadsFilename", forceLoadsFilename, "__none");
	configFile.addOptionOptional("elasticTensorFilename",elasticTensorFilename,"__none");
	configFile.addOptionOptional("planeNumber",&planeNumber,planeNumber);
	configFile.addOptionOptional("planesFilename",planesFilename,planesFilename);
	configFile.addOptionOptional("extraObjectsNumber",&extraObjectsNum,extraObjectsNum);

	configFile.addOptionOptional("addGravity", &addGravity, addGravity);
	configFile.addOptionOptional("g", &g, g);
	configFile.addOptionOptional("useRealTimeNormals", &useRealTimeNormals, 0);

	configFile.addOption("timeStep", &timeStep);
	configFile.addOptionOptional("frame_rate", &frame_rate, frame_rate);	
	configFile.addOptionOptional("totalFrame", &totalFrame, totalFrame);
	configFile.addOptionOptional("syncTimestepWithGraphics", &syncTimestepWithGraphics, syncTimestepWithGraphics);
	configFile.addOptionOptional("pauseSimulation", &pauseSimulation, pauseSimulation);
	
	// option for corotational linear FEM: if warp is disabled, one gets purely linear FEM
	configFile.addOptionOptional("corotationalLinearFEM_warp", &corotationalLinearFEM_warp, corotationalLinearFEM_warp);
	configFile.addOptionOptional("invertibleMaterial", invertibleMaterialString, invertibleMaterialString);
	configFile.addOptionOptional("implicitSolverMethod", implicitSolverMethod, "none"); // this is now obsolete, but preserved for backward compatibility, use "solver" below
	configFile.addOptionOptional("solver", solverMethod, "implicitNewmark");
	configFile.addOptionOptional("centralDifferencesTangentialDampingUpdateMode", &centralDifferencesTangentialDampingUpdateMode, centralDifferencesTangentialDampingUpdateMode);
	configFile.addOption("dampingMassCoef", &dampingMassCoef);
	configFile.addOption("dampingStiffnessCoef", &dampingStiffnessCoef);
	configFile.addOption("deformableObjectCompliance", &deformableObjectCompliance);
	configFile.addOption("baseFrequency", &baseFrequency);
	configFile.addOptionOptional("dampingLaplacianCoef", &dampingLaplacianCoef, dampingLaplacianCoef);
	configFile.addOptionOptional("newmarkBeta", &newmarkBeta, newmarkBeta);
	configFile.addOptionOptional("newmarkGamma", &newmarkGamma, newmarkGamma);	
	configFile.addOptionOptional("forceNeighborhoodSize", &forceNeighborhoodSize, forceNeighborhoodSize);
	configFile.addOptionOptional("maxIterations", &maxIterations, 1);
	configFile.addOptionOptional("epsilon", &epsilon, 1E-6);
	configFile.addOptionOptional("principalStretchThreshold", &principalStretchThreshold, -DBL_MAX);
	configFile.addOptionOptional("numInternalForceThreads", &numInternalForceThreads, 0);
	configFile.addOptionOptional("numSolverThreads", &numSolverThreads, 1);

	configFile.addOptionOptional("renderWireframe", &renderWireframe, 1);
	configFile.addOptionOptional("renderAxes", &renderAxes, renderAxes);
	configFile.addOptionOptional("extraSceneGeometry", extraSceneGeometryFilename, "__none");
	configFile.addOptionOptional("enableTextures", &enableTextures, enableTextures);
	configFile.addOptionOptional("backgroundColor", backgroundColorString, backgroundColorString);
	configFile.addOptionOptional("groundPlane", groundPlaneString, "__none");
	// parse the configuration file
	if (configFile.parseOptions((char*)configFilename.c_str()) != 0)
	{
		printf("Error parsing options.\n");
		exit(1);
	}
	// the config variables have now been loaded with their specified values

	// informatively print the variables (with assigned values) that were just parsed
	configFile.printOptions();
	// set the solver based on config file input
	solver = UNKNOWN;
	if (strcmp(implicitSolverMethod, "implicitNewmark") == 0)
		solver = IMPLICITNEWMARK;
	if (strcmp(implicitSolverMethod, "implicitBackwardEuler") == 0)
		solver = IMPLICITBACKWARDEULER;
	if (strcmp(solverMethod, "implicitNewmark") == 0)
		solver = IMPLICITNEWMARK;
	if (strcmp(solverMethod, "implicitBackwardEuler") == 0)
		solver = IMPLICITBACKWARDEULER;
	if (strcmp(solverMethod, "Euler") == 0)
		solver = EULER;
	if (strcmp(solverMethod, "symplecticEuler") == 0)
		solver = SYMPLECTICEULER;
	if (strcmp(solverMethod, "centralDifferences") == 0)
		solver = CENTRALDIFFERENCES;
	if (solver == UNKNOWN)
	{
		printf("Error: unknown implicit solver specified.\n");
		exit(1);
	}
}
//
//************************************************************GUI*********************************************
//// GLUI-related functions
void exit_buttonCallBack(int code)
{
	exit(0);
}
void initGLUI()
{
	// generate the UI, via the GLUI library
	glui = GLUI_Master.create_glui("Controls", 0, windowWidth + 52, 80);
	glui->add_button("Exit", 0, exit_buttonCallBack);
	glui->sync_live();
	glui->set_main_gfx_window( windowID );		
}
//***********************************************************************************************************

// main function
int main(int argc, char* argv[])
{
	int numFixedArgs = 2;
	if ( argc < numFixedArgs ) 
	{
		printf("Real-time deformable object simulator.\n");
		printf("Usage: %s [config file]\n", argv[0]);
		return 1;
	}
	// parse command line options
	char * configFilenameC = argv[1];
	opt_t opttable[] =
	{
		{ NULL, 0, NULL }
	};
	argv += (numFixedArgs-1);
	argc -= (numFixedArgs-1);
	int optup = getopts(argc,argv,opttable);
	if (optup != argc)
	{
		printf("Error parsing options. Error at option %s.\n",argv[optup]);
	}
	
	cout<<"Enter output mode:\n";
	cout<<"1.glut_windows\n";
	cout<<"2.output files\n";
	cin>>output_mode;
	test_case=0;
	cout<<"Enter test case number:\n";
	cout<<"0.gravity\n";
	cout<<"1.twist\n";
	cout<<"2.fall off on the slope\n";
	cout<<"3.\n";
	cout<<"4.flower swing in the wind.\n";
	cin>>test_case;
	printf("Starting application.\n");
	configFilename = string(configFilenameC);
	printf("Loading scene configuration from %s.\n", configFilename.c_str());
	initConfigurations(); // parse the config file
	if(output_mode==1)
	{
		initGLUT(argc, argv, windowTitleBase , windowWidth, windowHeight, &windowID);
		initGraphics(windowWidth, windowHeight); // more OpenGL initialization calls
		initGLUI(); // init the UI
	}	
	initSimulation(); // init the simulation
	initFunction(test_case);
	if(output_mode==2)
	{
		outputFilesLoop();
	}
	if(output_mode==1)
	{
		glutMainLoop(); // you have reached the point of no return..
	}	
	return 0;
}
void simulationFunction(int test_case_)
{
}
void initFunction(int test_case_)
{
	if(test_case_==1)
	{
		//the top of the bar is fixed, the bottom rotate
		Vec3d bottom_center(0,-2.55,0);
		double fixed_min_height=2.55;
		double bottom_angle=120;
		double bottom_radian=bottom_angle*2.0*PI/360.0;
		double current_radian=0.0;
		double bottom_min_height=2.55;
		double shrink_bottom_height=1;
		//find the minimal fixed height position along y direction
		for(unsigned int i=0;i<numFixedVertices*3;)
		{
			Vec3d vertex_pos=*volumetricMesh->getVertex(fixedDOFs[i]/3);
			if(vertex_pos[1]<fixed_min_height)
				fixed_min_height=vertex_pos[1];
			i=i+3;
		}
		//find the minimal position along y direction
		for(unsigned int i=0;i<volumetricMesh->getNumVertices();++i)
		{
			Vec3d vetex_pos=*volumetricMesh->getVertex(i);
			if(vetex_pos[1]<bottom_min_height)
				bottom_min_height=vetex_pos[1];
			//set initial velocity of each vertex
			velInitial[3*i+0]=velInitial[3*i+1]=velInitial[3*i+2]=0.0;
		}
		//set initial twist bar position
		for(unsigned int i=0;i<volumetricMesh->getNumVertices();++i)
		{				
			Vec3d vec_pos=*volumetricMesh->getVertex(i);
			//get the twist radian for each vertex
			current_radian=-bottom_radian*(vec_pos[1]-fixed_min_height)/(bottom_min_height-fixed_min_height);				
			initial_pos[3*i]=bottom_center[0]+(vec_pos[0]-bottom_center[0])*cos(current_radian)-(vec_pos[2]-bottom_center[2])*sin(current_radian);
			initial_pos[3*i+1]=vec_pos[1]/*+shrink_bottom_height*(vec_pos[1]-fixed_min_height)/(bottom_min_height-fixed_min_height)*/;
			initial_pos[3*i+2]=bottom_center[2]+(vec_pos[0]-bottom_center[0])*sin(current_radian)+(vec_pos[2]-bottom_center[2])*cos(current_radian);

			u[3*i]=initial_pos[3*i]-vec_pos[0];	
			u[3*i+1]=initial_pos[3*i+1]-vec_pos[1];
			u[3*i+2]=initial_pos[3*i+2]-vec_pos[2];			
		}
		//set the displacement of the fixed vertices
		for(unsigned int i=0;i<numFixedVertices;++i)
		{
			u[fixedDOFs[3*i+0]]=u[fixedDOFs[3*i+1]]=u[fixedDOFs[3*i+2]]=0.0;
		}
		integratorBase->SetState(u, velInitial);
	}
	else if(test_case==2)
	{
		Vec3d rotate_center(-1.25,-0.5,-0.5);
		double rotate_omega=-60*2*PI/360.0;
		for(unsigned int i=0;i<volumetricMesh->getNumVertices();++i)
		{
		//	//velInitial[3*i+1]=-10.0;
		//	/*velInitial[3*i+1]=rotate_omega*(vect_pos[2]-rotate_center[2]);
		//	velInitial[3*i+2]=(-1.0)*rotate_omega*(vect_pos[1]-rotate_center[1]);
		//	velInitial[3*i]=0.0;*/		
		//	/*u[3*i]=rotate_center[0]+(vect_pos[0]-rotate_center[0])*cos(rotate_omega)-(vect_pos[1]-rotate_center[1])*sin(rotate_omega)-vect_pos[0];
		//	u[3*i+1]=rotate_center[1]+(vect_pos[0]-rotate_center[0])*sin(rotate_omega)+(vect_pos[1]-rotate_center[1])*cos(rotate_omega)-vect_pos[1];
		//	u[3*i+2]=0.0;*/
		//	/*velInitial[3*i+0]=velInitial[3*i+1]=velInitial[3*i+2]=0.0;*/	

 

			Vec3d vect_pos=*volumetricMesh->getVertex(i);
			u[3*i]=u[3*i+1]=u[3*i+2]=0.0;
			velInitial[3*i+1]=-2.0+/*-10.0+*/rotate_omega*(vect_pos[0]-rotate_center[0]);
			velInitial[3*i]=(-1.0)*rotate_omega*(vect_pos[1]-rotate_center[1]);
			velInitial[3*i+2]=0.0;
		}
		integratorBase->SetState(u,velInitial);
	}
	else if(test_case==3)
	{
		Vec3d rotate_center(1.2,0,0);
		//Vec3d rotate_center(0,-0.5,0);
		double rotate_omega=-60*2*PI/360.0;
		for(unsigned int i=0;i<volumetricMesh->getNumVertices();++i)
		{
			//
			/*Vec3d vect_pos=*volumetricMesh->getVertex(i);
			u[3*i]=(vect_pos[0]-rotate_center[0])*cos(rotate_omega)-(vect_pos[1]-rotate_center[1])*sin(rotate_omega)-vect_pos[0];
			u[3*i+1]=(vect_pos[0]-rotate_center[0])*sin(rotate_omega)+(vect_pos[1]-rotate_center[1])*cos(rotate_omega)-vect_pos[1];
			u[3*i+2]=0.0;
			velInitial[3*i]=velInitial[3*i+1]=velInitial[3*i+2]=0.0;*/
			//totate
			Vec3d vect_pos=*volumetricMesh->getVertex(i);
			velInitial[3*i+1]=10.0/*rotate_omega*(vect_pos[2]-rotate_center[2])*/;
			velInitial[3*i+2]=0.0/*(-1.0)*rotate_omega*(vect_pos[1]-rotate_center[1])*/;
			velInitial[3*i]=0.0;			
			u[3*i]=u[3*i+1]=u[3*i+2]=0.0;
		}
		integratorBase->SetState(u,velInitial);
	}
	else if(test_case==4)
	{
		//flower swing in the wind
		for(unsigned int i=0;i<volumetricMesh->getNumVertices();++i)
		{
			f_ext[3*i]=-1.0;
			f_ext[3*i+1]=f_ext[3*i+2]=0.0;
			/*velInitial[3*i+1]=-10.0;
			velInitial[3*i+2]=0.0;
			velInitial[3*i]=0.0;			
			u[3*i]=u[3*i+1]=u[3*i+2]=0.0;*/
		}
		for(unsigned int i=0;i<numFixedVertices;++i)
		{
			f_ext[fixedDOFs[3*i+0]]=f_ext[fixedDOFs[3*i+1]]=f_ext[fixedDOFs[3*i+2]]=0.0;
		}
		//integratorBase->SetState(u,velInitial);
		integratorBaseSparse->SetExternalForces(f_ext);
		integratorBase->DoTimestep();
	}
}
void testStiffnessMatrix(void)
{
	//test K
	double *u1,*du,*f1,*f0,*multi_k_du,*error_value;
	u1=(double*)malloc(sizeof(double)*3*simulation_vertice_num);
	du=(double*)malloc(sizeof(double)*3*simulation_vertice_num);
	f1=(double*)malloc(sizeof(double)*3*simulation_vertice_num);
	f0=(double*)malloc(sizeof(double)*3*simulation_vertice_num);
	multi_k_du=(double*)malloc(sizeof(double)*3*simulation_vertice_num);
	error_value=(double*)malloc(sizeof(double)*3*simulation_vertice_num);
	SparseMatrix *K0,*K1;
	double f1_max,f0_max,multi_k_du_max;
	forceModel->GetTangentStiffnessMatrixTopology(&K0);
	forceModel->GetTangentStiffnessMatrixTopology(&K1);
	double max_error=0.0;		
	for(unsigned int l=0;l<simulation_vertice_num;++l)
	{
		Vec3d vert_pos=*volumetricMesh->getVertex(l);
		if(vert_pos[0]>0)
			u[3*l]=-0.1;
		else
			u[3*l]=0.1;
		if(vert_pos[1]>0)
			u[3*l+1]=0.1;
		else
			u[3*l+1]=-0.2;
		if(vert_pos[2]>0)
			u[3*l+2]=0.1;
		else
			u[3*l+2]=-0.3;
	}	
	double ran_numf=0.0;
	srand((unsigned)time(0));
	for(unsigned int l=0;l<3*simulation_vertice_num;++l)
	{
		du[l]=(1.0e-5)*(rand()/(double)(RAND_MAX));
		error_value[l]=0.0;
	}
	forceModel->GetForceAndMatrix(u,f0,K0);
	K0->MultiplyVector(du,multi_k_du);
	for(unsigned int l=0;l<3*simulation_vertice_num;++l)
	{
		u1[l]=u[l]+du[l];
	}
	forceModel->GetForceAndMatrix(u1,f1,K1);

	for(unsigned int l=0;l<3*simulation_vertice_num;++l)
	{
		error_value[l]=f1[l]-f0[l]-multi_k_du[l];
	}

	for(unsigned int l=0;l<3*simulation_vertice_num;++l)
	{
		if(fabs(f1[l])>1.0e-8)
		{
			if(fabs(error_value[l]/f1[l])>max_error);
			{
				max_error=fabs(error_value[l]/f1[l]);
				f1_max=f1[l];
				f0_max=f0[l];
				multi_k_du_max=multi_k_du[l];
			}
		}
		else
		{
			std::cout<<"!";
			//max_error=0.0;
		}
	}
	std::cout<<max_error<<"--f1="<<f1_max<<"--f0="<<f0_max<<"--multi_k_du_max="<<multi_k_du_max<<"----------------------------";
	delete [] u1;
	delete [] du;
	delete [] f0;
	delete [] f1;
	delete [] multi_k_du;
	delete [] error_value;
}
