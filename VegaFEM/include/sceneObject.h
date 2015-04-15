/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 1.1                               *
 *                                                                       *
 * "sceneObject" library , Copyright (C) 2007 CMU, 2009 MIT, 2012 USC    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Daniel Schroeder                         *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  A library to load and render an obj mesh (called here a "scene object").
  This library depends on the "objMesh" library. See also objMesh.h.
*/

/*
  Fei Zhu added rendering with custom colors function:
  *	  virtual void RenderWithCustomColors();
  *	  void setCustomColors(Vec3d color); 
  *	  void setCustomColors(std::vector<Vec3d> colors); 
  *	  void resetRenderMode();
  2012 PKU
  zhuf@graphics.pku.edu.cn
*/

#ifndef _SCENEOBJECT_H_
#define _SCENEOBJECT_H_

#ifdef WIN32
  #include <windows.h>
#endif

#include "macros.h"
#include "objMesh.h"
#include "objMeshRender.h"
#include "minivector.h"

class SceneObject
{
public:

  enum LightingModulationType {MODULATE, REPLACE};
  enum MipmapType {USEMIPMAP, NOMIPMAP};

  // create a static scene object, by loading it from an Alias Wavefront OBJ file
  SceneObject(char * filename);
  virtual ~SceneObject();

  // ==== render ====

  virtual void Render();
  virtual void RenderWithCustomColors();
  virtual void RenderVertices();
  virtual void RenderEdges(); 
  virtual void RenderNormals();
  virtual void RenderFacesAndEdges();

  virtual void RenderVertex(int vertex);
  virtual void RenderVertices(int numVertices, int * vertexList); // 0-indexed vertices
  virtual void RenderVertices_Selection();

  virtual void RenderShadow(double ground[4], double light[4]);
  void RenderEdgesInGroup(char * groupName); // renders only the edges in the given group

  // ==== custom colors ===
  void setCustomColors(Vec3d color); // constant color for each mesh vertex
  void setCustomColors(std::vector<Vec3d> colors); // specific color for every mesh vertex
  void resetRenderMode();

  // ==== display lists ====

  void BuildDisplayList();
  void PurgeDisplayList();

  // ==== mesh info and geometric queries ====

  inline int Getn() { return n; }
  inline int GetNumVertices() { return n; }
  inline int GetNumFaces() { return mesh->getNumFaces(); }
  inline ObjMesh * GetMesh() { return mesh; }

  // smallest ball radius that encloses the model, with the ball centered at the given centroid
  void ComputeMeshRadius(Vec3d & centroid, double * radius);
  // compute mesh centroid and smallest enclosing radius
  void ComputeMeshGeometricParameters(Vec3d * centroid, double * radius);
  // export mesh data
  void ExportMeshGeometry(int * numVertices, double ** vertices, int * numTriangles, int ** triangles);

  // finds the closest vertex using an exhaustive search
  // returns distance in "distance", if distance is not NULL
  // in this class, you can safely ignore the last parameter (keep it NULL)
  virtual int GetClosestVertex(Vec3d & queryPos, double * distance=NULL, double * auxVertexBuffer=NULL);

  // ==== texture mapping and materials ====

  // lightingModulation determines whether to multiply or replace the object color with the texture color
  // mipmap determines whether to use mipmapping for texture rendering
  int SetUpTextures(LightingModulationType lightingModulation=MODULATE, MipmapType mipmap=USEMIPMAP); // you must call this **after** OpenGL has been initialized!!!
  void EnableTextures();
  void DisableTextures();
  inline bool hasTextures() { return hasTextures_; }

  void SetMaterialAlpha(double alpha); // sets the material alpha value for all the materials

  // ==== normals ====

  void BuildFaceNormals();
  void BuildNeighboringStructure();  // must be called before the vertex-normal functions below

  // second parameter is treshold angle for hard edges:
  void BuildVertexNormals(double thresholdAngle=85.0); // assumes pre-existing face normals
  void BuildNormals(double thresholdAngle=85.0); // builds both face and vertex normals
  void BuildNormalsFancy(double thresholdAngle=85.0);  // rebuilds facet normals + calls vertex-per-triangle normal update

  void SetNormalsToFaceNormals();

  // ==== vertex labelling ====

  // shows all point labels
  void ShowPointLabels();
  // shows point labels for points from k to l
  void ShowPointLabels(int k, int l);
  void HighlightVertex(int i); // same as RenderVertex, except it always renders in green color and with point size 8.0

  // model matrix for the shadow
  static void SetShadowingModelviewMatrix(double ground[4], double light[4]);

protected:

  int n;
  unsigned int renderMode;

  ObjMesh * mesh;
  ObjMeshRender * meshRender;
  GLuint displayList;
  bool displayListExists;

  GLuint displayListEdges;
  bool displayListEdgesExists;

  void PrintBitmapString(float x,float y, float z, char* s);
  void PrintBitmapInteger(float x,float y, float z,long i);

  bool hasTextures_;
};

#endif
