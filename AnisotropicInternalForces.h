/*AnistropicInternalForce.h*/

#ifndef _AnisotropicInternalForces_H_
#define _AnisotropicInternalForces_H_

/*
 It contains four vertices in a tetMesh. To compute the internal force for each vertex,we use the equation:
 H=[f1 f2 f3]=-VePDm^-T, f4=-f1-f2-f3;
 We should notice that here we compute the internal force f in Mu'' + Du' + f = f_ext .
 So we actually compute H=VePDm^-T
 P is the firstPiolaKirchhoffStress, which is computed as P=FS
 S is the secondPiolaKirchhoffStress, which is denoted as S=C:E
 C is the elastic tensor which is loaded by an input file.
 E is the strain tensor, which is divided by the Green strain tensor and the Cauchy strain tensor
 for Cauchy strain tensor: e=0.5*(F+F^T)-I;
 for Green strain tensor: E=0.5*(F^T*F-I);
 F denotes the deformation gradient, which is computed as F=DsDm^-1
 Ds is the displacement of current configuration and Dm is the displacement of the initial configuration
 Ds=[x1-x4 x2-x4 x3-x4]     Dm=[X1-X4 X2-X4 X3-X4]
    [y1-y4 y2-y4 y3-y4]        [Y1-Y4 Y2-Y4 Y3-Y4]
	[z1-z4 z2-z4 z3-z4]        [Z1-Z4 Z2-Z4 Z3-Z4]
*/

#include <vector>
#include <string>
#include "sparseMatrix.h"
#include "volumetricMesh.h"
#include "tetMesh.h"
#include "mat3d.h"
#include "vec3d.h"


class AnisotropicInternalForces
{
public:
  AnisotropicInternalForces(VolumetricMesh * volumetricMesh, bool addGravity=false, double g=9.81);

  virtual ~AnisotropicInternalForces();

  void loadElasticTensorOnCoaseMesh(const std::string &elastic_tensor_file_name);  //load elastic_tensor file and store the value in the coarse_mesh_elastic_tensor
  Mat3d getInitialDisplacementMatrixOnEachElement(unsigned int ele_idx) const;
  //the displacement is 3n*1 vector, n is the vertex number
  std::vector<Mat3d> getCurrentDisplacementMatrixOnAllElements(const double *vert_displacement) const;
  Mat3d getDeformationGradient(const Mat3d init_dis, const Mat3d current_dis) const;
  double * getGreenStrainTensor(const Mat3d F) const;
  double * getCauchyStrainTensor(const Mat3d F) const;
  enum StrainType{CAUCHY_STRAIN, GREEN_STRAIN} strainType;
  Mat3d firstPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,AnisotropicInternalForces::StrainType strain_type) const;
  Mat3d secondPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,AnisotropicInternalForces::StrainType strain_type) const;
  virtual void ComputeForces(const double * vertexDisplacements, double * internalForces);
  void ResetVector(double * vec); // aux function
  // enables or disables the gravity (note: you can also set this in the constructor; use this routine to turn the gravity on/off during the simulation)
  void SetGravity(bool addGravity) { this->add_gravity_ = addGravity; InitGravity(); } // if addGravity is enabled, ComputeForces will subtract the gravity force from the internal forces (note: subtraction, not addition, is used because the internal forces are returned with the sign as described in the f_int(x) comment above)
  inline VolumetricMesh * GetVolumetricMesh() { return volumetric_mesh_; }
protected:
	void InitGravity(); // aux function
    int ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const;    //modified SVD for inversion handling
    // given a vector, find a unit vector that is orthogonal to it
    void FindOrthonormalVector(Vec3d & v, Vec3d & result) const;
protected:
	VolumetricMesh * volumetric_mesh_;
	double *** ele_elastic_tensor_vector_;
	std::vector<double> strain_tensor_;
	int vert_num_;
	int ele_num_;
	int ele_vert_num_;
	TetMesh * tet_mesh_;
	bool add_gravity_;
	double * gravity_force_;
	double g_;
};

#endif

