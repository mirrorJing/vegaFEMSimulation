/*AnistropicInternalForce.h*/

#ifndef _AnisotropicInternalForces_H_
#define _AnisotropicInternalForces_H_

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
  Mat3d getInistialDisplacementMatrixOnEachElement(unsigned int ele_idx) const;
 // Mat3d getCurrentDisplacementMatrixOnEachElement(unsigned int ele_idx,const double *vert_displacement) const; 
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

 // Vec3d * getWeightGradient(unsigned int ele_idx) const;//get weight gradient of each element
  //void getDisplacementDifferentiationMatrix(unsigned int ele_idx,double **dis_differentiation_vector);//get the displancement differentiation vector contains 6*12 entries

  // both vertex displacements and internal forces refer to the vertices of the simulation mesh
  // they must be (pre-allocated) vectors of length 3 * numVertices
  // the internal forces are returned with the sign corresponding to f_int(x) on the left side of the equation M * x'' + f_int(x) = f_ext
  // i.e., the computed internal forces are negatives of the actual physical internal forces acting on the material
  
  //SparseMatrix * ComputeStiffnessMatrixForEachElement(const unsigned int &ele_idx) const;//the stiffness matrix for each element is a 12*12 matrix
  
  //void ComputeStiffnessMatrix(SparseMatrix * sparseMatrix);
 // virtual double ComputeEnergy(const double * vertexDisplacements); // get the nonlinear elastic strain energy
protected:
	void InitGravity(); // aux function
protected:
	VolumetricMesh * volumetric_mesh_;
	double *** ele_elastic_tensor_vector_;
	//std::vector<SparseMatrix *> elastic_tensor_matrix_;
	std::vector<double> strain_tensor_;
	//std::vector<std::vector<double> > dis_differentiation_vector_;
	int vert_num_;
	int ele_num_;
	int ele_vert_num_;
	double * gravity_force_;
	TetMesh * tet_mesh_;
	bool add_gravity_;
	double g_;
};

#endif

