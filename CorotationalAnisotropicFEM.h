/*corotational anisotropic material*/

#ifndef _CorotationalAnisotropicFEM_H_
#define _CorotationalAnisotropicFEM_H_

#include <string>
#include "mat3d.h"
#include "sparseMatrix.h"
#include "tetMesh.h"
/*
   To compute the stiffness matrix for an corotational anistropic material. 
   The method to compute Stiffnessmatrix and internal force is same as anistropic material on each element.
   First,remove the rotation part of the current position. 
   Second, take the left small deformation part on the current position to compute the stiffness matrix and the internal force.

*/
class CorotationalAnisotropicFEM
{
public:
  CorotationalAnisotropicFEM(TetMesh * tetMesh,bool addGravity, double g);
  virtual ~CorotationalAnisotropicFEM();
  void GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology); //build stiffness matrix skeleton as 3n*3n
  inline void ResetStiffnessMatrix(SparseMatrix * sparseMatrix) {sparseMatrix->ResetToZero();} 
  virtual void ComputeForceAndStiffnessMatrix(double * vertexDisplacements, double * internalForces, SparseMatrix * stiffnessMatrix);
  //load elastic_tensor file and store the value in the coarse_mesh_elastic_tensor
  void loadElasticTensorOnCoaseMesh(const std::string &elastic_tensor_file_name); 
  inline TetMesh * GetTetMesh() { return tet_mesh_; }
  void SetGravity(bool addGravity) { this->add_gravity_ = addGravity; InitGravity(); } // if addGravity is enabled, ComputeForces will subtract the gravity force from the internal forces (note: subtraction, not addition, is used because the internal forces are returned with the sign as described in the f_int(x) comment above)
protected: 
  enum StrainType{CAUCHY_STRAIN, GREEN_STRAIN} strainType; 
  void WarpMatrix(double * K, double * R, double * RK, double * RKRT);
  void ComputeStiffnessMatrixOnEachElement(double * small_deformation_pos, unsigned int ele_idx, double * K_ele);
  void ComputeForceAndStiffnessMatrixOfSubmesh(double * current_ele_dis, double * internalForces, SparseMatrix * stiffnessMatrix, int elementLo, int elementHi);
  void ComputeForces(const double * vertexDisplacements, double * forces);
  Mat3d getInistialDisplacementMatrixOnEachElement(unsigned int ele_idx) const;//Dm for one element
  //compute Ds for all elements, "vert_displacement" is an array of vertex deformations, of length 3*n
  std::vector<Mat3d> getCurrentDisplacementMatrixOnAllElements(const double *vert_displacement) const;
  //input parameter is the current position on each element
  Mat3d getCurrentDisplacementMatrixOnEachElement(const double *current_ele_pos,unsigned int ele_idx) const;
  //F=dsDm^-1
  Mat3d getDeformationGradient(const Mat3d init_dis, const Mat3d current_dis) const;
  //---dF/dx_j^k=e_k*e_j^T*D_m^(-1),vert_idx denotes j, vert_idx_dim denotes k
  Mat3d getDeformationGradientDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim) const;
  double * getGreenStrainTensor(const Mat3d F) const;
  double * getCauchyStrainTensor(const Mat3d F) const;
  double * getGreenStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const;
  double * getCauchyStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const;
  Mat3d firstPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,CorotationalAnisotropicFEM::StrainType strain_type) const;
  Mat3d secondPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,CorotationalAnisotropicFEM::StrainType strain_type) const;
  // dP/dx=dF/dx*(C:E)+F(C:dE/dx);
  Mat3d firstPiolaKirchhoffStressDerivative(unsigned int ele_idx, const Mat3d F,unsigned int vert_idx, 
											 unsigned int vert_idx_dim,CorotationalAnisotropicFEM::StrainType strain_type) const; 
  void InitGravity(); // aux function
protected:
	TetMesh * tet_mesh_;
	double *** ele_elastic_tensor_vector_;
	double * undeformedPositions;
	int vert_num_;//vertices num
	int ele_num_;//elements num
	int ele_vert_num_;//vertices num in each element
	int ** row_;
	int ** column_;
	double * gravity_force_;
	bool add_gravity_;
	double g_;
};
#endif

