/*anisotropic material*/

#ifndef _AnisotropicStiffnessMatrix_H_
#define _AnisotropicStiffnessMatrix_H_

#include <string>
#include "AnistropicInternalForces.h"
#include "mat3d.h"
#include "sparseMatrix.h"
/*
   anisotropic material. Material properties are read from the tet mesh.
   1.F=DsDm^-1;
   2.F derivative:
   (1)dF/dx_j^k=e_k*e_j^T*D_m^(-1);  j:vertex index, k: the k-th coordinate of the derivative vertex----for j=1,2,3;k=1,2,3
   (2)dF/dx_4 is solved directly by ds for j=4,k=1,2,3;
   3.GreenStrainTensor:  E=0.5*(F^T*F-I);
   4.CauchyStrainTensor: e=0.5*(F+F^T)-I
   5.GreenStrainTensorDerivative: dE/dx=0.5*(dF^T*F+F^T*dF/dx);
   6.CauchyStrainTensorDerivative: de/dx=0.5*(dF/dx+dF^T/dx);
   7.
*/
class AnisotropicStiffnessMatrix
{
public:
  AnisotropicStiffnessMatrix(AnisotropicInternalForces * anistropicInternalForces);
  virtual ~AnisotropicStiffnessMatrix();
  void GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology); 
  inline void ResetStiffnessMatrix(SparseMatrix * sparseMatrix) {sparseMatrix->ResetToZero();} 
  inline VolumetricMesh * GetVolumetricMesh() { return volumetric_mesh_; }

  void loadElasticTensorOnCoaseMesh(const std::string &elastic_tensor_file_name);  //load elastic_tensor file and store the value in the coarse_mesh_elastic_tensor
  Mat3d getInistialDisplacementMatrixOnEachElement(unsigned int ele_idx) const;
  std::vector<Mat3d> getCurrentDisplacementMatrixOnAllElements(const double *vert_displacement) const;
  Mat3d getDeformationGradient(const Mat3d init_dis, const Mat3d current_dis) const;
  //---dF/dx_j^k=e_k*e_j^T*D_m^(-1)
  Mat3d getDeformationGradientDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim) const;
  double * getGreenStrainTensor(const Mat3d F) const;
  double * getCauchyStrainTensor(const Mat3d F) const;
  double * getGreenStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const;
  double * getCauchyStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const;
  enum StrainType{CAUCHY_STRAIN, GREEN_STRAIN} strainType;
  
  // dP/dx=dF/dx*CE+FC*dE/dx;
  Mat3d firstPiolaKirchhoffStressDerivative(unsigned int ele_idx, const Mat3d F,unsigned int vert_idx, 
											 unsigned int vert_idx_dim,AnisotropicStiffnessMatrix::StrainType strain_type) const;
  // evaluates the tangent stiffness matrix in the given deformation configuration
  // "vertexDisplacements" is an array of vertex deformations, of length 3*n, where n is the total number of mesh vertices
  virtual void ComputeStiffnessMatrix(double * vertexDisplacements, SparseMatrix * sparseMatrix);
  
protected:
	AnisotropicInternalForces * anistropic_Internal_Forces_;
	VolumetricMesh * volumetric_mesh_;
	TetMesh * tet_mesh_;
	double *** ele_elastic_tensor_vector_;
	int vert_num_;//vertices num
	int ele_num_;//elements num
	int ele_vert_num_;//vertices num in each element
	int ** row_;
	int ** column_;

	// adds a 3x3 block matrix corresponding to a derivative of force on vertex c wrt to vertex a
	// c is 0..7
	// a is 0..7
	inline void AddMatrix3x3Block(int c, int a, int element, Mat3d & matrix, SparseMatrix * sparseMatrix);
};

inline void AnisotropicStiffnessMatrix::AddMatrix3x3Block(int c, int a, int element, Mat3d & matrix, SparseMatrix * sparseMatrix)
{
	int * row = row_[element];
	int * column = column_[element];

	for(int k=0; k<3; k++)
		for(int l=0; l<3; l++)
			sparseMatrix->AddEntry(3*row[c]+k, 3*column[ele_vert_num_*c+a]+l, matrix[k][l]);
}
#endif

