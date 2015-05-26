/*anisotropic material*/

#ifndef _AnisotropicStiffnessMatrix_H_
#define _AnisotropicStiffnessMatrix_H_

#include <string>
#include "AnisotropicInternalForces.h"
#include "mat3d.h"
#include "sparseMatrix.h"
/*
   To compute the stiffness matrix for an anistropic material. The final matrix is a 3n*3n matrix.(n is the vert_num)
   For a tetMesh, we have to compute the derivative of each element on the four vertices in one element.
   for each internal force on the vertex, df/dx_j^k is a 3*3 matrix 
   df_i/dx_j^k denotes a derivative of force on vertex i wrt to vertex j on the k-th coordinate
   (j denotes the index of four vertices, k is the x,y,z coordinates)
   The computation equation is: dH/dx=[df1 df2 df3]/dx, df4/dx=-(df1/dx+df2/dx+df3/dx);
   dH/dx=Ve*(dP/dx)*Dm^-T;
   1.F=DsDm^-1;
   2.F derivative:dF/dx_j is a 3-order tensor, df/dx_j^k=(dDs/dx)*Dm^-1 is a 3*3 matrix
   (1)dF/dx_j^k=e_k*e_j^T*D_m^(-1); for j=1,2,3
   (2)dF/dx_4 is solved directly by ds for j=4
   3.GreenStrainTensor:  E=0.5*(F^T*F-I);
   4.CauchyStrainTensor: e=0.5*(F+F^T)-I
   5.GreenStrainTensorDerivative: dE/dx=0.5*(dF^T*F+F^T*dF/dx);
   6.CauchyStrainTensorDerivative: de/dx=0.5*(dF/dx+dF^T/dx);
	For dH/dx_j^k, the three columns denote [df1/dx_j^k df2/dx_j^k df3/dx_j^k]
	For k=1,2,3 we assemble the result and get df1/dx_j=[df1/dx_j^1 df1/dx_j^2 df1/dx_j^3] as an example
	We have to assemble the stiffness matrix for each df/dx_j
	Considering f=Kx.
	For each line of the global K,the row denotes f entry, the column denotes x entry.
*/
class AnisotropicStiffnessMatrix
{
public:
  AnisotropicStiffnessMatrix(AnisotropicInternalForces * anistropicInternalForces);
  virtual ~AnisotropicStiffnessMatrix();
  void GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology); //build stiffness matrix skeleton as 3n*3n
  inline void ResetStiffnessMatrix(SparseMatrix * sparseMatrix) {sparseMatrix->ResetToZero();} 
  inline VolumetricMesh * GetVolumetricMesh() { return volumetric_mesh_; }
  //load elastic_tensor file and store the value in the coarse_mesh_elastic_tensor
  void loadElasticTensorOnCoaseMesh(const std::string &elastic_tensor_file_name);  
  Mat3d getInistialDisplacementMatrixOnEachElement(unsigned int ele_idx) const;//Dm for one element
  //compute Ds for all elements, "vert_displacement" is an array of vertex deformations, of length 3*n
  std::vector<Mat3d> getCurrentDisplacementMatrixOnAllElements(const double *vert_displacement) const;
  //F=dsDm^-1
  Mat3d getDeformationGradient(const Mat3d init_dis, const Mat3d current_dis) const;
  //---dF/dx_j^k=e_k*e_j^T*D_m^(-1),vert_idx denotes j, vert_idx_dim denotes k
  Mat3d getDeformationGradientDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim) const;
  double * getGreenStrainTensor(const Mat3d F) const;
  double * getCauchyStrainTensor(const Mat3d F) const;
  double * getGreenStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const;
  double * getCauchyStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const;
  enum StrainType{CAUCHY_STRAIN, GREEN_STRAIN} strainType;
  
  // dP/dx=dF/dx*(C:E)+F(C:dE/dx);
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
	inline void AddMatrix3x3Block(int c, int a, int element, Mat3d & matrix, SparseMatrix * sparseMatrix);

    int ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const;    //modified SVD for inversion handling
    // given a vector, find a unit vector that is orthogonal to it
    void FindOrthonormalVector(Vec3d & v, Vec3d & result) const;
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

