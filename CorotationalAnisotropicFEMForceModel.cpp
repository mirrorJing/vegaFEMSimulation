/**
CorotationalAnisotropicFEMForceModel.cpp
@author Mirror
**/

#include "CorotationalAnisotropicFEMForceModel.h"
#include "tetMesh.h"

CorotationalAnisotropicFEMForceModel::CorotationalAnisotropicFEMForceModel(CorotationalAnisotropicFEM * corotationalAnisotropicFEM): corotational_anisotropic_fem_(corotationalAnisotropicFEM)
{
	r = 3 * corotational_anisotropic_fem_->GetTetMesh()->getNumVertices();
}

CorotationalAnisotropicFEMForceModel::~CorotationalAnisotropicFEMForceModel()
{
}

void CorotationalAnisotropicFEMForceModel::GetInternalForce(double * u, double * internalForces)
{
	corotational_anisotropic_fem_->ComputeForceAndStiffnessMatrix(u, internalForces, NULL);
}

void CorotationalAnisotropicFEMForceModel::GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix)
{
	corotational_anisotropic_fem_->GetStiffnessMatrixTopology(tangentStiffnessMatrix);
}

void CorotationalAnisotropicFEMForceModel::GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix)
{
	corotational_anisotropic_fem_->ComputeForceAndStiffnessMatrix(u, NULL, tangentStiffnessMatrix);
} 
void CorotationalAnisotropicFEMForceModel::GetForceAndMatrix(double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix)
{
	corotational_anisotropic_fem_->ComputeForceAndStiffnessMatrix(u, internalForces, tangentStiffnessMatrix);
}
