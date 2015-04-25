/**
anisotropicForceModel.cpp
@author Mirror
**/

#include "AnisotropicForceModel.h"
#include "tetMesh.h"
#include <iostream>

AnisotropicForceModel::AnisotropicForceModel(AnisotropicInternalForces * anisotropicInternalForces, AnisotropicStiffnessMatrix * anistropicStiffnessMatrix): anisotropic_internal_forces_(anisotropicInternalForces), anistropic_stiffness_matrix_(anistropicStiffnessMatrix)
{
	r = 3 * anisotropicInternalForces->GetVolumetricMesh()->getNumVertices();
	ownStiffnessMatrix = false;
	if (anistropicStiffnessMatrix == NULL)
	{
		anistropicStiffnessMatrix = new AnisotropicStiffnessMatrix(anisotropicInternalForces);
		ownStiffnessMatrix = true;
	}
}

AnisotropicForceModel::~AnisotropicForceModel()
{
	if (ownStiffnessMatrix)
		delete(anistropic_stiffness_matrix_);
}

void AnisotropicForceModel::GetInternalForce(double * u, double * internalForces)
{
	anisotropic_internal_forces_->ComputeForces(u, internalForces);

}

void AnisotropicForceModel::GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix)
{
	anistropic_stiffness_matrix_->GetStiffnessMatrixTopology(tangentStiffnessMatrix);
}

void AnisotropicForceModel::GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix)
{
	anistropic_stiffness_matrix_->ComputeStiffnessMatrix(u, tangentStiffnessMatrix);
	//tangentStiffnessMatrix->ResetToZero();
} 
