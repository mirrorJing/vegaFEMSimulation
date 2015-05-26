/*
	Force model corresponding to the anisotropic material.
*/

#ifndef _ANISOTROPICFORCEMODEL_H_
#define _ANISOTROPICFORCEMODEL_H_

#include "forceModel.h"
#include "AnisotropicInternalForces.h"
#include "AnisotropicStiffnessMatrix.h"

class AnisotropicForceModel : virtual public ForceModel
{
public:
  AnisotropicForceModel(AnisotropicInternalForces * anistropicInternalForces, AnisotropicStiffnessMatrix * anisotropicStiffnessMatrix=NULL);
  virtual ~AnisotropicForceModel(); 

  virtual void GetInternalForce(double * u, double * internalForces);
  virtual void GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix);
  virtual void GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix); 
protected:
	AnisotropicInternalForces * anisotropic_internal_forces_;
	AnisotropicStiffnessMatrix * anistropic_stiffness_matrix_;
	bool ownStiffnessMatrix;
};

#endif

