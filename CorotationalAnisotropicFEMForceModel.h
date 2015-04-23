/*
	Force model corresponding to the corotational anisotropic material.
*/

#ifndef _CorotationalAnisotropicFEMForceModel_H_
#define _CorotationalAnisotropicFEMForceModel_H_

#include "forceModel.h"
#include "CorotationalAnisotropicFEM.h"

class CorotationalAnisotropicFEMForceModel : virtual public ForceModel
{
public:
  CorotationalAnisotropicFEMForceModel(CorotationalAnisotropicFEM * corotationalAnisotropicFEM);
  virtual ~CorotationalAnisotropicFEMForceModel(); 

  virtual void GetInternalForce(double * u, double * internalForces);
  virtual void GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix);
  virtual void GetTangentStiffnessMatrix(double * u, SparseMatrix * tangentStiffnessMatrix); 
  virtual void GetForceAndMatrix(double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix);
protected:
	CorotationalAnisotropicFEM * corotational_anisotropic_fem_;
	bool ownStiffnessMatrix;
};

#endif

