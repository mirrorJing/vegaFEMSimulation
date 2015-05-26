/**
	AnistropicInternalForce.cpp
**/

#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "AnisotropicInternalForces.h"
#include "volumetricMesh.h"
#include "vec3d.h"
#include "mat3d.h"
#include "tetMesh.h"
#include "sparseMatrix.h"

using namespace std;
using std::fstream;
using std::vector;

AnisotropicInternalForces::AnisotropicInternalForces(VolumetricMesh * volumetricMesh,bool addGravity, double g):volumetric_mesh_(volumetricMesh), gravity_force_(NULL), add_gravity_(addGravity), g_(g)
{	
	tet_mesh_=(TetMesh*)volumetric_mesh_;
	vert_num_=volumetric_mesh_->getNumVertices();
	ele_num_=volumetric_mesh_->getNumElements();
	ele_vert_num_=volumetric_mesh_->getNumElementVertices();
	InitGravity();
}

AnisotropicInternalForces::~AnisotropicInternalForces() 
{
	free(gravity_force_);
	if(ele_elastic_tensor_vector_)
	{
		for(int i=0;i<ele_num_;++i)
			if(ele_elastic_tensor_vector_[i])
			{
				for(int j=0;j<6;++j)
					if(ele_elastic_tensor_vector_[i][j])					
						delete [] ele_elastic_tensor_vector_[i][j];
				delete [] ele_elastic_tensor_vector_[i];
			}
			delete [] ele_elastic_tensor_vector_;
	}
}

void AnisotropicInternalForces::InitGravity()
{
	if (add_gravity_ && (gravity_force_ == NULL))
	{
		gravity_force_ = (double*) malloc (sizeof(double) * 3 * vert_num_);
		volumetric_mesh_->computeGravity(gravity_force_, g_);
	}  
}

int AnisotropicInternalForces::ModifiedSVD(Mat3d & F, Mat3d & U, Vec3d & Fhat, Mat3d & V) const
{
    // The code handles the following necessary special situations (see the code below) :

    //---------------------------------------------------------
    // 1. det(V) == -1
    //    - simply multiply a column of V by -1
    //---------------------------------------------------------
    // 2. An entry of Fhat is near zero
    //---------------------------------------------------------
    // 3. Tet is inverted.
    //    - check if det(U) == -1
    //    - If yes, then negate the minimal element of Fhat
    //      and the corresponding column of U
    //---------------------------------------------------------

    double modifiedSVD_singularValue_eps = 1e-8;

    // form F^T F and do eigendecomposition
    Mat3d normalEq = trans(F) * F;
    Vec3d eigenValues;
    Vec3d eigenVectors[3];

    // note that normalEq is changed after calling eigen_sym
    eigen_sym(normalEq, eigenValues, eigenVectors);

    V.set(eigenVectors[0][0], eigenVectors[1][0], eigenVectors[2][0],
          eigenVectors[0][1], eigenVectors[1][1], eigenVectors[2][1],
          eigenVectors[0][2], eigenVectors[1][2], eigenVectors[2][2]);
    /*
      printf("--- original V ---\n");
      V.print();
      printf("--- eigenValues ---\n");
      printf("%G %G %G\n", eigenValues[0], eigenValues[1], eigenValues[2]);
    */

    // Handle situation:
    // 1. det(V) == -1
    //    - simply multiply a column of V by -1
    if (det(V) < 0.0)
    {
        // convert V into a rotation (multiply column 1 by -1)
        V[0][0] *= -1.0;
        V[1][0] *= -1.0;
        V[2][0] *= -1.0;
    }

    Fhat[0] = (eigenValues[0] > 0.0) ? sqrt(eigenValues[0]) : 0.0;
    Fhat[1] = (eigenValues[1] > 0.0) ? sqrt(eigenValues[1]) : 0.0;
    Fhat[2] = (eigenValues[2] > 0.0) ? sqrt(eigenValues[2]) : 0.0;

    //printf("--- Fhat ---\n");
    //printf("%G %G %G\n", Fhat[0][0], Fhat[1][1], Fhat[2][2]);

    // compute inverse of singular values
    // also check if singular values are close to zero
    Vec3d FhatInverse;
    FhatInverse[0] = (Fhat[0] > modifiedSVD_singularValue_eps) ? (1.0 / Fhat[0]) : 0.0;
    FhatInverse[1] = (Fhat[1] > modifiedSVD_singularValue_eps) ? (1.0 / Fhat[1]) : 0.0;
    FhatInverse[2] = (Fhat[2] > modifiedSVD_singularValue_eps) ? (1.0 / Fhat[2]) : 0.0;
  
    // compute U using the formula:
    // U = F * V * diag(FhatInverse)
    U = F * V;
    U.multiplyDiagRight(FhatInverse);

    // In theory, U is now orthonormal, U^T U = U U^T = I .. it may be a rotation or a reflection, depending on F.
    // But in practice, if singular values are small or zero, it may not be orthonormal, so we need to fix it.
    // Handle situation:
    // 2. An entry of Fhat is near zero
    // ---------------------------------------------------------

    /*
      printf("--- FhatInverse ---\n");
      FhatInverse.print();
      printf(" --- U ---\n");
      U.print();
    */
  
    if ((Fhat[0] < modifiedSVD_singularValue_eps) && (Fhat[1] < modifiedSVD_singularValue_eps) && (Fhat[2] < modifiedSVD_singularValue_eps))
    {
        // extreme case, material has collapsed almost to a point
        // see [Irving 04], p. 4
        U.set(1.0, 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, 1.0);
    }
    else 
    {
        int done = 0;
        for(int dim=0; dim<3; dim++)
        {
            int dimA = dim;
            int dimB = (dim + 1) % 3;
            int dimC = (dim + 2) % 3;
            if ((Fhat[dimB] < modifiedSVD_singularValue_eps) && (Fhat[dimC] < modifiedSVD_singularValue_eps))
            {
                // only the column dimA can be trusted, columns dimB and dimC correspond to tiny singular values
                Vec3d tmpVec1(U[0][dimA], U[1][dimA], U[2][dimA]); // column dimA
                Vec3d tmpVec2;
                FindOrthonormalVector(tmpVec1, tmpVec2);
                Vec3d tmpVec3 = norm(cross(tmpVec1, tmpVec2));
                U[0][dimB] = tmpVec2[0];
                U[1][dimB] = tmpVec2[1];
                U[2][dimB] = tmpVec2[2];
                U[0][dimC] = tmpVec3[0];
                U[1][dimC] = tmpVec3[1];
                U[2][dimC] = tmpVec3[2];
                if (det(U) < 0.0)
                {
                    U[0][dimB] *= -1.0;
                    U[1][dimB] *= -1.0;
                    U[2][dimB] *= -1.0;
                }
                done = 1;
                break; // out of for
            }
        }

        if (!done) 
        {
            for(int dim=0; dim<3; dim++)
            {
                int dimA = dim;
                int dimB = (dim + 1) % 3;
                int dimC = (dim + 2) % 3;

                if (Fhat[dimA] < modifiedSVD_singularValue_eps)
                {
                    // columns dimB and dimC are both good, but column dimA corresponds to a tiny singular value
                    Vec3d tmpVec1(U[0][dimB], U[1][dimB], U[2][dimB]); // column dimB
                    Vec3d tmpVec2(U[0][dimC], U[1][dimC], U[2][dimC]); // column dimC
                    Vec3d tmpVec3 = norm(cross(tmpVec1, tmpVec2));
                    U[0][dimA] = tmpVec3[0];
                    U[1][dimA] = tmpVec3[1];
                    U[2][dimA] = tmpVec3[2];
                    if (det(U) < 0.0)
                    {
                        U[0][dimA] *= -1.0;
                        U[1][dimA] *= -1.0;
                        U[2][dimA] *= -1.0;
                    }
                    done = 1;
                    break; // out of for
                }
            }
        }

        if (!done)
        {
            // Handle situation:
            // 3. Tet is inverted.
            //    - check if det(U) == -1
            //    - If yes, then negate the minimal element of Fhat
            //      and the corresponding column of U

            double detU = det(U);
            if (detU < 0.0)
            {
                // tet is inverted
                // find smallest singular value (they are all non-negative)
                int smallestSingularValueIndex = 0;
                for(int dim=1; dim<3; dim++)
                    if (Fhat[dim] < Fhat[smallestSingularValueIndex])
                        smallestSingularValueIndex = dim;

                // negate smallest singular value
                Fhat[smallestSingularValueIndex] *= -1.0;
                U[0][smallestSingularValueIndex] *= -1.0;
                U[1][smallestSingularValueIndex] *= -1.0;
                U[2][smallestSingularValueIndex] *= -1.0;
            }
        }
    }

    /*
      printf("U = \n");
      U.print();
      printf("Fhat = \n");
      Fhat.print();
      printf("V = \n");
      V.print();
    */

    return 0;
}

void AnisotropicInternalForces::FindOrthonormalVector(Vec3d & v, Vec3d & result) const
{
    // find smallest abs component of v
    int smallestIndex = 0;
    for(int dim=1; dim<3; dim++)
        if (fabs(v[dim]) < fabs(v[smallestIndex]))
            smallestIndex = dim;

    Vec3d axis(0.0, 0.0, 0.0);
    axis[smallestIndex] = 1.0;

    // this cross-product will be non-zero (as long as v is not zero)
    result = norm(cross(v, axis));
}

void AnisotropicInternalForces::loadElasticTensorOnCoaseMesh(const string &elastic_tensor_file_name)
{
	if(!volumetric_mesh_)
	{
		std::cout<<"volumetric mesh is not exist!\n";
		return;
	}
	fstream rfile;
	rfile.open(elastic_tensor_file_name,std::ios::in);
	if(!rfile)
	{
		std::cout<<"Can't open file "<<elastic_tensor_file_name<<"!\n";
		return;
	}
	unsigned int str_num=0;
	unsigned int line_num=0;
	unsigned int matrix_row=0;
	unsigned int matrix_col=0;	
	vector<vector<double> > elastic_tensor_vector;
	elastic_tensor_vector.resize(ele_num_);	
	unsigned int col_num=36;
	while((!rfile.eof())&&(rfile.peek()!=std::ifstream::traits_type::eof()))
	{
		double temp_value;
		rfile>>temp_value;
		str_num++;	
		elastic_tensor_vector[line_num].push_back(temp_value);
		if((str_num % col_num == 0)&&(str_num >= col_num))
		{
			line_num++;
		}
		if(str_num>=(ele_num_*col_num))
			break;	
	}
	rfile.close();
	//store the elastic_tensor for each element
	ele_elastic_tensor_vector_=(double***)malloc(sizeof(double)*ele_num_*6*6);
	for(unsigned int ele_idx = 0; ele_idx < ele_num_; ++ele_idx)
	{
		ele_elastic_tensor_vector_[ele_idx] = (double**)malloc(sizeof(double)*6*6);
		unsigned int num=0;
		for(unsigned int row_idx = 0; row_idx < 6; ++row_idx)
		{
			ele_elastic_tensor_vector_[ele_idx][row_idx] = (double*)malloc(sizeof(double)*6);
			for(unsigned int col_idx = 0; col_idx < 6; ++col_idx)
			{
				num=6*row_idx+col_idx;			
				ele_elastic_tensor_vector_[ele_idx][row_idx][col_idx]=elastic_tensor_vector[ele_idx][num];
			}				
		}
	}
	std::cout<<"------------------------------------------------";
}
Mat3d AnisotropicInternalForces::getInitialDisplacementMatrixOnEachElement(unsigned int ele_idx) const//Dm
{
	if(!volumetric_mesh_)
	{
		std::cout<<"Volumetric mesh is not exist!\n";
	}
	int * vert_idx;
	Vec3d * vert_pos;
	Vec3d * result_vec;
	vert_idx=(int*)malloc(sizeof(int)*ele_vert_num_);
	vert_pos=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
	result_vec=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
	
	for(unsigned int idx=0;idx < ele_vert_num_;++idx)
	{
		vert_idx[idx]=volumetric_mesh_->getVertexIndex(ele_idx,idx); //get global vertex index
		vert_pos[idx]=*volumetric_mesh_->getVertex(vert_idx[idx]);
	}
	for(unsigned int dim = 0; dim < 3; ++dim)
	{
		result_vec[dim][0]=vert_pos[0][dim]-vert_pos[3][dim];
		result_vec[dim][1]=vert_pos[1][dim]-vert_pos[3][dim];
		result_vec[dim][2]=vert_pos[2][dim]-vert_pos[3][dim];
	}
	Mat3d result_matrix(result_vec[0],result_vec[1],result_vec[2]);
	delete [] vert_idx;
	delete [] vert_pos;
	delete [] result_vec;
	return result_matrix;
}
vector<Mat3d> AnisotropicInternalForces::getCurrentDisplacementMatrixOnAllElements(const double *vert_displacement) const
{
	if(!volumetric_mesh_)
	{
		std::cout<<"Volumetric mesh is not exist!\n";
	}
	vector<Mat3d> result_matrix_vector(ele_num_);
	for(unsigned int ele_idx = 0; ele_idx < ele_num_; ++ele_idx)
	{
		int * vert_idx;
		Vec3d * vert_pos;
		Vec3d * origin_vert_pos;
		Vec3d * result_vec;
		vert_idx=(int*)malloc(sizeof(int)*ele_vert_num_);
		vert_pos=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
		origin_vert_pos=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
		result_vec=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
		
		for(unsigned int idx=0;idx < ele_vert_num_;++idx)
		{
			vert_idx[idx]=volumetric_mesh_->getVertexIndex(ele_idx,idx);
			origin_vert_pos[idx]=*volumetric_mesh_->getVertex(vert_idx[idx]);
			for(unsigned int i=0;i<3;++i)
			{
				vert_pos[idx][i]=vert_displacement[3*vert_idx[idx]+i]+origin_vert_pos[idx][i];
			}	
		}
		for(unsigned int dim = 0; dim < 3; ++dim)
		{
			result_vec[dim][0]=vert_pos[0][dim]-vert_pos[3][dim];
			result_vec[dim][1]=vert_pos[1][dim]-vert_pos[3][dim];
			result_vec[dim][2]=vert_pos[2][dim]-vert_pos[3][dim];
		}
		Mat3d result_matrix(result_vec[0],result_vec[1],result_vec[2]);
		result_matrix_vector[ele_idx]=result_matrix;
		delete [] vert_idx;
		delete [] vert_pos;
		delete [] origin_vert_pos;
		delete [] result_vec;
	}	
	return result_matrix_vector;

}
Mat3d AnisotropicInternalForces::getDeformationGradient(const Mat3d init_dis_matrix, const Mat3d current_dis_matrix) const
{
	//the initial F is identity
	Mat3d result(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	result=current_dis_matrix*inv(init_dis_matrix);

    //handle inversion
    Mat3d U,V;
    Vec3d Fhat;
    ModifiedSVD(result,U,Fhat,V);
    //clamphat if below the principle stretch threshold
    double principle_threshold = 0.6;
    for(unsigned int i = 0; i < 3 ; ++i)
        if(Fhat[i] < principle_threshold)
            Fhat[i] = principle_threshold;
    Mat3d Fhat_mat;
    for(unsigned int i = 0; i < 3; ++i)
        for(unsigned int j = 0; j < 3; ++j)
            Fhat_mat[i][j] = (i==j)?Fhat[i]:0;
    result = U*Fhat_mat*trans(V);
    return result;
}
double * AnisotropicInternalForces::getCauchyStrainTensor(const Mat3d F) const
{
	Mat3d identity_matrix(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	Mat3d result_matrix(0.0);
	result_matrix=0.5*(trans(F)+F)-identity_matrix;
	double * result_vector;
	result_vector=(double*)malloc(sizeof(double)*6);	
	for(unsigned int i=0;i<3;++i)
	{
		result_vector[i]=result_matrix[i][i];
	}
	result_vector[3]=2.0*result_matrix[1][2];
	result_vector[4]=2.0*result_matrix[0][2];
	result_vector[5]=2.0*result_matrix[0][1];
	return result_vector;
}
double * AnisotropicInternalForces::getGreenStrainTensor(const Mat3d F) const
{
	Mat3d identity_matrix(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	Mat3d result_matrix(0.0);
	result_matrix=trans(F)*F;
	result_matrix-=identity_matrix;
	result_matrix=0.5*result_matrix;
	double * result_vector;
	result_vector=(double*)malloc(sizeof(double)*6);	
	memset(result_vector, 0, sizeof(double) * 6);
	for(unsigned int i=0;i<3;++i)
	{
		result_vector[i]=result_matrix[i][i];		
	}
	result_vector[3]=2.0*result_matrix[1][2];
	result_vector[4]=2.0*result_matrix[0][2];
	result_vector[5]=2.0*result_matrix[0][1];
	return result_vector;
} 
Mat3d AnisotropicInternalForces::firstPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,AnisotropicInternalForces::StrainType strain_type) const
{
	return F*secondPiolaKirchhoffStress(ele_idx,F,strain_type);
}
Mat3d AnisotropicInternalForces::secondPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,AnisotropicInternalForces::StrainType strain_type) const
{
	double * result;
	result=(double*)malloc(sizeof(double)*6);
	memset(result, 0, sizeof(double) * 6);
	const double * strain_tensor_vector;
	if(strain_type==CAUCHY_STRAIN)
		strain_tensor_vector=getCauchyStrainTensor(F);
	else if(strain_type==GREEN_STRAIN)
		strain_tensor_vector=getGreenStrainTensor(F);
	for(unsigned int i=0;i<6;++i)
	{
		for(unsigned int j=0;j<6;++j)
		{
			result[i]+=ele_elastic_tensor_vector_[ele_idx][i][j]*strain_tensor_vector[j];
		}	
	}
	Mat3d result_matrix(0.0);
	for(unsigned int i=0;i<3;++i)
	{
		result_matrix[i][i]=result[i];
	}
	result_matrix[1][2]=result_matrix[2][1]=result[3];
	result_matrix[0][2]=result_matrix[2][0]=result[4];
	result_matrix[0][1]=result_matrix[1][0]=result[5];
	delete [] strain_tensor_vector;
	free(result);
	return result_matrix;
}

void AnisotropicInternalForces::ResetVector(double * vec)
{
	memset(vec, 0, sizeof(double) * 3 * vert_num_);
}

void AnisotropicInternalForces::ComputeForces(const double * vertexDisplacements, double * forces)
{
	ResetVector(forces);
	if(!tet_mesh_)
	{
		std::cout<<"Tet mesh is not exist!\n";
		return;
	}
	AnisotropicInternalForces::StrainType strain_type=GREEN_STRAIN;
	vector<Mat3d> current_dis_matrix;
	current_dis_matrix.resize(ele_num_);
	//set all the values of the current_dis_matrix[]=0
	for(unsigned int i=0;i<ele_num_;++i)
		current_dis_matrix[i].set(0.0);
	current_dis_matrix=getCurrentDisplacementMatrixOnAllElements(vertexDisplacements);
	//std::cout<<"forces:"<<"\n";
	for(unsigned int ele_idx = 0; ele_idx <ele_num_; ++ele_idx)
	{
		int vert_global_idx[4];
		//vert_global_idx=(int*)malloc(sizeof(int)*ele_vert_num_);
		for(unsigned int k=0;k<ele_vert_num_;++k)
		{
			vert_global_idx[k] = volumetric_mesh_->getVertexIndex(ele_idx,k);
		}
		double ele_volume=tet_mesh_->getTetVolume(tet_mesh_->getVertex(ele_idx,0),tet_mesh_->getVertex(ele_idx,1),tet_mesh_->getVertex(ele_idx,2),tet_mesh_->getVertex(ele_idx,3));		
		
		Mat3d init_dis_matrix=getInitialDisplacementMatrixOnEachElement(ele_idx);
		Mat3d F=getDeformationGradient(init_dis_matrix,current_dis_matrix[ele_idx]);	
		Mat3d force_matrix=ele_volume*firstPiolaKirchhoffStress(ele_idx,F,strain_type)*inv(trans(init_dis_matrix));	
		for(int i =0 ; i < ele_vert_num_; ++i)
		{
			for(int j=0;j<3;++j)//Dim
			{
				if(i == 3)//the force of the 4-th vertex is computed by the other three vertices of the same element
				{
					forces[3*vert_global_idx[i]+j]+=(-1.0)*(force_matrix[j][0]+force_matrix[j][1]+force_matrix[j][2]);
				}
				else
				{
					forces[3*vert_global_idx[i]+j]+=force_matrix[j][i];
				}
			}	
		//	std::cout<<forces[3*vert_global_idx[i]]<<","<<forces[3*vert_global_idx[i]]<<","<<forces[3*vert_global_idx[i]+1]<<","<<forces[3*vert_global_idx[i]+2]<<",";
		}
	//free(vert_global_idx);
		
	}
	//std::cout<<"\n";
//	std::cout<<"f[390]="<<forces[390]<<","<<forces[391]<<","<<forces[392]<<"\n";
	if (add_gravity_)
	{
		for(int i=0; i<3*vert_num_; i++)
		{
			forces[i] -= gravity_force_[i];			
		}
	}
//	std::cout<<"f+g[3591]="<<forces[3591]<<","<<forces[3592]<<","<<forces[3593]<<"\n";
}
