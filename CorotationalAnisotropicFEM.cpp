/**
	corotationalAnistropicFEM.cpp
**/

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "CorotationalAnisotropicFEM.h"
#include "polarDecomposition.h"
#include "tetMesh.h"
#include "volumetricMesh.h"
#include "vec3d.h"
#include "mat3d.h"

using std::string;
using std::vector;
using std::fstream;


CorotationalAnisotropicFEM::CorotationalAnisotropicFEM(TetMesh * tetMesh,bool addGravity, double g):tet_mesh_(tetMesh), gravity_force_(NULL), add_gravity_(addGravity), g_(g)	
{
	ele_num_ = tet_mesh_->getNumElements();
	vert_num_=tet_mesh_->getNumVertices();
	ele_vert_num_=tet_mesh_->getNumElementVertices();

	// build stiffness matrix skeleton 
	SparseMatrix * stiffnessMatrixTopology;
	GetStiffnessMatrixTopology(&stiffnessMatrixTopology);
	// build acceleration indices
	row_ = (int**) malloc (sizeof(int*) * ele_num_);
	column_ = (int**) malloc (sizeof(int*) * ele_num_);
	for (int el=0; el < ele_num_; el++)
	{
		row_[el] = (int*) malloc (sizeof(int) * ele_vert_num_);
		column_[el] = (int*) malloc (sizeof(int) * ele_vert_num_ * ele_vert_num_);
		for(int ver=0; ver<ele_vert_num_; ver++)
			row_[el][ver] = tet_mesh_->getVertexIndex(el, ver);
		// seek for value row[j] in list associated with row[i]
		for(int i=0; i<ele_vert_num_; i++)
			for(int j=0; j<ele_vert_num_; j++)
				column_[el][ele_vert_num_ * i + j] =
				stiffnessMatrixTopology->GetInverseIndex(3*row_[el][i],3*row_[el][j]) / 3;
	}
	undeformedPositions = (double*) malloc (sizeof(double) * 3 * vert_num_);
	for(int i=0; i < vert_num_; i++)
	{
		Vec3d * v = tetMesh->getVertex(i);
		for(int j=0; j<3; j++)
			undeformedPositions[3*i+j] = (*v)[j];
	}
	InitGravity();
	delete(stiffnessMatrixTopology);
}

CorotationalAnisotropicFEM::~CorotationalAnisotropicFEM() 
{
	free(undeformedPositions);
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
	for(int i=0; i<ele_num_; i++)
		free(row_[i]);
	free(row_);

	for(int i=0; i<ele_num_; i++)
		free(column_[i]);
	free(column_);
}

void CorotationalAnisotropicFEM::InitGravity()
{
	if (add_gravity_ && (gravity_force_ == NULL))
	{
		gravity_force_ = (double*) malloc (sizeof(double) * 3 * vert_num_);
		tet_mesh_->computeGravity(gravity_force_, g_);
	}  
}
void CorotationalAnisotropicFEM::GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology)
{
	int * vertices = (int*) malloc (sizeof(int) * ele_vert_num_);
	// build skeleton of sparseMatrix
	SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * vert_num_);
	for (int el=0; el < ele_num_; el++)
	{
		for(int ver=0; ver<ele_vert_num_; ver++)
			vertices[ver] = tet_mesh_->getVertexIndex(el, ver);

		for (int i=0; i<ele_vert_num_; i++)
			for (int j=0; j<ele_vert_num_; j++)
			{
				for(int k=0; k<3; k++)
					for(int l=0; l<3; l++)
					{
						emptyMatrix->AddEntry( 3*vertices[i]+k, 3*vertices[j]+l, 0.0 );						
					}
			}
	}
	*stiffnessMatrixTopology = new SparseMatrix(emptyMatrix);
	delete(emptyMatrix);
	free(vertices);
}

void CorotationalAnisotropicFEM::loadElasticTensorOnCoaseMesh(const string &elastic_tensor_file_name)
{
	if(!tet_mesh_)
	{
		std::cout<<"tet mesh is not exist!\n";
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
	//store the coarse_elastic_tensor for each element
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
}
Mat3d CorotationalAnisotropicFEM::getInistialDisplacementMatrixOnEachElement(unsigned int ele_idx) const
{//Dm
	if(!tet_mesh_)
	{
		std::cout<<"Tet mesh is not exist!\n";
	}
	int * vert_idx;
	Vec3d * vert_pos;
	Vec3d * result_vec;
	vert_idx=(int*)malloc(sizeof(int)*ele_vert_num_);
	vert_pos=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
	result_vec=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
	for(unsigned int idx=0;idx < ele_vert_num_;++idx)
	{
		vert_idx[idx]=tet_mesh_->getVertexIndex(ele_idx,idx);
		vert_pos[idx]=*tet_mesh_->getVertex(vert_idx[idx]);
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

Mat3d CorotationalAnisotropicFEM::getCurrentDisplacementMatrixOnEachElement(const double *current_ele_pos,unsigned int ele_idx) const
{
	if(!tet_mesh_)
	{
		std::cout<<"Tet mesh is not exist!\n";
	}
	Vec3d * result_vec;
	result_vec=(Vec3d*)malloc(sizeof(Vec3d)*ele_vert_num_);
	for(unsigned int dim = 0; dim < 3; ++dim)
	{
		result_vec[dim][0]=current_ele_pos[dim]-current_ele_pos[dim+9];
		result_vec[dim][1]=current_ele_pos[dim+3]-current_ele_pos[dim+9];
		result_vec[dim][2]=current_ele_pos[dim+6]-current_ele_pos[dim+9];
	}
	Mat3d result_matrix(result_vec[0],result_vec[1],result_vec[2]);
	return result_matrix;
}
vector<Mat3d> CorotationalAnisotropicFEM::getCurrentDisplacementMatrixOnAllElements(const double *vert_displacement) const
{
	if(!tet_mesh_)
	{
		std::cout<<"tet mesh is not exist!\n";
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
			vert_idx[idx]=tet_mesh_->getVertexIndex(ele_idx,idx);
			origin_vert_pos[idx]=*tet_mesh_->getVertex(vert_idx[idx]);
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
Mat3d CorotationalAnisotropicFEM::getDeformationGradient(const Mat3d init_dis_matrix, const Mat3d current_dis_matrix) const
{
	Mat3d result(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	result=current_dis_matrix*inv(init_dis_matrix);
	for(int i=0;i<3;++i)
	{
		for(int j=0;j<3;++j)
		{
			if(fabs(result[i][j])<1.0e-7)
			{
				result[i][j]=0;
			}
		}
	}
	if(fabs(det(result))<1.0e-7)
		return (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	else
		return result;
}

Mat3d CorotationalAnisotropicFEM::getDeformationGradientDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim) const
{
	//vert_idx denotes the j-th vertex, vert_idx_dim is the k-th coordinate of the j-th vertex
	//we get the result as dF_i/dx_j^k, which is the derivative force of vertex i to the vertex j on the coordinate k
	//identity vector definition
	if((vert_idx>4)||(vert_idx_dim>3))
	{
		std::cout<<"the vert_idx or the vert_idx_dim is out of range, they should be smaller than 3";
	}
	Mat3d result_matrix(0.0);
	vector<Vec3d> e_vector;
	e_vector.resize(3,Vec3d(0.0,0.0,0.0));
	e_vector[0][0]=1.0;
	e_vector[1][1]=1.0;
	e_vector[2][2]=1.0;
	Mat3d e_matrix(0.0);
	//for the j=1,2,3 we have dF/dx_j^k=e_k*e_j^T*Dm^-1
	for(unsigned int row_idx=0;row_idx<3;++row_idx)
	{
		for(unsigned int col_idx=0;col_idx<3;++col_idx)
		{
			e_matrix[row_idx][col_idx]=e_vector[vert_idx_dim][row_idx]*e_vector[vert_idx][col_idx];
			if(fabs(e_matrix[row_idx][col_idx])<1.0e-7)
				e_matrix[row_idx][col_idx]=0.0;
		}
	}
	Mat3d init_dis_matrix=getInistialDisplacementMatrixOnEachElement(ele_idx);
	Mat3d init_dis_matrix_inv=inv(init_dis_matrix);
	//compute dF/dx_4
	if(vert_idx==3)
	{
		Mat3d vert_cord_matrix(0.0);
		if(vert_idx_dim==0)
		{
			vert_cord_matrix[0][0]=vert_cord_matrix[0][1]=vert_cord_matrix[0][2]=-1.0;
		}
		else if(vert_idx_dim==1)
		{
			vert_cord_matrix[1][0]=vert_cord_matrix[1][1]=vert_cord_matrix[1][2]=-1.0;
		}
		else if(vert_idx_dim==2)
		{
			vert_cord_matrix[2][0]=vert_cord_matrix[2][1]=vert_cord_matrix[2][2]=-1.0;
		}
		else
		{
		}
		result_matrix=vert_cord_matrix*init_dis_matrix_inv;
	}
	else
	{
		//compute dF/dx_j^k,j=1,2,3;k=1,2,3
		result_matrix=e_matrix*init_dis_matrix_inv;
	}
	return result_matrix;
}

double * CorotationalAnisotropicFEM::getCauchyStrainTensor(const Mat3d F) const
{
	Mat3d identity_matrix(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	Mat3d result(0.0);
	result=0.5*(trans(F)+F)-identity_matrix;
	double * result_vector;
	result_vector=(double*)malloc(sizeof(double)*6);	
	for(unsigned int i=0;i<3;++i)
	{
		result_vector[i]=result[i][i];
	}
	result_vector[3]=2.0*result[1][2];
	result_vector[4]=2.0*result[0][2];
	result_vector[5]=2.0*result[0][1];
	for(int i=0;i<6;++i)
		if(fabs(result_vector[i])<1.0e-7)
			result_vector[i]=0;
	return result_vector;
}
double * CorotationalAnisotropicFEM::getGreenStrainTensor(const Mat3d F) const
{
	Mat3d identity_matrix(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	Mat3d temp_matrix(0.0);
	temp_matrix=trans(F)*F;
	temp_matrix-=identity_matrix;
	Mat3d result(0.0);
	result=0.5*temp_matrix;
	double * result_vector;
	result_vector=(double*)malloc(sizeof(double)*6);	
	memset(result_vector, 0, sizeof(double) * 6);
	for(unsigned int i=0;i<3;++i)
	{
		result_vector[i]=result[i][i];		
	}
	result_vector[3]=2.0*result[1][2];
	result_vector[4]=2.0*result[0][2];
	result_vector[5]=2.0*result[0][1];
	for(int i=0;i<6;++i)
		if(fabs(result_vector[i])<1.0e-7)
			result_vector[i]=0;
	return result_vector;
}

double * CorotationalAnisotropicFEM::getCauchyStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const
{
	double * result_vector;
	result_vector=(double*)malloc(sizeof(double)*6);
	const Mat3d F_derivative=getDeformationGradientDerivative(ele_idx,vert_idx,vert_idx_dim);
	Mat3d result_matrix(1.0,0,0,0,1.0,0,0,0,1.0);
	result_matrix=0.5*(trans(F_derivative)+F_derivative);
	for(unsigned int i=0;i<3;++i)
	{
		result_vector[i]=result_matrix[i][i];
	}
	result_vector[3]=2.0*result_matrix[1][2];
	result_vector[4]=2.0*result_matrix[0][2];
	result_vector[5]=2.0*result_matrix[0][1];
	for(int i=0;i<6;++i)
		if(fabs(result_vector[i])<1.0e-7)
			result_vector[i]=0;
	return result_vector;
}
double * CorotationalAnisotropicFEM::getGreenStrainTensorDerivative(unsigned int ele_idx,unsigned int vert_idx,unsigned int vert_idx_dim,const Mat3d F) const
{
	double * result_vector;
	result_vector=(double*)malloc(sizeof(double)*6);
	const Mat3d F_derivative=getDeformationGradientDerivative(ele_idx,vert_idx,vert_idx_dim);
	Mat3d result_matrix(1.0,0,0,0,1.0,0,0,0,1.0);
	result_matrix=0.5*trans(F_derivative)*F;
	result_matrix+=0.5*trans(F)*F_derivative;
	for(unsigned int i=0;i<3;++i)
	{
		result_vector[i]=result_matrix[i][i];
	}
	result_vector[3]=2.0*result_matrix[1][2];
	result_vector[4]=2.0*result_matrix[0][2];
	result_vector[5]=2.0*result_matrix[0][1];
	for(int i=0;i<6;++i)
		if(fabs(result_vector[i])<1.0e-7)
			result_vector[i]=0;
	return result_vector;
}


Mat3d CorotationalAnisotropicFEM::firstPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,CorotationalAnisotropicFEM::StrainType strain_type) const
{
	return F*secondPiolaKirchhoffStress(ele_idx,F,strain_type);
}
Mat3d CorotationalAnisotropicFEM::secondPiolaKirchhoffStress(unsigned int ele_idx,const Mat3d F,CorotationalAnisotropicFEM::StrainType strain_type) const
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
	result_matrix[1][2]=result_matrix[2][1]=0.5*result[3];
	result_matrix[0][2]=result_matrix[2][0]=0.5*result[4];
	result_matrix[0][1]=result_matrix[1][0]=0.5*result[5];
	delete [] strain_tensor_vector;
	free(result);
	return result_matrix;
}
Mat3d CorotationalAnisotropicFEM::firstPiolaKirchhoffStressDerivative(unsigned int ele_idx,const Mat3d F,
					unsigned int vert_idx, unsigned int vert_idx_dim,CorotationalAnisotropicFEM::StrainType strain_type) const
{
	//---dP/dx_j^k=dF/dx_j^k*(C:E)+F(C:dE/dx_j^k);
	//compute C:E and convert a 6*1 vector to a 3*3 matrix
	double * elastic_multiply_strain;
	elastic_multiply_strain=(double*)malloc(sizeof(double)*6);
	memset(elastic_multiply_strain, 0, sizeof(double) * 6);
	const double * strain_tensor_vector;	
	if(strain_type==CAUCHY_STRAIN)
		strain_tensor_vector=getCauchyStrainTensor(F);
	else if(strain_type==GREEN_STRAIN)
		strain_tensor_vector=getGreenStrainTensor(F);
	for(unsigned int i=0;i<6;++i)
	{
		for(unsigned int j=0;j<6;++j)
		{
			elastic_multiply_strain[i]+=ele_elastic_tensor_vector_[ele_idx][i][j]*strain_tensor_vector[j];
		}	
	}
	//convert to 3*3 matrix
	Mat3d elastic_multiply_strain_matrix(0.0);
	for(unsigned int i=0;i<3;++i)
	{
		elastic_multiply_strain_matrix[i][i]=elastic_multiply_strain[i];
	}
	elastic_multiply_strain_matrix[1][2]=elastic_multiply_strain_matrix[2][1]=0.5*elastic_multiply_strain[3];
	elastic_multiply_strain_matrix[0][2]=elastic_multiply_strain_matrix[2][0]=0.5*elastic_multiply_strain[4];
	elastic_multiply_strain_matrix[0][1]=elastic_multiply_strain_matrix[1][0]=0.5*elastic_multiply_strain[5];
	const Mat3d F_derivative=getDeformationGradientDerivative(ele_idx,vert_idx,vert_idx_dim);	
	//compute C*dE/dx_j^k and convert a 6*1 vector to a 3*3 matrix
	double * elastic_multiply_strain_derivative;
	elastic_multiply_strain_derivative=(double*)malloc(sizeof(double)*6);
	memset(elastic_multiply_strain_derivative, 0, sizeof(double) * 6);
	double * strain_tensor_derivative_vector;
	if(strain_type==CAUCHY_STRAIN)
		strain_tensor_derivative_vector=getCauchyStrainTensorDerivative(ele_idx,vert_idx,vert_idx_dim,F);
	else if(strain_type==GREEN_STRAIN)
		strain_tensor_derivative_vector=getGreenStrainTensorDerivative(ele_idx,vert_idx,vert_idx_dim,F);
	//convert to 3*3 matrix
	for(unsigned int i=0;i<6;++i)
	{
		for(unsigned int j=0;j<6;++j)
		{
			elastic_multiply_strain_derivative[i]+=ele_elastic_tensor_vector_[ele_idx][i][j]*strain_tensor_derivative_vector[j];
		}	
	}
	Mat3d elastic_multiply_strain_derivative_matrix(0.0);
	for(unsigned int i=0;i<3;++i)
	{
		elastic_multiply_strain_derivative_matrix[i][i]=elastic_multiply_strain_derivative[i];
	}
	elastic_multiply_strain_derivative_matrix[1][2]=elastic_multiply_strain_derivative_matrix[2][1]=0.5*elastic_multiply_strain_derivative[3];
	elastic_multiply_strain_derivative_matrix[0][2]=elastic_multiply_strain_derivative_matrix[2][0]=0.5*elastic_multiply_strain_derivative[4];
	elastic_multiply_strain_derivative_matrix[0][1]=elastic_multiply_strain_derivative_matrix[1][0]=0.5*elastic_multiply_strain_derivative[5];
	Mat3d result_matrix(0.0);
	result_matrix=F_derivative*elastic_multiply_strain_matrix;
	result_matrix+=F*elastic_multiply_strain_derivative_matrix;
	delete [] elastic_multiply_strain;
	delete [] elastic_multiply_strain_derivative;
	delete [] strain_tensor_vector;
	delete [] strain_tensor_derivative_vector;
	return result_matrix;
}

// compute RK = R * K and RKRT = R * K * R^T (block-wise)
// input: K, R
// output: RK, RKRT
void CorotationalAnisotropicFEM::WarpMatrix(double * K, double * R, double * RK, double * RKRT)
{
	memset(RK, 0, sizeof(double) * 144);
	memset(RKRT, 0, sizeof(double) * 144);
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<4; j++)
		{
			// RK = R * K
			for(int k=0; k<3; k++)
				for(int l=0; l<3; l++)
					for(int m=0; m<3; m++)
						RK[12 * (3 * i + k) + (3 * j + l)] += R[3 * k + m] * K[12 * (3 * i + m) + (3 * j + l)];

			// RKRT = RK * R^T
			for(int k=0; k<3; k++)
				for(int l=0; l<3; l++)
					for(int m=0; m<3; m++)
						RKRT[12 * (3 * i + k) + (3 * j + l)] += RK[12 * (3 * i + k) + (3 * j + m)] * R[3 * l + m];
		}
	}
}
void CorotationalAnisotropicFEM::ComputeStiffnessMatrixOnEachElement(double * small_deformation_pos, unsigned int ele_idx, double * K_ele)
{
	if(K_ele!=NULL)
	{
		memset(K_ele,0,sizeof(double)*144);
	}
	if(!tet_mesh_)
	{
		std::cout<<"Tet mesh is not exist!\n";
		return;
	}
	CorotationalAnisotropicFEM::StrainType strain_type=CAUCHY_STRAIN;
	Mat3d init_dis_matrix=getInistialDisplacementMatrixOnEachElement(ele_idx);
	//Mat3d current_dis_matrix=getCurrentDisplacementMatrixOnEachElement(current_ele_pos,ele_idx);	
	//std::cout<<"current_dis_matrix:"<<current_dis_matrix<<"\n";
	/*for(unsigned int i=0;i<12;++i)
	{
		std::cout<<"current_ele_dis:"<<current_ele_dis[i]<<",";
	}
	std::cout<<"\n";*/
	//Mat3d F=getDeformationGradient(init_dis_matrix,current_dis_matrix);
	//std::cout<<"F:"<<F<<"\n";
	//double F_vector[9];
	//for(unsigned int i=0;i<3;++i)
	//{
	//	for(unsigned int j=0;j<3;++j)
	//	{
	//		F_vector[3*i+j]=F[i][j];
	//		//std::cout<<"F_vector:"<<3*i+j<<"-"<<F[i][j]<<",";
	//	}
	//}

	//double R[9], S[9];
	//double det=PolarDecomposition::Compute(F_vector,R,S);
	//if(det<0)
	//{
	//	for(int i=0;i<9;++i)
	//		R[i]*=-1.0;
	//}
	/*double small_deformation_pos[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	for(unsigned int i=0;i<ele_vert_num_;++i)
	{
		for(unsigned int j=0;j<3;++j)
		{
			for(unsigned int k=0;k<3;++k)
			{
				small_deformation_pos[3*i+j]+=R[3*k+j]*current_ele_pos[3*i+k];
			}			
		}
	}*/
	/*std::cout<<"small_deformation_dis:\n";
	for(unsigned int i=0;i<12;++i)
		std::cout<<small_deformation_pos[i]<<",";
	std::cout<<"\n";*/
	Mat3d small_deformation_dis_matrix=getCurrentDisplacementMatrixOnEachElement(small_deformation_pos,ele_idx);
	Mat3d F_new=getDeformationGradient(init_dis_matrix,small_deformation_dis_matrix);
	double ele_volume=tet_mesh_->getTetVolume(tet_mesh_->getVertex(ele_idx,0),tet_mesh_->getVertex(ele_idx,1),tet_mesh_->getVertex(ele_idx,2),tet_mesh_->getVertex(ele_idx,3));
	int vert[4];
	//----dH/dx_j^k=d[f_0 f_1 f_2]/dx_j^k		
	for(unsigned int j=0;j<ele_vert_num_;++j)
	{
		//compute dH/dx_j^k
		vector<Mat3d> H_derivative_matrix(3);//stores dH/dx_j^0,dH/dx_j^1,dH/dx_j^2	
		H_derivative_matrix.clear();
		for(unsigned int k=0;k<3;++k)
		{
			H_derivative_matrix[k]=ele_volume*firstPiolaKirchhoffStressDerivative(ele_idx,F_new,j,k,strain_type)*inv(trans(init_dis_matrix));										
		}	
		vector<Mat3d> f_derivative_matrix(ele_vert_num_);//stores df0/dx_j,df1/dx_j,df2/dx_j,df3/dx_j	
		f_derivative_matrix.clear();
		f_derivative_matrix[j].set(0.0);
		for(unsigned int row_idx=0;row_idx<3;++row_idx)
		{
			for(unsigned int col_idx=0;col_idx<3;++col_idx)
			{
				f_derivative_matrix[0][row_idx][col_idx]=H_derivative_matrix[col_idx][row_idx][0];
				f_derivative_matrix[1][row_idx][col_idx]=H_derivative_matrix[col_idx][row_idx][1];
				f_derivative_matrix[2][row_idx][col_idx]=H_derivative_matrix[col_idx][row_idx][2];
			}					
		}
		f_derivative_matrix[3]=-(1.0)*(f_derivative_matrix[0]+f_derivative_matrix[1]+f_derivative_matrix[2]);
		for(unsigned int i=0;i<ele_vert_num_;++i)
		{
			for(unsigned int m=0;m<3;++m)
			{
				for(unsigned int n=0;n<3;++n)
				{
					K_ele[12*(3*i+m)+3*j+n]=f_derivative_matrix[i][m][n];
				}
			}
		}
	}
}

void CorotationalAnisotropicFEM::ComputeForceAndStiffnessMatrix(double * u, double * f, SparseMatrix * stiffnessMatrix)
{
	ComputeForceAndStiffnessMatrixOfSubmesh(u, f, stiffnessMatrix, 0, ele_num_);
}
void CorotationalAnisotropicFEM::ComputeForceAndStiffnessMatrixOfSubmesh(double * vertexDisplacements, double * f_int, SparseMatrix * stiffness_matrix, int start_ele, int end_ele)
{
	// clear stiffness matrix to zero
	if (stiffness_matrix != NULL)
		stiffness_matrix->ResetToZero();
	for(unsigned int ele_idx=start_ele;ele_idx < end_ele;++ele_idx)
	{
		double current_ele_pos[12];
		int vert_idx[4];
		Vec3d vert_pos[4];
		for(unsigned int idx=0;idx<ele_vert_num_;++idx)
		{
			vert_idx[idx]=tet_mesh_->getVertexIndex(ele_idx,idx);
			vert_pos[idx]=*tet_mesh_->getVertex(vert_idx[idx]);			
			for(unsigned int j=0;j<3;++j)
			{
				current_ele_pos[3*idx+j]=vertexDisplacements[3*vert_idx[idx]+j]+vert_pos[idx][j];
			}
		}		
		Mat3d init_dis_matrix=getInistialDisplacementMatrixOnEachElement(ele_idx);
		Mat3d current_dis_matrix=getCurrentDisplacementMatrixOnEachElement(current_ele_pos,ele_idx);
		Mat3d F=getDeformationGradient(init_dis_matrix,current_dis_matrix);
		//convert F from 3*3 matrix to 9 length vector 
		double F_vector[9];
		for(unsigned int i=0;i<3;++i)
			for(unsigned int j=0;j<3;++j)
			{
				F_vector[3*i+j]=F[i][j];
			}
		//F=RS
		double R[9], S[9];
		double det=PolarDecomposition::Compute(F_vector,R,S);
		if(det<0)
		{
			for(int i=0;i<9;++i)
				R[i]*=-1.0;
		}
		//remove the rotation of current pos to get a small deformation pos
		double small_deformation_pos[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		for(unsigned int i=0;i<ele_vert_num_;++i)
		{
			for(unsigned int j=0;j<3;++j)
			{
				for(unsigned int k=0;k<3;++k)
				{
					small_deformation_pos[3*i+j]+=R[3*k+j]*current_ele_pos[3*i+k];
				}			
			}
		}
		double K_ele[144];
		ComputeStiffnessMatrixOnEachElement(small_deformation_pos,ele_idx,K_ele);
		// RK = R * K
		// KElement = R * K * R^T
		double RK[144],KElement[144]; // row-major
		WarpMatrix(K_ele, R, RK, KElement);
		if (stiffness_matrix != NULL)
		{
			int * rowIndex = row_[ele_idx];
			int * columnIndex = column_[ele_idx];
			// add KElement to the global stiffness matrix
			for (int i=0; i<4; i++)
			{
				for (int j=0; j<4; j++)
				{
					for(int k=0; k<3; k++)
					{
						for(int l=0; l<3; l++)
						{
							stiffness_matrix->AddEntry(3 * rowIndex[i] + k, 3 * columnIndex[4 * i + j] + l, KElement[12 * (3 * i + k) + 3 * j + l]);
						}
					}
				}	
			}
		}
	}

	if (f_int != NULL)
	{	
		ComputeForces(vertexDisplacements,f_int);
	}	
}

void CorotationalAnisotropicFEM::ComputeForces(const double * vertexDisplacements, double * forces)
{
	//ResetVector(forces);
	if (forces != NULL)
		memset(forces, 0, sizeof(double) * 3 * vert_num_);
	if(!tet_mesh_)
	{
		std::cout<<"Tet mesh is not exist!\n";
		return;
	}
	CorotationalAnisotropicFEM::StrainType strain_type=CAUCHY_STRAIN;
	vector<Mat3d> current_dis_matrix;
	current_dis_matrix.resize(ele_num_);
	//set all the values of the current_dis_matrix[]=0
	for(unsigned int i=0;i<ele_num_;++i)
		current_dis_matrix[i].set(0.0);
	current_dis_matrix=getCurrentDisplacementMatrixOnAllElements(vertexDisplacements);
	for(unsigned int ele_idx = 0; ele_idx <ele_num_; ++ele_idx)
	{
		int vert_global_idx[4];
		Vec3d vert_pos[4];
		double current_ele_pos[12];
		//vert_global_idx=(int*)malloc(sizeof(int)*ele_vert_num_);
		for(unsigned int k=0;k<ele_vert_num_;++k)
		{
			vert_global_idx[k] = tet_mesh_->getVertexIndex(ele_idx,k);
			vert_pos[k]=*tet_mesh_->getVertex(vert_global_idx[k]);
			for(unsigned int j=0;j<3;++j)
				current_ele_pos[3*k+j]=vertexDisplacements[3*vert_global_idx[k]+j]+vert_pos[k][j];
		}
		double ele_volume=tet_mesh_->getTetVolume(tet_mesh_->getVertex(ele_idx,0),tet_mesh_->getVertex(ele_idx,1),tet_mesh_->getVertex(ele_idx,2),tet_mesh_->getVertex(ele_idx,3));		

		Mat3d init_dis_matrix=getInistialDisplacementMatrixOnEachElement(ele_idx);
		Mat3d F=getDeformationGradient(init_dis_matrix,current_dis_matrix[ele_idx]);	
		
		double F_vector[9];
		for(unsigned int i=0;i<3;++i)
		{
			for(unsigned int j=0;j<3;++j)
			{
				F_vector[3*i+j]=F[i][j];
			}
		}
		double R[9], S[9];
		double det=PolarDecomposition::Compute(F_vector,R,S);
		if(det<0)
		{
			for(int i=0;i<9;++i)
				R[i]*=-1.0;
		}
		//remove the rotation part of the current position 
		double small_deformation_pos[12]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		for(unsigned int i=0;i<4;++i)
		{
			for(unsigned int j=0;j<3;++j)
			{
				for(unsigned int k=0;k<3;++k)
				{
					small_deformation_pos[3*i+j]+=R[3*k+j]*current_ele_pos[3*i+k];
				}			
			}
		}
		//convert Rotation vector to 3*3 matrix
		Mat3d R_matrix;
		for(unsigned int i=0;i<3;++i)
		{
			for(unsigned int j=0;j<3;++j)
			{
				R_matrix[i][j]=R[3*i+j];
			}
		}
		Mat3d small_deformation_dis_matrix=getCurrentDisplacementMatrixOnEachElement(small_deformation_pos,ele_idx);
		//compute deformation gradient on small deformation
		Mat3d F_new=getDeformationGradient(init_dis_matrix,small_deformation_dis_matrix);
		//add the rotation part to the force matrix
		Mat3d force_matrix=ele_volume*R_matrix*firstPiolaKirchhoffStress(ele_idx,F_new,strain_type)*inv(trans(init_dis_matrix));	
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
		}

	}
	if (add_gravity_)
	{
		for(int i=0; i<3*vert_num_; i++)
		{
			forces[i] -= gravity_force_[i];			
		}
	}
}
