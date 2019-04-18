



#include<complex>
using namespace std;



#ifndef TOKAMAK_H_
#define TOKAMAK_H_

/*
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
*/

#include <mpi.h>
#include "parameter.h"
#define PI 3.1415926
#define eps 1e-16


extern double dphi;


/**********Math Function[0]************/

/*
static const int ncof=28;

const double cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
-1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
-1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};

double erfccheb(double z){
	int j;
	double t,ty,tmp,d=0.,dd=0.;
	if (z < 0.) 
		throw("erfccheb requires nonnegative argument");
	t = 2./(2.+z);
	ty = 4.*t - 2.;
	for (j=ncof-1;j>0;j--) {
		tmp = d;
		d = ty*d - dd + cof[j];
		dd = tmp;
	}
	return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
}



inline double erf(double x) {
	if (x >=0.) 
		return 1.0 - erfccheb(x);
	else 
		return erfccheb(-x) - 1.0;
}

*/


template <class T>
class particle_vector{
	template <class T1> friend particle_vector<T1> operator+(const particle_vector<T1>&A,const particle_vector<T1>& B);
	template <class T1> friend particle_vector<T1> operator+(const particle_vector<T1>&A,const T1 &B);
	template <class T1> friend particle_vector<T1> operator*(const T1 A,const particle_vector<T1> &B);

public:
	int dim;
	T p[3];
	particle_vector()
	{
		dim=3;
//		p = new T [dim];
		for(int d=0;d<dim;d++)
			p[d] = 0.0;
	}
	particle_vector(const particle_vector &rhs)
	{
		dim=3;
//		p = new T [dim];
		for(int d=0;d<dim;d++)
			p[d] = rhs.p[d];
	}
	~particle_vector()
	{
//		delete p;
	}
	T& operator[](int i)
	{
		return p[i];
	 }
	const T& operator[](int i)const
	{
		return p[i];
	}
	particle_vector& operator=(const T &rhs)
	{
		for(int d=0;d<dim;d++)
			p[d] = rhs;
		return *this;
	}
	particle_vector& operator+=(const particle_vector &rhs)
	{
		for(int d=0;d<dim;d++)
			p[d] += rhs.p[d];
		return *this;
	}
};



template <class T>
class scalar{
	template <class T1> friend scalar<T1> operator+(const scalar<T1>&A,const scalar<T1>& B);
	template <class T1> friend scalar<T1> operator*(const scalar<T1>& A,const scalar<T1>& B);
	template <class T1> friend scalar<T1> operator*(const scalar<T1>& A,const double& B);
public:
	int row,col,hei;
	T ***p;
    scalar(){}
	scalar(int r,int c,int h)
	{
		row = r;
		col = c;
		hei = h;
		p = new T **[row];
		for(int i=0;i<row;i++)
		{
			p[i] = new T *[col];
			for(int j=0;j<col;j++)
			  p[i][j] = new T[hei];
		}
	  for(int i=0;i<row;i++)
		for(int j=0;j<col;j++)
		  for(int k=0;k<hei;k++)
			p[i][j][k] = 0.0;
	}
	scalar(const scalar &rhs): row(rhs.row), col(rhs.col), hei(rhs.hei)
	{
		p = new T **[row];
		for(int i=0;i<row;i++)
		{
			p[i] = new T *[col];
			for(int j=0;j<col;j++)
				p[i][j] = new T[hei];
		}
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				for(int k=0;k<hei;k++)
					p[i][j][k] = rhs.p[i][j][k];
	}
    scalar<T> initia(int r,int c,int h)
	{
		row = r;
		col = c;
		hei = h;
		p = new T **[row];
		for(int i=0;i<row;i++)
		{
			p[i] = new T *[col];
			for(int j=0;j<col;j++)
			  p[i][j] = new T[hei];
		}
        for(int i=0;i<row;i++)
            for(int j=0;j<col;j++)
                for(int k=0;k<hei;k++)
                    p[i][j][k] = 0.0;
        return *this;
	}
	~scalar()
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				delete [] p[i][j];
		for(int i=0;i<row;i++)
			delete [] p[i];
		delete [] p;
	}
	scalar& operator=(const scalar &rhs)
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				for(int k=0;k<hei;k++)
					p[i][j][k] = rhs.p[i][j][k];
		return *this;
	}
	scalar& operator=(const double &rhs)
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				for(int k=0;k<hei;k++)
					p[i][j][k] = rhs;
		return *this;
	}
	scalar& operator+=(const scalar &rhs)
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				for(int k=0;k<hei;k++)
					p[i][j][k] += rhs.p[i][j][k];
		return *this;
	}
	scalar operator-()
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				for(int k=0;k<hei;k++)
					p[i][j][k] = -p[i][j][k];
		return *this;
	}
	inline T** operator[](const int i) //subscripting: pointer to row i
	{
		return p[i];
	}
	inline const T* const * operator[](const int i) const
	{
		return p[i];
	}
};


template <class T>
class vector{
	template <class T1> friend vector<T1> operator+(const vector<T1>&A,const vector<T1>& B);
	template <class T1> friend scalar<T1> operator*(const vector<T1>& A,const vector<T1>& B);
	template <class T1> friend vector<T1> operator*(const vector<T1>& A,const double& B);
public:
	int dim;
	int row,col,hei;
	T ****p;
    vector(){}
	vector(int r,int c,int h)
	{
		row = r;
		col = c;
		hei = h;
		dim = 3;
		p = new T ***[dim];
		for(int d=0;d<dim;d++)
		{
		p[d] = new T **[row];
		for(int i=0;i<row;i++)
		{
		  p[d][i] = new T *[col];
		  for(int j=0;j<col;j++)
			p[d][i][j] = new T[hei];
		}
		}
		for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
			  p[d][i][j][k] = 0.0;
	}
    vector<T> initia(int r,int c,int h)
	{
		row = r;
		col = c;
		hei = h;
		dim = 3;
		p = new T ***[dim];
		for(int d=0;d<dim;d++)
		{
		    p[d] = new T **[row];
		    for(int i=0;i<row;i++)
		    {
		      p[d][i] = new T *[col];
		      for(int j=0;j<col;j++)
                  p[d][i][j] = new T[hei];
		    }
		}
		for(int d=0;d<dim;d++)
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++)
                    for(int k=0;k<hei;k++)
                        p[d][i][j][k] = 0.0;
        return *this;
	}
	vector(const vector &rhs): row(rhs.row), col(rhs.col), hei(rhs.hei)
	{
		dim = 3;
		p = new T ***[dim];
		for(int d=0;d<dim;d++)
		{
		  p[d] = new T **[row];
		  for(int i=0;i<row;i++)
		  {
			p[d][i] = new T *[col];
			for(int j=0;j<col;j++)
			  p[d][i][j] = new T[hei];
		  }
		}
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
			  p[d][i][j][k] = rhs.p[d][i][j][k];
	}
	vector(const scalar<T> &d1,const scalar<T> &d2,const scalar<T> &d3): row(d1.row), col(d1.col), hei(d1.hei)
	{
		dim = 3;
		p = new T ***[dim];
		for(int d=0;d<dim;d++)
		{
		  p[d] = new T **[row];
		  for(int i=0;i<row;i++)
		  {
			p[d][i] = new T *[col];
			for(int j=0;j<col;j++)
			  p[d][i][j] = new T[hei];
		  }
		}
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
			{
				p[0][i][j][k] = d1.p[i][j][k];
				p[1][i][j][k] = d2.p[i][j][k];
				p[2][i][j][k] = d3.p[i][j][k];
			}
	}
	~vector()
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			delete [] p[d][i][j];
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  delete [] p[d][i];
	  for(int d=0;d<dim;d++)
		  delete [] p[d];
	  delete [] p;
	}
	vector& operator=(const vector &rhs)
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
			  p[d][i][j][k] = rhs.p[d][i][j][k];
	  return *this;
	}
	vector& operator=(const double &rhs)
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
			  p[d][i][j][k] = rhs;
	  return *this;
	}
	vector& operator+=(const vector &rhs)
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
			  p[d][i][j][k] += rhs.p[d][i][j][k];
	  return *this;
	}
	vector operator-()
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
			  p[d][i][j][k] = -p[d][i][j][k];
	  return *this;
	}
	inline T*** operator[](const int i) //subscripting: pointer to row i
	{
		return p[i];
	}


	inline const T* const * const * operator[](const int i) const //subscripting: pointer to row i
	{
		return p[i];
	}
	inline scalar<T> operator()(const int d)const
	{
		scalar<T> compoment(row,col,hei);
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			for(int k=0;k<hei;k++)
				compoment.p[i][j][k] = p[d][i][j][k];
		return compoment;
	}
	inline particle_vector<T> operator()(const int i,const int j,const int k)
	{
	  particle_vector<T> v;
		for(int d=0;d<dim;d++)
		  v.p[d] = p[d][i][j][k];
		return v;
	}

};


template <class T>
particle_vector<T> operator+(const particle_vector<T> &A,const particle_vector<T> &B)
{
        particle_vector<T> sum;
        for(int d=0;d<sum.dim;d++)
             sum[d] = A[d] + B[d];
        return sum;
}

template <class T>
particle_vector<T> operator+(const particle_vector<T>& A,const T& B)
{
	particle_vector<T> sum;
	for(int d=0;d<sum.dim;d++)
		sum[d] = A[d] + B;
	return sum;
}

template <class T>
particle_vector<T> operator-(const particle_vector<T>& A,const particle_vector<T>& B)
{
	particle_vector<T> sub;
	for(int d=0;d<sub.dim;d++)
		sub[d] = A[d] - B[d];
	return sub;
}

template <class T>
particle_vector<T> operator*(const T A,const particle_vector<T> &B)
{
        particle_vector<T> mul;
        for(int d=0;d<3;d++)
             mul.p[d] = A*B.p[d];
        return mul;
}

template <class T>
T operator*(const particle_vector<T> &A,const particle_vector<T> &B)
{
        T mul = 0;
        for(int d=0;d<3;d++)
             mul += A.p[d]*B.p[d];
        return mul;
}



template <class T>
particle_vector<T> operator/(const particle_vector<T> &A,const T &B)
{
        particle_vector<T> div;
        for(int d=0;d<3;d++)
             div.p[d] = A.p[d]/B;
        return div;
}




template <class T>
scalar<T> operator+(const scalar<T>& A,const scalar<T>& B)
{
	scalar<T> sum(A);
	for(int i=0;i<sum.row;i++)
		for(int j=0;j<sum.col;j++)
			for(int k=0;k<sum.hei;k++)
				sum[i][j][k] += B[i][j][k];
	return sum;
}

template <class T>
scalar<T> operator-(const scalar<T>& A,const scalar<T>& B)
{
	scalar<T> sub(A);
	for(int i=0;i<sub.row;i++)
		for(int j=0;j<sub.col;j++)
			for(int k=0;k<sub.hei;k++)
				sub[i][j][k] -= B[i][j][k];
	return sub;
}

template <class T>
scalar<T> operator*(const scalar<T>& A,const scalar<T>& B)
{
	scalar<T> mul(A);
	for(int i=0;i<mul.row;i++)
		for(int j=0;j<mul.col;j++)
			for(int k=0;k<mul.hei;k++)
				mul[i][j][k] *= B[i][j][k];
	return mul;
}


template <class T>
scalar<T> operator*(const T& A,const scalar<T>& B)
{
	scalar<T> mul(B);
	for(int i=0;i<mul.row;i++)
		for(int j=0;j<mul.col;j++)
			for(int k=0;k<mul.hei;k++)
				mul[i][j][k] *= A;
	return mul;
}


template <class T>
scalar<T> operator/(const T& B,const scalar<T>& A)
{
	scalar<T> div(A);
	for(int i=0;i<div.row;i++)
		for(int j=0;j<div.col;j++)
			for(int k=0;k<div.hei;k++)
				div[i][j][k] = B/div[i][j][k];
	return div;
}



template <class T>
vector<T> operator+(const vector<T>& A,const vector<T>& B)
{
	vector<T> sum(A);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<sum.row;i++)
			for(int j=0;j<sum.col;j++)
				for(int k=0;k<sum.hei;k++)
					sum[d][i][j][k] += B[d][i][j][k];
	return sum;
}

template <class T>
vector<T> operator-(const vector<T>& A,const vector<T>& B)
{
	vector<T> sub(A);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<sub.row;i++)
			for(int j=0;j<sub.col;j++)
				for(int k=0;k<sub.hei;k++)
					sub[d][i][j][k] -= B[d][i][j][k];
	return sub;
}

template <class T>
vector<T> operator*(const vector<T>& A,const double& B)
{
	vector<T> mul(A);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<mul.row;i++)
			for(int j=0;j<mul.col;j++)
				for(int k=0;k<mul.hei;k++)
					mul[d][i][j][k] *= B;
	return mul;
}

template <class T>
vector<T> operator*(const T& B,const vector<T>& A)
{
	vector<T> mul(A);
	for(int d=0;d<mul.dim;d++)
		for(int i=0;i<mul.row;i++)
			for(int j=0;j<mul.col;j++)
				for(int k=0;k<mul.hei;k++)
					mul[d][i][j][k] *= B;
	return mul;
}

template <class T>
vector<T> operator*(const scalar<T>& A,const vector<T>& B)
{
	vector<T> mul(B);
	for(int d=0;d<mul.dim;d++)
		for(int i=0;i<mul.row;i++)
			for(int j=0;j<mul.col;j++)
				for(int k=0;k<mul.hei;k++)
					mul[d][i][j][k] *= A[i][j][k];
	return mul;
}

template <class T>
scalar<T> operator*(const vector<T>& A,const vector<T>& B)
{
	scalar<T> mul(A.row,A.col,A.hei);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<mul.row;i++)
			for(int j=0;j<mul.col;j++)
				for(int k=0;k<mul.hei;k++)
					mul[i][j][k] += A[d][i][j][k]*B[d][i][j][k];
	return mul;
}


template <class T>
scalar<T> operator/(const scalar<T>& A,const scalar<T>& B)
{
	scalar<T> div(A.row,A.col,A.hei);
	for(int i=0;i<div.row;i++)
		for(int j=0;j<div.col;j++)
			for(int k=0;k<div.hei;k++)
				div[i][j][k] = A[i][j][k]/B[i][j][k];
	return div;
}


template <class T>
vector<T> operator/(const double& B,const vector<T>& A)
{
	vector<T> div(A);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<div.row;i++)
			for(int j=0;j<div.col;j++)
				for(int k=0;k<div.hei;k++)
					div[d][i][j][k] = B/div[d][i][j][k];
	return div;
}

template <class T>
vector<T> operator/(const vector<T>& A,const T& B)
{
	vector<T> div(A);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<div.row;i++)
			for(int j=0;j<div.col;j++)
				for(int k=0;k<div.hei;k++)
					div[d][i][j][k] = div[d][i][j][k]/B;
	return div;
}

template <class T>
vector<T> operator/(const vector<T>& A,const scalar<T>& B)
{
	vector<T> div(A);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<div.row;i++)
			for(int j=0;j<div.col;j++)
				for(int k=0;k<div.hei;k++)
					div[d][i][j][k] = div[d][i][j][k]/B[i][j][k];
	return div;
}





template <class T>
T* operator+(const T& A,const T& B)
{
        T *sum;
        sum =new T[3];
        for(int d=0;d<3;d++)
             sum[d] = A[d] + B[d];
        return sum;
}













template <class T>
scalar<T> ddR(const scalar<T> &f)
{
	scalar<T> res(f);
	int mx= f.row,my = f.col,mz=f.hei;
//center difference
	for(int i=1;i<mx-1;i++)
		for(int j=0;j<my;j++)
			for(int k=0;k<mz;k++)
				res[i][j][k] = (f[i+1][j][k]-f[i-1][j][k])/(2*dR);
//boudnary free BC dfdx=0
    for(int j=0;j<my;j++)
    	for(int k=0;k<mz;k++)
    	{
			res[0][j][k]    = (f[1][j][k]-f[0][j][k])/dR;
			res[mx-1][j][k] = (f[mx-1][j][k]-f[mx-2][j][k])/dR;
    	}
	return res;
}



template <class T>
scalar<T> ddZ(const scalar<T> &f)
{
	scalar<T> res(f);
	int mx= f.row,my = f.col,mz=f.hei;
//center difference
	for(int i=0;i<mx;i++)
		for(int j=1;j<my-1;j++)
			for(int k=0;k<mz;k++)
				res[i][j][k] = (f[i][j+1][k]-f[i][j-1][k])/(2*dZ);
//boudnary free BC dfdy=0
    for(int i=0;i<mx;i++)
    	for(int k=0;k<mz;k++)
    	{
			res[i][0][k]    = (f[i][1][k]-f[i][0][k])/dZ;
			res[i][my-1][k] = (f[i][my-1][k]-f[i][my-2][k])/dZ;
    	}
	return res;
}


template <class T>
scalar<T> ddphi(const scalar<T> &f)
{
	scalar<T> res(f);
	int mx= f.row,my = f.col,mz=f.hei;
//center difference
	for(int i=0;i<mx;i++)
		for(int j=0;j<my;j++)
			for(int k=1;k<mz-1;k++)
				res[i][j][k] = (f[i][j][k+1]-f[i][j][k-1])/(2*dphi)*1.0/(R0-a*a_b+i*dR);
//period boudnary BC
    for(int i=0;i<mx;i++)
    	for(int j=0;j<my;j++)
    	{
//			res[i][j][0]    = (f[i][j][1]-f[i][j][mz-1])/(2*dphi)*1.0/(R0-a*a_b+i*dR);
//			res[i][j][mz-1] = (f[i][j][0]-f[i][j][mz-2])/(2*dphi)*1.0/(R0-a*a_b+i*dR);
            res[i][j][0]    = (f[i][j][1]-f[i][j][mz-2])/(2*dphi)*1.0/(R0-a*a_b+i*dR);
			res[i][j][mz-1] = (f[i][j][1]-f[i][j][mz-2])/(2*dphi)*1.0/(R0-a*a_b+i*dR);
    	}
	return res;
}

scalar<double> ddphi(const scalar<double> &f,int myid,int numprocs)
{
	MPI::Status status;
	scalar<double> res(f);
	int mx= f.row,my = f.col,mz=f.hei;
//center difference
	int myleft,myright;
	
	if(myid==0)
		myleft = numprocs-1;
	else
		myleft = myid - 1;

	if(myid==numprocs-1)
		myright = 0;
	else
		myright = myid + 1;

		
	for(int i=0;i<mx;i++)
		for(int j=0;j<my;j++)
			for(int k=1;k<mz-1;k++)
				res[i][j][k] = (f[i][j][k+1]-f[i][j][k-1])/(2*dphi)*1.0/(R0-a*a_b+i*dR);
	double left,right;
	//period boudnary BC
	for(int i=0;i<mx;i++)
		for(int j=0;j<my;j++)
		{
			MPI::COMM_WORLD.Sendrecv(&f[i][j][mz-4],1,MPI::DOUBLE,myright,0,&left,1,MPI::DOUBLE,myleft,0,status);
			MPI::COMM_WORLD.Sendrecv(&f[i][j][3],1,MPI::DOUBLE,myleft,1,&right,1,MPI::DOUBLE,myright,1,status);
			res[i][j][0] = (f[i][j][1]-left)/(2*dphi)*1.0/(R0-a*a_b+i*dR);
			res[i][j][mz-1] = (right-f[i][j][mz-2])/(2*dphi)*1.0/(R0-a*a_b+i*dR);
		}
	return res;
}



template <class T>
vector<T> grad(const scalar<T> &f)
{
	vector<T> res(ddR(f),ddZ(f),ddphi(f));
	return res;
}

template <class T>
scalar<T> div(const vector<T> &f)
{
	scalar<T> R(f.row,f.col,f.hei);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				R[i][j][k] = R0-a*a_b+i*dR;
	scalar<T> res(f.row,f.col,f.hei);
	res = ddR(R*f(0))/R + ddZ(f(1)) + ddphi(f(2));
	return res;
}

template <class T>
vector<T> curl(const vector<T> &f)
{
	scalar<T> R(f.row,f.col,f.hei);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				R[i][j][k] = R0-a*a_b+i*dR;
	vector<T> res(ddZ(f(2))-ddphi(f(1)),ddphi(f(0))-ddR(R*f(2))/R,ddR(f(1))-ddZ(f(0)));
	return res;
}


template <class T>
scalar<T> abs(const vector<T> &f)
{
	scalar<T> res(f.row,f.col,f.hei);
	for(int i=0;i<f.row;i++)
		for(int j=0;j<f.col;j++)
			for(int k=0;k<f.hei;k++)
				res[i][j][k] = sqrt(pow(f[0][i][j][k],2.0) + pow(f[1][i][j][k],2.0) + pow(f[2][i][j][k],2.0));
	return res;
}

template <class T>
scalar<T> abs(const scalar<T> &f)
{
	scalar<T> res(f.row,f.col,f.hei);
	for(int i=0;i<f.row;i++)
		for(int j=0;j<f.col;j++)
			for(int k=0;k<f.hei;k++)
				res[i][j][k] = abs(f[i][j][k]);
	return res;
}


template <class T>
vector<T> times(const vector<T> &A,const vector<T> &B)
{
	vector<T> res(A(1)*B(2)-A(2)*B(1),A(2)*B(0)-A(0)*B(2),A(0)*B(1)-A(1)*B(0));
	return res;
}

template <class T>
particle_vector<T> times(const particle_vector<T> &A,const particle_vector<T> &B)
{
	particle_vector<T> res;
	res[0] = (A[1]*B[2]-A[2]*B[1]);
	res[1] = (A[2]*B[0]-A[0]*B[2]);
	res[2] = (A[0]*B[1]-A[1]*B[0]);
	return res;
}



vector<double> grad(const scalar<double> &f,int myid,int numprocs)
{
	vector<double> res(ddR(f),ddZ(f),ddphi(f,myid,numprocs));
	return res;
}

vector<double> grad_cos(const scalar<double> &f,int myid,int numprocs,const scalar<double> &f1)
{
	scalar<double> f2(f1);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				f2[i][j][k] = -n[0]*f1[i][j][k]*1.0/(R0-a*a_b+i*dR);
	vector<double> res(ddR(f),ddZ(f),f2);
	return res;
}

vector<double> grad_sin(const scalar<double> &f,int myid,int numprocs,const scalar<double> &f1)
{
	scalar<double> f2(f1);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				f2[i][j][k] = n[0]*f1[i][j][k]*1.0/(R0-a*a_b+i*dR);
	vector<double> res(ddR(f),ddZ(f),f2);
	return res;
}


vector<double> curl_cos(const vector<double> &f,int myid,int numprocs,const vector<double> &f1)
{
	scalar<double> R(f.row,f.col,f.hei);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				R[i][j][k] = R0-a*a_b+i*dR;
	vector<double> res(ddZ(f(2))-(-n[0]*f1(1)/R),(-n[0]*f1(0)/R)-ddR(R*f(2))/R,ddR(f(1))-ddZ(f(0)));
	return res;
}

vector<double> curl_sin(const vector<double> &f,int myid,int numprocs,const vector<double> &f1)
{
	scalar<double> R(f.row,f.col,f.hei);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				R[i][j][k] = R0-a*a_b+i*dR;
	vector<double> res(ddZ(f(2))-(n[0]*f1(1)/R),(n[0]*f1(0)/R)-ddR(R*f(2))/R,ddR(f(1))-ddZ(f(0)));
	return res;
}


scalar<double> div(const vector<double> &f,int myid,int numprocs)
{
	scalar<double> R(f.row,f.col,f.hei);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				R[i][j][k] = R0-a*a_b+i*dR;
	scalar<double> res(f.row,f.col,f.hei);
	res = ddR(R*f(0))/R + ddZ(f(1)) + ddphi(f(2),myid,numprocs);
	return res;
}


vector<double> curl(const vector<double> &f,int myid,int numprocs)
{
	scalar<double> R(f.row,f.col,f.hei);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				R[i][j][k] = R0-a*a_b+i*dR;
	vector<double> res(ddZ(f(2))-ddphi(f(1),myid,numprocs),ddphi(f(0),myid,numprocs)-ddR(R*f(2))/R,ddR(f(1))-ddZ(f(0)));
	return res;
}



//define particle's property
class Particle{
public:double v_par,v_per,w,f_over_g,g,mu,P_phi,K;
	   int id;
particle_vector<double> X;
  Particle(){
	  X[0] = R0+0.3;
	  X[1] = 0;
	  X[2] = 0;
	  w = 0.0;
  }
};


class Species{
public:double mass,charge;
	   string label;
	   Species(string name):label(name)
	   {
		   if(label == "electron")
		   {
			   mass = 1.0;
			   charge = -1.0;
	 	   }
		   else if(label == "ion")
		   {
			   mass = 1.0;
			   charge = 1.0;
			}
	   }
	   Species(string name,double m):label(name)
	   {
		   if(label == "electron")
		   {
			   mass = 1.0;
			   charge = -1.0;
		   }
		   else if(label == "ion")
		   {
			   mass = 1.0;
			   charge = 1.0;
		   }
		   else
		   {
			   mass = m;
			   charge = 1;
		   }
	   }
};


template <class T>
class scalar2d{
public:
	int row,col;
	T **p;
    scalar2d(){}
	scalar2d(int r,int c)
	{
		row = r;
		col = c;
		p = new T *[row];
		for(int i=0;i<row;i++)
			p[i] = new T[col];
	  for(int i=0;i<row;i++)
		for(int j=0;j<col;j++)
			p[i][j] = 0.0;
	}
    scalar2d<T> initia(int r,int c)
	{
		row = r;
		col = c;
		p = new T *[row];
		for(int i=0;i<row;i++)
			p[i] = new T[col];
	      for(int i=0;i<row;i++)
		    for(int j=0;j<col;j++)
			    p[i][j] = 0.0;
          return *this;
	}
	scalar2d(const scalar2d &rhs): row(rhs.row), col(rhs.col)
	{
		p = new T *[row];
		for(int i=0;i<row;i++)
			p[i] = new T[col];
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				p[i][j] = rhs.p[i][j];
	}
	~scalar2d()
	{
		for(int i=0;i<row;i++)
			delete [] p[i];
		delete [] p;
	}
	scalar2d& operator=(const scalar2d &rhs)
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				p[i][j] = rhs.p[i][j];
		return *this;
	}
	scalar2d& operator=(const double &rhs)
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				p[i][j] = rhs;
		return *this;
	}
	scalar2d& operator+=(const scalar2d &rhs)
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				p[i][j] += rhs.p[i][j];
		return *this;
	}
	scalar2d operator-()
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<col;j++)
				p[i][j] = -p[i][j];
		return *this;
	}
	inline T* operator[](const int i) //subscripting: pointer to row i
	{
		return p[i];
	}
    inline const T* operator[](const int i)const //subscripting: pointer to row i
	{
		return p[i];
	}
};


template <class T>
class vector2d{
public:
	int dim;
	int row,col;
	T ***p;
    vector2d(){}
	vector2d(int r,int c)
	{
		row = r;
		col = c;
		dim = 3;
		p = new T **[dim];
		for(int d=0;d<dim;d++)
		{
			p[d] = new T *[row];
			for(int i=0;i<row;i++)
				p[d][i] = new T[col];
		}
		for(int d=0;d<dim;d++)
			for(int i=0;i<row;i++)
				for(int j=0;j<col;j++)
					p[d][i][j] = 0.0;
	}
    vector2d<T> initia(int r,int c)
	{
		row = r;
		col = c;
		dim = 3;
		p = new T **[dim];
		for(int d=0;d<dim;d++)
		{
			p[d] = new T *[row];
			for(int i=0;i<row;i++)
				p[d][i] = new T[col];
		}
		for(int d=0;d<dim;d++)
			for(int i=0;i<row;i++)
				for(int j=0;j<col;j++)
					p[d][i][j] = 0.0;
        return *this;
	}
	vector2d(const vector2d &rhs): row(rhs.row), col(rhs.col)
	{
		dim = 3;
		p = new T **[dim];
		for(int d=0;d<dim;d++)
		{
			p[d] = new T *[row];
			for(int i=0;i<row;i++)
				p[d][i] = new T[col];
		}
		for(int d=0;d<dim;d++)
			for(int i=0;i<row;i++)
				for(int j=0;j<col;j++)
					p[d][i][j] = rhs.p[d][i][j];
	}
	vector2d(const scalar2d<T> &d1,const scalar2d<T> &d2,const scalar2d<T> &d3): row(d1.row), col(d1.col)
	{
		dim = 3;
		p = new T **[dim];
		for(int d=0;d<dim;d++)
		{
			p[d] = new T *[row];
			for(int i=0;i<row;i++)
				p[d][i] = new T[col];
		}
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			{
				p[0][i][j] = d1.p[i][j];
				p[1][i][j] = d2.p[i][j];
				p[2][i][j] = d3.p[i][j];
			}
	}
	~vector2d()
	{
		for(int d=0;d<dim;d++)
			for(int i=0;i<row;i++)
				delete [] p[d][i];
		for(int d=0;d<dim;d++)
			delete [] p[d];
		delete [] p;
	}
	vector2d& operator=(const vector2d &rhs)
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			  p[d][i][j] = rhs.p[d][i][j];
	  return *this;
	}
	vector2d& operator=(const double &rhs)
	{
		for(int d=0;d<dim;d++)
			for(int i=0;i<row;i++)
				for(int j=0;j<col;j++)
					p[d][i][j] = rhs;
		return *this;
	}
	vector2d& operator+=(const vector2d &rhs)
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			  p[d][i][j] += rhs.p[d][i][j];
	  return *this;
	}
	vector2d operator-()
	{
	  for(int d=0;d<dim;d++)
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			  p[d][i][j] = -p[d][i][j];
	  return *this;
	}
	inline T** operator[](const int i) //subscripting: pointer to row i
	{
		return p[i];
	}


	inline const T* const * operator[](const int i) const //subscripting: pointer to row i
	{
		return p[i];
	}
	inline scalar2d<T> operator()(const int d)const
	{
		scalar2d<T> compoment(row,col);
		for(int i=0;i<row;i++)
		  for(int j=0;j<col;j++)
			compoment.p[i][j] = p[d][i][j];
		return compoment;
	}
};

template <class T>
scalar2d<T> operator*(const vector2d<T>& A,const vector2d<T>& B)
{
	scalar2d<T> mul(A.row,A.col);
	for(int d=0;d<A.dim;d++)
		for(int i=0;i<mul.row;i++)
			for(int j=0;j<mul.col;j++)
				mul[i][j] += A[d][i][j]*B[d][i][j];
	return mul;
}

template <class T>
scalar2d<T> operator-(const scalar2d<T>& A,const scalar2d<T>& B)
{
	scalar2d<T> sub(A.row,A.col);
	for(int i=0;i<sub.row;i++)
		for(int j=0;j<sub.col;j++)
			sub[i][j] = A[i][j] - B[i][j];
	return sub;
}

template <class T>
scalar2d<T> operator-(const scalar2d<T>& A,const T& B)
{
	scalar2d<T> sub(A.row,A.col);
	for(int i=0;i<sub.row;i++)
		for(int j=0;j<sub.col;j++)
			sub[i][j] = A[i][j] - B;
	return sub;
}

template <class T>
scalar2d<T> operator-(const T& A,const scalar2d<T>& B)
{
	scalar2d<T> sub(B.row,B.col);
	for(int i=0;i<sub.row;i++)
		for(int j=0;j<sub.col;j++)
			sub[i][j] = A - B[i][j];
	return sub;
}

template <class T>
scalar2d<T> operator*(const scalar2d<T>& A,const scalar2d<T>& B)
{
	scalar2d<T> mul(A.row,A.col);
	for(int i=0;i<mul.row;i++)
		for(int j=0;j<mul.col;j++)
			mul[i][j] = A[i][j]*B[i][j];
	return mul;
}

template <class T>
scalar2d<T> operator/(const scalar2d<T>& A,const scalar2d<T>& B)
{
	scalar2d<T> div(A.row,A.col);
	for(int i=0;i<A.row;i++)
		for(int j=0;j<A.col;j++)
			div[i][j] = A[i][j]/B[i][j];
	return div;
}

template <class T>
scalar2d<T> operator/(const T A,const scalar2d<T>& B)
{
	scalar2d<T> div(B.row,B.col);
	for(int i=0;i<B.row;i++)
		for(int j=0;j<B.col;j++)
			div[i][j] = A/B[i][j];
	return div;
}

template <class T>
vector2d<T> operator/(const vector2d<T>& A,const scalar2d<T>& B)
{
	vector2d<T> div(A.row,A.col);
	for(int d=0;d<div.dim;d++)
		for(int i=0;i<div.row;i++)
			for(int j=0;j<div.col;j++)
				div[d][i][j] = A[d][i][j]/B[i][j];
	return div;
}


template <class T>
scalar2d<T> abs(const vector2d<T> &f)
{
	scalar2d<T> res(f.row,f.col);
	for(int i=0;i<f.row;i++)
		for(int j=0;j<f.col;j++)
			res[i][j] = sqrt(pow(f[0][i][j],2.0) + pow(f[1][i][j],2.0) + pow(f[2][i][j],2.0));
	return res;
}

template <class T>
scalar2d<T> ddR(const scalar2d<T> &f)
{
	scalar2d<T> res(f);
	int mx= f.row,my = f.col;
//center difference
	for(int i=1;i<mx-1;i++)
		for(int j=0;j<my;j++)
			res[i][j] = (f[i+1][j]-f[i-1][j])/(2*dR);
//boudnary free BC dfdx=0
    for(int j=0;j<my;j++)
   	{
		res[0][j]    = (f[1][j]-f[0][j])/dR;
		res[mx-1][j] = (f[mx-1][j]-f[mx-2][j])/dR;
   	}
	return res;
}

template <class T>
scalar2d<T> ddZ(const scalar2d<T> &f)
{
	scalar2d<T> res(f);
	int mx= f.row,my = f.col;
//center difference
	for(int i=0;i<mx;i++)
		for(int j=1;j<my-1;j++)
			res[i][j] = (f[i][j+1]-f[i][j-1])/(2*dZ);
//boudnary free BC dfdy=0
    for(int i=0;i<mx;i++)
   	{
		res[i][0]    = (f[i][1]-f[i][0])/dZ;
		res[i][my-1] = (f[i][my-1]-f[i][my-2])/dZ;
   	}
	return res;
}

template <class T>
vector2d<T> curl(const vector2d<T> &f)
{
	scalar2d<T> R(f.row,f.col);
	for(int i=0;i<f.row;i++)
	    for(int j=0;j<f.col;j++)
			R[i][j] = R0-a*a_b+i*dR;
	vector2d<T> res(ddZ(f(2))-0.0,0.0-ddR(R*f(2))/R,ddR(f(1))-ddZ(f(0)));
	return res;
}

template <class T>
vector2d<T> grad(const scalar2d<T> &f)
{
    scalar2d<T> zero(f.row,f.col);
    zero = 0.0;
	vector2d<T> res(ddR(f),ddZ(f),zero);
	return res;
}

template <class T>
vector2d<T> grad2d(const scalar2d<T> &f)
{
	vector2d<T> res(ddR(f),ddZ(f),f);
	return res;
}


template <class T>
vector<T> grad2d(const scalar<T> &f)
{
	vector<T> res(ddR(f),ddZ(f),f);
	return res;
}

template <class T>
vector<T> curl2d(const vector<T> &f)
{
	scalar<T> R(f.row,f.col,f.hei);
	for(int i=0;i<mR;i++)
	    for(int j=0;j<mZ;j++)
	    	for(int k=0;k<mphi;k++)
				R[i][j][k] = R0-a*a_b+i*dR;
	vector<T> res(ddZ(f(2))-f(1),f(0)-ddR(R*f(2))/R,ddR(f(1))-ddZ(f(0)));
	return res;
}

#endif
