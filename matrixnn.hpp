#ifndef __MATRIXNN_H__
#define __MATRIXNN_H__
#include "thread_help.hpp"

#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

class MatrixNN{
public:
  MatrixNN(const int &n){
		n_ = n;
		data_.resize(n_*n_);
		std::fill(data_.begin(),data_.end(),0.0);
	}
	
	MatrixNN(const MatrixNN &m){
		n_ = m.n_;
		std::copy(m.data_.begin(), m.data_.end(),data_.begin());
	}
	
	/// Set to identity matrix
	void Identity() {
		for(int r = 0; r < n_; r++){
			for(int c = 0; c < n_; c++){
				if(r==c)(*this)(r,c) = 1.0;
			}
		}	
	}
	/// Fill matrix with random data
	void FillRandom(){
		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(1.0,2.0);
		for(int i = 0; i < data_.size(); i++){
			data_[i] = distribution(generator);		
		}
	}
	
	/// Get size of matrix (returns NxN)
	int get_size()const{return data_.size();}
	
	/// Get number of rows/cols (matrix is symmetrical)
	int get_n()const{return n_;}

	/// Multiply matrices (single thread)
	/*MatrixNN SlowMult(const MatrixNN& m)const{
		MatrixNN tmp(m.get_n());
		const MatrixNN *a = this;
		for(int r = 0; r < n_; r++){
			for(int c = 0; c < n_; c++){
				for(int rm = 0; rm < n_; rm++){
					tmp(r,c) += (*this)(r,rm) * m(rm,c); 
				}
			}
		}
		return tmp;
	}*/

	/// Multiply matrices (concurrent threading)
	MatrixNN QuickMult(const MatrixNN& m)const{
		MatrixNN tmp(m.get_n());
		const MatrixNN *a = this;
		auto myFunc = [&](unsigned int i){
			int r = i/tmp.get_n();
			int c = i-(r*tmp.get_n());
			for(int rm = 0; rm < n_; rm++){
				tmp(r,c) += (*a)(r,rm) * m(rm,c);
			}
		};
		parallelFor(data_.size(),myFunc);
		return tmp;
	}
	
	double determ(){
		return QuickDet((*this),(*this).get_n());
	}

	void minor_mat(MatrixNN& b,const MatrixNN& a,int i,int n){
		int j,l,h=0,k=0;
		for(l=1;l<n;l++)
			for( j=0;j<n;j++){
				if(j == i)
					continue;
				b(h,k) = a(l,j);
				k++;
				if(k == (n-1)){
					h++;
					k=0;
			}
		}
	}

	//	calculate determinate of matrix
	double QuickDet(const MatrixNN& m,int n){
		int i;
		MatrixNN tmp(m.get_n());
		double sum;
		sum=0;
		if (n == 1)
			return m(0,0);
		else if(n == 2)
			return (m(0,0)*m(1,1)-m(0,1)*m(1,0));
		else
			for(i=0;i<n;i++){
				minor_mat(tmp,m,i,n);	// read function
				sum = (double) (sum+m(0,i)*pow(-1,i)*QuickDet(tmp,(n-1)));	// read function	// sum = determinte matrix

			}
		
		//parallelFor(data_.size(),myFunc);
		return sum;
	}

	MatrixNN QuickTrans()const{
		MatrixNN tmp((*this).get_n());
		auto myFunc = [&](unsigned int i){
			for (int j = 0; j < (*this).get_n(); ++j)
				{
					tmp(i,j) = (*this)(j,i);
				}
		};
		//printf("%ld\n",data_.size());
		parallelFor((long int)(*this).get_n(),myFunc);
		return tmp;
	}

	MatrixNN cofactor(const MatrixNN& a, int n,double determinte){
		MatrixNN tmp((*this).get_n());
		MatrixNN b((*this).get_n());
		MatrixNN c((*this).get_n());
		int l,h,m,k,i,j;
		for (h=0;h<n;h++){
			for (l=0;l<n;l++){
				m=0;
				k=0;
				for (i=0;i<n;i++){
					for (j=0;j<n;j++){
						if (i != h && j != l){
							b(m,k)=a(i,j);
							if (k<(n-2))
								k++;
							else{
								k=0;
								m++;
							}
						}
					}
				}
				c(h,l) = (double) pow(-1,(h+l))*QuickDet(b,(n-1));	// c = cofactor Matrix
			}
		}
		tmp = c.QuickTrans();
		return tmp;
	}

	MatrixNN inverse(){
		int n = (*this).get_n();
		MatrixNN tmp((*this).get_n());
		double det = (*this).determ();
		if(det == 0)
			printf("\nInverse of Entered Matrix is not possible\n");
		else if(n == 1)
			tmp(0,0) = 1;
		else

			tmp = cofactor((*this),n,det);

			for (int i=0;i<n;i++){
				for (int j=0;j<n;j++){
					tmp(i,j) = tmp(i,j)/det;					
				}
			}
			return tmp;
	}

	/// Multiply matrices (single thread)
	MatrixNN operator*(const MatrixNN& m)const{
		MatrixNN tmp(m.get_n());
		const MatrixNN *a = this;
		for(int r = 0; r < n_; r++){
			for(int c = 0; c < n_; c++){
				for(int rm = 0; rm < n_; rm++){
					tmp(r,c) += (*this)(r,rm) * m(rm,c); 
				}
			}
		}
		return tmp;
	}
	
	/// Equality operator compares size and cells
	bool operator==(const MatrixNN &rhs)const{
		if(data_.size() != rhs.get_size())return false;
		if(n_ != rhs.get_n())return false;
	
		for(int i = 0; i < data_.size(); i++){
			if(data_[i] != rhs.data_[i])return false;
		}	
		return true;
	}
	
	void operator=(const MatrixNN &m){
		n_ = m.n_;
		data_.resize(n_*n_);
		std::copy(m.data_.begin(), m.data_.end(),data_.begin());
	}
	
	double operator()(const int& row, const int& col)const{
		int index = row * n_ + col;
		return data_[index]; 
	}
	
	double &operator()(const int& row, const int& col){
		int index = row * n_ + col;
		return data_[index]; 
	}
	
	double &operator[](const int& index){
		return data_[index];
	}

	void print_mat()
	{
		for(int i=0;i<(*this).get_n();i++)
		{
			for(int j=0;j<(*this).get_n();j++)
			{
				printf("%f ",(*this)(i,j));
			}
			printf("\n");
		}
	}
private:
	std::vector<double> data_;
	int n_;
};

#endif
