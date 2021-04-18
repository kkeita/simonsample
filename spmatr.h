#ifndef SPMATR_H
#define SPMATR_H

#include <unordered_map>
#include <utility>
#include "mkl.h"  //For testing (timer function dsecnd() )
#include "mkl_spblas.h"
#include <iostream>  //For testing (print_matrix member function)
#include <iomanip>   // ''       ''
#include <cstring>
#include <stdint.h>

//Includes for sparsehash library, you should have the src directory inside the sparsehash package as an include directory in project properties
#ifdef _WIN32
#include <port.cc>
#include <dense_hash_map>
#else
#include <sparsehash/dense_hash_map>
#endif

using google::dense_hash_map;


/**
Function object for the hash function of SpMatr_hash
*/
struct hash_int_pair : std::unary_function<uint64_t, size_t>
{
public:
	std::hash<int> hash_fn; //Default hash function for int
	size_t custom_hash_fn(const unsigned int input) const;
	size_t operator()(const uint64_t index_pair) const; //index_pair passed as constant reference	
};

/**
Sparse matrix class built on a hash table for element insertion and access.
Operator += should be used during assembly to add a value to an existing element (it can also be used to insert a new element)
The assignment operator (ex: my_matrix(i,j) = 1.0) can be used to add a new element.
It can also be used to add a value to an existing element (ex: my_matrix(i,j) = my_matrix(i,j) + 1.0), but this is slower than using the += operator.
*/
class SpMatr_hash
{
public:	
	SpMatr_hash();  //Default constructor, generates empty matrix with 0 rows, 0 columns, only providing a default constructor...
	SpMatr_hash(int nb_rows, int nb_cols, int expected_nnz); //Constructor for an empty matrix, takes as arguments the matrix dimensions and the expected number of non-zeros 
															//(approximate or exact is fine, fixes hash table's nb. of buckets)
	void print_matrix(); //For testing purposes
	double& operator()(uint64_t i,uint64_t j); //Overloading of operator (int,int) to access matrix element
	int get_nnz();  //Get number of non-zeros in the matrix
	int get_nb_rows();  //Returns nb_rows
	int get_nb_cols();  //Returns nb_cols
	
	typedef std::unordered_map<uint64_t, double, hash_int_pair>::iterator iterator; //Define iterator type for the unordered_map val
	iterator begin();  //Return iterator corresponding to .begin() of unordered_map val
	iterator end();  //Same thing for .end()
	unsigned int first_int(iterator& iter);  //Returns first int corresponding to the key of the mapping in ordered_map pointed to by iterator argument
	unsigned int second_int(iterator& iter); //Same as above for second int
	//THIS IS NOT NEEDED SpMatr_CSR* convert_to_csr();  //Converts matrix to a SpMatr_CSR object and returns a pointer to it	
	std::unordered_map<uint64_t, double, hash_int_pair> val;  //Hash table for storing the matrix entries
private:
	int nb_rows;  //Number of rows in the matrix
	int nb_cols;  //Number of columns in the matrix
	
};

/**
Sparse matrix class very similar to SpMatr_hash but built using dense hash from Sparsehash library
*/
class SpMatr_hash_fast
{
public:	
	SpMatr_hash_fast();  //Default constructor, generates empty matrix with 0 rows, 0 columns, only providing a default constructor...
	SpMatr_hash_fast(int nb_rows, int nb_cols, int expected_nnz); //Constructor for an empty matrix, takes as arguments the matrix dimensions and the expected number of non-zeros 
															//(approximate or exact is fine, fixes hash table's nb. of buckets)
	void print_matrix(); //For testing purposes
	double& operator()(uint64_t i,uint64_t j); //Overloading of operator (int,int) to access matrix element
	int get_nnz();  //Get number of non-zeros in the matrix
	int get_nb_rows();  //Returns nb_rows
	int get_nb_cols();  //Returns nb_cols
	
	typedef dense_hash_map<uint64_t, double, hash_int_pair>::iterator iterator; //Define iterator type for the unordered_map val
	iterator begin();  //Return iterator corresponding to .begin() of unordered_map val
	iterator end();  //Same thing for .end()
	unsigned int first_int(iterator& iter);  //Returns first int corresponding to the key of the mapping in ordered_map pointed to by iterator argument
	unsigned int second_int(iterator& iter); //Same as above for second int
	//THIS IS NOT NEEDED SpMatr_CSR* convert_to_csr();  //Converts matrix to a SpMatr_CSR object and returns a pointer to it	
	dense_hash_map<uint64_t, double, hash_int_pair> val;  //Hash table for storing the matrix entries
private:
	int nb_rows;  //Number of rows in the matrix
	int nb_cols;  //Number of columns in the matrix
	
};

/**
CSR sparse matrix class, one-based indexing
*/
class SpMatr_CSR
{
public:
	SpMatr_CSR(SpMatr_hash& init_matrix, bool ordered); //Constructor that takes a one-indexed SpMatr_hash sparse matrix as input, ordered should be true (ordered) or false (not ordered)
	SpMatr_CSR(SpMatr_hash_fast& init_matrix, bool ordered); //Constructor that takes a one-indexed SpMatr_hash sparse matrix as input, ordered should be true (ordered) or false (not ordered)
	SpMatr_CSR();  //Default constructor (initializes all variables to 0 or NULL pointers)
	SpMatr_CSR(const SpMatr_CSR& original);  //Copy constructor
	~SpMatr_CSR();  //Destructor
	SpMatr_CSR& operator=(SpMatr_CSR assigned_matrix);  //Assignment operator (copy-swap style)

	void multiply_row(int row_index, double factor);  //Multiplies row determined by row_index (one-based indexing) by factor

	double frob_norm();  //Computes and returns Frobenius norm (sqrt of sum of all squared elements)

	//Some accessor methods
	int get_nnz();  //Get number of non-zeros in the matrix
	int get_nb_rows();  //Returns nb_rows
	int get_nb_cols();  //Returns nb_cols
	double* get_values_ptr();
	int* get_columns_ptr();
	int* get_rowst_index_ptr();

	void print_arrays(); //For testing purposes

	double* values;  //Array with the non-zero element values
	
private:
	int nb_rows;
	int nb_cols;
	int nnz; //Number of non-zeros in the matrix
	bool ordered;  //true if each row is ordered by column indices, false otherwise
	
	int* columns;  //Column index of each element
	int* rowst_index; //Array of size nb_rows containing the index in values of the first element (start) of each row	
};



void test_spmatr(); //For testing only, temporary

void test_spmatr_mult();  //For testing
void test_spmatr_spd_CG();  //Test MKL CG for symmetric positive-definite matrix
void test_spmatr_GMRES();  //Test MKL FGMRES - for preconditioner debugging


#endif