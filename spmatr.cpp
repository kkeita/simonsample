#include "spmatr.h"

using namespace std;

//SpMatr_hash definitions
size_t hash_int_pair::custom_hash_fn(const unsigned int input) const
{
	size_t hash_result;
	hash_result = ((input >> 16) ^ input) * 0x45d9f3b;
    hash_result = ((hash_result >> 16) ^ hash_result) * 0x45d9f3b;
    hash_result = (hash_result >> 16) ^ hash_result;
    return hash_result;
}

size_t hash_int_pair::operator()(const uint64_t index_pair) const
{
	//return hash_fn(index_pair.first) ^ hash_fn(index_pair.second);  //This has the undesired property that hash(a,a) = hash(b,b), for any a,b
	//return hash_fn(index_pair.first * 3) ^ hash_fn(index_pair.second);
	//return hash_fn(((int) (index_pair >> 32)) * 3) ^ hash_fn((int) (index_pair & 0x00000000FFFFFFFF));
	return custom_hash_fn(((unsigned int) (index_pair >> 32)) * 3) ^ custom_hash_fn((unsigned int) (index_pair & 0x00000000FFFFFFFF));	
}

//******************************************************
SpMatr_hash::SpMatr_hash() : nb_rows(0), nb_cols(0)
{
}

SpMatr_hash::SpMatr_hash(int nb_rows, int nb_cols, int expected_nnz) : nb_rows(nb_rows), nb_cols(nb_cols)
{	
	val.max_load_factor(1.0);
	//val.reserve(expected_nnz);
	val.rehash(expected_nnz);
}

void SpMatr_hash::print_matrix()
{
	for(SpMatr_hash::iterator iter = val.begin(); iter != val.end(); iter++)
	{
		//cout << "(" << (iter->first).first << ", " << (iter->first).second << ") = " << iter->second << "\n";
		cout << "(" << first_int(iter) << ", " << second_int(iter) << ") = " << iter->second << "\n";
	}
}

double& SpMatr_hash::operator()(uint64_t i,uint64_t j)
{
	return val[(i << 32) | j];
	//return val[pair<int,int>(i,j)];
}

int SpMatr_hash::get_nnz()
{
	return val.size();
}

int SpMatr_hash::get_nb_rows()
{
	return nb_rows;
}

int SpMatr_hash::get_nb_cols()
{
	return nb_cols;
}

SpMatr_hash::iterator SpMatr_hash::begin()
{
	return val.begin();
}

SpMatr_hash::iterator SpMatr_hash::end()
{
	return val.end();
}

unsigned int SpMatr_hash::first_int(iterator& iter)
{
	return (unsigned int) (iter->first >> 32);
}

unsigned int SpMatr_hash::second_int(iterator& iter)
{
	return (unsigned int) (iter->first & 0x00000000FFFFFFFF);
}
//***************************************************
SpMatr_hash_fast::SpMatr_hash_fast() : nb_rows(0), nb_cols(0)
{
}

SpMatr_hash_fast::SpMatr_hash_fast(int nb_rows, int nb_cols, int expected_nnz) : nb_rows(nb_rows), nb_cols(nb_cols)
{	
	val.set_empty_key(0xFFFFFFFFFFFFFFFF);
	val.set_deleted_key(0xFFFFFFFFFFFFFFFE);
	val.min_load_factor(0);
	val.max_load_factor(0.8);
	val.resize(expected_nnz*1.3);	
}

void SpMatr_hash_fast::print_matrix()
{
	for(SpMatr_hash_fast::iterator iter = val.begin(); iter != val.end(); iter++)
	{
		//cout << "(" << (iter->first).first << ", " << (iter->first).second << ") = " << iter->second << "\n";
		cout << "(" << first_int(iter) << ", " << second_int(iter) << ") = " << iter->second << "\n";
	}
}

double& SpMatr_hash_fast::operator()(uint64_t i,uint64_t j)
{	
	return val[(i << 32) | j];
	//return val[pair<int,int>(i,j)];
}

int SpMatr_hash_fast::get_nnz()
{
	return val.size();
}

int SpMatr_hash_fast::get_nb_rows()
{
	return nb_rows;
}

int SpMatr_hash_fast::get_nb_cols()
{
	return nb_cols;
}

SpMatr_hash_fast::iterator SpMatr_hash_fast::begin()
{
	return val.begin();
}

SpMatr_hash_fast::iterator SpMatr_hash_fast::end()
{
	return val.end();
}

unsigned int SpMatr_hash_fast::first_int(iterator& iter)
{
	return (unsigned int) (iter->first >> 32);
}

unsigned int SpMatr_hash_fast::second_int(iterator& iter)
{
	return (unsigned int) (iter->first & 0x00000000FFFFFFFF);
}

//***************************************************
//SpMatr_CSR definitions

SpMatr_CSR::SpMatr_CSR(SpMatr_hash& init_matrix, bool ordered) : ordered(ordered)
{	
	nb_rows = init_matrix.get_nb_rows();
	nb_cols = init_matrix.get_nb_cols();
	nnz = init_matrix.get_nnz();
	
	if(ordered == false) //Elements in each row will not necessarily be sorted by column index
	{
		//Allocate rowst_index array and set it to zero
		rowst_index = (int*) malloc(sizeof(int) * (nb_rows + 1)); //Size is nb_rows plus one dummy index at the end
		memset(rowst_index, 0, sizeof(int) * (nb_rows + 1));

		//Count number of non-zeros in each row and store it in rowst_index (temporarily)
		for(SpMatr_hash::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//rowst_index[(iter->first).first - 1] += 1;  //The -1 is because of one-based indexing
			rowst_index[init_matrix.first_int(iter) - 1] += 1;  //The -1 is because of one-based indexing
		}

		//Set rowst_index values to indices of values corresponding to the element after the last element in each row
		rowst_index[0] += 1;
		for(int i = 1; i < nb_rows; i++)
		{
			rowst_index[i] += rowst_index[i-1];
		}

		//Set last element of rowst_index (dummy index) to nnz + 1
		rowst_index[nb_rows] = nnz + 1;

		//Allocate memory for arrays values and columns
		values = (double*) malloc(sizeof(double) * nnz);
		columns = (int*) malloc(sizeof(int) * nnz);
		
		//Fill values and columns as we iterate through elements
		int element_index;
		int row_index;
		for(SpMatr_hash::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//row_index = (iter->first).first;
			row_index = init_matrix.first_int(iter);
			element_index = --(rowst_index[row_index - 1]);  // One-based index in values and columns arrays where element will be placed, corresponding rowst_index is decremented beforehand
			values[element_index-1] = iter->second;  // -1 because of one-based indexing
			//columns[element_index-1] = (iter->first).second;	
			columns[element_index-1] = init_matrix.second_int(iter);
		}		
	}
	else if(ordered == true) //Each row will be sorted by column indices
	{
		//First generate a column-linked sparse matrix representation (note that one-based indexing is used)

		//Allocate memory for matrix arrays 
		double* values_temp = (double*) malloc(sizeof(double) * nnz);  //Contains values of elements for the column-linked representation
		int* rows = (int*) malloc(sizeof(int) * nnz);  //rows contains the row index associated to each value element
		int* link = (int*) malloc(sizeof(int) * nnz);  //link arrays for column-based linked list (contains index number of next element in the list, or zero if its the last element in the list)	
		int* colst_index = (int*) malloc(sizeof(int) * nb_cols);  //colst_index array contains indices of the first element for each column in
																//arrays values_temp, row and link, or 0 if the column is empty.
		memset(values_temp, 0, sizeof(double) * nnz);  //zero-initialize
		memset(rows, 0, sizeof(int) * nnz);
		memset(link, 0, sizeof(int) * nnz);
		memset(colst_index, 0, sizeof(int) * nb_cols);		

		//Iterate through each matrix element and add it to column-linked matrix
		int col_index;
		int counter = 1;  //Will hold index in values_temp, link and rows where new elements have to be added
		for(SpMatr_hash::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//col_index = (iter->first).second;  //Column index of element to insert
			col_index = init_matrix.second_int(iter);  //Column index of element to insert
			values_temp[counter-1] = iter->second;  //-1 because of one-based indexing
			//rows[counter-1] = (iter->first).first;
			rows[counter-1] = init_matrix.first_int(iter);
			link[counter-1] = colst_index[col_index-1];  //Set link to previous first element in column list
			colst_index[col_index-1] = counter;  //Set first element in column list to the element just added
			counter++;
		}
		//End of column-linked sparse matrix generation

		//From the column-linked sparse matrix, generate the ordered CSR matrix		
		
		//Allocate rowst_index array and set it to zero
		rowst_index = (int*) malloc(sizeof(int) * (nb_rows + 1)); //Size is nb_rows plus one dummy index at the end
		memset(rowst_index, 0, sizeof(int) * (nb_rows + 1));
		
		//Count number of non-zeros in each row and store it in rowst_index (temporarily)
		for(SpMatr_hash::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//rowst_index[(iter->first).first - 1] += 1;  //The -1 is because of one-based indexing			
			rowst_index[init_matrix.first_int(iter) - 1] += 1;  //The -1 is because of one-based indexing
		}

		//Set rowst_index values to indices of values corresponding to the element after the last element in each row
		rowst_index[0] += 1;
		for(int i = 1; i < nb_rows; i++)
		{
			rowst_index[i] += rowst_index[i-1];
		}

		//Set last element of rowst_index (dummy index) to nnz + 1
		rowst_index[nb_rows] = nnz + 1;

		//Allocate memory for arrays values and columns
		values = (double*) malloc(sizeof(double) * nnz);
		columns = (int*) malloc(sizeof(int) * nnz);

		//Fill values and columns as we iterate through elements by columns, starting with the last column
		int element_index;
		int row_index;
		int next_link;
		for(int col_index = nb_cols; col_index > 0; col_index--)
		{
			next_link = colst_index[col_index-1];
			while(next_link != 0)
			{
				row_index = rows[next_link-1];
				element_index = --(rowst_index[row_index-1]);  // One-based index in values and columns arrays where element will be placed, corresponding rowst_index is decremented beforehand
				values[element_index-1] = values_temp[next_link-1];
				columns[element_index-1] = col_index;
				next_link = link[next_link-1];
			}
		}	

		//Free memory used to store column-linked sparse matrix
		free(values_temp);
		free(rows);
		free(link);
		free(colst_index);		
	}

}

SpMatr_CSR::SpMatr_CSR(SpMatr_hash_fast& init_matrix, bool ordered) : ordered(ordered)
{	
	nb_rows = init_matrix.get_nb_rows();
	nb_cols = init_matrix.get_nb_cols();
	nnz = init_matrix.get_nnz();
	
	if(ordered == false) //Elements in each row will not necessarily be sorted by column index
	{
		//Allocate rowst_index array and set it to zero
		rowst_index = (int*) malloc(sizeof(int) * (nb_rows + 1)); //Size is nb_rows plus one dummy index at the end
		memset(rowst_index, 0, sizeof(int) * (nb_rows + 1));

		//Count number of non-zeros in each row and store it in rowst_index (temporarily)
		for(SpMatr_hash_fast::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//rowst_index[(iter->first).first - 1] += 1;  //The -1 is because of one-based indexing
			rowst_index[init_matrix.first_int(iter) - 1] += 1;  //The -1 is because of one-based indexing
		}

		//Set rowst_index values to indices of values corresponding to the element after the last element in each row
		rowst_index[0] += 1;
		for(int i = 1; i < nb_rows; i++)
		{
			rowst_index[i] += rowst_index[i-1];
		}

		//Set last element of rowst_index (dummy index) to nnz + 1
		rowst_index[nb_rows] = nnz + 1;

		//Allocate memory for arrays values and columns
		values = (double*) malloc(sizeof(double) * nnz);
		columns = (int*) malloc(sizeof(int) * nnz);
		
		//Fill values and columns as we iterate through elements
		int element_index;
		int row_index;
		for(SpMatr_hash_fast::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//row_index = (iter->first).first;
			row_index = init_matrix.first_int(iter);
			element_index = --(rowst_index[row_index - 1]);  // One-based index in values and columns arrays where element will be placed, corresponding rowst_index is decremented beforehand
			values[element_index-1] = iter->second;  // -1 because of one-based indexing
			//columns[element_index-1] = (iter->first).second;	
			columns[element_index-1] = init_matrix.second_int(iter);
		}		
	}
	else if(ordered == true) //Each row will be sorted by column indices
	{
		//First generate a column-linked sparse matrix representation (note that one-based indexing is used)

		//Allocate memory for matrix arrays 
		double* values_temp = (double*) malloc(sizeof(double) * nnz);  //Contains values of elements for the column-linked representation
		int* rows = (int*) malloc(sizeof(int) * nnz);  //rows contains the row index associated to each value element
		int* link = (int*) malloc(sizeof(int) * nnz);  //link arrays for column-based linked list (contains index number of next element in the list, or zero if its the last element in the list)	
		int* colst_index = (int*) malloc(sizeof(int) * nb_cols);  //colst_index array contains indices of the first element for each column in
																//arrays values_temp, row and link, or 0 if the column is empty.
		memset(values_temp, 0, sizeof(double) * nnz);  //zero-initialize
		memset(rows, 0, sizeof(int) * nnz);
		memset(link, 0, sizeof(int) * nnz);
		memset(colst_index, 0, sizeof(int) * nb_cols);		

		//Iterate through each matrix element and add it to column-linked matrix
		int col_index;
		int counter = 1;  //Will hold index in values_temp, link and rows where new elements have to be added
		for(SpMatr_hash_fast::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//col_index = (iter->first).second;  //Column index of element to insert
			col_index = init_matrix.second_int(iter);  //Column index of element to insert
			values_temp[counter-1] = iter->second;  //-1 because of one-based indexing
			//rows[counter-1] = (iter->first).first;
			rows[counter-1] = init_matrix.first_int(iter);
			link[counter-1] = colst_index[col_index-1];  //Set link to previous first element in column list
			colst_index[col_index-1] = counter;  //Set first element in column list to the element just added
			counter++;
		}
		//End of column-linked sparse matrix generation

		//From the column-linked sparse matrix, generate the ordered CSR matrix		
		
		//Allocate rowst_index array and set it to zero
		rowst_index = (int*) malloc(sizeof(int) * (nb_rows + 1)); //Size is nb_rows plus one dummy index at the end
		memset(rowst_index, 0, sizeof(int) * (nb_rows + 1));
		
		//Count number of non-zeros in each row and store it in rowst_index (temporarily)
		for(SpMatr_hash_fast::iterator iter = init_matrix.begin(); iter != init_matrix.end(); iter++)
		{
			//rowst_index[(iter->first).first - 1] += 1;  //The -1 is because of one-based indexing			
			rowst_index[init_matrix.first_int(iter) - 1] += 1;  //The -1 is because of one-based indexing
		}

		//Set rowst_index values to indices of values corresponding to the element after the last element in each row
		rowst_index[0] += 1;
		for(int i = 1; i < nb_rows; i++)
		{
			rowst_index[i] += rowst_index[i-1];
		}

		//Set last element of rowst_index (dummy index) to nnz + 1
		rowst_index[nb_rows] = nnz + 1;

		//Allocate memory for arrays values and columns
		values = (double*) malloc(sizeof(double) * nnz);
		columns = (int*) malloc(sizeof(int) * nnz);

		//Fill values and columns as we iterate through elements by columns, starting with the last column
		int element_index;
		int row_index;
		int next_link;
		for(int col_index = nb_cols; col_index > 0; col_index--)
		{
			next_link = colst_index[col_index-1];
			while(next_link != 0)
			{
				row_index = rows[next_link-1];
				element_index = --(rowst_index[row_index-1]);  // One-based index in values and columns arrays where element will be placed, corresponding rowst_index is decremented beforehand
				values[element_index-1] = values_temp[next_link-1];
				columns[element_index-1] = col_index;
				next_link = link[next_link-1];
			}
		}	

		//Free memory used to store column-linked sparse matrix
		free(values_temp);
		free(rows);
		free(link);
		free(colst_index);		
	}

}

SpMatr_CSR::SpMatr_CSR() : nb_rows(0), nb_cols(0), nnz(0), ordered(true), values(NULL), columns(NULL), rowst_index(NULL)
{
}

SpMatr_CSR::SpMatr_CSR(const SpMatr_CSR& original) : nb_rows(original.nb_rows), nb_cols(original.nb_cols), nnz(original.nnz), ordered(original.ordered)
{
	//Allocate memory for arrays
	values = (double*) malloc(sizeof(double) * nnz);
	columns = (int*) malloc(sizeof(int) * nnz);
	rowst_index = (int*) malloc(sizeof(int) * (nb_rows + 1));

	//Copy array values
	memcpy(values, original.values, sizeof(double) * nnz);
	memcpy(columns, original.columns, sizeof(int) * nnz);
	memcpy(rowst_index, original.rowst_index, sizeof(int) * (nb_rows + 1));
}

SpMatr_CSR::~SpMatr_CSR()
{	
	free(values);
	free(columns);
	free(rowst_index);	
}

SpMatr_CSR& SpMatr_CSR::operator=(SpMatr_CSR assigned_matrix)  //assigned_matrix is passed by value (hence is copied using the copy constructor)
{
	//Now we just need to swap *this with assigned_matrix, the destructor will be called on assigned_matrix, therefore deleting the original *this
	swap(nb_rows, assigned_matrix.nb_rows);	
	swap(nb_cols, assigned_matrix.nb_cols);	
	swap(nnz, assigned_matrix.nnz);	
	swap(ordered, assigned_matrix.ordered);
	
	swap(values, assigned_matrix.values);  //Note that only pointers are swapped, not the whole array contents	
	swap(columns, assigned_matrix.columns);
	swap(rowst_index, assigned_matrix.rowst_index);
	return *this;
}

void SpMatr_CSR::multiply_row(int row_index, double factor)
{
	int start_index = rowst_index[row_index - 1] - 1;  //Index of row's first element in values array (-1's because of one-based indexing)
	int start_next_row_index = rowst_index[row_index] - 1;  //Index next row's first element in values array

	for(int i = start_index; i < start_next_row_index; i++)
	{
		values[i] *= factor;
	}
}

double SpMatr_CSR::frob_norm()
{
	double frob_norm = 0.0;

	for(int i = 0; i < nnz; i++)
	{
		frob_norm += pow(values[i],2);
	}

	frob_norm = sqrt(frob_norm);

	return frob_norm;
}

int SpMatr_CSR::get_nnz()
{
	return nnz;
}

int SpMatr_CSR::get_nb_rows()
{
	return nb_rows;
}

int SpMatr_CSR::get_nb_cols()
{
	return nb_cols;
}

double* SpMatr_CSR::get_values_ptr()
{
	return values;
}

int* SpMatr_CSR::get_columns_ptr()
{
	return columns;
}

int* SpMatr_CSR::get_rowst_index_ptr()
{
	return rowst_index;
}

void SpMatr_CSR::print_arrays()
{
	cout << setw(15) << "values: ";
	for(int i = 0; i < nnz; i++)
	{
		cout << values[i] << ", ";
	}
	cout << "\n";

	cout << setw(15) << "columns: ";
	for(int i = 0; i < nnz; i++)
	{
		cout << columns[i] << ", ";
	}
	cout << "\n";

	cout << setw(15) << "rowst_index: ";
	for(int i = 0; i < (nb_rows + 1); i++)
	{
		cout << rowst_index[i] << ", ";
	}
	cout << "\n";
}

//****************************************************
void test_spmatr()
{
	/*
	SpMatr_hash my_matrix(4,4,5);
	
	my_matrix(1,1) = 2.0;
	my_matrix(2,1) = 14.5;
	//my_matrix(1,1) = my_matrix(1,1) + 1.5;
	my_matrix(1,1) += 1.5;
	//my_matrix.print_matrix();
	*/

	//SpMatr_hash my_matrix2(5, 7, 9);
	SpMatr_hash_fast my_matrix2(5, 7, 9);

	my_matrix2(1,1) = 1;
	my_matrix2(1,2) += -1;
	my_matrix2(1,4) = -3;
	my_matrix2(2,2) = 5;
	my_matrix2(3,3) = 4;
	my_matrix2(3,4) = 6;
	my_matrix2(3,5) = 4;
	my_matrix2(4,4) = 7;
	my_matrix2(5,7) = -5;
	my_matrix2(5,7) += -3;
	

	//SpMatr_CSR my_CSR_matrix(my_matrix2, false);
	SpMatr_CSR my_CSR_matrix(my_matrix2, true);

	my_CSR_matrix.print_arrays();

	/*
	//Timing test to compare my_matrix(1,1) = my_matrix(1,1) + x with my_matrix(1,1 += x
	double cpu_time1, cpu_time2;
	my_matrix(1,1) = 0;

	cpu_time1 = dsecnd();
	for(int k = 0; k < 100000; k++)
	{
		//my_matrix(1,1) = my_matrix(1,1) + 1.0; //this takes 2.60 ms
		my_matrix(1,1) += 1.0;					//this takes 1.30 ms
	}
	cpu_time2 = dsecnd();
	cout << "1000 assembly operations in " << cpu_time2 - cpu_time1 << " s\n";
	// /Timing test
	*/
	
}

void test_spmatr_mult()
{
	SpMatr_hash_fast my_matrix2(5, 5, 13);

	my_matrix2(1,1) = 1;
	my_matrix2(1,3) += -1;
	my_matrix2(2,1) = -1;
	my_matrix2(2,2) = 1;
	my_matrix2(2,4) = -1;	
	my_matrix2(3,2) = 1;
	my_matrix2(3,3) = -2;
	my_matrix2(3,5) = 1;
	my_matrix2(4,3) = -1;
	my_matrix2(4,4) = 2;
	my_matrix2(4,5) = -1;
	my_matrix2(5,4) = -1;
	my_matrix2(5,5) += -3;

	SpMatr_CSR my_CSR_matrix(my_matrix2, true);
	my_CSR_matrix.print_arrays();

	double x[5] = {-1, 1, 0, 1, -1};

	printf("Array x:\n");
	for(int i = 0; i<5;i++)
	{
		printf("%f\n", x[i]);
	}

	double y[5];
	memset(y, 0, sizeof(double)*5);

	char transaction = 'n';
	MKL_INT m = 5;
	mkl_dcsrgemv(&transaction, &m, my_CSR_matrix.get_values_ptr(), my_CSR_matrix.get_rowst_index_ptr(), my_CSR_matrix.get_columns_ptr(), x, y);

	printf("Array y:\n");
	for(int i = 0; i<5;i++)
	{
		printf("%f\n", y[i]);
	}

}

void test_spmatr_spd_CG()
{
	SpMatr_hash_fast my_matrix(5, 5, 10);

	my_matrix(1,1) = 1;
	my_matrix(1,3) = 1;
	my_matrix(1,4) += 2;
	my_matrix(2,2) += 4;
	my_matrix(2,4) = 1;
	my_matrix(3,3) = 2;
	my_matrix(3,4) = 3;
	my_matrix(3,5) = 1;
	my_matrix(4,4) = 8;
	my_matrix(5,5) = 5;

	SpMatr_CSR my_CSR_matrix(my_matrix, true);

	//my_CSR_matrix.print_arrays();

	//Initialize CG solver
	double *CG_tmp;
	MKL_INT CG_rci_request;
	MKL_INT CG_itercount;
	MKL_INT CG_ipar[128];
	double CG_dpar[128];
	char uplo[1];
	double solution[5];  //Solution array
	double rhs[5] = {3, 5, 5, 13, -3};  //RHS

	int nb_dof = 5;
	CG_tmp = (double*) malloc(sizeof(double) * nb_dof * 4);	

	int CG_max_nb_iterations = 150;


	memset(solution, 0, sizeof(double)*5);  //Zero initial guess for solution
	//Use as initial guess 0
	uplo[0] = 'U';  //For matrix vector computation with mkl_dcsrsymv
	dcg_init(&nb_dof, solution, rhs, &CG_rci_request, CG_ipar, CG_dpar, CG_tmp);

	CG_ipar[4] = CG_max_nb_iterations;  //Max nb. of iterations
	CG_ipar[7] = 1;  //0: Do not check for max nb. of iterations reached, 1 (default): Stop when max nb. of iterations (ipar[4]) is reached
	CG_ipar[8] = 0;  //0 (default): Do not perform residual stopping test, 1: Perform residual stopping test
	//Residual stopping test: square norm of residual  <= relative_tolerance * square norm of initial residual + absolute tolerance
	CG_ipar[9] = 0;  //0: Do not request user-defined stopping test, 1 (default): Request user-defined stopping test
	CG_ipar[10] = 0;  //0 (default): Non-preconditioned version of CG, 1: preconditioned CG

	CG_dpar[0] = 1e-6;  //Relative tolerance (default: 1e-6)
	CG_dpar[1] = 0;  //Absolute tolerance (default: 0)

	//For debugging
	dcg_check(&nb_dof, solution, rhs, &CG_rci_request, CG_ipar, CG_dpar, CG_tmp);
	if(CG_rci_request != 0)
	{
		cerr << "Inconsistent data/parameters for MKL CG solver\n";
		exit(0);
	}
	
cg_rci:
	dcg(&nb_dof, solution, rhs, &CG_rci_request, CG_ipar, CG_dpar, CG_tmp);

	if(CG_rci_request == 0)  //Solution was found and stopping criterion satisfied (or max nb. of iterations reached)
		goto cg_getsln;

	if(CG_rci_request == -1)
	{
		//Max nb of iterations reached but convergence test not satisfied
		printf("Max nb. of iterations reached, failed to satisfy convergence criterion.\n");
		goto cg_getsln;
	}

	if(CG_rci_request == 1)  //Need to perform matrix-vector product
	{
		//Verbose for test:
		dcg_get(&nb_dof, solution, rhs, &CG_rci_request, CG_ipar, CG_dpar, CG_tmp, &CG_itercount);
		printf("Current iteration #: %d\n", CG_itercount);
		printf("Achieved norm of residual: %e\n", CG_dpar[4]);

		// /Verbose for test

		mkl_dcsrsymv(uplo, &nb_dof, my_CSR_matrix.get_values_ptr(), my_CSR_matrix.get_rowst_index_ptr(), my_CSR_matrix.get_columns_ptr(), CG_tmp, &CG_tmp[nb_dof]);
		goto cg_rci;
	}

cg_getsln:
	dcg_get(&nb_dof, solution, rhs, &CG_rci_request, CG_ipar, CG_dpar, CG_tmp, &CG_itercount);
	
	printf ("\nThe CG preconditioner system has been solved \n");
	printf ("Number of iterations: %d\n", CG_itercount);
	printf("Norm of Initial residual: %e\n", CG_dpar[2]);
	printf("Target norm of residual: %e\n", CG_dpar[3]);
	printf("Achieved norm of residual: %e\n", CG_dpar[4]);

	printf("Solution vector = [%f %f %f %f %f]", solution[0], solution[1], solution[2], solution[3], solution[4]);


	free(CG_tmp);
	exit(0);
}

void test_spmatr_GMRES()
{
	SpMatr_hash_fast my_matrix(5, 5, 13);

	my_matrix(1,1) = 1;
	my_matrix(1,3) += -1;
	my_matrix(2,1) = -1;
	my_matrix(2,2) = 1;
	my_matrix(2,4) = -1;	
	my_matrix(3,2) = 1;
	my_matrix(3,3) = -2;
	my_matrix(3,5) = 1;
	my_matrix(4,3) = -1;
	my_matrix(4,4) = 2;
	my_matrix(4,5) = -1;
	my_matrix(5,4) = -1;
	my_matrix(5,5) += -3;

	
	SpMatr_CSR my_CSR_matrix(my_matrix, true);

	//my_CSR_matrix.print_arrays();

	//Initialize CG solver
	double *FGMRES_tmp;
	MKL_INT FGMRES_rci_request;
	MKL_INT FGMRES_itercount;
	MKL_INT FGMRES_ipar[128];
	double FGMRES_dpar[128];
	
	double solution[5];  //Solution array
	//double rhs[5] = {-1, 1, 0, 3, 2};  //RHS, Simon example
	double rhs[5] = {-1,1,0,3,2};
	double residual[5];
	
	int nb_dof = 5;
	MKL_INT n = nb_dof;
	int input_index, output_index;

	char transaction = 'n';  //Doing a normal matrix-vector multiplication
		
	int FGMRES_max_nb_iterations = 150;
	FGMRES_tmp = (double*) malloc(sizeof(double) * ((2 * FGMRES_max_nb_iterations + 1) * nb_dof + FGMRES_max_nb_iterations * (FGMRES_max_nb_iterations + 9) / 2 + 1) ); 
	
	memset(solution, 0, sizeof(double)*5);  //Zero initial guess for solution
	//Use as initial guess 0
	
	dfgmres_init(&n, solution, rhs, &FGMRES_rci_request, FGMRES_ipar, FGMRES_dpar, FGMRES_tmp);

	FGMRES_ipar[4] = FGMRES_max_nb_iterations;  //Max nb. of iterations
	FGMRES_ipar[7] = 1;  //0: Do not check for max nb. of iterations reached, 1 (default): Stop when max nb. of iterations (ipar[4]) is reached
	FGMRES_ipar[8] = 0;  //0 (default): Do not perform residual stopping test, 1: Perform residual stopping test
	//Residual stopping test: square norm of residual  <= relative_tolerance * square norm of initial residual + absolute tolerance
	FGMRES_ipar[9] = 0;  //0: Do not request user-defined stopping test, 1 (default): Request user-defined stopping test
	FGMRES_ipar[10] = 0;  //0 (default): Non-preconditioned version of FGMRES, 1: preconditioned FGMRES
	FGMRES_ipar[11] = 1;
	FGMRES_ipar[12] = 0;  //0: Writes computed solution to x, -1: solution array is not updated
	//FGMRES_ipar[14] = FGMRES_max_nb_iterations;  //Used to specify the number of non-restarted iterations, we don't use restarted iterations so this is max_nb_iterations
	
	FGMRES_dpar[0] = 1e-6;  //Relative tolerance (default: 1e-6)
	FGMRES_dpar[1] = 0;  //Absolute tolerance (default: 0)

	//For debugging
	dfgmres_check(&n, solution, rhs, &FGMRES_rci_request, FGMRES_ipar, FGMRES_dpar, FGMRES_tmp);
	if(FGMRES_rci_request != 0)
	{
		cerr << "Inconsistent data/parameters for MKL FGMRES solver\n";
		exit(0);
	}

fgmres_rci:
	dfgmres(&n, solution, rhs, &FGMRES_rci_request, FGMRES_ipar, FGMRES_dpar, FGMRES_tmp);
	printf("Just called dfmgres, got rci_request = %d\n", FGMRES_rci_request);
	if(FGMRES_rci_request == 0)  //Solution was found and stopping criterion satisfied (or max nb. of iterations reached)
		goto fgmres_getsln;

	if(FGMRES_rci_request == -1)
	{
		//Max nb of iterations reached but convergence test not satisfied
		printf("Max nb. of iterations reached, failed to satisfy convergence criterion.\n");
		goto fgmres_getsln;
	}

	if(FGMRES_rci_request == 1)  //Need to perform matrix-vector product
	{
		//Verbose for test:
		//dfgmres_get(&n, solution, rhs, &FGMRES_rci_request, FGMRES_ipar, FGMRES_dpar, FGMRES_tmp, &FGMRES_itercount);
		printf("Current iteration #: %d\n", FGMRES_itercount);
		printf("Achieved norm of residual: %e\n", FGMRES_dpar[4]);

		// /Verbose for test
		input_index = FGMRES_ipar[21] - 1;
		output_index = FGMRES_ipar[22] - 1;
		mkl_dcsrgemv(&transaction, &n, my_CSR_matrix.get_values_ptr(), my_CSR_matrix.get_rowst_index_ptr(), my_CSR_matrix.get_columns_ptr(), &FGMRES_tmp[input_index], &FGMRES_tmp[output_index]);

		//printf("Matrix multiplication input: %f %f %f %f %f\n", FGMRES_tmp[input_index], FGMRES_tmp[input_index+1], FGMRES_tmp[input_index+2], FGMRES_tmp[input_index+3], FGMRES_tmp[input_index+4]);
		//printf("Matrix multiplication output: %f %f %f %f %f\n", FGMRES_tmp[output_index], FGMRES_tmp[output_index+1], FGMRES_tmp[output_index+2], FGMRES_tmp[output_index+3], FGMRES_tmp[output_index+4]);
		goto fgmres_rci;
	}

	
fgmres_getsln:
	dfgmres_get(&n, solution, rhs, &FGMRES_rci_request, FGMRES_ipar, FGMRES_dpar, FGMRES_tmp, &FGMRES_itercount);

	printf ("\nThe system has been solved \n");
	printf ("Number of iterations: %d\n", FGMRES_itercount);
	printf("Norm of Initial residual: %e\n", FGMRES_dpar[2]);
	printf("Target norm of residual: %e\n", FGMRES_dpar[3]);
	printf("Achieved norm of residual: %e\n", FGMRES_dpar[4]);

	printf("Solution vector = [%f %f %f %f %f]\n", solution[0], solution[1], solution[2], solution[3], solution[4]);

	mkl_dcsrgemv(&transaction, &n, my_CSR_matrix.get_values_ptr(), my_CSR_matrix.get_rowst_index_ptr(), my_CSR_matrix.get_columns_ptr(), solution, residual);
	for(int i = 0; i < n; i++)
	{
		residual[i] += -rhs[i];
	}
	printf("Norm of RHS = %e\n", cblas_dnrm2(n, rhs, 1));
	printf("Computed residual = %e\n", cblas_dnrm2(n, residual, 1));

	free(FGMRES_tmp);
	exit(0);
}