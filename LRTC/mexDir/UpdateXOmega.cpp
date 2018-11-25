/*
* MATLAB MEX function. Faster version of X update in Omega for noisy case...
* WARNING : the function modifies first input argument. The function should not be called if there is a 
* shared copy of the matrix given as first argument in the matlab workspace !!!!!!
*/


#include "mex.h"
#include "matrix.h"
#include <vector>

/*
* Input:
- 0. matrix X to update.
- 1. matrix M of known data (in Omega)
- 2. Matrix Omega : logical Matrix (false or 0 = unknown elt of M).
- 3. lambda
- 4. Column subset : vector of indices of the column subset for the update of X tolerates noise in Omega.
* Output: (no output, just update X)
*/
void checkInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nlhs > 0 )
		mexErrMsgTxt("No output argument required.");

	if (nrhs > 5)
		mexErrMsgTxt("Incorrect number of input arguments (must take 5).");

	//Check Matrices
		//Dimensions
	mwSize n0 = mxGetNumberOfDimensions(prhs[0]);
	const mwSize* dims0 = mxGetDimensions(prhs[0]);
	mwSize n1 = mxGetNumberOfDimensions(prhs[1]);
	const mwSize* dims1 = mxGetDimensions(prhs[1]);
	mwSize n2 = mxGetNumberOfDimensions(prhs[2]);
	const mwSize* dims2 = mxGetDimensions(prhs[2]);

	if(n0 != n1 || n0 != n2)
		mexErrMsgTxt("Input Matrices must have the same number of dimensions.");
	for(int i=0; i<n0; ++i){
		if (dims0[i] != dims1[i] || dims0[i] != dims2[i]) mexErrMsgTxt("Input Matrices must have the same size.");
	}
		// types
	if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
		mexErrMsgTxt("Matrix in first argument must be of type double.");
	if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
		mexErrMsgTxt("Matrix in second argument must be of type double.");
	if (mxGetClassID(prhs[2]) != mxLOGICAL_CLASS)
		mexErrMsgTxt("Mask Matrix in third argument must be of class logical.");
	
	// Check lambda
	if (mxGetClassID(prhs[3]) != mxDOUBLE_CLASS || mxGetNumberOfElements(prhs[3]) != 1)
		mexErrMsgTxt("lambda parameter should be a scalar double.");

	// Check the vector of column indices
		// number of dimensions
	if (mxGetM(prhs[4]) > 1 && mxGetN(prhs[4]) > 1)
		mexErrMsgTxt("Fourth input argument must be a vector.");
		// type : double (only accept doubles for simplicity)
	if (mxGetClassID(prhs[4]) != mxDOUBLE_CLASS)
		mexErrMsgTxt("Vector in fourth argument must be of class double.");

	return;
}



/*
* mex function
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//Inputs:
	checkInputs(nlhs, plhs, nrhs, prhs);

	mwSize ndims = mxGetNumberOfDimensions(prhs[0]);
	const mwSize* dims = mxGetDimensions(prhs[0]);
	double *X = (double*)mxGetData(prhs[0]);
	double *M = (double*)mxGetData(prhs[1]);
	bool *Omega = (bool*) mxGetData(prhs[2]);

	double lambda = mxGetScalar(prhs[3]);

	size_t nb_CIdx = 0;
	double *cIndices = NULL;
	nb_CIdx = mxGetNumberOfElements(prhs[4]);
	cIndices = (double*)mxGetData(prhs[4]);


//Initialize variables:
	double sumSqFull = 0;
	double sumSqSub = 0;
	int i, j, k, ijk = 0;
	double sqErrElt;
	char isInSubset;

	mwSize num_i=dims[0], num_j=dims[1], num_k=0;
	
	if (ndims > 2) {
		for (int dim = 2; dim < ndims; ++dim)
			num_k += dims[dim];
	}
	else num_k = 1;

	//Initialize table of columns in the subset
	std::vector<bool> inSubsetCol(num_j, false);
	int icIdx;
	for (j = 0; j < nb_CIdx; ++j) {
		icIdx = (int)cIndices[j];
		if (icIdx <= 0 || icIdx>num_j || cIndices[j] != icIdx) mexErrMsgTxt("Column indices must have integer values in a range corresponding to 2nd dimension of the matrices.");
		inSubsetCol[icIdx - 1] = true;
	}

//Main Loop
	for (k = 0; k < num_k; ++k)
	for (j = 0; j < num_j; ++j)
	for (i = 0; i < num_i; ++i)
	{
		if (Omega[ijk]) {
			if (inSubsetCol[j]){
				X[ijk] = (X[ijk] + lambda*M[ijk]) / (1 + lambda);
			}
			else {
				X[ijk] = M[ijk];
			}
		}
		++ijk;
	}
	
	//if (ndims > 3) delete[] inSubsetCol;

	return;

}


