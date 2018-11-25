/*
* MATLAB MEX function for computing squared quadratic errors between two matrices.
* -Can take a mask as input to ignore some elements.
* -Can also take column indices to ignore the other columns. (both the sq errors on the full matrix and column subset are computed in one step and returned.)
*/

#include "mex.h"
#include "matrix.h"
#include <vector>

/*
* Input:
- 0. First input Matrix.
- 1. Second input Matrix (both matrices must have the same dimensions).
- 2. (Optional) Mask : logical Matrix (false=ignored element).
- 3. (Optional) Column subset : vector of indices of the column subset for which a quadratic error must be measured. (the Mask specified in argument#2 also applies to the column subset).
* Output:
- quadratic error (ignoring masked elements) on all the columns.
- quadratic error (ignoring masked elements) on the column subset (same as previous output if no column subset is given).
*/
void checkInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nlhs > 2 )
		mexErrMsgTxt("Too many output arguments (only 2 required).");

	if (nlhs < 1)
		mexErrMsgTxt("No output argument.");

	if (nrhs > 4)
		mexErrMsgTxt("Incorrect number of input arguments (can take at most 4).");

	//Check data matrices
		//Dimensions
	mwSize n0 = mxGetNumberOfDimensions(prhs[0]);
	const mwSize* dims0 = mxGetDimensions(prhs[0]);
	mwSize n1 = mxGetNumberOfDimensions(prhs[1]);
	const mwSize* dims1 = mxGetDimensions(prhs[1]);
	if(n0!=n1)
		mexErrMsgTxt("Input Matrices must have the same number of dimensions.");
	for(int i=0; i<n0; ++i){
		if (dims0[i] != dims1[i]) mexErrMsgTxt("Input Matrices must have the same size.");
	}
		// type
	if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
		mexErrMsgTxt("Matrix in first argument must be of type single.");
	if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
		mexErrMsgTxt("Matrix in second argument must be of type single.");

	//Check Mask Matrix
	if (nrhs > 2) {
		mwSize n2 = mxGetNumberOfDimensions(prhs[2]);
		const mwSize* dims2 = mxGetDimensions(prhs[2]);
		if (n0 != n2)
			mexErrMsgTxt("Mask Matrix must have the same number of dimensions as data matrices.");
		for (int i = 0; i<n0; ++i) {
			if (dims0[i] != dims2[i]) mexErrMsgTxt("Mask Matrix must have the same size as data matrices.");
		}
		//type : logical
		if (mxGetClassID(prhs[2]) != mxLOGICAL_CLASS)
			mexErrMsgTxt("Mask Matrix must be of class logical.");
	}
	
	// Check the vector of column indices
		// number of dimensions
	if (nrhs > 3) {
		if (mxGetM(prhs[3]) > 1 && mxGetN(prhs[3]) > 1)
			mexErrMsgTxt("Fourth input argument must be a vector.");
		// type : double (only accept doubles for simplicity)
		if (mxGetClassID(prhs[3]) != mxUINT32_CLASS)
			mexErrMsgTxt("Vector in fourth argument must be of class uint32.");
	}

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
	float *m1 = (float*)mxGetData(prhs[0]);
	float *m2 = (float*)mxGetData(prhs[1]);

	bool *mask = NULL;
	if (nrhs > 2) {
		mask = (bool*) mxGetData(prhs[2]);
	}

	size_t nb_CIdx = 0;
	unsigned int *cIndices = NULL;
	if (nrhs > 3) {
		nb_CIdx = mxGetNumberOfElements(prhs[3]);
		cIndices = (unsigned int*)mxGetData(prhs[3]);
	}


//Initialize variables:
	float sumSqFull = 0;
	float sumSqSub = 0;
	int i, j, k, ijk = 0;
	float sqErrElt;
	char isInSubset;

	mwSize num_i=dims[0], num_j=dims[1], num_k=0;
	
	if (ndims > 2) {
		for (int dim = 2; dim < ndims; ++dim)
			num_k += dims[dim];
	}
	else num_k = 1;

	//Initialize table of columns in the subset
	std::vector<bool> inSubsetCol(num_j, false);
	if (nrhs > 3) {
		//bool* inSubsetCol = new bool[num_j]();
		int icIdx;
		for (j = 0; j < nb_CIdx; ++j) {
			icIdx =  cIndices[j];
			if( icIdx <= 0 || icIdx>num_j || cIndices[j]!=icIdx ) mexErrMsgTxt("Column indices must have integer values in a range corresponding to 2nd dimension of the matrices.");
			inSubsetCol[icIdx - 1] = true;
		}
	}
	else
	{
		std::fill(inSubsetCol.begin(), inSubsetCol.end(), true);
	}
	

//Main Loop
	if (nrhs > 2)
	{
		for (k = 0; k < num_k; ++k)
		for (j = 0; j < num_j; ++j)
		for (i = 0; i < num_i; ++i)
		{
			if (mask[ijk]) {
				sqErrElt = (m1[ijk] - m2[ijk]) * (m1[ijk] - m2[ijk]);
				sumSqFull += sqErrElt;
				if (inSubsetCol[j]) sumSqSub += sqErrElt;
			}
			++ijk;
		}
	}
	else // simple version : no mask / no subset of columns defined.
	{
		for (k = 0; k < num_k; ++k)
		for (j = 0; j < num_j; ++j)
		for (i = 0; i < num_i; ++i)
		{
			sqErrElt = (m1[ijk] - m2[ijk]) * (m1[ijk] - m2[ijk]);
			sumSqFull += sqErrElt;
			++ijk;
		}
		sumSqSub = sumSqFull;
	}
		
	//if (ndims > 3) delete[] inSubsetCol;

//Output
	//out = mxCreateDoubleMatrix(1, 1, mxREAL);
	//*mxGetPr(out) = sumSqFull;
	plhs[0] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	float * data = (float *)mxGetData(plhs[0]);
	data[0] = sumSqFull;
	plhs[1] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
	data = (float *)mxGetData(plhs[1]);
	data[0] = sumSqFull;
//	plhs[0] = mxCreateDoubleScalar(sumSqFull);
//	plhs[1] = mxCreateDoubleScalar(sumSqSub);

	return;

}


