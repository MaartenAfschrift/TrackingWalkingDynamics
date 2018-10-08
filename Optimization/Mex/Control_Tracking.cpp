//-Start header -----------------------------------------------------------
#include <OpenSim/OpenSim.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>  
#include <windows.h>
#include "matrix.h"
#include "mex.h"
using namespace OpenSim;
using namespace SimTK;
using namespace std;

// ........HEADER........................

// ........FUNCTIONS........................

//-End header -------------------------------------------------------------

//-Start MEX enterance ----------------------------------------------------
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// outIndex. Determines to which pointer in lhs to output which result. Starts at 1 since xDot is always an output (1st output)
	int outIndex = 1;
	// maximum expected values
	int max_states = 50;
	int max_pts = 900;
	int num_options = 4;
	int max_markers = 30;


	// Buffer sizes for shared memory (using maximum expected values)
	DWORD BUF_SIZE_STATES = sizeof(double)*max_states*max_pts;
	DWORD BUF_SIZE_TIME = sizeof(double)*max_pts;
	DWORD BUF_SIZE_BOOL = sizeof(bool);
	DWORD BUF_SIZE_INT = sizeof(int);
	DWORD BUF_SIZE_CHAR = sizeof(char)*max_markers;
	DWORD BUF_SIZE_BODIES = sizeof(double)*max_markers;
	DWORD BUF_SIZE_MarkerPos = sizeof(double)*max_markers*max_pts*3;

	// Names that provide key to unlocking shared memory
	TCHAR szNameq[] = TEXT("Global\\q");
	TCHAR szNameqdot[] = TEXT("Global\\qdot");
	TCHAR szNameqddot[] = TEXT("Global\\qddot");
	TCHAR szNameGRF_F_D[] = TEXT("Global\\GRF_F_D");
	TCHAR szNameGRF_P_D[] = TEXT("Global\\GRF_P_D");
	TCHAR szNameGRF_F_S[] = TEXT("Global\\GRF_F_S");
	TCHAR szNameGRF_P_S[] = TEXT("Global\\GRF_P_S");
	TCHAR szNamePT_coordB[] = TEXT("Global\\PT_coordB");
	TCHAR szNamePT_bodynames[] = TEXT("Global\\PT_bodynames");
	TCHAR szNameNumCoords[] = TEXT("Global\\NumCoords");
	TCHAR szNameNMarkers[] = TEXT("Global\\NMarkers");
	TCHAR szNameResults1[] = TEXT("Global\\Results1");
	TCHAR szNameResults2[] = TEXT("Global\\Results2");

	TCHAR szNameq_R[] = TEXT("Global\\q_R");
	TCHAR szNameqdot_R[] = TEXT("Global\\qdot_R");
	TCHAR szNameqddot_R[] = TEXT("Global\\qddot_R");
	TCHAR szNameGRF_F_D_R[] = TEXT("Global\\GRF_F_D_R");
	TCHAR szNameGRF_P_D_R[] = TEXT("Global\\GRF_P_D_R");
	TCHAR szNameGRF_F_S_R[] = TEXT("Global\\GRF_F_S_R");
	TCHAR szNameGRF_P_S_R[] = TEXT("Global\\GRF_P_S_R");
	TCHAR szNamePT_coordB_R[] = TEXT("Global\\PT_coordB_R");
	TCHAR szNamePT_bodynames_R[] = TEXT("Global\\PT_bodynames_R");
	TCHAR szNameNumCoords_R[] = TEXT("Global\\NumCoords_R");
	TCHAR szNameNMarkers_R[] = TEXT("Global\\NMarkers_R");
	TCHAR szNameResults3[] = TEXT("Global\\Results3");
	TCHAR szNameResults4[] = TEXT("Global\\Results4");

	TCHAR szNameKeepOpen[] = TEXT("Global\\KeepOpen");
	TCHAR szNameRunJob[] = TEXT("Global\\RunJob");
	TCHAR szNameNumPts[] = TEXT("Global\\NumPts");
	TCHAR szNameResReady[] = TEXT("Global\\ResReady");

	if (mxGetNumberOfElements(prhs[0]) == 1) // close program if only passing a single number
	{
		mexPrintf("TestMessage_close");
		bool* pBuf_KeepOpen;
		HANDLE hMapFile_KeepOpen; // keep program open
		hMapFile_KeepOpen = OpenFileMapping(FILE_MAP_ALL_ACCESS, FALSE, szNameKeepOpen);
		if (hMapFile_KeepOpen == NULL)
		{
			mexErrMsgTxt("Could not create file mapping object.\n");
		}
		pBuf_KeepOpen = (bool*)MapViewOfFile(hMapFile_KeepOpen, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
		if (pBuf_KeepOpen == NULL)
		{
			CloseHandle(hMapFile_KeepOpen);
			mexErrMsgTxt("Could not map view of file.\n");
		}
		// end other program
		*pBuf_KeepOpen = false;
		// close
		UnmapViewOfFile(pBuf_KeepOpen);
		CloseHandle(hMapFile_KeepOpen);
	}
	else
	{

		// mex inputs
		//mexPrintf("TestMessage");
		char *osimFileName_ptr = mxArrayToString(prhs[0]);
		//mexPrintf("TestMessage 2");

		//std::string osimFileName = std::string(&osimFileName_ptr[0]);
		//Model osimModel = Model(osimFileName);
		
		//mexErrMsgTxt("!!!Input Error!!!\n");

		// get the q - qdot and qddot
		double *q = mxGetPr(prhs[1]);
		const int *dimPosiz = mxGetDimensions(prhs[1]);
		int Ncoord = *(dimPosiz + 1);
		int numPts = *(dimPosiz + 0);
		double *qdot = mxGetPr(prhs[2]);
		double *qddot = mxGetPr(prhs[3]);

		// get the markenames and bodies
		double *PT_coordB = mxGetPr(prhs[4]);
		const int *dimPT_coordB = mxGetDimensions(prhs[4]);
		int NMarkers = *(dimPT_coordB + 1);
		double *PT_bodynames = mxGetPr(prhs[5]);

		// ge the ground reaction forces
		double *GRF_F_D = mxGetPr(prhs[6]);
		double *GRF_P_D = mxGetPr(prhs[7]);
		double *GRF_F_S = mxGetPr(prhs[8]);
		double *GRF_P_S = mxGetPr(prhs[9]);

		// get input information of the second model
		double *q_R = mxGetPr(prhs[10]);
		const int *dimPosiz_R = mxGetDimensions(prhs[10]);
		int Ncoord_R = *(dimPosiz_R + 1);
		double *qdot_R = mxGetPr(prhs[11]);
		double *qddot_R = mxGetPr(prhs[12]);		
		double *PT_coordB_R = mxGetPr(prhs[13]);
		const int *dimPT_coordB_R = mxGetDimensions(prhs[13]);
		int NMarkers_R = *(dimPT_coordB_R + 1);
		double *PT_bodynames_R = mxGetPr(prhs[14]);
		double *GRF_F_D_R = mxGetPr(prhs[15]);
		double *GRF_P_D_R = mxGetPr(prhs[16]);
		double *GRF_F_S_R= mxGetPr(prhs[17]);
		double *GRF_P_S_R = mxGetPr(prhs[18]);



		//mexPrintf("TestMessage 3");
		// mex outputs (only first output)
		plhs[0] = mxCreateDoubleMatrix(numPts, Ncoord, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(numPts * 3, NMarkers, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(numPts, Ncoord_R, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(numPts * 3, NMarkers_R, mxREAL);

		double *results1 = mxGetPr(plhs[0]);
		double *results2 = mxGetPr(plhs[1]);
		double *results3 = mxGetPr(plhs[2]);
		double *results4 = mxGetPr(plhs[3]);

		// mapfiles
		HANDLE hMapFile_q;					// data on which to operate
		HANDLE hMapFile_qdot;				// data on which to operate
		HANDLE hMapFile_qddot;				// data on which to operate
		HANDLE hMapFile_GRF_F_D;			// data on which to operate
		HANDLE hMapFile_GRF_P_D;			// data on which to operate
		HANDLE hMapFile_GRF_F_S;			// data on which to operate
		HANDLE hMapFile_GRF_P_S;			// data on which to operate
		HANDLE hMapFile_PT_coordB;			// data on which to operate
		HANDLE hMapFile_PT_bodynames;		// data on which to operate
		HANDLE hMapFile_NumCoords;			// number of coordinates in model
		HANDLE hMapFile_NMarkers;			// number of markers on the bodies
		HANDLE hMapFile_Results1;			// data to get to matlab
		HANDLE hMapFile_Results2;			// data to get to matlab

		HANDLE hMapFile_q_R;				// data on which to operate
		HANDLE hMapFile_qdot_R;				// data on which to operate
		HANDLE hMapFile_qddot_R;			// data on which to operate
		HANDLE hMapFile_GRF_F_D_R;			// data on which to operate
		HANDLE hMapFile_GRF_P_D_R;			// data on which to operate
		HANDLE hMapFile_GRF_F_S_R;			// data on which to operate
		HANDLE hMapFile_GRF_P_S_R;			// data on which to operate
		HANDLE hMapFile_PT_coordB_R;		// data on which to operate
		HANDLE hMapFile_PT_bodynames_R;		// data on which to operate
		HANDLE hMapFile_NumCoords_R;		// number of coordinates in model
		HANDLE hMapFile_NMarkers_R;			// number of markers on the bodies
		HANDLE hMapFile_Results3;			// data to get to matlab
		HANDLE hMapFile_Results4;			// data to get to matlab

		HANDLE hMapFile_KeepOpen;			// keep program open
		HANDLE hMapFile_RunJob;				// run the job
		HANDLE hMapFile_NumPts;				// number of data points (rows)
		HANDLE hMapFile_ResReady;			// results are ready
		//mexPrintf("TestMessage 5");

		bool* pBuf_KeepOpen;
		bool* pBuf_RunJob;		
		int* pBuf_NumPts;
		bool* pBuf_ResReady;

		int* pBuf_NumCoords;
		int* pBuf_NMarkers;
		double* pBuf_q;
		double* pBuf_qdot;
		double* pBuf_qddot;
		double* pBuf_GRF_F_D;
		double* pBuf_GRF_P_D;
		double* pBuf_GRF_F_S;
		double* pBuf_GRF_P_S;
		double*	pBuf_PT_coordB;
		double* pBuf_PT_bodynames;		
		double* pBuf_Results1;
		double* pBuf_Results2;

		int* pBuf_NumCoords_R;
		int* pBuf_NMarkers_R;
		double* pBuf_q_R;
		double* pBuf_qdot_R;
		double* pBuf_qddot_R;
		double* pBuf_GRF_F_D_R;
		double* pBuf_GRF_P_D_R;
		double* pBuf_GRF_F_S_R;
		double* pBuf_GRF_P_S_R;
		double*	pBuf_PT_coordB_R;
		double* pBuf_PT_bodynames_R;
		double* pBuf_Results3;
		double* pBuf_Results4;

		hMapFile_q = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameq);
		hMapFile_qdot = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameqdot);
		hMapFile_qddot = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameqddot);
		hMapFile_GRF_F_D = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_F_D);
		hMapFile_GRF_P_D = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_P_D);
		hMapFile_GRF_F_S = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_F_S);
		hMapFile_GRF_P_S = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_P_S);
		hMapFile_NumCoords = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_INT, szNameNumCoords);
		hMapFile_NMarkers = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_INT, szNameNMarkers);
		hMapFile_Results1 = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameResults1);
		hMapFile_Results2 = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_MarkerPos, szNameResults2);
		hMapFile_PT_bodynames = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BODIES, szNamePT_bodynames);
		hMapFile_PT_coordB = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNamePT_coordB);

		hMapFile_q_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameq_R);
		hMapFile_qdot_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameqdot_R);
		hMapFile_qddot_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameqddot_R);
		hMapFile_GRF_F_D_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_F_D_R);
		hMapFile_GRF_P_D_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_P_D_R);
		hMapFile_GRF_F_S_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_F_S_R);
		hMapFile_GRF_P_S_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameGRF_P_S_R);
		hMapFile_NumCoords_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_INT, szNameNumCoords_R);
		hMapFile_NMarkers_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_INT, szNameNMarkers_R);
		hMapFile_Results3 = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNameResults3);
		hMapFile_Results4 = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_MarkerPos, szNameResults4);
		hMapFile_PT_bodynames_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BODIES, szNamePT_bodynames_R);
		hMapFile_PT_coordB_R = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_STATES, szNamePT_coordB_R);

		hMapFile_KeepOpen = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameKeepOpen);
		hMapFile_RunJob = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameRunJob);		
		hMapFile_NumPts = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_INT, szNameNumPts);		
		hMapFile_ResReady = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameResReady);
		//mexPrintf("TestMessage 7");

		//if (hMapFile_Data1 == NULL || hMapFile_RunJob == NULL || hMapFile_ResReady == NULL || hMapFile_Results1 == NULL || ((hMapFile_FMo == NULL || hMapFile_lMo == NULL) && (nrhs == 9))) // I don't have all File mappings here
	   //{
		  //mexErrMsgTxt("Could not create file mapping object.\n");
	   //}

		pBuf_q = (double*)MapViewOfFile(hMapFile_q, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_qdot = (double*)MapViewOfFile(hMapFile_qdot, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_qddot = (double*)MapViewOfFile(hMapFile_qddot, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_F_D = (double*)MapViewOfFile(hMapFile_GRF_F_D, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_P_D = (double*)MapViewOfFile(hMapFile_GRF_P_D, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_F_S = (double*)MapViewOfFile(hMapFile_GRF_F_S, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_P_S = (double*)MapViewOfFile(hMapFile_GRF_P_S, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_NumCoords = (int*)MapViewOfFile(hMapFile_NumCoords, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_INT);
		pBuf_NMarkers = (int*)MapViewOfFile(hMapFile_NMarkers, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_INT);
		pBuf_Results1 = (double*)MapViewOfFile(hMapFile_Results1, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_Results2 = (double*)MapViewOfFile(hMapFile_Results2, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_MarkerPos);
		pBuf_PT_bodynames = (double*)MapViewOfFile(hMapFile_PT_bodynames, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BODIES);
		pBuf_PT_coordB = (double*)MapViewOfFile(hMapFile_PT_coordB, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);

		pBuf_q_R = (double*)MapViewOfFile(hMapFile_q_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_qdot_R = (double*)MapViewOfFile(hMapFile_qdot_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_qddot_R = (double*)MapViewOfFile(hMapFile_qddot_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_F_D_R = (double*)MapViewOfFile(hMapFile_GRF_F_D_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_P_D_R = (double*)MapViewOfFile(hMapFile_GRF_P_D_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_F_S_R = (double*)MapViewOfFile(hMapFile_GRF_F_S_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_GRF_P_S_R = (double*)MapViewOfFile(hMapFile_GRF_P_S_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_NumCoords_R = (int*)MapViewOfFile(hMapFile_NumCoords_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_INT);
		pBuf_NMarkers_R = (int*)MapViewOfFile(hMapFile_NMarkers_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_INT);
		pBuf_Results3 = (double*)MapViewOfFile(hMapFile_Results3, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);
		pBuf_Results4 = (double*)MapViewOfFile(hMapFile_Results4, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_MarkerPos);
		pBuf_PT_bodynames_R = (double*)MapViewOfFile(hMapFile_PT_bodynames_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BODIES);
		pBuf_PT_coordB_R = (double*)MapViewOfFile(hMapFile_PT_coordB_R, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_STATES);

		pBuf_KeepOpen = (bool*)MapViewOfFile(hMapFile_KeepOpen, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
		pBuf_RunJob = (bool*)MapViewOfFile(hMapFile_RunJob, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
		pBuf_ResReady = (bool*)MapViewOfFile(hMapFile_ResReady, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
		pBuf_NumPts = (int*)MapViewOfFile(hMapFile_NumPts, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_INT);

		
		// Assign all the variables to the shared memory
		for (int i = 0; i < numPts; ++i) // for all rows of data i
		{
			for (int j = 0; j < Ncoord; j++) // for all column of data state j
			{
				pBuf_q[i*Ncoord + j] = q[i + j*numPts];
				pBuf_qdot[i*Ncoord + j] = qdot[i + j*numPts];
				pBuf_qddot[i*Ncoord + j] = qddot[i + j*numPts];
			}
		}
		for (int i = 0; i < numPts; ++i) // for all rows of data i
		{
			for (int col = 0; col < 3; ++col) {
				pBuf_GRF_F_D[i*3 + col] = GRF_F_D[i + col*numPts];
				pBuf_GRF_P_D[i*3 + col] = GRF_P_D[i + col*numPts];
				pBuf_GRF_F_S[i*3 + col] = GRF_F_S[i + col*numPts];
				pBuf_GRF_P_S[i*3 + col] = GRF_P_S[i + col*numPts];
			}

		}
		for (int i = 0; i < NMarkers; i++)
		{
			pBuf_PT_bodynames[i] = PT_bodynames[i];
			for (int j = 0; j < 3; j++)
			{
				pBuf_PT_coordB[j + 3 * i] = PT_coordB[j + 3 * i];
			}
		}
		*pBuf_NumPts = numPts;
		*pBuf_NMarkers = NMarkers;
		*pBuf_NumCoords = Ncoord;


		// assign variables of the right-side model the shared memory
		for (int i = 0; i < numPts; ++i) // for all rows of data i
		{
			for (int j = 0; j < Ncoord_R; j++) // for all column of data state j
			{
				pBuf_q_R[i*Ncoord_R + j] = q_R[i + j*numPts];
				pBuf_qdot_R[i*Ncoord_R + j] = qdot_R[i + j*numPts];
				pBuf_qddot_R[i*Ncoord_R + j] = qddot_R[i + j*numPts];
			}
		}
		for (int i = 0; i < numPts; ++i) // for all rows of data i
		{
			for (int col = 0; col < 3; ++col) {
				pBuf_GRF_F_D_R[i * 3 + col] = GRF_F_D_R[i + col*numPts];
				pBuf_GRF_P_D_R[i * 3 + col] = GRF_P_D_R[i + col*numPts];
				pBuf_GRF_F_S_R[i * 3 + col] = GRF_F_S_R[i + col*numPts];
				pBuf_GRF_P_S_R[i * 3 + col] = GRF_P_S_R[i + col*numPts];
			}

		}
		for (int i = 0; i < NMarkers_R; i++)
		{
			pBuf_PT_bodynames_R[i] = PT_bodynames_R[i];
			for (int j = 0; j < 3; j++)
			{
				pBuf_PT_coordB_R[j + 3 * i] = PT_coordB_R[j + 3 * i];
			}
		}
		*pBuf_NMarkers_R = NMarkers_R;
		*pBuf_NumCoords_R = Ncoord_R;


		// Set flag RunJob to true
		*pBuf_RunJob = true;

		// wait for results
		bool printWait = true;
		while (*pBuf_ResReady == false)
		{
			if (printWait == true)
			{
				//mexPrintf("Waiting for results... \n");
				mexPrintf("");
			}
			printWait = false;
		}

		// get the resulst from the buffer
		for (int i = 0; i<numPts; ++i) // for all rows of data i
		{
			for (int j = 0; j<Ncoord; j++) // for all column of data state j
			{
				results1[i + j*numPts]= pBuf_Results1[i*Ncoord + j];
			}
			for (int j = 0; j<Ncoord_R; j++) // for all column of data state j
			{
				results3[i + j*numPts] = pBuf_Results3[i*Ncoord_R + j];
			}
			for (int j = 0; j<NMarkers; j++) // for all column of data state j
			{
				results2[i + j*numPts * 3] = pBuf_Results2[i*NMarkers * 3 + j];
				results2[i + j*numPts * 3 + numPts] = pBuf_Results2[i*NMarkers * 3 + NMarkers + j];
				results2[i + j*numPts * 3 + 2 * numPts] = pBuf_Results2[i*NMarkers * 3 + 2 * NMarkers + j];
			}
			for (int j = 0; j<NMarkers_R; j++) // for all column of data state j
			{
				results4[i + j*numPts * 3] = pBuf_Results4[i*NMarkers_R * 3 + j];
				results4[i + j*numPts * 3 + numPts] = pBuf_Results4[i*NMarkers_R * 3 + NMarkers_R + j];
				results4[i + j*numPts * 3 + 2 * numPts] = pBuf_Results4[i*NMarkers_R * 3 + 2 * NMarkers_R + j];
			}
		}
		
		*pBuf_ResReady = false;

		//mexPrintf("%d %d %d",*pBuf_NumCoords,*pBuf_NumMarkers,*pBuf_NumForces);
		//close
		UnmapViewOfFile(pBuf_q);
		UnmapViewOfFile(pBuf_qdot);
		UnmapViewOfFile(pBuf_qddot);
		UnmapViewOfFile(pBuf_GRF_F_D);
		UnmapViewOfFile(pBuf_GRF_P_D);
		UnmapViewOfFile(pBuf_GRF_F_S);
		UnmapViewOfFile(pBuf_GRF_P_S);
		UnmapViewOfFile(pBuf_PT_coordB);
		UnmapViewOfFile(pBuf_PT_bodynames);
		UnmapViewOfFile(pBuf_NumCoords);
		UnmapViewOfFile(pBuf_NMarkers);
		UnmapViewOfFile(pBuf_Results1);
		UnmapViewOfFile(pBuf_Results2);

		UnmapViewOfFile(pBuf_q_R);
		UnmapViewOfFile(pBuf_qdot_R);
		UnmapViewOfFile(pBuf_qddot_R);
		UnmapViewOfFile(pBuf_GRF_F_D_R);
		UnmapViewOfFile(pBuf_GRF_P_D_R);
		UnmapViewOfFile(pBuf_GRF_F_S_R);
		UnmapViewOfFile(pBuf_GRF_P_S_R);
		UnmapViewOfFile(pBuf_PT_coordB_R);
		UnmapViewOfFile(pBuf_PT_bodynames_R);
		UnmapViewOfFile(pBuf_NumCoords_R);
		UnmapViewOfFile(pBuf_NMarkers_R);
		UnmapViewOfFile(pBuf_Results3);
		UnmapViewOfFile(pBuf_Results4);

		UnmapViewOfFile(pBuf_KeepOpen);
		UnmapViewOfFile(pBuf_RunJob);		
		UnmapViewOfFile(pBuf_NumPts);

		CloseHandle(hMapFile_q);
		CloseHandle(hMapFile_qdot);
		CloseHandle(hMapFile_qddot);
		CloseHandle(hMapFile_GRF_F_D);
		CloseHandle(hMapFile_GRF_P_D);
		CloseHandle(hMapFile_GRF_F_S);
		CloseHandle(hMapFile_GRF_P_S);
		CloseHandle(hMapFile_PT_coordB);
		CloseHandle(hMapFile_PT_bodynames);
		CloseHandle(hMapFile_NumCoords);
		CloseHandle(hMapFile_NMarkers);
		CloseHandle(hMapFile_Results1);
		CloseHandle(hMapFile_Results2);

		CloseHandle(hMapFile_q_R);
		CloseHandle(hMapFile_qdot_R);
		CloseHandle(hMapFile_qddot_R);
		CloseHandle(hMapFile_GRF_F_D_R);
		CloseHandle(hMapFile_GRF_P_D_R);
		CloseHandle(hMapFile_GRF_F_S_R);
		CloseHandle(hMapFile_GRF_P_S_R);
		CloseHandle(hMapFile_PT_coordB_R);
		CloseHandle(hMapFile_PT_bodynames_R);
		CloseHandle(hMapFile_NumCoords_R);
		CloseHandle(hMapFile_NMarkers_R);
		CloseHandle(hMapFile_Results3);
		CloseHandle(hMapFile_Results4);

		CloseHandle(hMapFile_KeepOpen);
		CloseHandle(hMapFile_RunJob);
		CloseHandle(hMapFile_NumPts);

	}
}
//-End MEX enterance ------------------------------------------------------

//-Start Code (FUNCTIONS) -------------------------------------------------------------

//-End Code (FUNCTIONS)---------------------------------------------------------------
