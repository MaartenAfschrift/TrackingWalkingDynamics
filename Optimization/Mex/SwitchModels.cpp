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

	DWORD BUF_SIZE_BOOL = sizeof(bool);


	// Names that provide key to unlocking shared memory	
	TCHAR szNameSwitchModel[] = TEXT("Global\\SwitchModels");	
	//TCHAR szNameResReady[] = TEXT("Global\\ResReady");
	HANDLE hMapFile_SwitchModel;				// switch models
	//HANDLE hMapFile_ResReady;			// results are ready


	//mexPrintf("TestMessage 5");

	bool* pBuf_SwitchModel;
	//bool* pBuf_ResReady;

	hMapFile_SwitchModel = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameSwitchModel);
	//hMapFile_ResReady = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameResReady);
	
	pBuf_SwitchModel = (bool*)MapViewOfFile(hMapFile_SwitchModel, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
	//pBuf_ResReady = (bool*)MapViewOfFile(hMapFile_ResReady, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);

	// Set flag switch the models
	*pBuf_SwitchModel = true;

	// wait for results
	/*
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

	*pBuf_ResReady		= false;	
	*pBuf_SwitchModel   = false;		//Is already in MPI program
	*/

	//mexPrintf("%d %d %d",*pBuf_NumCoords,*pBuf_NumMarkers,*pBuf_NumForces);
	//close
	UnmapViewOfFile(pBuf_SwitchModel);
	CloseHandle(pBuf_SwitchModel);
}
//-End MEX enterance ------------------------------------------------------

//-Start Code (FUNCTIONS) -------------------------------------------------------------

//-End Code (FUNCTIONS)---------------------------------------------------------------
