#include <OpenSim/OpenSim.h>
#include <windows.h> // windows libraries for shared memory
#include <stdio.h>
#include <conio.h>
#include <tchar.h>
#include <iostream> // input-output libraries e.g. cout, cin
#include <mpi.h> // MPI libraries
#include <math.h> // Math libraries
#include <InverseDynamicsSolver.h>
#include <string>

// size allocated to shared memory
#define MAX_STATES 50 
#define MAX_PTS 900
#define NUM_OPTIONS 4
#define MAX_MARKERS 100

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
TCHAR szNameTime[] = TEXT("Global\\Time");
TCHAR szNameKeepOpen[] = TEXT("Global\\KeepOpen");
TCHAR szNameRunJob[] = TEXT("Global\\RunJob");
TCHAR szNameNumPts[] = TEXT("Global\\NumPts");
TCHAR szNameNumCoords[] = TEXT("Global\\NumCoords");
TCHAR szNameNMarkers[] = TEXT("Global\\NMarkers");
TCHAR szNameResults1[] = TEXT("Global\\Results1");
TCHAR szNameResults2[] = TEXT("Global\\Results2");
TCHAR szNameResReady[] = TEXT("Global\\ResReady");


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

using namespace OpenSim;
using namespace SimTK;
using namespace std; // included here for using cout easily

static Model osimModel;
static Model osimModel_R;

typedef vector<int> intvec;


void TestFunction(int nFr, int nCoords, double* jointTorques, double* PTinGout, int CdimPT_coordB, SimTK::State& sss, const SimTK::MultibodySystem* mbs, const SimbodyEngine& SBE)
{
	for (int i = 0; i < nFr; i++)
	{
		for (int c = 0; c < CdimPT_coordB; c++)
		{
			PTinGout[i*CdimPT_coordB * 3 + c] = 0;
			PTinGout[i*CdimPT_coordB * 3 + CdimPT_coordB + c] = 0;
			PTinGout[i*CdimPT_coordB * 3 + CdimPT_coordB * 2 + c] = 0;
		}
		for (int j = 0; j < nCoords; j++)
		{
			jointTorques[i*nCoords + j] = 0;
		}

	}
}

void GetID_Torques(int nFr, int nCoords, double* const q, double* const qdot, double* const qddot, double* const GRF_F_D, double* const GRF_P_D,
	double* const GRF_F_S, double* const GRF_P_S, double* jointTorques, double* PTinGout, int CdimPT_coordB, double* PT_coordB, double* PT_bodynames,
	OpenSim::Model* ModelInput, SimTK::State sss, const SimTK::MultibodySystem* mbs, const SimbodyEngine& SBE)
{
	Vec3 &PTinG = Vec3(0.0);
	// create vector with kinematic information

	SimTK::Vector qvector = Vector(nCoords, 0.0);//SimTK::Vector qvector = sss.getQ();
	SimTK::Vector qdotvector = Vector(nCoords, 0.0);
	SimTK::Vector qddotvector = Vector(nCoords, 0.0);
	int nb = ModelInput->getNumBodies();
	SimTK::Vec3 g = ModelInput->getGravity();
	
	// create vector with forces
	const SimTK::SimbodyMatterSubsystem& SMS = ModelInput->getMatterSubsystem();
	SimTK::MobilizedBody mobod = SMS.getGround();
	SimTK::Vector_<SimTK::SpatialVec> totCorForce_bodyFrame; //coriolis and gyroscopic wrenches of each body expressed in body origin
	SimTK::Vector_<SimTK::SpatialVec> gravForces_bodyFrame; //gravitational forces in the body frame
	SimTK::Vector_<SimTK::SpatialVec> GRF_D_bodyFrame; //modeled GRF at heel in body frame
	SimTK::Vector_<SimTK::SpatialVec> GRF_S_bodyFrame; //modeled GRF at heel in body frame
	SimTK::Vector MobilityForces_bodyFrame; //moobility forces always 0
	totCorForce_bodyFrame.resize(nb);
	gravForces_bodyFrame.resize(nb);
	GRF_D_bodyFrame.resize(nb);
	GRF_S_bodyFrame.resize(nb);
	MobilityForces_bodyFrame.resize(nCoords);
	SimTK::Vector_<SimTK::SpatialVec> totF_frame(SMS.getNumBodies(), SimTK::SpatialVec());
	SimTK::Vector totF_system(nCoords, 0.0);
	SimTK::Vector m_udot(nCoords, 0.0);
	const BodySet& BS = ModelInput->getBodySet();
	SimTK::Vector jointTorquesVect;


	SimTK::Vec3 F_GRF_D;
	SimTK::Vec3 F_GRF_S;
	SimTK::Vec3 P_GRF_D;
	SimTK::Vec3 P_GRF_S;


	// loop over the collocation points
	for (int i = 0; i < nFr; i++) 
	{
		GRF_S_bodyFrame.setToZero();
		GRF_D_bodyFrame.setToZero();
		MobilityForces_bodyFrame.setToZero();
		totCorForce_bodyFrame.setToZero();
		gravForces_bodyFrame.setToZero();
		jointTorquesVect.setToZero();
		// Read GRF forces at the current frame. If they are not given, it just supply zeros
		for (int col = 0; col < 3; ++col) {
			if (true) {

				F_GRF_D[col] = GRF_F_D[i*3 + col];
				F_GRF_S[col] = GRF_F_S[i*3 + col];
				P_GRF_D[col] = GRF_P_D[i*3 + col];
				P_GRF_S[col] = GRF_P_S[i*3 + col];
			}
			else {
				F_GRF_D[col] = 0;
				F_GRF_S[col] = 0;
				P_GRF_D[col] = 0;
				P_GRF_S[col] = 0;
			}
		}
		for (int j = 0; j < nCoords; j++)
		{
			qvector[j] = q[i*nCoords + j];
			qdotvector[j] = qdot[i*nCoords + j];
			qddotvector[j] = qddot[i*nCoords + j];

		}
		
		sss.setQ(qvector);
		sss.setU(qdotvector);
		ModelInput->assemble(sss);
		ModelInput->getMultibodySystem().realize(sss, SimTK::Stage::Velocity);
		// Loop over the markers
		for (int c = 0; c < CdimPT_coordB; c++)
		{
			Vec3 &PTinB = Vec3(PT_coordB[0 + 3 * c], PT_coordB[1 + 3 * c], PT_coordB[2 + 3 * c]);
			SBE.getPosition(sss, BS.get(PT_bodynames[c]-1), PTinB, PTinG);
			PTinGout[i*CdimPT_coordB*3 + c] = PTinG.get(0);
			PTinGout[i*CdimPT_coordB*3 + CdimPT_coordB + c] = PTinG.get(1);
			PTinGout[i*CdimPT_coordB*3 + CdimPT_coordB*2 + c] = PTinG.get(2);
		}
				// Get the joint torque at each collocation point
		for (SimTK::MobilizedBodyIndex mbi(1); mbi < SMS.getNumBodies(); ++mbi)
		{
			mobod = SMS.getMobilizedBody(mbi);
			SMS.addInStationForce(sss, mbi, mobod.getBodyMassCenterStation(sss), mobod.getBodyMass(sss)*g, gravForces_bodyFrame);
			if (mbi == BS.getIndex("pelvis"))
			{
				SMS.addInStationForce(sss, mbi, P_GRF_D, F_GRF_D, GRF_D_bodyFrame);
				SMS.addInBodyTorque(sss, mbi, F_GRF_S, GRF_S_bodyFrame);
			}
		}
		totF_frame = GRF_D_bodyFrame + GRF_S_bodyFrame + gravForces_bodyFrame;
		SMS.calcResidualForceIgnoringConstraints(sss, MobilityForces_bodyFrame, totF_frame, qddotvector, jointTorquesVect);
		for (int j = 0; j < nCoords; j++)
		{
			jointTorques[i*nCoords + j] = jointTorquesVect[j];
		}
	}
};



// Function that computes number of contiguous elements and spacing to scatter and/or gather
//	--> I think this function divides the distributes the frames over the different cores
void countsAndDispls(int nthreads, int nRows, int nCols, int* scounts, int* displs)
{
	int splitRowsNominal = floor((double)nRows / nthreads);
	int nRowsLeft = nRows % nthreads;
	int addition = 1;
	for (int i = 0; i < nthreads; i++)
	{
		// addition to add extra row to process, when negative, no need to add extra row
		if (nRowsLeft > 0)
			addition = 1;
		else
			addition = 0;

		scounts[i] = (splitRowsNominal + addition)*nCols;
		nRowsLeft--;
		// Creating displs array by taking previous displacement and addind previous number of counts
		if (i == 0)
		{
			displs[i] = 0;
		}
		else
		{
			displs[i] = displs[i - 1] + scounts[i - 1];
		}
	}
}


int _tmain(int argc, char *argv[])

//****************************************************************************80
{

	DWORD BUF_SIZE_STATES = sizeof(double)*MAX_STATES*MAX_PTS;
	DWORD BUF_SIZE_TIME = sizeof(double)*MAX_PTS;
	DWORD BUF_SIZE_BOOL = sizeof(bool);
	DWORD BUF_SIZE_INT = sizeof(int);
	DWORD BUF_SIZE_CHAR = sizeof(char)*MAX_MARKERS;
	DWORD BUF_SIZE_BODIES = sizeof(int)*MAX_MARKERS;
	DWORD BUF_SIZE_MarkerPos = sizeof(double)*MAX_MARKERS*MAX_PTS*3;

	// default way to create an mpi program

	// read input information
	string modelName;
	modelName = argv[1];
	osimModel = Model(modelName);
	int NumCoords = osimModel.getNumCoordinates(); //(int)malloc(BUF_SIZE_INT);
	SimTK::State sss = osimModel.initSystem();
	const SimTK::MultibodySystem* mbs = &osimModel.getMultibodySystem();
	const SimbodyEngine& SBE = osimModel.getSimbodyEngine();

	string modelName_R;
	modelName_R = argv[2];
	osimModel_R = Model(modelName_R);
	int NumCoords_R = osimModel_R.getNumCoordinates(); //(int)malloc(BUF_SIZE_INT);
	SimTK::State sss_R = osimModel_R.initSystem();
	const SimTK::MultibodySystem* mbs_R = &osimModel_R.getMultibodySystem();
	const SimbodyEngine& SBE_R = osimModel_R.getSimbodyEngine();

	// shared memory handles
	HANDLE hMapFile_q;				// data on which to operate
	HANDLE hMapFile_qdot;			// data on which to operate
	HANDLE hMapFile_qddot;			// data on which to operate
	HANDLE hMapFile_GRF_F_D;		// data on which to operate
	HANDLE hMapFile_GRF_P_D;		// data on which to operate
	HANDLE hMapFile_GRF_F_S;		// data on which to operate
	HANDLE hMapFile_GRF_P_S;		// data on which to operate
	HANDLE hMapFile_PT_coordB;		// data on which to operate
	HANDLE hMapFile_PT_bodynames;		// data on which to operate
	HANDLE hMapFile_KeepOpen; // keep program open
	HANDLE hMapFile_RunJob; // run the job
	HANDLE hMapFile_NumCoords; // number of coordinates in model
	HANDLE hMapFile_NumPts; // number of data points (rows)
	HANDLE hMapFile_NMarkers; // number of markers on the bodies
	HANDLE hMapFile_Results1; // data to get to matlab
	HANDLE hMapFile_Results2; // data to get to matlab
	HANDLE hMapFile_ResReady; // results are ready


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



	bool* pBuf_KeepOpen;
	bool* pBuf_RunJob;
	bool* pBuf_ResReady;
	int* pBuf_NumPts;

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


	int tid, nthreads; // processor id and number of threads
	bool first_time = true; // Initializing first_time running program
	int remoteGoOrQuit = 0;
	int numPts = 193;
	int NMarkers = 23;
	double* PT_coordB = (double*)malloc(BUF_SIZE_STATES);
	double* PT_bodynames = (double*)malloc(BUF_SIZE_BODIES);

	int NMarkers_R = 19;
	double* PT_coordB_R = (double*)malloc(BUF_SIZE_STATES);
	double* PT_bodynames_R = (double*)malloc(BUF_SIZE_BODIES);

	// Start mpi program
	MPI_Init(&argc, &argv);

	// Request a thread id, sometimes called a "rank" from the MPI master process, which has rank or tid == 0
	MPI_Comm_rank(MPI_COMM_WORLD, &tid);

	// Request number of thread or processes launched by MPI, this should be NCPUs-1
	MPI_Comm_size(MPI_COMM_WORLD, &nthreads);

	// initializing pointers for data transfer
	int* displs_q= new int[nthreads];
	int* scounts_q = new int[nthreads];
	int* displs_GRF = new int[nthreads];
	int* scounts_GRF = new int[nthreads];
	int* displs_Out2 = new int[nthreads];
	int* scounts_Out2 = new int[nthreads];

	int* displs_q_R = new int[nthreads];
	int* scounts_q_R = new int[nthreads];
	int* displs_Out2_R = new int[nthreads];
	int* scounts_Out2_R = new int[nthreads];

	if (tid == 0)
	{
		if (first_time)
		{

			// here we have to create file mapping objects
			//		--> directly mapping between files and adress => apperently the I/O goes faster in this way
			hMapFile_KeepOpen = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameKeepOpen);
			hMapFile_RunJob = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameRunJob);
			hMapFile_NumPts = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_INT, szNameNumPts);
			hMapFile_ResReady = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BOOL, szNameResReady);

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
			hMapFile_PT_bodynames= CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, BUF_SIZE_BODIES, szNamePT_bodynames);
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

			// Assign shared memory to pointers

			pBuf_ResReady = (bool*)MapViewOfFile(hMapFile_ResReady, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
			pBuf_KeepOpen = (bool*)MapViewOfFile(hMapFile_KeepOpen, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
			pBuf_RunJob = (bool*)MapViewOfFile(hMapFile_RunJob, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BOOL);
			pBuf_NumPts = (int*)MapViewOfFile(hMapFile_NumPts, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_INT);

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
			pBuf_PT_bodynames= (double*)MapViewOfFile(hMapFile_PT_bodynames, FILE_MAP_ALL_ACCESS, 0, 0, BUF_SIZE_BODIES);
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


			// Check for Null pointers in shared memory
			if (pBuf_q == NULL || pBuf_RunJob == NULL || pBuf_KeepOpen == NULL)
			{
				_tprintf(TEXT("Could not map view of file (%d).\n"), GetLastError());
				CloseHandle(hMapFile_q);
				CloseHandle(hMapFile_KeepOpen);
				CloseHandle(hMapFile_RunJob);
				CloseHandle(hMapFile_NumCoords);
				CloseHandle(hMapFile_NumPts);
				CloseHandle(hMapFile_NMarkers);
				CloseHandle(hMapFile_Results1);
				CloseHandle(hMapFile_Results2);
				return 1;
			}
			first_time = false;
			*pBuf_NumCoords = NumCoords;
			*pBuf_NumCoords_R = NumCoords_R;
		}
		*pBuf_KeepOpen = true;
		*pBuf_RunJob = false;


		// Start while loop that runs continuously
		while (*pBuf_KeepOpen)
		{
			// Check if we need to run a job
			if (*pBuf_RunJob)
			{
				remoteGoOrQuit = 1; // tell process to go
				MPI_Bcast(&remoteGoOrQuit, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast GO information to run job
				{
					numPts = *pBuf_NumPts;
					NumCoords = *pBuf_NumCoords;
					NMarkers = *pBuf_NMarkers;
					PT_coordB = pBuf_PT_coordB;
					PT_bodynames = pBuf_PT_bodynames;

					NumCoords_R = *pBuf_NumCoords_R;
					NMarkers_R = *pBuf_NMarkers_R;
					PT_coordB_R = pBuf_PT_coordB_R;
					PT_bodynames_R = pBuf_PT_bodynames_R;

					// Calculate ncounts_Data1 and displs_Data1
					countsAndDispls(nthreads, numPts, NumCoords, scounts_q, displs_q);			// this is for q
					countsAndDispls(nthreads, numPts, 3, scounts_GRF, displs_GRF);				// this is for q
					countsAndDispls(nthreads, numPts, NMarkers*3, scounts_Out2, displs_Out2);
					
					countsAndDispls(nthreads, numPts, NumCoords_R, scounts_q_R, displs_q_R);			// this is for q
					countsAndDispls(nthreads, numPts, NMarkers_R * 3, scounts_Out2_R, displs_Out2_R);
					
					// broadcast information 
					MPI_Bcast(&numPts, 1, MPI_INT, 0, MPI_COMM_WORLD);
					
					MPI_Bcast(&NumCoords, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(&NMarkers, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(PT_bodynames, NMarkers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Bcast(PT_coordB, 3 * NMarkers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					
					MPI_Bcast(&NumCoords_R, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(&NMarkers_R, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(PT_bodynames_R, NMarkers_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Bcast(PT_coordB_R, 3 * NMarkers_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);

					// allocate memory for the variables local to process 0
					double* qLocal = new double[numPts*NumCoords];
					double* qdotLocal = new double[numPts*NumCoords];
					double* qddotLocal = new double[numPts*NumCoords];
					double* GRF_F_DLocal = new double[numPts * 3];
					double* GRF_P_DLocal = new double[numPts * 3];
					double* GRF_F_SLocal = new double[numPts * 3];
					double* GRF_P_SLocal = new double[numPts * 3];
					double* JointTorquesLocal = new double[numPts*NumCoords];
					double* JointTorquesGlobal = new double[numPts*NumCoords];
					double* PTinGoutLocal = new double[3 * NMarkers * numPts];
					double* PTinGoutGlobal = new double[3 * NMarkers* numPts];

					double* qLocal_R = new double[numPts*NumCoords_R];
					double* qdotLocal_R = new double[numPts*NumCoords_R];
					double* qddotLocal_R = new double[numPts*NumCoords_R];
					double* GRF_F_DLocal_R = new double[numPts * 3];
					double* GRF_P_DLocal_R = new double[numPts * 3];
					double* GRF_F_SLocal_R = new double[numPts * 3];
					double* GRF_P_SLocal_R = new double[numPts * 3];
					double* JointTorquesLocal_R = new double[numPts*NumCoords_R];
					double* JointTorquesGlobal_R = new double[numPts*NumCoords_R];
					double* PTinGoutLocal_R = new double[3 * NMarkers_R * numPts];
					double* PTinGoutGlobal_R = new double[3 * NMarkers_R* numPts];

					// Scatter all variables along the processes
					MPI_Scatterv(pBuf_q, scounts_q, displs_q, MPI_DOUBLE, qLocal, scounts_q[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_qdot, scounts_q, displs_q, MPI_DOUBLE, qdotLocal, scounts_q[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_qddot, scounts_q, displs_q, MPI_DOUBLE, qddotLocal, scounts_q[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_F_D, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_DLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_P_D, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_DLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_F_S, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_SLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_P_S, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_SLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);

					MPI_Scatterv(pBuf_q_R, scounts_q_R, displs_q_R, MPI_DOUBLE, qLocal_R, scounts_q_R[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_qdot_R, scounts_q_R, displs_q_R, MPI_DOUBLE, qdotLocal_R, scounts_q_R[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_qddot_R, scounts_q_R, displs_q_R, MPI_DOUBLE, qddotLocal_R, scounts_q_R[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_F_D_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_DLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_P_D_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_DLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_F_S_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_SLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Scatterv(pBuf_GRF_P_S_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_SLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);

					// perform jobs
					
					GetID_Torques(scounts_q[tid] / NumCoords, NumCoords, qLocal, qdotLocal, qddotLocal, GRF_F_DLocal, GRF_P_DLocal,
						GRF_F_SLocal, GRF_P_SLocal, JointTorquesLocal, PTinGoutLocal, NMarkers, PT_coordB, PT_bodynames, &osimModel,
						sss, mbs, SBE);
					GetID_Torques(scounts_q_R[tid] / NumCoords_R, NumCoords_R, qLocal_R, qdotLocal_R, qddotLocal_R, GRF_F_DLocal_R, GRF_P_DLocal_R,
						GRF_F_SLocal_R, GRF_P_SLocal_R, JointTorquesLocal_R, PTinGoutLocal_R, NMarkers_R, PT_coordB_R, PT_bodynames_R, &osimModel_R,
						sss_R, mbs_R, SBE_R);

					//TestFunction(scounts_q[tid] / NumCoords, NumCoords, JointTorquesLocal, PTinGoutLocal, NMarkers,sss,mbs,SBE);
					//TestFunction(scounts_q_R[tid] / NumCoords_R, NumCoords_R, JointTorquesLocal_R, PTinGoutLocal_R, NMarkers_R,sss_R, mbs_R, SBE_R);

					// get the results
					MPI_Gatherv(JointTorquesLocal, scounts_q[tid], MPI_DOUBLE,
						JointTorquesGlobal, scounts_q, displs_q, MPI_DOUBLE, 0, MPI_COMM_WORLD);

					MPI_Gatherv(PTinGoutLocal, scounts_Out2[tid], MPI_DOUBLE,
						PTinGoutGlobal, scounts_Out2, displs_Out2, MPI_DOUBLE, 0, MPI_COMM_WORLD);		// we still have a problem here => compute number of output arguments

					MPI_Gatherv(JointTorquesLocal_R, scounts_q_R[tid], MPI_DOUBLE,
						JointTorquesGlobal_R, scounts_q_R, displs_q_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);

					MPI_Gatherv(PTinGoutLocal_R, scounts_Out2_R[tid], MPI_DOUBLE,
						PTinGoutGlobal_R, scounts_Out2_R, displs_Out2_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);		// we still have a problem here => compute number of output arguments
					// Put results into shared memory buffer
					CopyMemory((PVOID)pBuf_Results1, JointTorquesGlobal, (numPts*NumCoords * sizeof(double)));
					CopyMemory((PVOID)pBuf_Results2, PTinGoutGlobal, (3 * NMarkers * numPts * sizeof(double)));
					CopyMemory((PVOID)pBuf_Results3, JointTorquesGlobal_R, (numPts*NumCoords_R * sizeof(double)));
					CopyMemory((PVOID)pBuf_Results4, PTinGoutGlobal_R, (3 * NMarkers_R * numPts * sizeof(double)));
					// De-allocating memory
					//delete[] timeLocal;
					delete[] qLocal;
					delete[] qdotLocal;
					delete[] qddotLocal;
					delete[] GRF_F_DLocal;
					delete[] GRF_P_DLocal;
					delete[] GRF_F_SLocal;
					delete[] GRF_P_SLocal;
					delete[] JointTorquesLocal;
					delete[] JointTorquesGlobal;
					delete[] PTinGoutLocal;
					delete[] PTinGoutGlobal;

					delete[] qLocal_R;
					delete[] qdotLocal_R;
					delete[] qddotLocal_R;
					delete[] GRF_F_DLocal_R;
					delete[] GRF_P_DLocal_R;
					delete[] GRF_F_SLocal_R;
					delete[] GRF_P_SLocal_R;
					delete[] JointTorquesLocal_R;
					delete[] JointTorquesGlobal_R;
					delete[] PTinGoutLocal_R;
					delete[] PTinGoutGlobal_R;

					
				}
				// Adjust flow control variables
				*pBuf_RunJob = false;
				*pBuf_ResReady = true;
			}
			// Otherwise we do nothing, and loop until (*pBuf_KeepOpen == false)
		}

		// closing instructions have been sent => close the program properly
		remoteGoOrQuit = 0; // tell process to QUIT
		MPI_Bcast(&remoteGoOrQuit, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast QUIT information to exit all processes and finish


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
		



	}		// end of tid == 0 
	else // for all other processes not root
	{
		while (1)
		{
			// blocking receive. Waiting for signal to either go or quit
			MPI_Bcast(&remoteGoOrQuit, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (remoteGoOrQuit == 1)
			{
				// First time setup for process (e.g. load model)
				if (first_time)
				{
					first_time = false;
				}
				// recieve information
				MPI_Bcast(&numPts, 1, MPI_INT, 0, MPI_COMM_WORLD);
				
				MPI_Bcast(&NumCoords, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&NMarkers, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(PT_bodynames, NMarkers, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(PT_coordB, 3 * NMarkers, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				MPI_Bcast(&NumCoords_R, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&NMarkers_R, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(PT_bodynames_R, NMarkers_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(PT_coordB_R, 3 * NMarkers_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				// allocate memory to local variables
				double* qLocal = new double[numPts*NumCoords];
				double* qdotLocal = new double[numPts*NumCoords];
				double* qddotLocal = new double[numPts*NumCoords];
				double* GRF_F_DLocal = new double[numPts * 3];
				double* GRF_P_DLocal = new double[numPts * 3];
				double* GRF_F_SLocal = new double[numPts * 3];
				double* GRF_P_SLocal = new double[numPts * 3];
				double* JointTorquesLocal = new double[numPts*NumCoords];
				double* JointTorquesGlobal = new double[numPts*NumCoords];
				double* PTinGoutLocal = new double[3 * NMarkers * numPts];
				double* PTinGoutGlobal = new double[3 * NMarkers* numPts];
		
				double* qLocal_R = new double[numPts*NumCoords_R];
				double* qdotLocal_R = new double[numPts*NumCoords_R];
				double* qddotLocal_R = new double[numPts*NumCoords_R];
				double* GRF_F_DLocal_R = new double[numPts * 3];
				double* GRF_P_DLocal_R = new double[numPts * 3];
				double* GRF_F_SLocal_R = new double[numPts * 3];
				double* GRF_P_SLocal_R = new double[numPts * 3];
				double* JointTorquesLocal_R = new double[numPts*NumCoords_R];
				double* JointTorquesGlobal_R = new double[numPts*NumCoords_R];
				double* PTinGoutLocal_R = new double[3 * NMarkers_R * numPts];
				double* PTinGoutGlobal_R = new double[3 * NMarkers_R* numPts];
				// Receive xlocal or state array local version
				pBuf_q = NULL; // to initialize it to something and not get warning. (Initialization not needed since ignored in receiving call.)
				pBuf_qdot = NULL;
				pBuf_qddot = NULL;
				pBuf_GRF_F_D = NULL;
				pBuf_GRF_P_D = NULL;
				pBuf_GRF_F_S = NULL;
				pBuf_GRF_P_S = NULL;
				
				pBuf_q_R = NULL; // to initialize it to something and not get warning. (Initialization not needed since ignored in receiving call.)
				pBuf_qdot_R = NULL;
				pBuf_qddot_R = NULL;
				pBuf_GRF_F_D_R = NULL;
				pBuf_GRF_P_D_R = NULL;
				pBuf_GRF_F_S_R = NULL;
				pBuf_GRF_P_S_R = NULL;

				countsAndDispls(nthreads, numPts, NumCoords, scounts_q, displs_q);			// this is for q
				countsAndDispls(nthreads, numPts, 3, scounts_GRF, displs_GRF);				// this is for q
				countsAndDispls(nthreads, numPts, NMarkers * 3, scounts_Out2, displs_Out2);
				countsAndDispls(nthreads, numPts, NumCoords_R, scounts_q_R, displs_q_R);			// this is for q
				countsAndDispls(nthreads, numPts, NMarkers_R * 3, scounts_Out2_R, displs_Out2_R);


				// recieve information from torques, marker location and time => we will also need information about all the other variables 
				MPI_Scatterv(pBuf_q, scounts_q, displs_q, MPI_DOUBLE, qLocal, scounts_q[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_qdot, scounts_q, displs_q, MPI_DOUBLE, qdotLocal, scounts_q[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_qddot, scounts_q, displs_q, MPI_DOUBLE, qddotLocal, scounts_q[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_F_D, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_DLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_P_D, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_DLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_F_S, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_SLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_P_S, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_SLocal, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);

				MPI_Scatterv(pBuf_q_R, scounts_q_R, displs_q_R, MPI_DOUBLE, qLocal_R, scounts_q_R[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_qdot_R, scounts_q_R, displs_q_R, MPI_DOUBLE, qdotLocal_R, scounts_q_R[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_qddot_R, scounts_q_R, displs_q_R, MPI_DOUBLE, qddotLocal_R, scounts_q_R[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_F_D_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_DLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_P_D_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_DLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_F_S_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_F_SLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Scatterv(pBuf_GRF_P_S_R, scounts_GRF, displs_GRF, MPI_DOUBLE, GRF_P_SLocal_R, scounts_GRF[tid], MPI_DOUBLE, 0, MPI_COMM_WORLD);
				
				
				
				GetID_Torques(scounts_q[tid] / NumCoords, NumCoords, qLocal, qdotLocal, qddotLocal, GRF_F_DLocal, GRF_P_DLocal,
					GRF_F_SLocal, GRF_P_SLocal, JointTorquesLocal, PTinGoutLocal, NMarkers, PT_coordB, PT_bodynames, &osimModel,
					sss, mbs, SBE);
				GetID_Torques(scounts_q_R[tid] / NumCoords_R, NumCoords_R, qLocal_R, qdotLocal_R, qddotLocal_R, GRF_F_DLocal_R, GRF_P_DLocal_R,
					GRF_F_SLocal_R, GRF_P_SLocal_R, JointTorquesLocal_R, PTinGoutLocal_R, NMarkers_R, PT_coordB_R, PT_bodynames_R, &osimModel_R,
					sss_R, mbs_R, SBE_R);
				
				//TestFunction(scounts_q[tid] / NumCoords, NumCoords, JointTorquesLocal, PTinGoutLocal, NMarkers, sss, mbs, SBE);
				//TestFunction(scounts_q_R[tid] / NumCoords_R, NumCoords_R, JointTorquesLocal_R, PTinGoutLocal_R, NMarkers_R, sss_R, mbs_R, SBE_R);

				// Send results back to process 0
				MPI_Gatherv(JointTorquesLocal, scounts_q[tid], MPI_DOUBLE,
					JointTorquesGlobal, scounts_q, displs_q, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				MPI_Gatherv(PTinGoutLocal, scounts_Out2[tid], MPI_DOUBLE,
					PTinGoutGlobal, scounts_Out2, displs_Out2, MPI_DOUBLE, 0, MPI_COMM_WORLD);		// we still have a problem here => compute number of output arguments

				MPI_Gatherv(JointTorquesLocal_R, scounts_q_R[tid], MPI_DOUBLE,
					JointTorquesGlobal_R, scounts_q_R, displs_q_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);

				MPI_Gatherv(PTinGoutLocal_R, scounts_Out2_R[tid], MPI_DOUBLE,
					PTinGoutGlobal_R, scounts_Out2_R, displs_Out2_R, MPI_DOUBLE, 0, MPI_COMM_WORLD);		// we still have a problem here => compute number of output arguments

				// De-allocating memory
				delete[] qLocal;
				delete[] qdotLocal;
				delete[] qddotLocal;
				delete[] GRF_F_DLocal;
				delete[] GRF_P_DLocal;
				delete[] GRF_F_SLocal;
				delete[] GRF_P_SLocal;
				delete[] JointTorquesLocal;
				delete[] JointTorquesGlobal;
				delete[] PTinGoutLocal;
				delete[] PTinGoutGlobal;

				delete[] qLocal_R;
				delete[] qdotLocal_R;
				delete[] qddotLocal_R;
				delete[] GRF_F_DLocal_R;
				delete[] GRF_P_DLocal_R;
				delete[] GRF_F_SLocal_R;
				delete[] GRF_P_SLocal_R;
				delete[] JointTorquesLocal_R;
				delete[] JointTorquesGlobal_R;
				delete[] PTinGoutLocal_R;
				delete[] PTinGoutGlobal_R;


			}
			else if (remoteGoOrQuit == 0)
			{
				break;
			}
		}
	}


	// Cleaning up memory
	delete[] displs_q;
	delete[] scounts_q;
	delete[] displs_GRF;
	delete[] scounts_GRF;
	delete[] displs_Out2;
	delete[] scounts_Out2;

	delete[] displs_q_R;
	delete[] scounts_q_R;
	delete[] displs_Out2_R;
	delete[] scounts_Out2_R;

	free (PT_coordB);
	free (PT_bodynames);
	free (PT_coordB_R);
	free (PT_bodynames_R);

	// closing MPI
	MPI_Finalize();


	return 0;
};