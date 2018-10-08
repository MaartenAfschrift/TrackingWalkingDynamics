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

int main()
{
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// INPUT:
//		  0 - ModelName
//		  1 - angles
//		  2 - angular velocity
//		  3 - accelerations
//		  4 - Tracing points for Bodies
//		  5 - Names of Bodies for tracking
//		  6 - GRF right
//		  7 - Pisition GRF right
//		  8 - GRF left
//		  9 - position GRF left
//		 10 - Platform acceleration
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Load the model
	char *osimFileName_ptr;
	osimFileName_ptr = mxArrayToString(prhs[0]);
	std::string osimFileName = std::string(&osimFileName_ptr[0]);
	Model osimModel = Model(osimFileName);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Get the coordinates
	double *posiz;
	const int *dimPosiz;
	int Cposiz, Rposiz;
	posiz = mxGetPr(prhs[1]);
	dimPosiz = mxGetDimensions(prhs[1]);
	Cposiz = *(dimPosiz + 1);
	Rposiz = *(dimPosiz + 0);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Get the velocities
	double *veloc;
	veloc = mxGetPr(prhs[2]);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	get the accelerations
	double *accel;
	accel = mxGetPr(prhs[3]);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Get the coordinates for the tracking bodies
	double *PT_coordB;
	const int *dimPT_coordB;
	int CdimPT_coordB, RdimPT_coordB;
	PT_coordB = mxGetPr(prhs[4]);
	dimPT_coordB = mxGetDimensions(prhs[4]);
	CdimPT_coordB = *(dimPT_coordB + 1);
	RdimPT_coordB = *(dimPT_coordB + 0);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Get the body names to track
	Array<std::string> PT_bodynames;
	char *PT_corpi_ptr;
	for (int c = 0; c < CdimPT_coordB; c++)
	{
		PT_corpi_ptr = mxArrayToString(mxGetCell(prhs[5], c));
		PT_bodynames.append(PT_corpi_ptr);
		//mexPrintf("%s\n", PT_bodynames.get(c));
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Get the ground reaction forces
	double *GRF_F_D;
	GRF_F_D = mxGetPr(prhs[6]);
	double *GRF_P_D;
	GRF_P_D = mxGetPr(prhs[7]);
	double *GRF_F_S;
	GRF_F_S = mxGetPr(prhs[8]);
	double *GRF_P_S;
	GRF_P_S = mxGetPr(prhs[9]);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Get the platform acceleration
	double *Platform_Acc;
	Platform_Acc = mxGetPr(prhs[10]);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Create output matrix
	plhs[1] = mxCreateDoubleMatrix(Rposiz * 3, CdimPT_coordB, mxREAL);
	double *PTinGout = mxGetPr(plhs[1]);
	Vec3 &PTinG = Vec3(0.0);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	OpenSim
	//		---		Get variables
	SimTK::State &sss = osimModel.initSystem();
	const SimTK::MultibodySystem* mbs = &osimModel.getMultibodySystem();
	const SimbodyEngine& SBE = osimModel.getSimbodyEngine();
	SimTK::Vector posizvector = Vector(Cposiz, 0.0);//SimTK::Vector posizvector = sss.getQ();
	SimTK::Vector velocvector = Vector(Cposiz, 0.0);
	SimTK::Vector accelvector = Vector(Cposiz, 0.0);
	int nb = osimModel.getNumBodies();
	SimTK::Vec3 g = osimModel.getGravity();
	SimTK::Array_<SimTK::Vector> jointTorques(Rposiz, SimTK::Vector(Cposiz, 0.0));
	SimTK::Vector q(Cposiz, 0.0);
	SimTK::Vector u(Cposiz, 0.0);
	SimTK::Vector udot(Cposiz, 0.0);
	const SimTK::SimbodyMatterSubsystem& SMS = osimModel.getMatterSubsystem();
	SimTK::MobilizedBody mobod = SMS.getGround();;
	SimTK::Vector_<SimTK::SpatialVec> totCorForce_bodyFrame; //coriolis and gyroscopic wrenches of each body expressed in body origin
	SimTK::Vector_<SimTK::SpatialVec> gravForces_bodyFrame; //gravitational forces in the body frame
	SimTK::Vector_<SimTK::SpatialVec> GRF_D_bodyFrame; //modeled GRF at heel in body frame
	SimTK::Vector_<SimTK::SpatialVec> GRF_S_bodyFrame; //modeled GRF at heel in body frame
	SimTK::Vector MobilityForces_bodyFrame; //moobility forces always 0
	totCorForce_bodyFrame.resize(nb);
	gravForces_bodyFrame.resize(nb);
	GRF_D_bodyFrame.resize(nb);
	GRF_S_bodyFrame.resize(nb);
	MobilityForces_bodyFrame.resize(Cposiz);
	SimTK::Vector_<SimTK::SpatialVec> totF_frame(SMS.getNumBodies(), SimTK::SpatialVec());
	SimTK::Vector totF_system(Cposiz, 0.0);
	SimTK::Vector m_udot(Cposiz, 0.0);
	const BodySet& BS = osimModel.getBodySet();


	SimTK::Vector_<SimTK::SpatialVec> AccelerationForce_bodyFrame;
	AccelerationForce_bodyFrame.resize(nb);
	//mexPrintf("Value: %f\n", Platform_Acc[0]);

	// Get the joint torques at each collocation point
	for (int i = 0; i < Rposiz; i++) 
	{
		GRF_S_bodyFrame.setToZero();
		GRF_D_bodyFrame.setToZero();
		MobilityForces_bodyFrame.setToZero();
		totCorForce_bodyFrame.setToZero();
		gravForces_bodyFrame.setToZero();
		SimTK::Vec3 F_GRF_D;
		SimTK::Vec3 F_GRF_S;
		SimTK::Vec3 P_GRF_D;
		SimTK::Vec3 P_GRF_S;
		SimTK::Vec3 Pacc_sel;	
		AccelerationForce_bodyFrame.setToZero();
		// Read GRF forces at the current frame. If they are not given, it just supply zeros
		for (int col = 0; col < 3; ++col) {
			if (true) {
				
				F_GRF_D[col] = GRF_F_D[i + col*Rposiz];
				F_GRF_S[col] = GRF_F_S[i + col*Rposiz];
				P_GRF_D[col] = GRF_P_D[i + col*Rposiz];
				P_GRF_S[col] = GRF_P_S[i + col*Rposiz];
			}
			else {
				F_GRF_D[col] = 0;
				F_GRF_S[col] = 0;
				P_GRF_D[col] = 0;
				P_GRF_S[col] = 0;
			}			
		}

		// add the perturbation force
		Pacc_sel[0] = 0;
		Pacc_sel[1] = 0;
		Pacc_sel[2] = Platform_Acc[i];


		for (int j = 0; j < Cposiz; j++)
		{
			posizvector[j] = posiz[i + j*Rposiz];
			velocvector[j] = veloc[i + j*Rposiz];
			accelvector[j] = accel[i + j*Rposiz];
		}
		sss.setQ(posizvector);
		sss.setU(velocvector);
		mbs->realize(sss, SimTK::Stage::Dynamics);
		// For all collocation points
		for (int c = 0; c < CdimPT_coordB; c++)
		{
			Vec3 &PTinB = Vec3(PT_coordB[0 + 3 * c], PT_coordB[1 + 3 * c], PT_coordB[2 + 3 * c]);
			SBE.getPosition(sss, BS.get(PT_bodynames.get(c)), PTinB, PTinG);
			PTinGout[i + c*Rposiz * 3] = PTinG.get(0);
			PTinGout[i + Rposiz + c*Rposiz * 3] = PTinG.get(1);
			PTinGout[i + Rposiz * 2 + c*Rposiz * 3] = PTinG.get(2);
		}
		// Loop over all the bodies
		for (SimTK::MobilizedBodyIndex mbi(1); mbi < SMS.getNumBodies(); ++mbi)
		{
			mobod = SMS.getMobilizedBody(mbi);
			SMS.addInStationForce(sss, mbi, mobod.getBodyMassCenterStation(sss), mobod.getBodyMass(sss)*g, gravForces_bodyFrame);	// Maarten: Here we will have to add the perturbation force (at the location of the center of mass)
			if (mbi == BS.getIndex("pelvis"))
			{
				SMS.addInStationForce(sss, mbi, P_GRF_D, F_GRF_D, GRF_D_bodyFrame);
				SMS.addInBodyTorque(sss, mbi, F_GRF_S, GRF_S_bodyFrame);
			}
			SMS.addInStationForce(sss, mbi, mobod.getBodyMassCenterStation(sss), mobod.getBodyMass(sss) * Pacc_sel, AccelerationForce_bodyFrame); // remark => Do we have to convert this to world frame first ?
			//mexPrintf(AccelerationForce_bodyFrame);
			//mexPrintf("Value: %f\n", Pacc_sel[2]);
		}
		totF_frame = GRF_D_bodyFrame + GRF_S_bodyFrame + gravForces_bodyFrame + AccelerationForce_bodyFrame;
		SMS.calcResidualForceIgnoringConstraints(sss, MobilityForces_bodyFrame, totF_frame, accelvector, jointTorques[i]);
	}
	//mexEvalString("drawnow;");
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	OUTPUT
	// 0 - matrice ID
	// 1 - coordinate in Ground dei punti da seguire -> definito sopra
	plhs[0] = mxCreateDoubleMatrix(Rposiz, Cposiz, mxREAL);
	double*qqqID = mxGetPr(plhs[0]);
	for (int r = 0; r < Rposiz; r++)
	{
		for (int c = 0; c < Cposiz; c++)
		{
			qqqID[r + c*(Rposiz)] = jointTorques[r][c];
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return;
}