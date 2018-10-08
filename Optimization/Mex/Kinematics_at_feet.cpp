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
//		  0 - nome modello
//		  1 - dati posizione coordinate
//		  2 - dati velocità coordinate
//		  3 - nomi corpi di apparteneza dei punti da seguire
//		  4 - coordinate in Body dei punti da seguire
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//	0 - APRI MODELLO
	char *osimFileName_ptr;
	osimFileName_ptr = mxArrayToString(prhs[0]);
	std::string osimFileName = std::string(&osimFileName_ptr[0]);
	Model osimModel = Model(osimFileName);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	I 1 - IMPORTA POSIZIONE COORDINATE
	double *posiz;
	const int *dimPosiz;
	int Cposiz, Rposiz;
	posiz = mxGetPr(prhs[1]);
	dimPosiz = mxGetDimensions(prhs[1]);
	Cposiz = *(dimPosiz + 1);
	Rposiz = *(dimPosiz + 0);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	I 2 - IMPORTA VELOCITÀ COORDINATE
	double *veloc;
	veloc = mxGetPr(prhs[2]);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	I 4 - IMPORTA NOMI CORPI DI APPARTENENZA DEI PUNTI DA SEGUIRE
	double *CoordPunti;
	const int *dimPunti;
	int Npunti;
	CoordPunti = mxGetPr(prhs[4]);
	dimPunti = mxGetDimensions(prhs[4]);
	Npunti = *(dimPunti + 1);
	//mexPrintf("%i\n", Npunti);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	I 3 - IMPORTA NOMI CORPI DI APPARTENENZA DEI PUNTI DA SEGUIRE
	Array<std::string> PT_bodynames;
	char *PT_corpi_ptr;
	for (int c = 0; c < Npunti; c++)
	{
		PT_corpi_ptr = mxArrayToString(mxGetCell(prhs[3], c));
		PT_bodynames.append(PT_corpi_ptr);
		//mexPrintf("%s\n", PT_bodynames.get(c));
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	O 0 - CREA IL PUNTATORE ALL'OUTPUT DELLA POSIZIONE DEI PUNTI
	plhs[0] = mxCreateDoubleMatrix(Rposiz * 3, Npunti, mxREAL);
	double *PTinGout = mxGetPr(plhs[0]);
	Vec3 &PTinG = Vec3(0.0);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	O 1 - CREA IL PUNTATORE ALL'OUTPUT DELLA VELOCITÀ DEI PUNTI
	plhs[1] = mxCreateDoubleMatrix(Rposiz * 3, Npunti, mxREAL);
	double *PTVELinGout = mxGetPr(plhs[1]);
	Vec3 &PTVELinG = Vec3(0.0);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	O 2 - CREA IL PUNTATORE ALL'OUTPUT DEGLI ANGOLI DEI PUNTI
	plhs[2] = mxCreateDoubleMatrix(Rposiz * 3, Npunti, mxREAL);
	double *PTANGinGout = mxGetPr(plhs[2]);
	double DirCos[3][3];
	double rE[3];
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	O 3 - CREA IL PUNTATORE ALL'OUTPUT DELLA VELOCITÀ ANGOLARE DEI PUNTI
	plhs[3] = mxCreateDoubleMatrix(Rposiz * 3, Npunti, mxREAL);
	double *PTANGVELinGout = mxGetPr(plhs[3]);
	Vec3 &AngVel = Vec3(0.0);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	USA OPENSIM
	// ¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤  //
	//		---		QUESTA PARTE USATA PER DEFINIRE LE VARIABILI CHE SERVONO DOPO
	SimTK::State &sss = osimModel.initSystem();
	const SimTK::MultibodySystem* mbs = &osimModel.getMultibodySystem();
	const SimbodyEngine& SBE = osimModel.getSimbodyEngine();
	SimTK::Vector posizvector = Vector(Cposiz, 0.0);//SimTK::Vector posizvector = sss.getQ();
	SimTK::Vector velocvector = Vector(Cposiz, 0.0);
	int nb = osimModel.getNumBodies();
	SimTK::Vector q(Cposiz, 0.0);
	SimTK::Vector u(Cposiz, 0.0);
	// ¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤  //
	// CICLO SUI COLLOCATION POINTS
	for (int i = 0; i < Rposiz; i++) // Rposiz sono le righe della matrice delle coordinate che corrispondono ai varii istanti di tempo
	{
		for (int j = 0; j < Cposiz; j++)
		{
			posizvector[j] = posiz[i + j*Rposiz];
			velocvector[j] = veloc[i + j*Rposiz];
		}
		sss.setQ(posizvector);
		sss.setU(velocvector);
		mbs->realize(sss, SimTK::Stage::Velocity);
		// CICLO SUI PUNTI DA CALCOLARE 
		for (int c = 0; c < Npunti; c++)
		{
			Vec3 TMPpt = Vec3(CoordPunti[0 + c * 3], CoordPunti[1 + c * 3], CoordPunti[2 + c * 3]);
			SBE.getPosition(sss, osimModel.getBodySet().get(PT_bodynames.get(c)), TMPpt, PTinG);
			SBE.getVelocity(sss, osimModel.getBodySet().get(PT_bodynames.get(c)), TMPpt, PTVELinG);
			SBE.getDirectionCosines(sss, osimModel.getBodySet().get(PT_bodynames.get(c)), DirCos);
			SBE.getAngularVelocity(sss, osimModel.getBodySet().get(PT_bodynames.get(c)), AngVel);
			//SBE.getAngularVelocityBodyLocal(sss, osimModel.getBodySet().get(PT_bodynames.get(c)), AngVel);
			
			/*
			for (int www = 0; www < 3; www++)
			{
				mexPrintf("%f\t", CoordPunti[www + c * 3]);
			}
			mexPrintf("\n");
			*/

			SBE.convertDirectionCosinesToAngles(DirCos, &rE[0], &rE[1], &rE[2]);
			PTinGout[i + c*Rposiz * 3] = PTinG.get(0);
			PTinGout[i + Rposiz + c*Rposiz * 3] = PTinG.get(1);
			PTinGout[i + Rposiz * 2 + c*Rposiz * 3] = PTinG.get(2);
			PTVELinGout[i + c*Rposiz * 3] = PTVELinG.get(0);
			PTVELinGout[i + Rposiz + c*Rposiz * 3] = PTVELinG.get(1);
			PTVELinGout[i + Rposiz * 2 + c*Rposiz * 3] = PTVELinG.get(2);
			PTANGinGout[i + c*Rposiz * 3] = rE[0];
			PTANGinGout[i + Rposiz + c*Rposiz * 3] = rE[1];
			PTANGinGout[i + Rposiz * 2 + c*Rposiz * 3] = rE[2];
			PTANGVELinGout[i + c*Rposiz * 3] = AngVel.get(0);
			PTANGVELinGout[i + Rposiz + c*Rposiz * 3] = AngVel.get(1);
			PTANGVELinGout[i + Rposiz * 2 + c*Rposiz * 3] = AngVel.get(2);
		}
	}
}