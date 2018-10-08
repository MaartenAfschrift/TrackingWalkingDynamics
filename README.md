# Tracking simulation

Tracking simulation software used in publication:

Afschrift et al.; Modulation of gluteus medius activity reflects the potential of the muscle to meet the mechanical demands during perturbed walking, Scientific Reports 2018

https://www.nature.com/articles/s41598-018-30139-9

## Install instruction

Required programs:

- Matlab (tested on R2015b)
- GPOPS II
- OpenMPI
- Opensim (tested for v 3.3)

### Install OpenMPI program:

MPI was used to parralize the simulation of the skeletal dynamics (i.e. decrease total time simulation). The state at the different collocation points are distributed between multiple cores to solve inverse dynamics.

*  compile Control_Tracking.cpp as a mex file, this files communicates between cpp and matlab:
  * Matlab command: compilamex(pwd,'Control_Tracking.cpp','C:\OpenSim33');

- I had to change some files in the opensim installation: C:\OpenSim33\sdk\include\OpenSim\Common, DebugUtilities.h and LoadOpenSimLibrary.h). You can find these files in the folder:~\ModelSoft/OpenSim_Installation

*  built Opensim-mpi libraries: Currently I added the executable, so you don't have to build the cpp project. Add cleaned up full projec tin the next version. You can find the source code in the folder: ~/Optimization/SolveID_OpenMPI/OpenSimMPI_muscle.cpp:
  *  Install open-mpi for windows ( *http://www.microsoft.com/en-us/download/details.aspx?id=41634*)
  *  Install it to C:/msmpi
  *  Open cmakelist and change path MPI_HEADERS and MPI_LIBS_DIR (depending on where you installen msmpi)
  *  cmake to folder mpi_install
  *  built somewhere

* open command line, go to Release folder and run mpi program with following command:
  * *mpiexec -np 2 opensimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim*           
  * (When running the examples, the program will ask you to run this command)    

### Compile mex files

You will have to compile some mex files that are used to communicate between matlab and opensim. You can find these files in the folder ~/Optimization/Mex and use the function compilamex to compile all the files.

* Control_Tracking.cpp          (Solving inverse dynamics with OpenMPI)
* Kinematics_at_feet.cpp      (Computation kinematics )
* SwitchModels.cpp              (Reload models in OpenMPI program)
* TrovaIDePath_MD               (Used for solving inverse dynamics without OpenMPI)

(example of matlab command: compilamex(pwd,'Control_Tracking.cpp','C:\OpenSim33') )

  ### Example

The file Example Tracking contains the code for the tracking simulation of the gait2392 model with a locked mtp joint (subtalar joint is open).





​                    

   