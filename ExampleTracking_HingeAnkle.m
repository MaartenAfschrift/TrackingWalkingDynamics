%% Example Tracking Simulation


% Tracking of motion capture of skeletal model with ground contact model
% using direct collocation (Afschrift et al. 2018,
% https://doi.org/10.1038/s41598-018-30139-9)

% Maarten Afschrift, 8 October 2018


%% 0) Input information

% Path to inputs
InstallFolder   = 'C:\Users\u0088756\Documents\Phd\software\TrackingSim_Simtk';     % Installation folder (everything is relative to this path)
ExampleFolder   = fullfile(InstallFolder,'Data\ExampleHinge');                   % folder with example data
S.ModelFile     = fullfile(ExampleFolder,'OsimModel.osim');                % OpenSim Model (.osim)
S.IK_File       = fullfile(ExampleFolder,'KS.mot');                                 % inverse kinematics solution (.mot)
S.GRFFile       = fullfile(ExampleFolder,'GRF_data.mot');                           % ground reaction forces (.mot)
S.trcFile       = fullfile(ExampleFolder,'MarkerFile.trc');                         % marker coordinates (.trc)
S.ID_File       = fullfile(ExampleFolder,'InverseDynamics.sto');                    % inverse dynamics solution(.sto)
S.OutPath       = fullfile(ExampleFolder,'TrackingResults_hinge');                        % folder with output

S.ifPause       = 0;                                                                % (1) codes ask to run MPI programs

% Flow control
S.Run_GuessTracking         = 1;
S.Run_Tracking              = 1;
S.Run_ComputeMomentArms     = 1;
S.Run_SolveMuscleRedundancy = 1;

% compute moment arms (polynomial fitting)
S.do_muscle_analysis=1;
S.do_polynomial_fitting=1;
S.MuscleFunctionsGeneric=1;

% redundancy solver
S.if_RunMuscleAnalaysis =   1;          % Run muscle analysis
S.LeftModel             =   1;          % Compute for left side model
S.RightModel            =   1;          % Compute for right side model        
S.RunAdigator           =   1;          % Adigator installation
S.a_min                 = 0.02;         % minimal muscle activity
S.AdigatorPath = 'C:\Users\u0088756\Documents\Software\adigator';   % installation location adigator

% settings for the mesh refinment
S.nMeshIt               = 0;
S.meshTolerance         = 1e-6;
S.NMeshPoints           = 60;
S.NCollit               = 7;
S.damping               = 5;
    


%% 1) Run first version of tracking simulation         (TrackingSimulation)

% based on Lorenzo Pitto & Maarten Afschrift their code to convert the .osim model to a left and
% right leg model and run tracking simulation on a rough mesh.

if S.Run_GuessTracking
    MPI_path=fullfile(InstallFolder,'Optimization','SolveID_OpenMPI');
    MexPath=fullfile(InstallFolder,'Optimization','Mex');   addpath(genpath(MexPath));
    addpath(genpath(InstallFolder));
    addpath(genpath(fullfile(InstallFolder,'Data','ModelTemplates')));    
    OutName='Opt_Guess';    
    if S.ifPause
        disp(' Open Models in command line as administrator and hit ENTER when finished');
        disp([' cd ' MPI_path]);
        disp(' mpiexec -np 2 opensimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim');
        pause();
    end    
    [output] = TrackingSim_ID_ImprovedContact(S,MPI_path,OutName);
    close all;
end

%% 2) Improve mesh of tracking simulation           (TrackingSimulation)

% refine the collocation mesh based on the previous solution
if S.Run_Tracking        
    
    % path information
    MPI_path=fullfile(InstallFolder,'Optimization','SolveID_OpenMPI');
    MexPath=fullfile(InstallFolder,'Optimization','Mex');   addpath(genpath(MexPath));
    addpath(genpath(InstallFolder));
    addpath(genpath(fullfile(InstallFolder,'Data','ModelTemplates')));
    
    if S.ifPause
        disp(' Open Models in command line as administrator and hit ENTER when finished');
        disp([' cd ' MPI_path]);
        disp(' mpiexec -np 12 opensimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim');       % Runs on 12 cores
        pause();
    end
    
    InName='Opt_Guess';
    ContactProp.K=5;                % Stiffness Contact model
    ContactProp.c=1;                % scale for friction foces
    ContactProp.vt=0.1;             % damping
    ContactProp.u=1;                % ? [look this one up]
    OutName='TrackingSim';
    ImproveMesh(S,MPI_path,InName,OutName,ContactProp);
end



%% 3) Compute Moment arms functions                 (FitMuscleFunctions)

% Describe muscle geometry using polynomials (needed for forward interation, not necessary for tracking simulation)

if S.Run_ComputeMomentArms
    addpath(genpath(fullfile(MainPath,'Software','FitMuscleFunctions')));
    for s=S.sVect
        S.SubjName          = PathInfo.SubjNames{s};
        S.OutFolderName     = PathInfo.SubjNames{s};
        S.SubjNumber        = s;
        S.do_muscle_analysis=1;
        GetMusclePolynomials_subtalar(S,MainPath);           % polynomial fitting
        MuscleFunctions_Left_subtalar(S,MainPath);           % export symbolic equations
        MuscleFunctions_Right_subtalar(S,MainPath);          % export symbolic equations
    end
end


%% 4) Solve muscle redundancy                       (SolveMuscleRedundancy)
if S.Run_SolveMuscleRedundancy
    MatlabData          = fullfile(MainPath,'Data');
    Adigator_Run        = fullfile(S.AdigatorPath,'startupadigator.m');
    if S.RunAdigator
        run(Adigator_Run);
    end    
    TrackingResults='TrackingSim';
    MuscleRedundancy_EMG_LimitForceDamped_subtalar(S,MainPath,MatlabData,TrackingResults);
end






