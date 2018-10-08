%% Example Tracking Simulation


% Tracking of motion capture of skeletal model with ground contact model
% using direct collocation (Afschrift et al. 2018,
% https://doi.org/10.1038/s41598-018-30139-9)

% Maarten Afschrift, 8 October 2018


%% 0) Input information

% Path to inputs
InstallFolder   = 'C:\Users\u0088756\Documents\Phd\software\TrackingSim_Simtk';     % Installation folder (everything is relative to this path)
ExampleFolder   = fullfile(InstallFolder,'Data\ExampleSubtalar');                   % folder with example data
S.ModelFile     = fullfile(ExampleFolder,'OsimModel_subtalar.osim');                % OpenSim Model (.osim)
S.IK_File       = fullfile(ExampleFolder,'OsimModel_subtalar.osim');                % inverse kinematics solution (.mot)
S.GRFFile       = fullfile(ExampleFolder,'OsimModel_subtalar.osim');                % ground reaction forces (.mot)
S.trcFile       = fullfile(ExampleFolder,'OsimModel_subtalar.osim');                % marker coordinates (.trc)
S.ID_File       = fullfile(ExampleFolder,'OsimModel_subtalar.osim');                % inverse dynamics solution(.sto)
S.OutPath       = fullfile(ExampleFolder,'TrackingResults');                        % folder with output

S.ifPause       = 1;                                                                % (1) codes ask to run MPI programs

%% 1) Run first version of tracking simulation         (TrackingSimulation)

% based on Lorenzo Pitto & Maarten Afschrift their code to convert the .osim model to a left and
% right leg model and run tracking simulation on a rough mesh.

if S.Run_GuessTracking
    MPI_path=fullfile(InstallFolder,'Optimization','SolveID_OpenMPI');
    MexPath=fullfile(InstallFolder,'Optimization','Mex');   addpath(genpath(MexPath));
    addpath(genpath(InstallFolder));
    addpath(genpath(fullfile(InstallFolder,'Data','ModelTemplates')));
    
    OutName='Opt_Guess_subtalar';
    disp(' Open Models in command line as administrator and hit ENTER when finished');
    disp([' cd ' MPI_path]);
    disp(' mpiexec -np 2 opensimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim');
    
    if S.ifPause
        pause();        
    end
    
    [output] = TrackingSim_ID_ImprovedContact_Subtalar(S,MainPath,MPI_path,OutName);   
    close all;
end

%% 3) Improve mesh of tracking simulation           (TrackingSimulation)

% refine the collocation mesh based on the previous solution

if S.Run_Tracking
    MPI_path=fullfile(MainPath,'Executables\opensimmpi_muscle\build2\Release');
    MexPath=fullfile(MainPath,'Software','MexFiles');   addpath(genpath(MexPath));
    addpath(genpath(fullfile(MainPath,'Software','TrackingSimulation')));
    disp(' Open Models in command line as administrator and hit ENTER when finished');
    disp([' cd ' MPI_path]);
    disp(' mpiexec -np 12 opensimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim');
    if S.ifPause
        %pause();
        S.ifPause_copy=S.ifPause;		S.ifPause=0;	% models switch automatically
    end
    for s=S.sVect
        S.SubjName          = PathInfo.SubjNames{s};
        S.OutFolderName     = PathInfo.SubjNames{s};
        S.SubjNumber        = s;
        InName='Opt_Guess_subtalar';
        ContactProp.K=5;
        ContactProp.c=1;
        ContactProp.vt=0.1;
        ContactProp.u=1;
        OutName='TrackingSim_subtalar';
        ImproveMesh_TempPelvis_Contact_interpolate_subtalar(S,MainPath,MPI_path,InName,OutName,ContactProp);
    end
end



