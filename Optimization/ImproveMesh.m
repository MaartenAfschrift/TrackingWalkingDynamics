function [] = ImproveMesh(S,MPI_path,InName,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



% default Contact Properties
ContactProp.K=100;
ContactProp.c=1;
ContactProp.vt=0.1;
ContactProp.u=1;

% load varargin information
OutName='Tracking_I';
if ~isempty(varargin)
   OutName=varargin{1};   
end
if length(varargin)>1
    ContactProp=varargin{2};
end

ContinuousFunction='DoubleModel_AD_CONTINUOUS_CtStiff';      % default continuous function
if length(varargin)>2
   ContinuousFunction=varargin{3}; 
end


%% Path information
OutFolder       = S.OutPath;
ResultsName     = [InName '.mat'];

%% Setup the optimal control problem
global III IIIarr BoolFirst
IIIarr=0;   III=0;

OptRes=load(fullfile(OutFolder,ResultsName));
setup=OptRes.output.result.setup;
r=OptRes.output.result.solution;

t_or=r.phase.time;
x_or=r.phase.state;
u_or=r.phase.control;

NcollNew=S.NMeshPoints*S.NCollit;
tNew=linspace(t_or(1),t_or(end),NcollNew);
x_new=interp1(t_or,x_or,tNew);
u_new=interp1(t_or,u_or,tNew);

% update initial guess
setup.guess.parameter=OptRes.output.result.solution.parameter;
setup.guess.phase.time=tNew';
setup.guess.phase.control=u_new;
setup.guess.phase.state=x_new;
setup.guess.phase.integral=OptRes.output.result.solution.phase.integral;

% update mesh tolerance
setup.mesh.phase.colpoints=ones(1,S.NMeshPoints)*S.NCollit;
setup.mesh.phase.fraction=(1/(S.NMeshPoints))*ones(1,S.NMeshPoints);
setup.mesh.tolerance=S.meshTolerance;
setup.mesh.maxiterations=S.nMeshIt;
setup.mesh.method='hp-LiuRao';
setup.bounds.phase.path.lower(13:end)=-10^-10;
setup.bounds.phase.path.upper(13:end)=10^-10;
setup.mesh.colpointsmax=5;

% change contact parameter bounds
disp('decreased stiffness used for bounds');
setup.bounds = rmfield(setup.bounds,'parameter');
setup.guess = rmfield(setup.guess,'parameter');

Ncontrols=length(setup.guess.phase.control(1,:));
Nstates=length(setup.guess.phase.state(1,:));

auxdata=OptRes.output.result.setup.auxdata;
auxdata.LEFTMODELNAME=[OutFolder,'\LeftSideModel.osim'];
auxdata.RIGHTMODELNAME=[OutFolder,'\RightSideModel.osim'];
auxdata.ContactProp=ContactProp;
if Ncontrols~=Nstates/2+13+6      % Add damping at joints
    %% Compute values of damping    
    time=OptRes.output.result.solution.phase.time;
    nfr=length(time);
    [Q_L_Xcpp,Q_R_Xcpp,U_L_Xcpp ,U_R_Xcpp,A_L_Xcpp,A_R_Xcpp,Fp,Mp,Q_L,Q_R,U_L,U_R] = Get_KinematicsInfo_TrackingSim(OptRes);
    ZERO_F=zeros(nfr,3);
    [T_cpp_L,PT_posG_cpp_L]=...
        trovaIDePATH_MD(auxdata.LEFTMODELNAME,Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,auxdata.PT_marker_L,...
        auxdata.PT_markerbody_L,-Fp,ZERO_F,-Mp,ZERO_F);
    [T_cpp_R,PT_posG_cpp_R]=...
        trovaIDePATH_MD(auxdata.RIGHTMODELNAME,Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,auxdata.PT_marker_R,...
        auxdata.PT_markerbody_R,Fp,ZERO_F,Mp,ZERO_F);
    
    qdotL=U_L_Xcpp(:,8:15);
    qdotR=U_R_Xcpp(:,8:12);
    
    TL=T_cpp_L(:,8:15);
    TR=T_cpp_R(:,8:12);
    
    damping=S.damping;
    figure();plot(TL(:,1)); hold on;  plot(qdotL(:,1).*damping,'--r');
    disp(['Damping values = ' num2str(damping)]);
    %% Add ID torques to controls and path constraint
    Ndof_damping=13;
    Tmax=500;
    setup.bounds.phase.path.lower=[setup.bounds.phase.path.lower zeros(1,Ndof_damping)];
    setup.bounds.phase.path.upper=[setup.bounds.phase.path.upper zeros(1,Ndof_damping)];
    setup.bounds.phase.control.lower=[setup.bounds.phase.control.lower -ones(1,Ndof_damping)];
    setup.bounds.phase.control.upper=[setup.bounds.phase.control.upper ones(1,Ndof_damping)];
    
    T_cpp_L_int=interp1(t_or,T_cpp_L,tNew);
    T_cpp_R_int=interp1(t_or,T_cpp_R,tNew);
    U_L_Xcppint=interp1(t_or,U_L_Xcpp,tNew);
    U_R_Xcppint=interp1(t_or,U_R_Xcpp,tNew);%     

    ID_guess=[T_cpp_L_int(:,8:15) T_cpp_R_int(:,8:12)];
    Damping_guess=[U_L_Xcppint(:,8:15) U_R_Xcppint(:,8:12)].*damping;
    ID_damped_guess=ID_guess + Damping_guess;
    setup.guess.phase.control=[setup.guess.phase.control ID_damped_guess./Tmax];
    
    auxdata.Ndof_damping=Ndof_damping;
    auxdata.Tmax=Tmax;
    auxdata.damping=damping;
end

%% Run optimization problem
setup.functions.continuous=ContinuousFunction;
setup.functions.endpoint='DoubleModel_AD_ENDPOINT';

LEFTMODELNAME=[OutFolder,'\LeftSideModel.osim'];
RIGHTMODELNAME=[OutFolder,'\RightSideModel.osim'];
copyfile(LEFTMODELNAME,[MPI_path '\LeftSideModel.osim']);
copyfile(RIGHTMODELNAME,[MPI_path '\RightSideModel.osim']);setup.auxdata=auxdata;

% change contact model properties
% auxdata.ContactInfo.VectorRadius=[0.05 0.02 0.02 0.018];

disp(' Open Models in command line as administrator and hit ENTER when finished');
disp([' cd ' MPI_path]);
disp(' mpiexec -np 2 opensimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim');
copyfile(LEFTMODELNAME,[MPI_path '\LeftSideModel.osim']);
copyfile(RIGHTMODELNAME,[MPI_path '\RightSideModel.osim']);setup.auxdata=auxdata;
SwitchModels();
pause(1);       % pause 1 second to make sure that models are loaded before starting optimization
if S.ifPause
    pause();    
end
% last adjustments
% setup.scales.method='none';

% run the optimization
BoolFirst=true;     tic;    figure();
output = gpops2(setup);

% save results
save(fullfile(OutFolder,[OutName '.mat']),'output','auxdata');
end

