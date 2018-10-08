%% Test switch models script
addpath(genpath('C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Software'));

PT_marker_R=ones(3,1);
IndBodies_R=2;
PT_marker_L=ones(3,1);
IndBodies_L=2;

Q_L_Xcpp=zeros(100,16);   U_L_Xcpp=zeros(100,16);   A_L_Xcpp=zeros(100,16);
nfr=length(Q_L_Xcpp(:,1));
ZERO_F=zeros(nfr,3);
Q_R_Xcpp=zeros(100,13);   U_R_Xcpp=zeros(100,13);   A_R_Xcpp=zeros(100,13);
Pacc=zeros(100,1);

MPI_path='C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Executables\IDMPI_BodyKin_perturb\build\Release';
disp('Run MPI program and hit ENTER');
disp(['cd ' MPI_path]);
disp('mpiexec -np 1 GetIDTorques_Perturb.exe LeftSideModel.osim RightSideModel.osim');
pause();

[T_cpp_L,~,T_cpp_R,~,L_PelvisQ,L_PelvisPos,R_PelvisQ,R_PelvisPos]=...
    Control_Tracking_BodyKin_Perturb('notneeded',Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,PT_marker_L,IndBodies_L,ZERO_F,ZERO_F,...
    Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,PT_marker_R,IndBodies_R,ZERO_F,ZERO_F,Pacc);

ModelsPath='C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Data\temp_debug\TestSwitchModels';
copyfile(fullfile(ModelsPath,'Models2','LeftSideModel.osim'),fullfile(MPI_path,'LeftSideModel.osim'));
copyfile(fullfile(ModelsPath,'Models2','RightSideModel.osim'),fullfile(MPI_path,'RightSideModel.osim'));
SwitchModels();

[T_cpp_L2,~,T_cpp_R,~,L_PelvisQ,L_PelvisPos,R_PelvisQ,R_PelvisPos]=...
    Control_Tracking_BodyKin_Perturb('notneeded',Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,PT_marker_L,IndBodies_L,ZERO_F,ZERO_F,...
    Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,PT_marker_R,IndBodies_R,ZERO_F,ZERO_F,Pacc);

ModelsPath='C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Data\temp_debug\TestSwitchModels';
copyfile(fullfile(ModelsPath,'Models1','LeftSideModel.osim'),fullfile(MPI_path,'LeftSideModel.osim'));
copyfile(fullfile(ModelsPath,'Models1','RightSideModel.osim'),fullfile(MPI_path,'RightSideModel.osim'));
SwitchModels();


[T_cpp_L3,~,T_cpp_R,~,L_PelvisQ,L_PelvisPos,R_PelvisQ,R_PelvisPos]=...
    Control_Tracking_BodyKin_Perturb('notneeded',Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,PT_marker_L,IndBodies_L,ZERO_F,ZERO_F,...
    Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,PT_marker_R,IndBodies_R,ZERO_F,ZERO_F,Pacc);

figure();bar([T_cpp_L(1,:); T_cpp_L2(1,:); T_cpp_L3(1,:)]');


% ok => switch model function werkt !!