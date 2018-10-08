%% Test switch models script
addpath(genpath('C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Software'));

PT_marker_R=ones(3,1);
IndBodies_R=2;
PT_marker_L=ones(3,1);
IndBodies_L=2;
nfr=100;

Q_L_Xcpp=zeros(nfr,16);   U_L_Xcpp=zeros(nfr,16);   A_L_Xcpp=zeros(nfr,16);
nfr=length(Q_L_Xcpp(:,1));
ZERO_F=zeros(nfr,3);
Q_R_Xcpp=zeros(nfr,13);   U_R_Xcpp=zeros(nfr,13);   A_R_Xcpp=zeros(nfr,13);
Pacc=zeros(nfr,1);

MPI_path='C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Executables\opensimmpi_muscle\build2\Release';
disp('Run MPI program and hit ENTER');
disp(['cd ' MPI_path]);
disp('mpiexec -np 1 OpenSimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim');
% pause();
tic
for i=1:100
    [T_cpp_L,PT_posG_cpp_L,T_cpp_R,PT_posG_cpp_R] =...
        Control_Tracking('notneeded',...
        Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,PT_marker_L,IndBodies_L,ZERO_F,ZERO_F,ZERO_F,ZERO_F,...
        Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,PT_marker_R,IndBodies_R,ZERO_F,ZERO_F,ZERO_F,ZERO_F);
end
toc

ModelsPath='C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Data\temp_debug\TestSwitchModels';
copyfile(fullfile(ModelsPath,'Models2','LeftSideModel.osim'),fullfile(MPI_path,'LeftSideModel.osim'));
copyfile(fullfile(ModelsPath,'Models2','RightSideModel.osim'),fullfile(MPI_path,'RightSideModel.osim'));
% SwitchModels();
pause(1);
tic
for i=1:100
    [T_cpp_L2,PT_posG_cpp_L,T_cpp_R2,PT_posG_cpp_R] =...
        Control_Tracking('notneeded',...
        Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,PT_marker_L,IndBodies_L,ZERO_F,ZERO_F,ZERO_F,ZERO_F,...
        Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,PT_marker_R,IndBodies_R,ZERO_F,ZERO_F,ZERO_F,ZERO_F);
end
toc
ModelsPath='C:\Users\u0088756\Documents\documenten\software\SoftwareProjects\DC_TrackingSimulation\Data\temp_debug\TestSwitchModels';
copyfile(fullfile(ModelsPath,'Models1','LeftSideModel.osim'),fullfile(MPI_path,'LeftSideModel.osim'));
copyfile(fullfile(ModelsPath,'Models1','RightSideModel.osim'),fullfile(MPI_path,'RightSideModel.osim'));
% SwitchModels();
pause(1);
tic
for i=1:100
    [T_cpp_L3,PT_posG_cpp_L,T_cpp_R3,PT_posG_cpp_R] =...
        Control_Tracking('notneeded',...
        Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,PT_marker_L,IndBodies_L,ZERO_F,ZERO_F,ZERO_F,ZERO_F,...
        Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,PT_marker_R,IndBodies_R,ZERO_F,ZERO_F,ZERO_F,ZERO_F);
end
toc
figure();
subplot(1,2,1)
bar([T_cpp_L(1,:); T_cpp_L2(1,:); T_cpp_L3(1,:)]');
subplot(1,2,2)
bar([T_cpp_R(1,:); T_cpp_R2(1,:); T_cpp_R3(1,:)]');

