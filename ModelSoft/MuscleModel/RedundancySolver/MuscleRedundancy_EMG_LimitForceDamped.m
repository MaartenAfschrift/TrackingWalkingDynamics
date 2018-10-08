function [] = MuscleRedundancy_EMG_LimitForceDamped(S,EMGs,MatlabData,OutName,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
ExtName=[];
if ~isempty(varargin)
    ExtName=varargin{1};
end

bool_CreateAdigatorFiles=true;
if length(varargin)>1
	bool_CreateAdigatorFiles=varargin{2};
end
if length(varargin)>2
	S.a_min=varargin{3};
end

% new added: EMGs is a structure with the EMG information
%   EMGs.RefL =>  Left leg muscle activity of refernece walking interpolated between heelstrikes
%   EMGs.RefR =>  Right leg muscle activity of refernece walking interpolated between heelstrikes


%% Path info
OutPath     = S.OutPath;
OptRes      = load(fullfile(OutPath,OutName));
OptRes.info.CARTELLAOUTPUT=OutPath;
MuscleInfoPathLeft  = fullfile(S.OutPath,'Left');
MuscleInfoPathRight = fullfile(S.OutPath,'Right');

S.KsName_Left   =fullfile(OutPath,'KS_Left_TrackinSim.mot');
S.KsName_Right  =fullfile(OutPath,'KS_Right_TrackinSim.mot');
S.IdName_Left   =fullfile(OutPath,'ID_Left_TrackingSim.sto');
S.IdName_Right  =fullfile(OutPath,'ID_Right_TrackingSim.sto');

%% Add paths to matlab
addpath(genpath(fullfile(OutPath)));
if S.MuscleFunctionsGeneric == 0
    addpath(genpath(MuscleInfoPathLeft));      % functions to compute moment arms
    addpath(genpath(MuscleInfoPathRight));     % functions to compute moment arms
end

import org.opensim.modeling.*

%% get the ID and KS solutions
LmodelPath=fullfile(OutPath,'LeftSideModel.osim');
RmodelPath=fullfile(OutPath,'RightSideModel.osim');

auxdata=OptRes.output.result.setup.auxdata;
solution=OptRes.output.result.solution;
state=solution.phase.state;
time=solution.phase.time;

Ncoord_L=auxdata.Ncoord_L;          Ncoord_R=auxdata.Ncoord_R;
IcOK_L=auxdata.IND_coord_OK_L;      IcNO_L=auxdata.IND_coord_NO_L;
IcOK_R=auxdata.IND_coord_OK_R;      IcNO_R=auxdata.IND_coord_NO_R;

Q_L=state(:,1:Ncoord_L);
Q_R=state(:,Ncoord_L+(1:Ncoord_R));
U_L=state(:,Ncoord_L+Ncoord_R+(1:Ncoord_L));
U_R=state(:,Ncoord_L+Ncoord_R+Ncoord_L+(1:Ncoord_R));

Q_L_Xcpp(:,IcOK_L)=Q_L;
Q_R_Xcpp(:,IcOK_R)=Q_R;
U_L_Xcpp(:,IcOK_L)=U_L;
U_R_Xcpp(:,IcOK_R)=U_R;
Q_L_Xcpp(:,IcNO_L)=0;
Q_R_Xcpp(:,IcNO_R)=0;
U_L_Xcpp(:,IcNO_L)=0;
U_R_Xcpp(:,IcNO_R)=0;

A_L=solution.phase.control(:,1:Ncoord_L);                     % left model qddot
A_R=solution.phase.control(:,Ncoord_L+(1:Ncoord_R));          % right model qddot
Fp      =solution.phase.control(:,Ncoord_L+Ncoord_R+(1:3));   % force pelvis
Mp      =solution.phase.control(:,Ncoord_L+Ncoord_R+3+(1:3)); % torques pelvis

A_L_Xcpp(:,IcOK_L)=A_L;
A_R_Xcpp(:,IcOK_R)=A_R;
A_L_Xcpp(:,IcNO_L)=0;
A_R_Xcpp(:,IcNO_R)=0;

% Get ID torques with cpp function used in the optimization
ZERO_F=zeros(length(time),3);
[T_cpp_L,PT_posG_cpp_L]=...
    trovaIDePATH_MD(LmodelPath,Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,auxdata.PT_marker_L,auxdata.PT_markerbody_L,-Fp,ZERO_F,-Mp,ZERO_F);
[T_cpp_R,PT_posG_cpp_R]=...
    trovaIDePATH_MD(RmodelPath,Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,auxdata.PT_marker_R,auxdata.PT_markerbody_R,Fp,ZERO_F,Mp,ZERO_F);

ID_headersL=cell(size(auxdata.NomiCoord_L));        ID_headersR=cell(size(auxdata.NomiCoord_R));
for i=1:length(auxdata.NomiCoord_L)
    ID_headersL{i}=[auxdata.NomiCoord_L{i} '_moment'];
end
for i=1:length(auxdata.NomiCoord_R)
    ID_headersR{i}=[auxdata.NomiCoord_R{i} '_moment'];
end

input.auxdata.ID_initial=T_cpp_L;

%% set heelstrike information
S.t_HS=solution.phase.time(1);         
S.t_HS2=solution.phase.time(end);
% remark ! important to track information betwen heelstrike, otherwise we run into problems here

%% Add EMG constraints
GML_ref=EMGs.RefL;
GMR_ref=EMGs.RefR;
GML={'glut_med1_l','glut_med2_l','glut_med3_l'};
GMR={'glut_med1_r','glut_med2_r','glut_med3_r'};

ML=Model(LmodelPath);    LMuscleNames=[];
for i=1:ML.getMuscles().getSize;
    LMuscleNames{i}= char(ML.getMuscles().get(i-1).getName());
end
MR=Model(RmodelPath);   RMuscleNames=[];
for i=1:MR.getMuscles().getSize;
    RMuscleNames{i}= char(MR.getMuscles().get(i-1).getName());
end

% interpolate to the correct time stamp
dt_stride=S.t_HS2-S.t_HS;
time_ref=linspace(S.t_HS-dt_stride,S.t_HS2+dt_stride,300);
GML_int=interp1(time_ref,repmat(GML_ref,1,3)',time);
pp_GML=spline(time,GML_int);
GMR_int=interp1(time_ref,repmat(GMR_ref,1,3)',time);
pp_GMR=spline(time,GMR_int);

IndsGML=zeros(1,3); IndsGMR=zeros(1,3);
for i=1:3
    IndsGML(i)=find(strcmp(LMuscleNames,GML{i}));
    IndsGMR(i)=find(strcmp(RMuscleNames,GMR{i}));
end

%% Export files
generateMotFile([time T_cpp_L],['time' ID_headersL],S.IdName_Left);
generateMotFile([time T_cpp_R],['time' ID_headersR],S.IdName_Right);
stampaIKfromGPOPS(S.KsName_Left,Q_L_Xcpp,auxdata.NomiCoord_L,time);
stampaIKfromGPOPS(S.KsName_Right,Q_R_Xcpp,auxdata.NomiCoord_R,time);

%% Get the muscle information
MuscleAnalysisPathL=fullfile(OutPath,'MuscleAnalysisModelLeft'); if ~isdir(MuscleAnalysisPathL); mkdir(MuscleAnalysisPathL); end;
if S.if_RunMuscleAnalaysis; OpenSim_Muscle_Analysis(S.KsName_Left,LmodelPath,MuscleAnalysisPathL,[time(1) time(end)]); end;
MuscleAnalysisPathR=fullfile(OutPath,'MuscleAnalysisModelRight'); if ~isdir(MuscleAnalysisPathR); mkdir(MuscleAnalysisPathR); end;
if S.if_RunMuscleAnalaysis; OpenSim_Muscle_Analysis(S.KsName_Right,RmodelPath,MuscleAnalysisPathR,[time(1) time(end)]); end;

%% Get Input information Left Model
if S.LeftModel
    import org.opensim.modeling.*
    ModelL=Model(LmodelPath);
    for i=1:ModelL.getMuscles().getSize();
        Misc.MuscleNames_Input{i}=char(ModelL.getMuscles().get(i-1).getName());
    end
    Misc.MuscleAnalysisPath=MuscleAnalysisPathL;
    if S.SubtalarDof
        Misc.DofNames_Input={'subtalar_angle_l','ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};
    else
        Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};
    end
    Misc.time=[time(1) time(end)];
    [DatStore] = getMuscleInfo(S.KsName_Left,S.IdName_Left,Misc);
    DatStore.ID_initial=T_cpp_L; DatStore.NumCol_initial=length(T_cpp_L(:,1));
    [DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=ReadMuscleParameters(LmodelPath,DatStore.MuscleNames);
    
    Settings.ScaleMuscleF=1.5;
    DatStore.Fiso=DatStore.Fiso.*Settings.ScaleMuscleF;
    
    %% Add coordinate limit force to simulation
    IndExp=2+1;
    k=[-12 15 11.03 -11.33];    theta=[-2.4 0.13];   c=0.025;
    q=Q_L_Xcpp(:,9); q_dot=U_L_Xcpp(:,9);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    TKneeL=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-TKneeL;
    
    %% Add possive force to the ankle joint
    k=[-4 15 4 -15];    theta=[-0.74 0.52]; c=0.025;
    Ind=8;  IndExp=1+1;
    q=Q_L_Xcpp(:,Ind); q_dot=U_L_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    TAnkleL=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-TAnkleL;
    
    %% Add passive force to the hip flexion
    Ind=10;     IndExp=3+1;
    k=[5 5.05 3 -10];   theta=[-0.47 1.81]; c=0.025;
    q=Q_L_Xcpp(:,Ind); q_dot=U_L_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    Thf=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-Thf;
    
    %% Add passive force to the hip adduction
    Ind=11; IndExp=4+1;
    k=[-5 15 5 -15];   theta=[-0.5235987755982988 0.2618]; c=0.025;
    q=Q_L_Xcpp(:,Ind); q_dot=U_L_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    Tadd=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-Tadd;
    
    %% Add passive force to the hip rotation
    Ind=12; IndExp=5+1;
    k=[-5 15 5 -15];   theta=[-0.5235987755982988 0.5235987755982988]; c=0.025;
    q=Q_L_Xcpp(:,Ind); q_dot=U_L_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    Trot=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-Trot;
    
    %% get LMT and mometn arm info form polynomials -- fast
    [nfr, ndof]=size(Q_L_Xcpp);
    load(fullfile(MuscleInfoPathLeft,'MuscleInfoL.mat'));
    dM_Store=zeros(nfr,ndof,length(MuscleInfoL.muscle));  	LMT_Store=zeros(nfr,length(MuscleInfoL.muscle));
    VMT_Store=zeros(nfr,length(MuscleInfoL.muscle));
    Q7=Q_L_Xcpp(:,7);
    Q8=Q_L_Xcpp(:,8);       Q9=Q_L_Xcpp(:,9);       Q10=Q_L_Xcpp(:,10);
    Q11=Q_L_Xcpp(:,11);     Q12=Q_L_Xcpp(:,12);     Q13=Q_L_Xcpp(:,13);
    Q14=Q_L_Xcpp(:,14);     Q15=Q_L_Xcpp(:,15);     x1=ones(nfr,1);
    x0=zeros(nfr,1);
    if S.SubtalarDof
        MatdiffAll=Get_Mat_Q_Left(Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,x0,x1);
        MatAll=Get_Mat_Left(Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,x1);
    else
        MatdiffAll=Get_Mat_Q_Left(Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,x0,x1);
        MatAll=Get_Mat_Left(Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,x1);
    end
    for m=1:length(MuscleInfoL.muscle)
        
        % get the matrices with kinematics info
        index_dof_crossing=MuscleInfoL.muscle(m).index_dof_crossing;
        nr_dof_crossing=MuscleInfoL.muscle(m).nr_dof_crossing;
        [mat,diff_mat_q] = getMatrices_MomentArmComp(MuscleInfoL,nfr,m,MatAll,MatdiffAll);
        nr_coeffs=length(MuscleInfoL.muscle(m).coeff);
        
        % evaluate the polynomials
        coeff=MuscleInfoL.muscle(m).coeff;
        dM_recon = zeros(nfr, nr_dof_crossing);
        for dof_nr = 1:nr_dof_crossing
            dM_recon(:,dof_nr) = (-squeeze(diff_mat_q(:,:,dof_nr)))*coeff;
        end
        LMT_recon=mat*coeff;
        dM_Store(:,index_dof_crossing,m)=dM_recon;
        LMT_Store(:,m)=LMT_recon;
        dM_sel=squeeze(dM_Store(:,index_dof_crossing,m));
        VMT_Store(:,m)=sum(-dM_sel.*U_L_Xcpp(:,index_dof_crossing),2);     % d(LMT)/dt = d(LMT)/dqi * dqi/dt + d(LMT)/dqj * dqj/dt
    end
    DatStore.dM=squeeze(dM_Store(:,8:12,:));
    DatStore.LMT=LMT_Store;
    DatStore.VMT=VMT_Store;
    DatStore.time=time;
    
    %% Solve static optimization
    DatStore = SolveStaticOptimization_IPOPT(DatStore);
    
    %% Add EMG constraint information
    DatStore.pp_GM=pp_GML;
    DatStore.IndsMuscle=IndsGML;
    
    %% Add lower bound EMG informatin
    if isfield(S,'a_min') && ~isempty(S.a_min);
        Misc.ActLower = S.a_min;
    end
    %% Solve the optimal control problem
    Misc.Atendon=ones(1,length(Misc.MuscleNames_Input)).*35;
    Misc.Mesh_Frequency=100;
    [output] = SolveMRed_Function_Va_EMG_Damped(DatStore,Misc,OptRes,bool_CreateAdigatorFiles);
    save(fullfile(OutPath,[OutName  ExtName 'LeftModel_vA_EMG_damped.mat']),'output','DatStore','Misc','OptRes');
    clear DatStore Misc
end

%% Get Input information Right Model
if S.RightModel
    import org.opensim.modeling.*
    ModelR=Model(RmodelPath);
    for i=1:ModelR.getMuscles().getSize();
        Misc.MuscleNames_Input{i}=char(ModelR.getMuscles().get(i-1).getName());
    end
    Misc.MuscleAnalysisPath=MuscleAnalysisPathR;
    if S.SubtalarDof
        Misc.DofNames_Input={'subtalar_angle_r','ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};
    else
        Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};
    end
    Misc.time=[time(1) time(end)];
    [DatStore] = getMuscleInfo(S.KsName_Right,S.IdName_Right,Misc);
    DatStore.ID_initial=T_cpp_R; DatStore.NumCol_initial=length(T_cpp_R(:,1));
    [DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=ReadMuscleParameters(RmodelPath,DatStore.MuscleNames);
    
    Settings.ScaleMuscleF=1.5;
    DatStore.Fiso=DatStore.Fiso.*Settings.ScaleMuscleF;
    
    %% Add coordinate limit force to knee
    k=[-12 15 11.03 -11.33];    theta=[-2.4 0.13];   c=0.025;
    IndExp=2+1;
    q=Q_R_Xcpp(:,9); q_dot=U_R_Xcpp(:,9);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    TKneeR=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-TKneeR; 
    
    %% Add possive force to the ankle joint
    k=[-4 15 4 -15];    theta=[-0.74 0.52]; c=0.025;
    Ind=8;  IndExp=1+1;
    q=Q_R_Xcpp(:,Ind); q_dot=U_R_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    TAnkleR=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-TAnkleR;
    
    %% Add passive force to the hip flexion
    Ind=10; IndExp=3+1;
    k=[5 5.05 3 -10];   theta=[-0.47 1.81]; c=0.025;
    q=Q_R_Xcpp(:,Ind); q_dot=U_R_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    Thf=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-Thf;
    
    %% Add passive force to the hip adduction
    Ind=11; IndExp=4+1;
    k=[-5 15 5 -15];   theta=[-0.5235987755982988 0.2618]; c=0.025;
    q=Q_R_Xcpp(:,Ind); q_dot=U_R_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    Tadd=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-Tadd;
    
    %% Add passive force to the hip rotation
    Ind=12; IndExp=5+1;
    k=[-5 15 5 -15];   theta=[-0.5235987755982988 0.5235987755982988]; c=0.025;    
    q=Q_R_Xcpp(:,Ind); q_dot=U_R_Xcpp(:,Ind);
    TP=k(1)*exp(k(2)*(q-theta(2)))+k(3)*exp(k(4)*(q-theta(1)));
    TD=-c*q_dot;    Trot=TP+TD;
    DatStore.T_exp(:,IndExp)=DatStore.T_exp(:,IndExp)-Trot;    
    
    %% get LMT and mometn arm info form polynomials -- fast
    [nfr, ndof]=size(Q_R_Xcpp);
    load(fullfile(MuscleInfoPathRight,'MuscleInfoR.mat'));
    dM_Store=zeros(nfr,ndof,length(MuscleInfoR.muscle));    LMT_Store=zeros(nfr,length(MuscleInfoR.muscle));
    VMT_Store=zeros(nfr,length(MuscleInfoR.muscle));
    Q7=Q_R_Xcpp(:,7);
    Q8=Q_R_Xcpp(:,8);       Q9=Q_R_Xcpp(:,9);       Q10=Q_R_Xcpp(:,10);
    Q11=Q_R_Xcpp(:,11);     Q12=Q_R_Xcpp(:,12);     x1=ones(nfr,1);
    x0=zeros(nfr,1);
    if S.SubtalarDof
        MatdiffAll=Get_Mat_Q_Right(Q7,Q8,Q9,Q10,Q11,Q12,x0,x1);
        MatAll=Get_Mat_Right(Q7,Q8,Q9,Q10,Q11,Q12,x1);
    else    
        MatdiffAll=Get_Mat_Q_Right(Q8,Q9,Q10,Q11,Q12,x0,x1);
        MatAll=Get_Mat_Right(Q8,Q9,Q10,Q11,Q12,x1);
    end
    for m=1:length(MuscleInfoR.muscle)
        
        % get the matrices with kinematics info
        index_dof_crossing=MuscleInfoR.muscle(m).index_dof_crossing;
        nr_dof_crossing=MuscleInfoR.muscle(m).nr_dof_crossing;
        [mat,diff_mat_q] = getMatrices_MomentArmComp(MuscleInfoR,nfr,m,MatAll,MatdiffAll);
        nr_coeffs=length(MuscleInfoR.muscle(m).coeff);
        
        % evaluate the polynomials
        coeff=MuscleInfoR.muscle(m).coeff;
        dM_recon = zeros(nfr, nr_dof_crossing);
        for dof_nr = 1:nr_dof_crossing
            dM_recon(:,dof_nr) = (-squeeze(diff_mat_q(:,:,dof_nr)))*coeff;
        end
        LMT_recon=mat*coeff;
        dM_Store(:,index_dof_crossing,m)=dM_recon;
        LMT_Store(:,m)=LMT_recon;
        dM_sel=squeeze(dM_Store(:,index_dof_crossing,m));
        VMT_Store(:,m)=sum(-dM_sel.*U_R_Xcpp(:,index_dof_crossing),2);     % d(LMT)/dt = d(LMT)/dqi * dqi/dt + d(LMT)/dqj * dqj/dt
        
    end
    if S.SubtalarDof
        DatStore.dM=squeeze(dM_Store(:,7:12,:));
    else
        DatStore.dM=squeeze(dM_Store(:,8:12,:));
    end
    DatStore.LMT=LMT_Store;
    DatStore.VMT=VMT_Store;
    DatStore.time=time;
    
    %% Solve static optimization
    DatStore = SolveStaticOptimization_IPOPT(DatStore);
    
    %% Add EMG constraint information
    DatStore.pp_GM=pp_GMR;
    DatStore.IndsMuscle=IndsGMR;
    
    %% Solve the optimal control problem
    Misc.Atendon=ones(1,length(Misc.MuscleNames_Input)).*35;
    Misc.Mesh_Frequency=100;
    addpath(pwd);
    [output] = SolveMRed_Function_Va_EMG_Damped(DatStore,Misc,OptRes,bool_CreateAdigatorFiles);
    %     save(fullfile(OutPath,'Results_MuscleRedundancyRightModel_vA_EMG.mat'),'output','DatStore','Misc','OptRes');
    save(fullfile(OutPath,[OutName ExtName 'RightModel_vA_EMG_damped.mat']),'output','DatStore','Misc','OptRes');
    clear DatStore Misc OptRes
    
    
end
if S.MuscleFunctionsGeneric == 0
    rmpath(genpath(MuscleInfoPathLeft));      % functions to compute moment arms
    rmpath(genpath(MuscleInfoPathRight));     % functions to compute moment arms
end
end

