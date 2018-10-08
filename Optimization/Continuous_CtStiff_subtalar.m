function OUT=Continuous_CtStiff_subtalar(IN)
global BoolFirst
persistent NCollPoints jtToTrack  TRC GRF

% copia roba varia
Ncoord_L=IN.auxdata.Ncoord_L;
Ncoord_R=IN.auxdata.Ncoord_R;
time=IN.phase.time;
N=length(time);
ZERO_F=zeros(N,3);
IcOK_L=IN.auxdata.IND_coord_OK_L;
IcNO_L=IN.auxdata.IND_coord_NO_L;
IcOK_R=IN.auxdata.IND_coord_OK_R;
IcNO_R=IN.auxdata.IND_coord_NO_R;
SCALA=IN.auxdata.SCALA;


% evaluate splines
if BoolFirst || N~=NCollPoints
   jtToTrack=ppval(IN.auxdata.ppID,time)';
   TRC=[ppval(IN.auxdata.ppTRC_L,time)' ppval(IN.auxdata.ppTRC_R,time)'];
   GRF=[ppval(IN.auxdata.ppGRF_F_L,time)' ppval(IN.auxdata.ppGRF_F_R,time)'];
   NCollPoints=N;
   BoolFirst=false;
end

% Copy the states
Q_L=IN.phase.state(:,1:Ncoord_L);
Q_R=IN.phase.state(:,Ncoord_L+(1:Ncoord_R));
U_L=IN.phase.state(:,Ncoord_L+Ncoord_R+(1:Ncoord_L));
U_R=IN.phase.state(:,Ncoord_L+Ncoord_R+Ncoord_L+(1:Ncoord_R));
if SCALA
    Q_L=Q_L*IN.auxdata.scaling.Q(1:Ncoord_L,1:Ncoord_L);
    Q_R=Q_R*IN.auxdata.scaling.Q(Ncoord_L+1:end,Ncoord_L+1:end);
    U_L=U_L*IN.auxdata.scaling.U(1:Ncoord_L,1:Ncoord_L);
    U_R=U_R*IN.auxdata.scaling.U(Ncoord_L+1:end,Ncoord_L+1:end);
end

% Create input files for cpp function
Q_L_Xcpp(:,IcOK_L)=Q_L;
Q_R_Xcpp(:,IcOK_R)=Q_R;
U_L_Xcpp(:,IcOK_L)=U_L;
U_R_Xcpp(:,IcOK_R)=U_R;
Q_L_Xcpp(:,IcNO_L)=0;
Q_R_Xcpp(:,IcNO_R)=0;
U_L_Xcpp(:,IcNO_L)=0;
U_R_Xcpp(:,IcNO_R)=0;

% copy the controls and create input for cpp
A_L=IN.phase.control(:,1:Ncoord_L);                     % left model qddot
A_R=IN.phase.control(:,Ncoord_L+(1:Ncoord_R));          % right model qddot
Fp      =IN.phase.control(:,Ncoord_L+Ncoord_R+(1:3));   % force pelvis
Mp      =IN.phase.control(:,Ncoord_L+Ncoord_R+3+(1:3)); % torques pelvis
Tc_u  =IN.phase.control(:,Ncoord_L+Ncoord_R+7:Ncoord_L+Ncoord_R+6+IN.auxdata.Ndof_damping); % torques pelvis
Tc    = Tc_u.*IN.auxdata.Tmax;

if SCALA
    A_L=A_L*IN.auxdata.scaling.A(1:Ncoord_L,1:Ncoord_L);
    A_R=A_R*IN.auxdata.scaling.A(Ncoord_L+1:end,Ncoord_L+1:end);
    Fp=Fp*IN.auxdata.scaling.Fp;
    Mp=Mp*IN.auxdata.scaling.Mp;
end
A_L_Xcpp(:,IcOK_L)=A_L;
A_R_Xcpp(:,IcOK_R)=A_R;
A_L_Xcpp(:,IcNO_L)=0;
A_R_Xcpp(:,IcNO_R)=0;

% get contact model properties
Stiffness   =   IN.auxdata.ContactProp.K;
u           =   IN.auxdata.ContactProp.u;
c           =   IN.auxdata.ContactProp.c;
vt          =   IN.auxdata.ContactProp.vt;

%% Compute ground reaction force in matlab
PV_PT_L=contactPointsQeU_GPOPS_Vect_Contact(Q_L(:,1:3),Q_L(:,4:6),U_L(:,1:3),U_L(:,4:6),IN.auxdata.ContactInfo);
[F_L,M_L]=POSVELaar2GRF_GPOPS_Contact2(PV_PT_L,IN.auxdata.ContactInfo,Stiffness,u,c,vt);
M_L=-M_L;
PV_PT_R=contactPointsQeU_GPOPS_Vect_Contact(Q_R(:,1:3),Q_R(:,4:6),U_R(:,1:3),U_R(:,4:6),IN.auxdata.ContactInfo);
[F_R,M_R]=POSVELaar2GRF_GPOPS_Contact2(PV_PT_R,IN.auxdata.ContactInfo,Stiffness,u,c,vt);
M_R=-M_R;


%% Solve ID problem in OpenSim
[T_cpp_L,PT_posG_cpp_L,T_cpp_R,PT_posG_cpp_R] =...
    Control_Tracking('notneeded',...
    Q_L_Xcpp,U_L_Xcpp,A_L_Xcpp,IN.auxdata.PT_marker_L,IN.auxdata.IndsBodies,-Fp,ZERO_F,-Mp,ZERO_F,...
    Q_R_Xcpp,U_R_Xcpp,A_R_Xcpp,IN.auxdata.PT_marker_R,IN.auxdata.IndsBodies_R,Fp,ZERO_F,Mp,ZERO_F);

[T_ccp_L_inG,T_ccp_R_inG ] = ComputeGRFConstraint_vect(T_cpp_L,T_cpp_R,Q_L,Q_R);


nCpt=3;
nMark_L=length(IN.auxdata.PT_markerbody_L)-nCpt;
nMark_R=length(IN.auxdata.PT_markerbody_R)-nCpt;
PT_pelv_L=reshape(PT_posG_cpp_L(:,1:nCpt),size(PT_posG_cpp_L,1)/3,nCpt*3);
PT_pelv_R=reshape(PT_posG_cpp_R(:,1:nCpt),size(PT_posG_cpp_R,1)/3,nCpt*3);
TRCcpp_L=reshape(PT_posG_cpp_L(:,nCpt+1:end),size(PT_posG_cpp_L,1)/3,nMark_L*3);
TRCcpp_R=reshape(PT_posG_cpp_R(:,nCpt+1:end),size(PT_posG_cpp_R,1)/3,nMark_R*3);

%% damping constraint

Damping_Torques=[U_L_Xcpp(:,7:15) U_R_Xcpp(:,7:12)].*IN.auxdata.damping;
T_ID_damped=Tc-Damping_Torques;
ID_dat=[T_cpp_L(:,7:15) T_cpp_R(:,7:12)];
Damping_Constraint=T_ID_damped-ID_dat;


%% Evaluate objective function (tracking)

% deltaTrackGRF=[F_L-ppval(IN.auxdata.ppGRF_F_L,time)' F_R-ppval(IN.auxdata.ppGRF_F_R,time)'];
deltaF_L=T_cpp_L(:,4:6)-F_L;
deltaF_R=T_cpp_R(:,4:6)-F_R;
deltaM_L=T_ccp_L_inG-M_L;
deltaM_R=T_ccp_R_inG-M_R;
deltaPT=PT_pelv_L-PT_pelv_R;

% track joint torques
% jtToTrack=ppval(IN.auxdata.ppID,time)';
deltaTrackID=[jtToTrack(:,IN.auxdata.ppIDindL(2,:)) jtToTrack(:,IN.auxdata.ppIDindR(2,:))]- ...
    [T_cpp_L(:,IN.auxdata.ppIDindL(1,:))   T_cpp_R(:,IN.auxdata.ppIDindR(1,:))  ];

% track marker position
% TRC=[ppval(IN.auxdata.ppTRC_L,time)' ppval(IN.auxdata.ppTRC_R,time)'];
deltaTrackTRC=[(TRC(:,1:nMark_L*3)-TRCcpp_L)*IN.auxdata.MarkerWeightsL (TRC(:,nMark_L*3+1:end)-TRCcpp_R)*IN.auxdata.MarkerWeightsR];
deltaTrackTRC(:,2:3:end)=deltaTrackTRC(:,2:3:end)*0.5;     % less weight on tracking vertical position

% track GRF
% GRF=[ppval(IN.auxdata.ppGRF_F_L,time)' ppval(IN.auxdata.ppGRF_F_R,time)'];
deltaTracGRF=GRF-[T_cpp_L(:,4:6) T_cpp_R(:,4:6)];

%%
if SCALA
    OUT.dynamics=[[U_L U_R]/IN.auxdata.scaling.U [A_L A_R]/IN.auxdata.scaling.A];
    OUT.path=[deltaF_L/IN.auxdata.scaling.F_L deltaM_L/IN.auxdata.scaling.M_L deltaF_R/IN.auxdata.scaling.F_R...
        deltaM_R/IN.auxdata.scaling.M_R deltaPT/blkdiag(IN.auxdata.scaling.P000,IN.auxdata.scaling.P010,IN.auxdata.scaling.P001)];
    NNN=10000000;
else
    NNN=100000;
    OUT.dynamics=[U_L U_R A_L A_R];
    OUT.path=[deltaF_L deltaM_L deltaF_R deltaM_R deltaPT Damping_Constraint];
end


q1=sum(deltaTrackTRC.^2,2)*10/NNN;
q2=sum(deltaTracGRF.^2,2)/200/NNN;
q3=sum(deltaTrackID.^2,2)/10/NNN;

q1=q1*300;     % increase weight on tracking markers
q=q1/10+q2+q3;
q_acc=sum([A_L A_R].^2,2) .* 0.00000001;
q_grfAdd=(U_R(:,4).^4.*F_R(:,1).^2).*10^-6;

OUT.integrand=q+ q_acc + q_grfAdd;