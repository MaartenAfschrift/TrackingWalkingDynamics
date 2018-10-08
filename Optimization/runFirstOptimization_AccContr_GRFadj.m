function [varargout]=runFirstOptimization_AccContr_GRFadj(varargin)
% reorganize input parameters
if rem(nargin,2)
    error('<<<<< Wrong number of input')
end
InputName=varargin(1:2:end);
InputValue=varargin(2:2:end);
% mesh interval
tmpI=strcmpi(InputName,'NMeshIntervals');
if ~any(tmpI)
    NMeshIntervals=10;
else
    NMeshIntervals=InputValue{tmpI};
end
% scaling bound
tmpI=strcmpi(InputName,'ScalingBound');
if ~any(tmpI)
    SCALINGBOUND=1.5;
else
    SCALINGBOUND=InputValue{tmpI};
end
% collocation points
tmpI=strcmpi(InputName,'NColPoints');
if ~any(tmpI)
    NCP=4;
else
    NCP=InputValue{tmpI};
end
% model names
tmpI=strcmpi(InputName,'LeftModelName');
if ~any(tmpI)
    error('<<<<< No Left Side Model Specified')
else
    LEFTMODELNAME=InputValue{tmpI};
end
tmpI=strcmpi(InputName,'RightModelName');
if ~any(tmpI)
    error('<<<<< No Right Side Model Specified')
else
    RIGHTMODELNAME=InputValue{tmpI};
end
% GRF data
tmpI=strcmpi(InputName,'GRFdata');
if ~any(tmpI)
    error('<<<<< No GRF Data Specified')
else
    TxGRF=InputValue{tmpI}.TxGRF;
    GRFinput.head=InputValue{tmpI}.head;
    GRFdata=InputValue{tmpI}.GRFdata;
end
% time
tmpI=strcmpi(InputName,'Time');
if ~any(tmpI)
    error('<<<<< No Time Vector Specified')
else
    T=InputValue{tmpI};
end
% marker info
tmpI=strcmpi(InputName,'MarkerInfo');
if ~any(tmpI)
    error('<<<<< No Marker Info Left Side Specified')
else
    Markerinfo_L=InputValue{tmpI}.Markerinfo_L;
    Markerinfo_R=InputValue{tmpI}.Markerinfo_R;
    MarkerTime_L=InputValue{tmpI}.MarkerTime_L;
    MarkerTime_R=InputValue{tmpI}.MarkerTime_R;
end
% output from tuning GRF
tmpI=strcmpi(InputName,'OUTtuneGRF');
if ~any(tmpI)
    error('<<<<< No GRF Tuning Data Specified')
else
    OUT_tuneGRF=InputValue{tmpI};
end
% handles to continuous and endpoint functions
tmpI=strcmpi(InputName,'HandleContinuous');
if ~any(tmpI)
    error('<<<<< No Continuous Function Specified')
else
    HandleCONTINUOUS=InputValue{tmpI};
end
tmpI=strcmpi(InputName,'HandleEndpoint');
if ~any(tmpI)
    error('<<<<< No Endpoint Function Specified')
else
    HandleENDPOINT=InputValue{tmpI};
end
% output initial guess ?
tmpI=strcmpi(InputName,'OUTinitGuessL');
if ~any(tmpI)
    error('<<<<< No Left Side Init Guess Data Specified')
else
    OUT_L=InputValue{tmpI};
end
tmpI=strcmpi(InputName,'OUTinitGuessR');
if ~any(tmpI)
    error('<<<<< No Right Side Init Guess Data Specified')
else
    OUT_R=InputValue{tmpI};
end
% index for guess
tmpI=strcmpi(InputName,'LineIndexForGuess');
if ~any(tmpI)
    LineIndexForGuess=1:length(T);
else
    LineIndexForGuess=InputValue{tmpI};
end
% data about states and control for initial guess
tmpI=strcmpi(InputName,'StatesControls');
if ~any(tmpI)
    error('<<<<< No Info About States And Controls Specified')
else
    Q_L=InputValue{tmpI}.Q_L;
    Q_R=InputValue{tmpI}.Q_R;
    U_L=InputValue{tmpI}.U_L;
    U_R=InputValue{tmpI}.U_R;
    A_L=InputValue{tmpI}.A_L;
    A_R=InputValue{tmpI}.A_R;
end
% forces between lest and right
tmpI=strcmpi(InputName,'ForcePelvis');
if ~any(tmpI)
    error('<<<<< No Pelvis Force Data Specified')
else
    Fp=InputValue{tmpI};
end
% torques between lest and right
tmpI=strcmpi(InputName,'TorquePelvis');
if ~any(tmpI)
    error('<<<<< No Pelvis Torque Data Specified')
else
    Mp=InputValue{tmpI};
end
% auxdata
tmpI=strcmpi(InputName,'Auxdata');
if any(tmpI)
    auxdata=InputValue{tmpI};
end
% stampa
tmpI=strcmpi(InputName,'stampa');
if ~any(tmpI)
    STAMPA=true;
else
    STAMPA=InputValue{tmpI};
end
% scale model variables manually
tmpI=strcmpi(InputName,'scala');
if ~any(tmpI)
    SCALA=false;
else
    SCALA=InputValue{tmpI};
end


%% definisci problem bounds
% stati         - [Q_l Q_r U_l U_r A_l A_r]
% controlli     - [A_l A_r Fp Mp]

% defines the coordinates range of motion to be used as boundary
[Range_L,NomiCoord_L]=definisciJointRange(LEFTMODELNAME);
Range_L(11,:)=[0 5/180*pi];
[Range_R,NomiCoord_R]=definisciJointRange(RIGHTMODELNAME);
Range_R(11,:)=[0 5/180*pi];

% take the mtp angle out of the coordinates to be controlled
IND_coord_NO_L=strcmpi(NomiCoord_L,'mtp_angle_l');
IND_coord_OK_L=~IND_coord_NO_L;
IND_coord_OK_L=find(IND_coord_OK_L);
IND_coord_NO_L=find(IND_coord_NO_L);
IND_coord_NO_R=strcmpi(NomiCoord_R,'mtp_angle_r');
IND_coord_OK_R=~IND_coord_NO_R;
IND_coord_OK_R=find(IND_coord_OK_R);
IND_coord_NO_R=find(IND_coord_NO_R);
auxdata.IND_coord_OK_L=IND_coord_OK_L;
auxdata.IND_coord_NO_L=IND_coord_NO_L;
auxdata.IND_coord_OK_R=IND_coord_OK_R;
auxdata.IND_coord_NO_R=IND_coord_NO_R;
Ncoord_L=length(IND_coord_OK_L);
Ncoord_R=length(IND_coord_OK_R);
Q_L=Q_L(:,IND_coord_OK_L);
Q_R=Q_R(:,IND_coord_OK_R);
U_L=U_L(:,IND_coord_OK_L);
U_R=U_R(:,IND_coord_OK_R);
A_L=A_L(:,IND_coord_OK_L);
A_R=A_R(:,IND_coord_OK_R);
auxdata.Ncoord_L=Ncoord_L;
auxdata.Ncoord_R=Ncoord_R;

del_store=[];
for i=1:length(IND_coord_NO_L)
    del_store=[del_store find(auxdata.ppIDindL(1,:)==IND_coord_NO_L(i))];
end
auxdata.ppIDindL(:,del_store)=[];
del_store=[];
for i=1:length(IND_coord_NO_R)
    del_store=[del_store find(auxdata.ppIDindR(1,:)==IND_coord_NO_R(i))];
end
auxdata.ppIDindR(:,del_store)=[];
% data to track
auxdata.ppP000=spline(TxGRF,OUT_L.Pelv000');
auxdata.ppQ_R=spline(T(LineIndexForGuess),Q_R');
auxdata.ppQ_L=spline(T(LineIndexForGuess),Q_L');
auxdata.ppGRF_F_R=spline(GRFdata(:,1),GRFdata(:,8:10)');%spline(GRFdata(:,1),GRFdata(:,9)');%auxdata.ppGRF_F_R=spline(TxGRF,OUT_tuneGRF.GRFopt(:,7:9)');
auxdata.ppGRF_F_L=spline(GRFdata(:,1),GRFdata(:,2:4)');%spline(GRFdata(:,1),GRFdata(:,3)');%auxdata.ppGRF_F_L=spline(TxGRF,OUT_tuneGRF.GRFopt(:,1:3)');
auxdata.ppTinG_L=spline(TxGRF,OUT_L.T_ccp_inG_FINAL');
auxdata.ppTinG_R=spline(TxGRF,OUT_R.T_ccp_inG_FINAL');
auxdata.scaling.Q=diag(max(abs([Q_L Q_R]),[],1));
auxdata.scaling.U=diag(max(abs([U_L U_R]),[],1));
auxdata.scaling.A=diag(max(abs([A_L A_R]),[],1));
auxdata.scaling.Fp=[800 0 0 ; 0 800 0 ; 0 0 300];%diag(max(abs(Fp),[],1));
auxdata.scaling.Mp=[200 0 0 ; 0 400 0 ; 0 0 200];%diag(max(abs(Mp),[],1));

% stati
PosMIN=[Range_L(IND_coord_OK_L,1)' Range_R(IND_coord_OK_R,1)'];
PosMAX=[Range_L(IND_coord_OK_L,2)' Range_R(IND_coord_OK_R,2)'];

VelMIN=-ones(1,Ncoord_L+Ncoord_R)*SCALINGBOUND*auxdata.scaling.U;
VelMAX=ones(1,Ncoord_L+Ncoord_R)*SCALINGBOUND*auxdata.scaling.U;
if SCALA
    PosMIN=PosMIN/auxdata.scaling.Q;
    PosMAX=PosMAX/auxdata.scaling.Q;
    VelMIN=VelMIN/auxdata.scaling.U;
    VelMAX=VelMAX/auxdata.scaling.U;
end
QboundL=[PosMIN VelMIN];
QboundU=[PosMAX VelMAX];

% controlli
Fp_MIN=[-1200 -1200 -700];
Fp_MAX=[1200 1200 700];
Mp_MIN=[-500 -500 -500];
Mp_MAX=[500 500 500];
AccMIN=-ones(1,Ncoord_L+Ncoord_R)*SCALINGBOUND*auxdata.scaling.A;
AccMAX=ones(1,Ncoord_L+Ncoord_R)*SCALINGBOUND*auxdata.scaling.A;
if SCALA
    Fp_MIN=Fp_MIN/auxdata.scaling.Fp;
    Fp_MAX=Fp_MAX/auxdata.scaling.Fp;
    Mp_MIN=Mp_MIN/auxdata.scaling.Mp;
    Mp_MAX=Mp_MAX/auxdata.scaling.Mp;
    AccMIN=AccMIN/auxdata.scaling.A;
    AccMAX=AccMAX/auxdata.scaling.A;
end

% Bounds
T0=T(1);
TF=T(end);
bounds.phase.initialtime.lower=T0;
bounds.phase.initialtime.upper=T0;
bounds.phase.finaltime.lower=TF;
bounds.phase.finaltime.upper=TF;
deltaT=(TF-T0)/NMeshIntervals;
bounds.phase.control.lower=[AccMIN Fp_MIN Mp_MIN ];
bounds.phase.control.upper=[AccMAX Fp_MAX Mp_MAX ];
bounds.phase.initialstate.lower=QboundL;
bounds.phase.initialstate.upper=QboundU;


if SCALA
    % scale joint kinematics at first frame
    auxdata.scaling.QL=diag(max(abs(Q_L),[],1));
    auxdata.scaling.QR=diag(max(abs(Q_R),[],1));
    Q_L_scaled=Q_L/auxdata.scaling.QL;
    Q_R_scaled=Q_R/auxdata.scaling.QR;
    dq_firstFrame_L=[0.03 0.1 0.03]/auxdata.scaling.QL(4:6,4:6);
    dq_firstFrame_R=[0.03 0.1 0.03]/auxdata.scaling.QR(4:6,4:6);
    bounds.phase.initialstate.lower([4:6 (4:6)+Ncoord_L])=[Q_L_scaled(1,4:6) Q_R_scaled(1,4:6)]-[dq_firstFrame_L    dq_firstFrame_R];
    bounds.phase.initialstate.upper([4:6 (4:6)+Ncoord_L])=[Q_L_scaled(1,4:6) Q_R_scaled(1,4:6)]+[dq_firstFrame_L    dq_firstFrame_R];
else
    bounds.phase.initialstate.lower([4:6 (4:6)+Ncoord_L])=[Q_L(1,4:6) Q_R(1,4:6)]-[0.03 0.1 0.03   0.03 0.1 0.03];
    bounds.phase.initialstate.upper([4:6 (4:6)+Ncoord_L])=[Q_L(1,4:6) Q_R(1,4:6)]+[0.03 0.1 0.03   0.03 0.1 0.03];
end
bounds.phase.state.lower=QboundL;
bounds.phase.state.upper=QboundU;
if SCALA
    bounds.phase.state.lower([4:6 (4:6)+Ncoord_L])=min([Q_L_scaled(:,4:6) Q_R_scaled(:,4:6)],[],1)-[dq_firstFrame_L    dq_firstFrame_R];
    bounds.phase.state.upper([4:6 (4:6)+Ncoord_L])=max([Q_L_scaled(:,4:6) Q_R_scaled(:,4:6)],[],1)+[dq_firstFrame_L    dq_firstFrame_R];
else
    bounds.phase.state.lower([4:6 (4:6)+Ncoord_L])=min([Q_L(:,4:6) Q_R(:,4:6)],[],1)-[0.03 0.15 0.03    0.03 0.15 0.03];
    bounds.phase.state.upper([4:6 (4:6)+Ncoord_L])=max([Q_L(:,4:6) Q_R(:,4:6)],[],1)+[0.03 0.03 0.03    0.03 0.03 0.03];
end
bounds.phase.finalstate.lower=QboundL;
bounds.phase.finalstate.upper=QboundU;
if SCALA
    bounds.phase.finalstate.lower([4:6 (4:6)+Ncoord_L])=[Q_L_scaled(end,4:6) Q_R_scaled(end,4:6)]-[dq_firstFrame_L    dq_firstFrame_R];
    bounds.phase.finalstate.upper([4:6 (4:6)+Ncoord_L])=[Q_L_scaled(end,4:6) Q_R_scaled(end,4:6)]+[dq_firstFrame_L    dq_firstFrame_R];
else
    bounds.phase.finalstate.lower([4:6 (4:6)+Ncoord_L])=[Q_L(end,4:6) Q_R(end,4:6)]-[0.03 0.15 0.03    0.03 0.15 0.03];
    bounds.phase.finalstate.upper([4:6 (4:6)+Ncoord_L])=[Q_L(end,4:6) Q_R(end,4:6)]+[0.03 0.03 0.03    0.03 0.03 0.03];
end

bounds.phase.initialstate.lower=bounds.phase.initialstate.lower;
bounds.phase.initialstate.upper=bounds.phase.initialstate.upper;
bounds.phase.finalstate.lower=bounds.phase.finalstate.lower;
bounds.phase.finalstate.upper=bounds.phase.finalstate.upper;
bounds.phase.state.lower=bounds.phase.state.lower;
bounds.phase.state.upper=bounds.phase.state.upper;

% defines weights to be assigned to the tracking of markers
% markers on the feet have higher weight
PPWW=2.5;
tmpMW=ones(1,length(Markerinfo_L));
indMW=strncmpi(Markerinfo_L(2,:),'calcn_',length('calcn_'));
tmpMW(indMW)=PPWW;
tmp2MW=[];
for q=1:length(tmpMW)
    tmp2MW=[tmp2MW ones(1,3)*tmpMW(q)];
end
auxdata.MarkerWeightsL=diag(tmp2MW);
tmpMW=ones(1,length(Markerinfo_R));
indMW=strncmpi(Markerinfo_R(2,:),'calcn_',length('calcn_'));
tmpMW(indMW)=PPWW;
tmp2MW=[];
for q=1:length(tmpMW)
    tmp2MW=[tmp2MW ones(1,3)*tmpMW(q)];
end
auxdata.MarkerWeightsR=diag(tmp2MW);

% Integral bounds
bounds.phase.integral.lower = 0;
auxdata.scaling.F_L=diag([1 1 1])*1;
auxdata.scaling.F_R=diag([1 1 1])*1;
auxdata.scaling.M_L=diag([1 1 1])*0.1;
auxdata.scaling.M_R=diag([1 1 1])*0.1;
auxdata.scaling.P000=diag([1 1 1])*1e-6;
auxdata.scaling.P010=diag([1 1 1])*1e-6;
auxdata.scaling.P001=diag([1 1 1])*1e-6;
PathBound_F_L=zeros(1,3);%ones(1,3)*auxdata.scaling.F_L;
PathBound_M_L=zeros(1,3);%ones(1,3)*auxdata.scaling.M_L;
PathBound_F_R=zeros(1,3);%ones(1,3)*auxdata.scaling.F_R;
PathBound_M_R=zeros(1,3);%ones(1,3)*auxdata.scaling.M_R;
PathBound_P000=ones(1,3)*auxdata.scaling.P000;
PathBound_P010=ones(1,3)*auxdata.scaling.P010;
PathBound_P001=ones(1,3)*auxdata.scaling.P001;
if SCALA
    PathBound_F_L=PathBound_F_L/auxdata.scaling.F_L;
    PathBound_M_L=PathBound_M_L/auxdata.scaling.M_L;
    PathBound_F_R=PathBound_F_R/auxdata.scaling.F_R;
    PathBound_M_R=PathBound_M_R/auxdata.scaling.M_R;
    PathBound_P000=PathBound_P000/auxdata.scaling.P000;
    PathBound_P010=PathBound_P010/auxdata.scaling.P010;
    PathBound_P001=PathBound_P001/auxdata.scaling.P001;
end
bounds.phase.path.lower=-[PathBound_F_L PathBound_M_L PathBound_F_R PathBound_M_R PathBound_P000 PathBound_P010 PathBound_P001];
bounds.phase.path.upper= [PathBound_F_L PathBound_M_L PathBound_F_R PathBound_M_R PathBound_P000 PathBound_P010 PathBound_P001];

% Eventgroup  -> impose periodicity
Ncoord=length(bounds.phase.state.lower);
bounds.eventgroup.lower=-ones(1,Ncoord-6)*0.01;
bounds.eventgroup.upper=ones(1,Ncoord-6)*0.01;

% Static parameters
bounds.parameter.lower=[100];       % stiffness lower
bounds.parameter.upper=[250];       % stiffness upper
auxdata.scalingoffsetparameter=mean([bounds.parameter.lower;bounds.parameter.upper],1);
auxdata.scalingparameter=diag(max([bounds.parameter.lower;bounds.parameter.upper]-[auxdata.scalingoffsetparameter;auxdata.scalingoffsetparameter],[],1));
if SCALA
    bounds.parameter.lower=(bounds.parameter.lower-auxdata.scalingoffsetparameter)/auxdata.scalingparameter;
    bounds.parameter.upper=(bounds.parameter.upper-auxdata.scalingoffsetparameter)/auxdata.scalingparameter;
end
% Initial guess
guess.phase.time(:,1)=T(LineIndexForGuess);
guess.phase.control=[[A_L A_R] Fp Mp];
guess.phase.state=[[Q_L Q_R] [U_L U_R]];
guess.phase.integral=1e-5;
guess.parameter=[120];
bounds.phase.integral.upper = 1;
bounds.phase.integral.lower = 0;
if SCALA
    guess.phase.control=guess.phase.control/blkdiag(auxdata.scaling.A,auxdata.scaling.Fp,auxdata.scaling.Mp);
    guess.phase.state=guess.phase.state/blkdiag(auxdata.scaling.Q,auxdata.scaling.U);
    guess.parameter=(guess.parameter-auxdata.scalingoffsetparameter)/auxdata.scalingparameter;
end

auxdata.LEFTMODELNAME=LEFTMODELNAME;
auxdata.RIGHTMODELNAME=RIGHTMODELNAME;
auxdata.Ncoord_L=Ncoord_L;
auxdata.Ncoord_R=Ncoord_R;
auxdata.PT_posB_pelv=[0 0 0
    1 0 0
    0 1 0]'/10;
auxdata.PT_bodynames={'pelvis','pelvis','pelvis'};
% add marker positions to the points to be tracked
auxdata.PT_marker_L=auxdata.PT_posB_pelv;
auxdata.PT_markerbody_L=auxdata.PT_bodynames;
TRC_L=[];
for q=1:size(Markerinfo_L,2)
    auxdata.PT_markerbody_L{end+1}=Markerinfo_L{2,q};
    auxdata.PT_marker_L=[auxdata.PT_marker_L Markerinfo_L{3,q}];
    TRC_L=[TRC_L Markerinfo_L{4,q}];
end
auxdata.PT_marker_R=auxdata.PT_posB_pelv;
auxdata.PT_markerbody_R=auxdata.PT_bodynames;
TRC_R=[];
for q=1:size(Markerinfo_R,2)
    auxdata.PT_markerbody_R{end+1}=Markerinfo_R{2,q};
    auxdata.PT_marker_R=[auxdata.PT_marker_R Markerinfo_R{3,q}];
    TRC_R=[TRC_R Markerinfo_R{4,q}];
end

%% adjust trc file => track KS solution and not raw markerdata

Q_L_Xcpp(:,IND_coord_OK_L)=Q_L;
Q_R_Xcpp(:,IND_coord_OK_R)=Q_R;
U_L_Xcpp(:,IND_coord_OK_L)=U_L;
U_R_Xcpp(:,IND_coord_OK_R)=U_R;
Q_L_Xcpp(:,IND_coord_NO_L)=0;
Q_R_Xcpp(:,IND_coord_NO_R)=0;
U_L_Xcpp(:,IND_coord_NO_L)=0;
U_R_Xcpp(:,IND_coord_NO_R)=0;
ZERO_F=zeros(length(Q_L_Xcpp),3);

[T_cpp_L,PT_posG_cpp_L]=...
    trovaIDePATH_MD(auxdata.LEFTMODELNAME,Q_L_Xcpp,U_L_Xcpp,U_L_Xcpp,auxdata.PT_marker_L,auxdata.PT_markerbody_L,-Fp,ZERO_F,-Mp,ZERO_F);
[T_cpp_R,PT_posG_cpp_R]=...
    trovaIDePATH_MD(auxdata.RIGHTMODELNAME,Q_R_Xcpp,U_R_Xcpp,U_R_Xcpp,auxdata.PT_marker_R,auxdata.PT_markerbody_R,Fp,ZERO_F,Mp,ZERO_F);

nCpt=3;
nMark_L=length(auxdata.PT_markerbody_L)-nCpt;
nMark_R=length(auxdata.PT_markerbody_R)-nCpt;
TRCcpp_L=reshape(PT_posG_cpp_L(:,nCpt+1:end),size(PT_posG_cpp_L,1)/3,nMark_L*3);
TRCcpp_R=reshape(PT_posG_cpp_R(:,nCpt+1:end),size(PT_posG_cpp_R,1)/3,nMark_R*3);

Markers=[TRCcpp_L TRCcpp_R];
MLabels=[Markerinfo_L(1,:) Markerinfo_R(1,:)];
trcfile=fullfile(pwd,'KSMarker.trc');
Rate=1./nanmean(diff(T));
Frames=1:length(T);
Time=T;
Units='m';
writeMarkersToTRC(trcfile, Markers, MLabels, Rate, Frames', Time', Units);
%% spline trc file
% figure();plot(T,TRCcpp_L,'r'); hold on;
% plot(MarkerTime_L,TRC_L,'--b');
auxdata.ppTRC_L=spline(T,TRCcpp_L');
auxdata.ppTRC_R=spline(T,TRCcpp_R');
% auxdata.ppTRC_L=spline(MarkerTime_L,TRC_L');
% auxdata.ppTRC_R=spline(MarkerTime_R,TRC_R');
auxdata.ContactInfo=OUT_tuneGRF;
auxdata.NomiCoord_L=NomiCoord_L;
auxdata.NomiCoord_R=NomiCoord_R;
auxdata.GRFhead=GRFinput.head;
auxdata.SCALA=SCALA;
%% get index for the bodynames
import org.opensim.modeling.*
m=Model(LEFTMODELNAME);
bodies=m.getBodySet();
for i=0:bodies.getSize()-1
    BodyNames{i+1}=char(bodies.get(i).getName());
end
IndsBodies=zeros(1,length(auxdata.PT_markerbody_L));
for i=1:length(auxdata.PT_markerbody_L)
    IndsBodies(i)=find(strcmp(BodyNames,auxdata.PT_markerbody_L{i}));
end
auxdata.IndsBodies=IndsBodies;

m_R=Model(RIGHTMODELNAME);
bodies_R=m_R.getBodySet();
for i=0:bodies_R.getSize()-1
    BodyNames_R{i+1}=char(bodies_R.get(i).getName());
end
IndsBodies_R=zeros(1,length(auxdata.PT_markerbody_R));
for i=1:length(auxdata.PT_markerbody_R)
    IndsBodies_R(i)=find(strcmp(BodyNames_R,auxdata.PT_markerbody_R{i}));
end
auxdata.IndsBodies_R=IndsBodies_R;

%% GPOPS setup
setup.name='AAA';
setup.auxdata=auxdata;
setup.bounds=bounds;
setup.guess=guess;
setup.nlp.solver='ipopt';
% setup.nlp.ipoptoptions.linear_solver='mumps';
setup.nlp.ipoptoptions.linear_solver='ma57';
setup.derivatives.supplier='sparseFD';
setup.derivatives.dependencies = 'sparseNaN';
setup.derivatives.derivativelevel='first';
setup.nlp.ipoptoptions.tolerance=1e-5;			% I believe that this should be higher !
setup.nlp.ipoptoptions.maxiterations=2000;
% setup.derivatives.supplier = 'adigator'; % se si passa per opensim non dovrebbe essere possibile
if SCALA
    setup.scales.method='none';
else
    setup.scales.method='automatic-guessUpdate';
end
% setup.mesh.method='hp-DarbyRao';
% setup.mesh.method='hp-LiuRao-Legendre';
setup.mesh.method='hp-LiuRao';
setup.mesh.tolerance=1e-4;
setup.mesh.maxiterations=0;
setup.mesh.colpointsmin=NCP;
setup.mesh.colpointsmax=20;
setup.derivatives.stepsize1 = 1e-6;
% setup.method='RPM-integration';
setup.method='RPM-differentiation';
setup.displaylevel=2;
setup.mesh.phase.colpoints=NCP*ones(1,NMeshIntervals);
setup.mesh.phase.fraction=(1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous=HandleCONTINUOUS;
setup.functions.endpoint=HandleENDPOINT;


%% solve optimization with GPOPS
% close all
figure
subplot(1,2,1)
plot(TxGRF,Q_L(:,4:6))
subplot(1,2,2)
plot(TxGRF,Q_R(:,4:6))
tic
varargout{1}=gpops2(setup);
toc
