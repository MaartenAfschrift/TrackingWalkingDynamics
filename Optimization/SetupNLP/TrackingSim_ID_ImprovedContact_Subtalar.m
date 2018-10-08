function [output] = TrackingSim_ID_ImprovedContact_Subtalar(S,MainPath,MPI_path,OutName,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

global III IIIarr BoolFirst
IIIarr=0;
BoolFirst =true;

ContinuousFunctionName='DoubleModel_AD_CONTINUOUS_mpi_Contact';
if ~isempty(varargin)
    ContinuousFunctionName=varargin{1};
end


%% Path Information

% To DO: remove this in the figure, always use S structue
OutPath     = S.OutPath;
ModelFile   = S.ModelFile;
KS_File     = S.IK_File;
GRFFile     = S.GRFFile;
trcFile     = S.trcFile;
ID_File     = S.ID_File;

% Check if we are still tracking the correct ID resutls !!
%% Translate to italian (Lorenzo Pitto his code)
cd(MainPath);
if ~isdir(OutPath); mkdir(OutPath); end
LefModelName=[OutPath,'\LeftSideModel.osim'];
RightModelName=[OutPath,'\RightSideModel.osim'];

SINGLEDOUBLE='d';
%% compiles the mex files

%% split model
% splitModel(NOMEMODELLO,CARTELLADATI,CARTELLAOUTPUT);
[L,R,~,Markerinfo_L,Markerinfo_R,Markerinfo_Sing]=invertedModel_reverseJoint(ModelFile,OutPath);
%% import marker trajectories to be used in the first optimization
[Markerinfo_L,MarkerTime_L]=importMarkerTRC_m(trcFile,Markerinfo_L);
[Markerinfo_R,MarkerTime_R]=importMarkerTRC_m(trcFile,Markerinfo_R);
%% generates the IK files for the split models
[OUT_split]=splitModelKS_m(ModelFile,KS_File,GRFFile,OutPath,ID_File);
GRFdata=OUT_split.GRFdata;
GRFhead=OUT_split.GRFhead;
T=OUT_split.T;
POS=OUT_split.POS;
VEL=OUT_split.VEL;
ANG=OUT_split.ANG;
ANGVEL=OUT_split.ANGVEL;
Q=OUT_split.Q;
U=OUT_split.U;
%% position of contact bodies with respect to ground


[OUT_tuneGRF]=tuneContact_Adj2(S,OutPath);


%% MEX INPUT
KS=importdata(KS_File);
TxGRF=linspace(KS.data(1,1),KS.data(end,1),100);
GRF_exp=importdata(GRFFile);
[~, nc]=size(GRF_exp.data);
GRFopt=zeros(100,nc-1);
for i=2:length(GRF_exp.data(1,:))
    GRFopt(:,i-1)=interp1(GRF_exp.data(:,1),GRF_exp.data(:,i),TxGRF);
end
GRFinput.T=TxGRF;
GRFinput.data=[TxGRF' GRFopt];
GRFinput.head=strsplit(GRF_exp.textdata{end});
[OUT_L]=DoubleModelInitailGuess(LefModelName,[OutPath,'\INVERSEMODELKINEMATICS_L.mot'],GRFinput);
[OUT_R]=DoubleModelInitailGuess(RightModelName,[OutPath,'\INVERSEMODELKINEMATICS_R.mot'],GRFinput);
printGRFmotFile(GRFinput.data,GRFinput.head,[OutPath,'\tunedGRF.mot'])
%% define splines for state and controls initial guess

ppIDindL=[];
ppIDindR=[];
ppIDind=[];
for q=1:length(OUT_L.Qhead)
    tmpind=find(strncmp(OUT_split.IDhead,OUT_L.Qhead{q},length(OUT_L.Qhead{q})));
    if ~isempty(tmpind)
        ppIDindL=[ppIDindL [q;tmpind]];
    end
end
for q=1:length(OUT_R.Qhead)
    tmpind=find(strncmp(OUT_split.IDhead,OUT_R.Qhead{q},length(OUT_R.Qhead{q})));
    if ~isempty(tmpind)
        ppIDindR=[ppIDindR [q;tmpind]];
    end
end
OUT_split.Qhead(strcmp(OUT_split.Qhead,'time'))=[];
for q=1:length(OUT_split.Qhead)
    tmpind=find(strncmp(OUT_split.IDhead,OUT_split.Qhead{q},length(OUT_split.Qhead{q})));
    if ~isempty(tmpind)
        ppIDind=[ppIDind [q;tmpind]];
    end
end
auxdata.ppID=OUT_split.ppID;
auxdata.ppIDindL=ppIDindL;
auxdata.ppIDindR=ppIDindR;
auxdata.ppIDind=ppIDind;
% index to definet the time points to be used in the initial guess
T=TxGRF;
LineIndexForGuess=1:length(T);
% -=- Model Single -=-
% Q
ppQ=spline(OUT_split.T,OUT_split.Q');
Q=ppval(ppQ,TxGRF(LineIndexForGuess))';
% U
ppU=fnder(ppQ,1);
U=ppval(ppU,TxGRF(LineIndexForGuess))';
% A
ppA=fnder(ppQ,2);
A=ppval(ppA,TxGRF(LineIndexForGuess))';
% J
ppJ=fnder(ppQ,3);
J=ppval(ppJ,TxGRF(LineIndexForGuess))';
% F_R
ppF_R=spline(TxGRF,OUT_R.F_ccp_inG_FINAL');
F_R=ppval(ppF_R,TxGRF(LineIndexForGuess))';
% M_R
ppM_R=spline(TxGRF,OUT_R.T_ccp_inG_FINAL');
M_R=ppval(ppM_R,TxGRF(LineIndexForGuess))';
% F_L
ppF_L=spline(TxGRF,OUT_L.F_ccp_inG_FINAL');
F_L=ppval(ppF_L,TxGRF(LineIndexForGuess))';
% M_L
ppM_L=spline(TxGRF,OUT_L.T_ccp_inG_FINAL');
M_L=ppval(ppM_L,TxGRF(LineIndexForGuess))';
% -=- Model L -=-
% Q
Q_L=ppval(OUT_L.ppQ,T(LineIndexForGuess))';
% U
U_L=ppval(OUT_L.ppU,T(LineIndexForGuess))';
% A
A_L=ppval(OUT_L.ppA,T(LineIndexForGuess))';
% J
J_L=ppval(OUT_L.ppJ,T(LineIndexForGuess))';
% Fp
ppFp_L=spline(T,OUT_L.Fp');
Fp_L=ppval(ppFp_L,T(LineIndexForGuess))';
% Mp
ppMp_L=spline(T,OUT_L.Mp');
Mp_L=ppval(ppMp_L,T(LineIndexForGuess))';
% -=- Model R -=-
% Q
Q_R=ppval(OUT_R.ppQ,T(LineIndexForGuess))';
% U
U_R=ppval(OUT_R.ppU,T(LineIndexForGuess))';
% A
A_R=ppval(OUT_R.ppA,T(LineIndexForGuess))';
% J
J_R=ppval(OUT_R.ppJ,T(LineIndexForGuess))';
% Fp
ppFp_R=spline(T,OUT_R.Fp');
Fp_R=ppval(ppFp_R,T(LineIndexForGuess))';
% Mp
ppMp_R=spline(T,OUT_R.Mp');
Mp_R=ppval(ppMp_R,T(LineIndexForGuess))';
% for the right side is the difference with the GRF, for the left side FR+res
Fp=(Fp_R-Fp_L)/2;
Mp=(Mp_R-Mp_L)/2;
FpRES=(Fp_R+Fp_L)/2;
MpRES=(Mp_R+Mp_L)/2;


%% Construct the optimal control problem
if strcmpi(SINGLEDOUBLE,'d')
    GRDdataXopt.GRFdata=GRFdata;
    GRDdataXopt.TxGRF=TxGRF;
    GRDdataXopt.head=GRFinput.head;
    StCo.Q_L=Q_L;
    StCo.Q_R=Q_R;
    StCo.U_L=U_L;
    StCo.U_R=U_R;
    StCo.A_L=A_L;
    StCo.A_R=A_R;
    MArkerInfoXopt.MarkerTime_L=MarkerTime_L;
    MArkerInfoXopt.MarkerTime_R=MarkerTime_R;
    MArkerInfoXopt.Markerinfo_L=Markerinfo_L;
    MArkerInfoXopt.Markerinfo_R=Markerinfo_R;
    IIIarr=0;
    III=0;
    
    disp(' Open Models in command line as administrator and hit ENTER when finished');
    disp([' cd ' MPI_path]);
    disp(' mpiexec -np 12 opensimMPI_muscle.exe LeftSideModel.osim RightSideModel.osim');
    copyfile(LEFTMODELNAME,[MPI_path '\LeftSideModel.osim']);
    copyfile(RIGHTMODELNAME,[MPI_path '\RightSideModel.osim']);
%     SwitchModels();
    pause(1);       % pause 1 second to make sure that models are loaded before starting optimization
    if S.ifPause
        pause();    
    end
    output=runFirstOptimization_AccContr_GRFadj('NMeshIntervals',40,'ScalingBound',1.5,...
        'NColPoints',4,'LeftModelName',LEFTMODELNAME,'scala',false,...
        'RightModelName',RIGHTMODELNAME,'GRFdata',GRDdataXopt,...
        'Time',T,'MarkerInfo',MArkerInfoXopt,'OUTtuneGRF',OUT_tuneGRF,...
        'HandleContinuous',ContinuousFunctionName,'HandleEndpoint','DoubleModel_AD_ENDPOINT',...
        'OUTinitGuessL',OUT_L,'OUTinitGuessR',OUT_R,...
        'LineIndexForGuess',LineIndexForGuess,'StatesControls',StCo,...
        'ForcePelvis',Fp,'TorquePelvis',Mp,'auxdata',auxdata,'stampa',true);
    
end
save(fullfile(OutPath,[OutName '.mat']),'output');


end

