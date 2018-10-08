function [OUT]=splitModelKS_m(NOMEMODELLO,NOMETRIAL,NOMEGRF,CARTELLAOUTPUT,NOMEID)
% load inverse kinematics data
fid=fopen(NOMETRIAL,'r');
[IKhead,IKdata]=readosim(fid);
T=IKdata(:,1);
fclose(fid);
% load GRF experimental measures
fid=fopen(NOMEGRF,'r');
[GRFhead,GRFdata]=readosim(fid);
fclose(fid);
% load Inverse Dynamics
fid=fopen(NOMEID,'r');
[IDhead,IDdata]=readosim(fid);
IDhead=IDhead(2:end);
fclose(fid);
IDtime=IDdata(:,1);
ppID=spline(IDtime,IDdata(:,2:end)');
% filter kinematics
frqFilt=15;
frqKin=1/mean(diff(T));
[b,a]=butter(6,frqFilt/(frqKin/2),'low');
IKdataFilt=[T filtfilt(b,a,IKdata(:,2:end))];
% compute velocity
splPos=spline(T,IKdataFilt(:,2:end)');
splVel=fnder(splPos,1);
% splAcc=fnder(splPos,2);
% splJer=fnder(splPos,3);
% built matrix for mex functions
nomi_corpi={'calcn_l','calcn_r'};           % nomi_corpi={'pelvis};      zeros(3,1)
Q=ppval(splPos,T)';
indDEG=true(1,size(Q,2));
indDEG([4 5 6])=false;
U=ppval(splVel,T)';
Q(:,indDEG)=Q(:,indDEG)/180*pi;
U(:,indDEG)=U(:,indDEG)/180*pi;
[POSarr,VELarr,ANGarr,ANGVELarr]=Kinematics_at_feet(NOMEMODELLO,Q,U,nomi_corpi,zeros(3,2));
POS.L=reshape(POSarr(:,1),length(POSarr(:,1))/3,3);
POS.R=reshape(POSarr(:,2),length(POSarr(:,2))/3,3);
VEL.L=reshape(VELarr(:,1),length(VELarr(:,1))/3,3);
VEL.R=reshape(VELarr(:,2),length(VELarr(:,2))/3,3);
ANG.L=reshape(ANGarr(:,1),length(ANGarr(:,1))/3,3);
ANG.R=reshape(ANGarr(:,2),length(ANGarr(:,2))/3,3);
ANGVEL.L=reshape(ANGVELarr(:,1),length(ANGVELarr(:,1))/3,3);
ANGVEL.R=reshape(ANGVELarr(:,2),length(ANGVELarr(:,2))/3,3);
%% write kinematics for the invreted models
stampaIKdoublemodel([CARTELLAOUTPUT,'\INVERSEMODELKINEMATICS_L.mot'],IKdataFilt,IKhead,ANG.L,POS.L);
stampaIKdoublemodel([CARTELLAOUTPUT,'\INVERSEMODELKINEMATICS_R.mot'],IKdataFilt,IKhead,ANG.R,POS.R);
%%
OUT.GRFhead=GRFhead;
OUT.GRFdata=GRFdata;
OUT.T=T;
OUT.POS=POS;
OUT.VEL=VEL;
OUT.ANG=ANG;
OUT.ANGVEL=ANGVEL;
OUT.Q=Q;
OUT.U=U;
OUT.Qhead=IKhead;
OUT.IDtime=IDtime;
OUT.ppID=ppID;
OUT.IDhead=IDhead;
% OUT.=;
% OUT.=;
% OUT.=;