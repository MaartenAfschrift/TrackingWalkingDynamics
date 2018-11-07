function [Q_L_Xcpp,Q_R_Xcpp,U_L_Xcpp ,U_R_Xcpp,A_L_Xcpp,A_R_Xcpp,Fp,Mp,Q_L,Q_R,U_L,U_R] = Get_KinematicsInfo_TrackingSim(OptRes,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
bool_scaled=0;
if ~isempty(varargin)
    bool_scaled=varargin{1};
end



auxdata=OptRes.output.result.setup.auxdata;
solution=OptRes.output.result.solution;
state=solution.phase.state;

t=solution.phase.time;
N=length(t);

if bool_scaled
    solution.phase.state=solution.phase.state.*repmat(auxdata.scales.Kx,N,1);
    solution.phase.control=solution.phase.control.*repmat(auxdata.scales.Ku,N,1);
    solution.parameter=solution.parameter.*auxdata.scales.Kparam;
end

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

end

