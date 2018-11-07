function [T_ccp_L_inG_fast,T_ccp_R_inG_fast ] = ComputeGRFConstraint_vect(T_cpp_L,T_cpp_R,Q_L,Q_R)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% This function converts the residuals forces and moments from ID to the world frame

testoutL= GRF_Constraint_Sym(T_cpp_L(:,1),T_cpp_L(:,2),T_cpp_L(:,3),Q_L(:,1),Q_L(:,2));		% function constructed with matlab symbolic toolbox (ToDo: add original function)
testoutR= GRF_Constraint_Sym(T_cpp_R(:,1),T_cpp_R(:,2),T_cpp_R(:,3),Q_R(:,1),Q_R(:,2));		% function constructed with matlab symbolic toolbox (ToDo: add original function)
nfr=length(testoutL)/3;
T_ccp_L_inG_fast=[testoutL(1:nfr) testoutL(nfr+1:nfr*2)  testoutL(nfr*2+1:nfr*3)];
T_ccp_R_inG_fast=[testoutR(1:nfr) testoutR(nfr+1:nfr*2)  testoutR(nfr*2+1:nfr*3)];

end

