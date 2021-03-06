function t_out = Euler2Mat_vector_tot(X,Y,Z,px,py,pz)
%EULER2MAT_VECTOR_TOT
%    T_OUT = EULER2MAT_VECTOR_TOT(X,Y,Z,PX,PY,PZ)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    04-Aug-2017 12:06:50


% Function to convert euler angles to rotation matrices constructed with the matlab symbolic toolbox to work with vectors (instead of for loop over all the collocation points)


t2 = cos(Y);
t3 = sin(Z);
t4 = cos(Z);
t5 = sin(Y);
t6 = cos(X);
t7 = sin(X);
t_out = [pz.*t5+px.*t2.*t4-py.*t2.*t3;px.*(t3.*t6+t4.*t5.*t7)+py.*(t4.*t6-t3.*t5.*t7)-pz.*t2.*t7;px.*(t3.*t7-t4.*t5.*t6)+py.*(t4.*t7+t3.*t5.*t6)+pz.*t2.*t6];
