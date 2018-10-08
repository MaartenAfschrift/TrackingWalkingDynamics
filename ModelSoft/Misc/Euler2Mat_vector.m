function OUT = Euler2Mat_vector(X,Y,Z)
%EULER2MAT_VECTOR
%    OUT = EULER2MAT_VECTOR(X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    04-Aug-2017 12:00:57

t2 = cos(Y);
t3 = sin(Z);
t4 = cos(Z);
t5 = sin(Y);
t6 = cos(X);
t7 = sin(X);
OUT = reshape([t2.*t4,t3.*t6+t4.*t5.*t7,t3.*t7-t4.*t5.*t6,-t2.*t3,t4.*t6-t3.*t5.*t7,t4.*t7+t3.*t5.*t6,t5,-t2.*t7,t2.*t6],[3,3]);
