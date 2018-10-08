function [OUT]=tuneContact_Adj2(S,ModelFolder)
% DOOPTIMIZE=PARAMS.DOOPTIMIZE;

% settings
MatrixLoc=[ 0.03 0.01 0;
    0.175  -0.015   -0.03;
    0.175   -0.015  0.04;
    0.175   -0.015  0.005];
VectorRadius    = [0.035   0.015    0.015   0.012];
S.stiffness     = 400;
S.dissipation   = 1;

% add contact spheres to left model
ModelName=[ModelFolder,'\LeftSideModel.osim'];       BodyName='calcn_l';
AddContactSpheres(S,ModelName,BodyName,MatrixLoc,VectorRadius);

% add contact spheres to right model
ModelName=[ModelFolder,'\RightSideModel.osim'];       BodyName='calcn_r';
AddContactSpheres(S,ModelName,BodyName,MatrixLoc,VectorRadius);

% save output
OUT.MatrixLoc       = MatrixLoc;
OUT.VectorRadius    = VectorRadius;
OUT.BodyNameL       = 'calcn_l';
OUT.BodyNameR       = 'calcn_r';
