%% [mat]=Euler2mat(r)
% mat:  rotation matrix - Rx*Ry*Rz
% r:    [rx, ry, rz] rad

function [OUT]=Euler2mat(r,ordine)
if nargin==1
    X=r(1);
    Y=r(2);
    Z=r(3);
    
    Rx=[1 0       0
        0 cos(X) -sin(X)
        0 sin(X)  cos(X)];
    
    Ry=[cos(Y) 0 sin(Y)
        0      1 0
        -sin(Y) 0 cos(Y)];
    
    Rz=[cos(Z) -sin(Z) 0
        sin(Z)  cos(Z) 0
        0       0      1];
    
    OUT=Rx*Ry*Rz;
elseif nargin==2
%     switch ordine
%         case {'zxy','ZXY'}
%             X=r(2);
%             Y=r(3);
%             Z=r(1);
%         case {'zyx','ZYX'}
%             X=r(3);
%             Y=r(2);
%             Z=r(1);
%         case{'xyz','XYZ'}
            X=r(1);
            Y=r(2);
            Z=r(3);
%     end
    Rx=[1 0       0
        0 cos(X) -sin(X)
        0 sin(X)  cos(X)];
    
    Ry=[cos(Y) 0 sin(Y)
        0      1 0
        -sin(Y) 0 cos(Y)];
    
    Rz=[cos(Z) -sin(Z) 0
        sin(Z)  cos(Z) 0
        0       0      1];
    switch ordine
        case 'zxy'
            OUT=(Rz*Rx)*Ry;
        case 'ZXY'
            OUT=Ry*(Rx*Rz);
        case 'zyx'
            OUT=(Rz*Ry)*Rx;
        case 'ZYX'
            OUT=Rx*(Ry*Rz);
        case 'xyz'
             OUT=(Rx*Ry)*Rz;
        case 'XYZ'
             OUT=Rz*(Ry*Rx);
    end
end