function [F]=GRF_xSim2(y,v,params,vt)


% Compute the ground reaction forces based on indentation and indentation velocity

us 				= params.u;
E1 				= params.k;        % stiffness
c  				= params.c;
radius_heel_r 	= params.y0;

pdot_heel_r=-v(:,2);
vs_heel_r=v(:,[1 3]);
p_heel_r=-(y(:,2)-radius_heel_r); 			%y(:,2);%-(y(:,2)-radius_heel_r);           % indentation
E=E1*(1/2)^(3/2);                           % adjusted stiffness ?
k_heel_r =(4/3)*sqrt(radius_heel_r ).*E;    % stiffness adjusted again ?
III=    p_heel_r>0;                         % III => when there is contact
F1=(10000000*k_heel_r.*p_heel_r .^3)+p_heel_r/10000;
F2=p_heel_r/10000;
    
fn_heel_r(III,1)=F1(III);                   % only F1 if there is contact
fn_heel_r(~III,1)=F2(~III);                 % otherwise equation F2
fn_heel_r=fn_heel_r.*(1+c.*pdot_heel_r);            % add demping
fSx=-fn_heel_r.*tanh(vs_heel_r(:,1)./vt).*us;   % compute x force from normal force and x-velocity)
fSz=-fn_heel_r.*tanh(vs_heel_r(:,2)./vt).*us;   % compute z force from normal force and z-velocity) 
    
F=[fSx fn_heel_r fSz];
end