function phaseout = musdynContinous_FtildeState_GM_old(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);
Scale_GM        = input.phase.parameter(:,1);

% Get controls
e   = input.phase.control(:,1:NMuscles);
aT  = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde  = 10*input.phase.control(:,NMuscles+Ndof+1:end);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:end);

% PATH CONSTRAINTS
% Hill-equilibrium constraint
[Hilldiff,F] = ForceEquilibrium_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,...
	splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam,input.auxdata.Atendon);

% Moments constraint
Topt = 150;
Tdiff = zeros(numColPoints,Ndof);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    Tdiff(:,dof) =  (T_exp-T_sim);
end

phaseout.path = [Tdiff Hilldiff];

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct(m),tauDeact(m),input.auxdata.b);
end

% Contraction dynamics is implicit
phaseout.dynamics = [dadt dFtilde];

% compute constraint on GM activity
SimAct1=a(:,input.auxdata.IndsMuscle(1));
SimAct2=a(:,input.auxdata.IndsMuscle(2));
SimAct3=a(:,input.auxdata.IndsMuscle(3));
MeasAct=splinestruct.EMG.*Scale_GM;
GM_Error=(SimAct1-MeasAct).^2+(SimAct2-MeasAct).^2+(SimAct3-MeasAct).^2;

% OBJECTIVE FUNCTION
w1 = 1000;
w2 = 50;
phaseout.integrand = sum(e.^2,2)+ w1.*sum(aT.^2,2)+w2.*GM_Error;





