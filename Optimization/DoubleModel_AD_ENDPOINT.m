function OUT=DoubleModel_AD_ENDPOINT(IN)
global III 

q = IN.phase.integral;
OUT.objective = q(1);

% symmetry constraint between initial and final state
tmp=IN.phase.initialstate(:,1:(IN.auxdata.Ncoord_L+IN.auxdata.Ncoord_R)*2)-IN.phase.finalstate(:,1:(IN.auxdata.Ncoord_L+IN.auxdata.Ncoord_R)*2);
tmp(:,[4:6 (4:6)+IN.auxdata.Ncoord_L])=[];
tmpF=IN.phase.finalstate(:,1:(IN.auxdata.Ncoord_L+IN.auxdata.Ncoord_R)*2)-IN.phase.finalstate(:,1:(IN.auxdata.Ncoord_L+IN.auxdata.Ncoord_R)*2);
tmpF(:,[4:6 (4:6)+IN.auxdata.Ncoord_L])=[];
OUT.eventgroup.event=tmp-tmpF;
III=q;