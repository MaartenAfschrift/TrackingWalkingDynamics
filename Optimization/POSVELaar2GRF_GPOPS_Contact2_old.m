function [R,M,qR,Mq,pCont]=POSVELaar2GRF_GPOPS_Contact2_old(PV_PT,ContactInfo,Stiffness,u,c,vt)

% compute ground contact forces and moments

L=size(PV_PT(1).PiG,1);
nPT=length(ContactInfo.VectorRadius);
R=zeros(L,3);
M=zeros(L,3);
Mglob=zeros(L,3);
for q=1:nPT
    pCont(q).PinG_wrt000=PV_PT(q).PinG_wrt000;
    pCont(q).PiG=PV_PT(q).PiG;
    pCont(q).ViG=PV_PT(q).ViG;
    PT(q).params.k=Stiffness;
    PT(q).params.y0=ContactInfo.VectorRadius(q);
    PT(q).params.u=u;
    PT(q).params.c=c;
    pCont(q).GRF=GRF_xSim2(pCont(q).PiG,pCont(q).ViG,PT(q).params,vt);
    R=R+pCont(q).GRF;
    M=M+cross(pCont(q).GRF,pCont(q).PinG_wrt000);
    Mglob=Mglob+cross(pCont(q).GRF,pCont(q).PiG);
end



if nargout > 2
    nR=sum(R.^2,2)*[1 1 1];
    dq=1./nR.*cross(R,Mglob);
    qR=dq-R.*(dq(:,2)./(R(:,2))*[1 1 1]);
    qR(nR==0)=0;
    Mq=zeros(L,3);
    for q=1:nPT
        Mq=Mq+cross(pCont(q).GRF,pCont(q).PiG+qR);
    end
    qR=-qR;
end
end