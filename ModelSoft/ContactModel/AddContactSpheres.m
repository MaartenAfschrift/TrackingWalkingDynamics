function [OUTname] = AddContactSpheres(S,ModelName,BodyName,MatrixLoc,VectorRadius)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


import org.opensim.modeling.*
M=Model(ModelName);
CGS=M.updContactGeometrySet();
PLAT=ContactHalfSpace();
PLAT.setName('Platform');
PLAT.setBodyName(M.getBodySet.get(0).getName());
PLAT.setLocation(Vec3(0,0,0));
PLAT.setOrientation(Vec3(0,0,-pi/2));
PLAT.set_display_preference(1);
CGS.cloneAndAppend(PLAT);


% createVector with contact spheres

Side=BodyName(end);
for i=1:length(VectorRadius);
    % add ContactSpheres
    
    F=['ContactElemnt'  num2str(i) '_' Side];
    N=['ContactGeom'    num2str(i) '_' Side];
    tmpCS.(F)=ContactSphere();
    tmpCS.(F).setName(N);
    tmpCS.(F).setBodyName(BodyName);
    tmpCS.(F).setRadius(VectorRadius(i));
    tmpCS.(F).setLocation(Vec3(MatrixLoc(i,1), MatrixLoc(i,2), MatrixLoc(i,3)));
    tmpCS.(F).setOrientation(Vec3(0,0,0));
    tmpCS.(F).set_display_preference(4);
    CG.(F)=ContactGeometry.safeDownCast(tmpCS.(F));
    HCF.(F)=HuntCrossleyForce();
    HCF.(F).setName(['Forza_',N]);
    HCF.(F).setStiffness(S.stiffness);
    HCF.(F).setDissipation(S.dissipation);
    HCF.(F).setStaticFriction(1);
    HCF.(F).setDynamicFriction(1);
    HCF.(F).setViscousFriction(1);
    HCF.(F).addGeometry('Platform');
    HCF.(F).addGeometry(N);
    FOR.(F)=Force.safeDownCast(HCF.(F));
    M.addForce(FOR.(F));
    CGS.cloneAndAppend(CG.(F));
end

OUTname=[ModelName(1:end-5) '_contact.osim'];
M.print(OUTname);
M.disownAllComponents();

end

