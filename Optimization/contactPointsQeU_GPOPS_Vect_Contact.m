% assumes that the joint of the foottoground is in (0 0 0)
function [PV_PT]=contactPointsQeU_GPOPS_Vect_Contact(ANG,POS,ANGVEL,VEL,ContactInfo)

% get indentation and indentation velocity of the contact spheres

NF=size(ANG,1);
for q=1:length(ContactInfo.VectorRadius)
    PiB=ContactInfo.MatrixLoc(q,:);
    PinG_wrtB000_temp=Euler2Mat_vector_tot(ANG(:,1),ANG(:,2),ANG(:,3),PiB(1),PiB(2),PiB(3));
    PinG_wrtB000=[PinG_wrtB000_temp(1:NF) PinG_wrtB000_temp(NF+1:NF*2) PinG_wrtB000_temp(NF*2+1:NF*3)];
    PiG=PinG_wrtB000+POS;
    ViG=(cross(ANGVEL,PinG_wrtB000))+VEL; 
    PV_PT(q).PiG=PiG;
    PV_PT(q).PinG_wrt000=PinG_wrtB000;
    PV_PT(q).ViG=ViG;
end