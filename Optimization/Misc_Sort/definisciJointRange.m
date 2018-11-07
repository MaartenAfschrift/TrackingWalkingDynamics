function [Range,DofName]=definisciJointRange(ModelName)
import org.opensim.modeling.*
m=Model(ModelName);
CS=m.getCoordinateSet();
Q=CS.getSize();
for q=1:Q
    C=CS.get(q-1);
    DofName{q}=C.getName().toCharArray()';
    Range(q,:)=sort([C.getRangeMax() C.getRangeMin()]);
end
m.disownAllComponents();