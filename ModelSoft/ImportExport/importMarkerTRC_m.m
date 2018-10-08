function [MarkerinfoOUT,T]=importMarkerTRC_m(NOMETRC,Markerinfo)
fid=fopen(NOMETRC);
[TRChead,TRCdata]=readTRC_PredictiveSimulations(fid);
fclose(fid);
NumMark=1;
IndsDel=[];
for q=1:size(Markerinfo,2)
    I=strcmp(Markerinfo{1,q},TRChead);
    if ~isempty(I)
        MarkerinfoOUT{1,NumMark}=Markerinfo{1,q};
        MarkerinfoOUT{2,NumMark}=Markerinfo{2,q};
        MarkerinfoOUT{3,NumMark}=Markerinfo{3,q};
        I=find(I)-1;
        if ~isempty(I)
            MarkerinfoOUT{4,NumMark}=TRCdata(:,1+(1:3)+3*(I-1));
        else
            IndsDel=[IndsDel q];
        end
        NumMark=NumMark+1;
    end
end
if ~isempty(IndsDel)
    MarkerinfoOUT(:,IndsDel)=[];
end

stampa=true;
figure,hold all,axis equal
if stampa
    for q=1:size(MarkerinfoOUT,2)
        line(MarkerinfoOUT{4,q}(:,1),MarkerinfoOUT{4,q}(:,2),MarkerinfoOUT{4,q}(:,3),'LineWidth',2)
    end
end
T=TRCdata(:,1);