function [head,data]=readTRC_PredictiveSimulations(fid)
tmp='';
tmpSTR{1}='';
while ~strcmp(tmpSTR{1},'Frame#')
    tmp=fgets(fid);
    tmpSTR=strread(tmp,'%s');
end
head(1,:)=tmpSTR;
tmp=fgets(fid);
tmpSTR=strread(tmp,'%s');
XYXead=tmpSTR;
tmp=fgets(fid);
while isempty(str2num(tmp))
    tmp=fgets(fid);
end
data=[];
while ischar(tmp)
    data=[data;str2num(tmp)];
    tmp=fgets(fid);
end

frames=data(:,1);
data(:,1)=[];
head(1)=[];
