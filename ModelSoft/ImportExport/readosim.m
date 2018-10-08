%% [colhead data]=readosim(fid)
% legge un file output di osim estarendo i dati e i nomi delle colonne,
% prende in input il fid del file da leggere

%%
function [colhead data]=readosim(fid)
tmp=fgets(fid);
while strcmp(strtrim(tmp),'endheader')==0
    tmp=fgets(fid);
end
colhead=strread(fgets(fid),'%s')';
while isempty(colhead)
    colhead=strread(fgets(fid),'%s')';
end
data=[];
tmp=fgets(fid);
while ischar(tmp)
    data=[data;str2num(tmp)];%#ok
    tmp=fgets(fid);
end