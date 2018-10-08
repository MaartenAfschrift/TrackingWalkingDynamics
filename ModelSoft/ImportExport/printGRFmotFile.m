function []=printGRFmotFile(Data,Nomi,NOMEfiLE)
[~,q,w]=fileparts(NOMEfiLE);
fid=fopen(NOMEfiLE,'w');
fprintf(fid,[q,w,'\n']);
fprintf(fid,'version=1\n');
fprintf(fid,'%s',['nRows=',num2str(size(Data,1))]);
fprintf(fid,'\n');
fprintf(fid,'%s',['nColumns=',num2str(size(Data,2))]);
fprintf(fid,'\n');
fprintf(fid,'inDegrees=yes\n\n');
fprintf(fid,'endheader\n');
fprintf(fid,'time');
for q=2:length(Nomi)
    fprintf(fid,'\t');
    fprintf(fid,Nomi{q});
end
fprintf(fid,'\n');
for r=1:size(Data,1)
    for c=1:size(Data,2)
        fprintf(fid,'%f',Data(r,c));
        fprintf(fid,'\t');
    end
    fprintf(fid,'\n');
end
fclose(fid);