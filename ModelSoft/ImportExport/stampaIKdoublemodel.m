function []=stampaIKdoublemodel(NOMEfiLE,DATIvecchi,HEADvecchi,DATInuoviANG,DATInuoviPOS)
namesNEWIK=HEADvecchi;
namesNEWIK{2}='first_body_tilt';
namesNEWIK{3}='first_body_list';
namesNEWIK{4}='first_body_rotation';
namesNEWIK{5}='first_body_tx';
namesNEWIK{6}='first_body_ty';
namesNEWIK{7}='first_body_tz';
dataNEWIK=DATIvecchi;
dataNEWIK(:,2:4)=DATInuoviANG(:,[3 1 2])*180/pi;
dataNEWIK(:,5:7)=DATInuoviPOS;
fid=fopen(NOMEfiLE,'w');
fprintf(fid,'Coordinates\n');
fprintf(fid,'version=1\n');
fprintf(fid,'%s',['nRows=',num2str(size(DATIvecchi,1))]);
fprintf(fid,'\n');
fprintf(fid,'%s',['nColumns=',num2str(size(DATIvecchi,2))]);
fprintf(fid,'\n');
fprintf(fid,'inDegrees=yes\n\n');
fprintf(fid,'Units are S.I. units (second, meters, Newtons, ...)\nAngles are in degrees.\n\n');
fprintf(fid,'endheader\n');
fprintf(fid,'time');
for q=2:length(namesNEWIK)
    fprintf(fid,'\t');
    fprintf(fid,namesNEWIK{q});
end
fprintf(fid,'\n');
for r=1:size(dataNEWIK,1)
   for c=1:size(dataNEWIK,2)
       fprintf(fid,'%f',dataNEWIK(r,c));
       fprintf(fid,'\t');
   end
   fprintf(fid,'\n');
end
fclose(fid);