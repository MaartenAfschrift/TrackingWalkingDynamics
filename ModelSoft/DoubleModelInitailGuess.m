function [OUT]=DoubleModelInitailGuess(NOMEMODELLO,IKFILE,GRF)
% //		  0 - nome modello
% //		  1 - dati posizione coordinate
% //		  2 - dati velocità coordinate
% //		  3 - dati accelerazione coordinate
% //		  4 - coordinate in Body dei punti da seguire
% //		  5 - nomi corpi di apparteneza dei punti da seguire
% //		  6 - GRF forza destra
% //		  7 - GRF posizione destra
% //		  8 - GRF forza sinistra
% //		  9 - GRF posizione sinistra
% //		 10 - tempo
% ---- in questo caso le forze di input non le metterei nei piedi ma nel punto in comune tra i due bacini, che poi devono essere uguali per
% i due modelli, le GRF le calcolo applicate nel punto di clcn2grnd e le impongo uguali al risultato di ID
% nella funzione di costo minimizzare il jerk, o le accelerazioni
% per gli istanti iniziali e finali costringere magari la x del bacino, e lasciare con un po' di margine le altre coordinate, visto che
% tanto per ora è solo sul single stance
STAMPA=true;
import org.opensim.modeling.*
% process input
T=GRF.T;
GRFtimeind=strcmp(GRF.head,'time');
GRFtimeOrig=GRF.data(:,GRFtimeind);
GRFdata=GRF.data(:,~GRFtimeind);
GRFhead=GRF.head(~GRFtimeind);
GRFdata=interp1(GRFtimeOrig,GRFdata,T,'spline');
% lato
lato=strrep(IKFILE,'.mot','');
lato=lato(end);
% read the IK file for one model (L)
fid=fopen(IKFILE,'r');
[IKhead,IKdata]=readosim(fid);
fclose(fid);
% open the model
modello=Model(NOMEMODELLO);
% imports the coordinates from the model and find which one are the first-body-to-ground translations
CS=modello.getCoordinateSet();
IndDeg=true(1,CS.getSize());
IndCS=zeros(1,CS.getSize());
for q=1:CS.getSize()
    NomiCS{q}=CS.get(q-1).getName().toCharArray()';
    IndCS(q)=find(strcmpi(IKhead,NomiCS{q}));
    if strcmpi(NomiCS{q}(end-2:end),'_tx')||strcmpi(NomiCS{q}(end-2:end),'_ty')||strcmpi(NomiCS{q}(end-2:end),'_tz')
        IndDeg(q)=false;
    end
end
% interpolate and generate spline from coordinates
Qhead=IKhead(IndCS);
Q=IKdata(:,IndCS);
Q=interp1(IKdata(:,1),Q,T,'spline');
Q(:,IndDeg)=Q(:,IndDeg)/180*pi;
ppQ=spline(T,Q');
ppU=fnder(ppQ,1);
ppA=fnder(ppQ,2);
ppJ=fnder(ppQ,3);
figure,subplot(3,1,1),plot(T,ppval(ppQ,T)),subplot(3,1,2),plot(T,ppval(ppU,T)),subplot(3,1,3),plot(T,ppval(ppA,T))
% derives spline for velocities and accelerations
Qpos=ppval(ppQ,T)';
Qvel=ppval(ppU,T)';
Qacc=ppval(ppA,T)';
PT_bodynames={'pelvis'};
PT_posB=[0 0 0]';
% computes ID with no additional forces and 
% GRF_F_D=zeros(length(T),3);
% GRF_P_D=zeros(length(T),3);
% GRF_F_S=zeros(length(T),3);
% GRF_P_S=zeros(length(T),3);
ZERO_F=zeros(length(T),3);
[T_cpp,~]=trovaIDePATH_MD(NOMEMODELLO,Qpos,Qvel,Qacc,PT_posB,PT_bodynames,ZERO_F,ZERO_F,ZERO_F,ZERO_F);
% PT_posG_cpp=reshape(PT_posG_cpp,length(T),3);
% interpolates the GRF data (imported before)
% GRFmat=[[GRFdata(:,1);GRFdata(:,2);GRFdata(:,3)]/100;GRFdata(:,4);GRFdata(:,5);GRFdata(:,6);GRFdata(:,13);GRFdata(:,14);GRFdata(:,15)];


% stem name force
if strcmpi(lato,'r')
    stem='1_ground_';
elseif strcmpi(lato,'l')
    stem='ground_';
end 
% B2G txyz force from GRF file
IndGRFforce=strncmpi(GRFhead,[stem,'force_v'],length([stem,'force_v']));
ForzaRexp=GRFdata(:,IndGRFforce);
% B2G txyz force from ID
ForzaRcpp=T_cpp(:,4:6);
% center of pressure from GRF file, world frame
IndGRFpoint=strncmpi(GRFhead,[stem,'force_p'],length([stem,'force_p']));
PuntoRexp=GRFdata(:,IndGRFpoint);
% position of first body 000, location of B2G
% PuntoRcpp=zeros(size(PuntoRexp));
PuntoRcppinG=Qpos(:,4:6);
% torque from GRF file, world frame
IndGRFmoment=strncmpi(GRFhead,[stem,'torque'],length([stem,'torque']));
Momentoexp=GRFdata(:,IndGRFmoment);
% B2G torque, from ID
% Momentocpp=T_cpp(:,1:3);
% difference from the force application point from GRF and B000
DeltaPunto=PuntoRcppinG-PuntoRexp;
if STAMPA
%     figure
%     subplot(2,1,1)
%     plot(PuntoRcppinG),hold all,plot(PuntoRexp,'DisplayName','PuntoRexp')
%     subplot(2,1,2)
%     plot(DeltaPunto)
%     figure
%     subplot(2,1,1)
%     plot(MomentoNuovo)
%     figure
%     subplot(2,2,1),hold all,xlim([TxGRF([1 end])])
%     plot(TxGRF,GRFdata(:,1),'color','r','LineStyle','-','LineWidth',2)
%     plot(TxGRF,T_cpp(:,4),'color','r','LineStyle','--','LineWidth',2)
%     plot(TxGRF,GRFdata(:,2),'color','g','LineStyle','-','LineWidth',2)
%     plot(TxGRF,T_cpp(:,5),'color','g','LineStyle','--','LineWidth',2)
%     plot(TxGRF,GRFdata(:,3),'color','b','LineStyle','-','LineWidth',2)
%     plot(TxGRF,T_cpp(:,6),'color','b','LineStyle','--','LineWidth',2)
%     subplot(2,2,3),hold all,xlim([TxGRF([1 end])])
%     plot(TxGRF,GRFdata(:,1)-T_cpp(:,4),'color','r','LineStyle','-','LineWidth',2)
%     plot(TxGRF,GRFdata(:,2)-T_cpp(:,5),'color','g','LineStyle','-','LineWidth',2)
%     plot(TxGRF,GRFdata(:,3)-T_cpp(:,6),'color','b','LineStyle','-','LineWidth',2)
%     subplot(2,2,2),hold all,xlim([TxGRF([1 end])])
%     plot(TxGRF,GRFdata(:,13),'color','r','LineStyle','-','LineWidth',2)
%     plot(TxGRF,T_cpp(:,1),'color','r','LineStyle','--','LineWidth',2)
%     plot(TxGRF,GRFdata(:,14),'color','g','LineStyle','-','LineWidth',2)
%     plot(TxGRF,T_cpp(:,2),'color','g','LineStyle','--','LineWidth',2)
%     plot(TxGRF,GRFdata(:,15),'color','b','LineStyle','-','LineWidth',2)
%     plot(TxGRF,T_cpp(:,3),'color','b','LineStyle','--','LineWidth',2)
%     plot(TxGRF,MomentoNuovo,'color','m','LineStyle','-','LineWidth',2)
%     subplot(2,2,4),hold all,xlim([TxGRF([1 end])])
%     plot(TxGRF,GRFdata(:,13)-T_cpp(:,1),'color','r','LineStyle','-','LineWidth',2)
%     plot(TxGRF,GRFdata(:,14)-T_cpp(:,2),'color','g','LineStyle','-','LineWidth',2)
%     plot(TxGRF,GRFdata(:,15)-T_cpp(:,3),'color','b','LineStyle','-','LineWidth',2)
end
% computes the moment after the translation of GRF to B000
MomentoNuovo=Momentoexp+cross(ForzaRexp,DeltaPunto);
%% define the forces acting on the pelvis (as for initial guess)
% generates the reference systems to compute the projections of the moments in the ground frame, to check if the ID function works fine
% Mx from I
% My from Rx
% Mz from Rx*Ry
for q=1:length(T)
    MxID_X(:,:,q)=Euler2mat([0,0,0],'xyz');
    MxID_Y(:,:,q)=Euler2mat([Qpos(q,1),0,0],'xyz');
    MxID_Z(:,:,q)=Euler2mat([Qpos(q,1),Qpos(q,2),0],'xyz');
end
% the force that has to be applied to the pelvis is the difference between the measured GRF and the compted force from ID
ResiduiPelvisFor=-GRFdata(:,IndGRFforce)+ForzaRcpp;
% B000 forces applying the forces at the pelvis
% position in world frame of Pelvis000
% T_cpp_2 here is in non orthogonal axis
[T_cpp_2,PT_posG_cpp_2]=trovaIDePATH_MD(NOMEMODELLO,Qpos,Qvel,Qacc,PT_posB,PT_bodynames,ResiduiPelvisFor,ZERO_F,ZERO_F,ZERO_F);
PT_posG_cpp_2=reshape(PT_posG_cpp_2,length(T),3);
% moment arm as the difference from P000 and B000, world frame
Braccio=PT_posG_cpp_2-Qpos(:,4:6);
% torque from the application of the forces in P000
% MomentoDaBraccio=cross(ResiduiPelvisFor,Braccio);
% computes the residuals to be applied to P as the difference between the moment after the application of the force in P000 and the moment
% obtained translating the GRF to B000
% transforms the moment from ID into the world frame
for q=1:length(T)
    MAT_MOM(:,:,q)=[MxID_X(:,1,q)'
                    MxID_Y(:,2,q)'
                    MxID_Z(:,3,q)'];
    T_ccp_2_inG(q,:)=MAT_MOM(:,:,q)\T_cpp_2(q,1:3)';
end
ResiduiPelvisMom=-MomentoNuovo+T_ccp_2_inG;
% defines 3 points to be trached on P to define position and orientation
PT_posB_2=[0 0 0
           1 0 0 
           0 1 0]';
PT_bodynames_2={'pelvis','pelvis','pelvis'};
% ID appling force and torque at P and trachiong the 3 points
[T_cpp_3,PT_posG_cpp_2]=...
    trovaIDePATH_MD(NOMEMODELLO,Qpos,Qvel,Qacc,PT_posB_2,PT_bodynames_2,ResiduiPelvisFor,ZERO_F,ResiduiPelvisMom,ZERO_F);
% extract the orientation and position of P
for q=1:size(PT_posG_cpp_2,2)
    PT_posG_cpp_2RES(:,:,q)=reshape(PT_posG_cpp_2(:,q),length(T),3);
end
for q=1:length(T)
    MATpelv(:,1,q)=(PT_posG_cpp_2RES(q,:,2)-PT_posG_cpp_2RES(q,:,1))/norm(PT_posG_cpp_2RES(q,:,2)-PT_posG_cpp_2RES(q,:,1));
    MATpelv(:,2,q)=(PT_posG_cpp_2RES(q,:,3)-PT_posG_cpp_2RES(q,:,1))/norm(PT_posG_cpp_2RES(q,:,3)-PT_posG_cpp_2RES(q,:,1));
    MATpelv(:,3,q)=(cross(MATpelv(:,1,q),MATpelv(:,2,q)))/norm(cross(MATpelv(:,1,q),MATpelv(:,2,q)));
%     ROTpelv(q,:)=mat2Euler(MATpelv(:,:,q));
end
% transforms the moment from ID into the world frame
for q=1:length(T)
    T_ccp_3_inG(q,:)=MAT_MOM(:,:,q)\T_cpp_3(q,1:3)';
end
for q=0:modello.getCoordinateSet().getSize()-1
    modello.getCoordinateSet().get(q)
end
for q=0:modello.getJointSet().getSize()-1
    modello.getJointSet().get(q)
end
if STAMPA
    figure('name','vedi nello spazio i punti di applicazione e vettori di F e M'),hold all,axis equal
    scatter3(PT_posG_cpp_2RES(:,1),PT_posG_cpp_2RES(:,2),PT_posG_cpp_2RES(:,3),'b.')
    scatter3(Qpos(:,4),Qpos(:,5),Qpos(:,6),'r.')
    for q=1:size(Braccio,1)
        line([PT_posG_cpp_2RES(q,1,1) PT_posG_cpp_2RES(q,1,1)-Braccio(q,1)/norm(Braccio(q,:))/4],...
            [PT_posG_cpp_2RES(q,2,1) PT_posG_cpp_2RES(q,2,1)-Braccio(q,2)/norm(Braccio(q,:))/4],...
            [PT_posG_cpp_2RES(q,3,1) PT_posG_cpp_2RES(q,3,1)-Braccio(q,3)/norm(Braccio(q,:))/4])
        line([PT_posG_cpp_2RES(q,1,1) PT_posG_cpp_2RES(q,1,1)+T_cpp_3(q,4)/norm(T_cpp_3(q,4:6))/4],...
            [PT_posG_cpp_2RES(q,2,1) PT_posG_cpp_2RES(q,2,1)+T_cpp_3(q,5)/norm(T_cpp_3(q,4:6))/4],...
            [PT_posG_cpp_2RES(q,3,1) PT_posG_cpp_2RES(q,3,1)+T_cpp_3(q,6)/norm(T_cpp_3(q,4:6))/4],'color','g')
    end
    figure('name','compara con F e M da GRF')
    subplot(1,2,1)
    plot(ForzaRexp,'color','m','LineWidth',2,'LineStyle','-'),hold all
    plot(T_cpp(:,4:6),'DisplayName','T_cpp')
    scatter(1:length(T_cpp_3),T_cpp_3(:,4),'k')
    scatter(1:length(T_cpp_3),T_cpp_3(:,5),'k')
    scatter(1:length(T_cpp_3),T_cpp_3(:,6),'k')
    subplot(1,2,2),hold all
    plot(MomentoNuovo,'color','m','LineWidth',2,'LineStyle','-')
    scatter(1:length(T_ccp_3_inG),T_ccp_3_inG(:,1),'k')
    scatter(1:length(T_ccp_3_inG),T_ccp_3_inG(:,2),'k')
    scatter(1:length(T_ccp_3_inG),T_ccp_3_inG(:,3),'k')
end
%%
OUT.Fp=ResiduiPelvisFor;
OUT.Mp=ResiduiPelvisMom;
OUT.MATpelv=MATpelv;
OUT.T_ccp_inG_FINAL=T_ccp_3_inG;
OUT.F_ccp_inG_FINAL=T_cpp_3(:,4:6);
OUT.ppQ=ppQ;
OUT.ppU=ppU;
OUT.ppA=ppA;
OUT.ppJ=ppJ;
OUT.Qhead=Qhead;
OUT.Pelv000=reshape(PT_posG_cpp_2(:,1),length(PT_posG_cpp_2(:,1))/3,3);