% conviene creare due modelli con due bacini, l'idea era di fare la divisione del bacino in due parti, ma magari è più semplice dimezzare le
% proprietà inerziali. Così il centro di massa rimane lo stesso, e i calcoli sono più semplici

% Lorenzo Pitto

function [varargout]=invertedModel_reverseJoint(varargin)
NOMEMODELLO=varargin{1};
CARTELLAOUTPUT=varargin{2};
% names of the body for building the new chains
BodiesR{1}={'calcn_r','talus_r','tibia_r','femur_r','pelvis'};
BodiesR{2}={'toes_r'};
BodiesL{1}={'calcn_l','talus_l','tibia_l','femur_l','pelvis','torso'};
BodiesL{2}={'toes_l'};
% calls the function
NOMEOUT_L=[CARTELLAOUTPUT,'\LeftSideModel.osim'];
NOMEOUT_R=[CARTELLAOUTPUT,'\RightSideModel.osim'];
[MarkerInfoL]=createChain(BodiesL,'Model_Left_Side',NOMEMODELLO,NOMEOUT_L);
[MarkerInfoR,MarkerInfoSingle]=createChain(BodiesR,'Model_Right_Side',NOMEMODELLO,NOMEOUT_R);
% output
varargout{1}=NOMEOUT_L;
varargout{2}=NOMEOUT_R;
varargout{3}=[];%MOD_C;
varargout{4}=MarkerInfoL;
varargout{5}=MarkerInfoR;
varargout{6}=MarkerInfoSingle;
end

function [varargout]=createChain(BN,nome,NOMEMODELLO,NOMEOUT)
import org.opensim.modeling.*
oldModello=Model(NOMEMODELLO);
TemplateModelPath=which('TemplateJointToGround.osim');
TMPJOINTMOD=Model(TemplateModelPath);
TMPJOINT=CustomJoint.safeDownCast(TMPJOINTMOD.getBodySet().get(1).getJoint());
OUT=Model();
OUT.setName(nome);
GRND=OUT.getGroundBody();
BS=oldModello.getBodySet();
oldI=ArrayDouble(0.0,9);
oldCoM=Vec3(3,0.0);
oldScFac=Vec3(3,0.0);
for q=1:length(BN)
    for w=1:length(BN{q})
        if q==1
            if w==1
                oldC=BS.get(BN{q}{w});
                oldC.getInertia(oldI);
                oldM=oldC.getMass();
                oldC.getScaleFactors(oldScFac);
                oldC.getMassCenter(oldCoM);
                newC=Body();
                newC.scale(oldScFac);
                newC.setName(BN{q}{w});
                newJ=TMPJOINT;
                newJ.setName([oldC.getName().toCharArray()','_toGround']);
                newJ.set_parent_body(GRND.getName());
            else
                newP=OUT.getBodySet().get(BN{q}{w-1});
                oldP=BS.get(newP.getName());
                newC=Body();
                newC.setName(BN{q}{w});
                oldC=BS.get(newC.getName());
                oldC.getScaleFactors(oldScFac);
                newC.scale(oldScFac);
                if strcmp(BN{q}{w},'torso')
                    oldJ=CustomJoint.safeDownCast(BS.get(BN{q}{w}).getJoint());
                else
                    oldJ=CustomJoint.safeDownCast(oldP.getJoint());
                end
                oldCS=oldJ.getCoordinateSet();
                oldST=oldJ.getSpatialTransform();
                oldOIP=oldJ.get_orientation_in_parent();
                oldOIC=oldJ.get_orientation();
                oldLIP=oldJ.get_location_in_parent();
                oldLIC=oldJ.get_location();
                newJ=CustomJoint();
                newJ.setName(oldJ.getName());
                newJ.set_parent_body(oldP.getName());
                newJ.set_CoordinateSet(oldCS);
                newJ.set_SpatialTransform(oldST);
                oldC.getInertia(oldI);
                oldM=oldC.getMass();
                oldC.getMassCenter(oldCoM);
                if ~strcmp(BN{q}{w},'torso')
                    newJ.set_location(oldLIP);
                    newJ.set_location_in_parent(oldLIC);
                    newJ.set_orientation_in_parent(oldOIC);
                    newJ.set_orientation(oldOIP);
                    newJ.set_reverse(true);
                else
                    newJ.set_location(oldLIC);
                    newJ.set_location_in_parent(oldLIP);
                    newJ.set_orientation_in_parent(oldOIP);
                    newJ.set_orientation(oldOIC);
                end
            end
        else % if q==1
            if w==1
                oldP=BS.get(BN{q}{w});
                newC=Body();
                newC.setName(BN{q}{w});
                oldC=BS.get(newC.getName());
                oldC.getScaleFactors(oldScFac);
                newC.scale(oldScFac);
                oldJ=CustomJoint.safeDownCast(oldP.getJoint);
                oldCS=oldJ.getCoordinateSet();
                oldST=oldJ.getSpatialTransform();
                newJ=CustomJoint();
                newJ.setName(oldJ.getName());
                newJ.set_parent_body(oldJ.getParentBody().getName());
                newJ.set_CoordinateSet(oldCS);
                newJ.set_SpatialTransform(oldST);
                oldOIP=oldJ.get_orientation_in_parent();
                oldOIC=oldJ.get_orientation();
                oldLIP=oldJ.get_location_in_parent();
                oldLIC=oldJ.get_location();
                oldC.getInertia(oldI);
                oldM=oldC.getMass();
                oldC.getMassCenter(oldCoM);
                newJ.set_location(oldLIC);
                newJ.set_location_in_parent(oldLIP);
                newJ.set_orientation_in_parent(oldOIP);
                newJ.set_orientation(oldOIC);
            else
            end
        end
        [oldM,oldI]=splitBody(BN{q}{w},oldM,oldI);
        newC.setMass(oldM);
        newC.setInertia(oldI);
        newC.setMassCenter(oldCoM);
        oldDisp=oldC.getDisplayer();
        for GScnt=0:oldDisp.getGeometrySet().getSize()-1;
            newC.addDisplayGeometry(oldDisp.getGeometryFileName(GScnt));
        end
        newC.setJoint(newJ);
        OUT.addBody(newC);
    end
end
% aggiungi i muscoli
MUSC=oldModello.getMuscles();
nMUSC=MUSC.getSize();
for m=1:nMUSC
    M=MUSC.get(m-1).clone();
    GP=M.getGeometryPath;
    PPS=GP.getPathPointSet();
    nPT=PPS.getSize();
    BodOK=[];
    for p=1:nPT
        BNM=PPS.get(p-1).getBodyName().toCharArray()'; 
        tmpBodOK=false;
        for q=1:length(BN)
            tmpBodOK=tmpBodOK||any(strcmp(BN{q},BNM));
        end
        BodOK=[BodOK tmpBodOK];
    end
    if all(BodOK)
%         M.getName
        OUT.getMuscles().adoptAndAppend(M);
        OUT.addForce(M);
    end
end
% aggiungi i marker e esporta i nomi e le traiettorie dei marker per ciascun modello
MSold=oldModello.getMarkerSet();
MSnew=OUT.getMarkerSet();
NumMark=0;
for m=1:MSold.getSize()
    Mpos=Vec3(0);
    M=MSold.get(m-1).clone();
    M.getOffset(Mpos);
    mbody=M.getBodyName().toCharArray()';
    BodOK=false;
    
    SingleMarkerInfo{1,m}=M.getName().toCharArray()';
    SingleMarkerInfo{2,m}=mbody;
    SingleMarkerInfo{3,m}=[Mpos.get(0);Mpos.get(1);Mpos.get(2)];
    
    for q=1:OUT.getBodySet().getSize()
        BodOK=BodOK || strcmp(OUT.getBodySet().get(q-1).getName().toCharArray()',mbody);
    end
    if BodOK
        Mnew=Marker();
        Mnew.setBodyName(OUT.getBodySet().get(mbody).getName());
        Mnew.setOffset(Mpos);
        Mnew.setName(M.getName());
%         Mnew.connectMarkerToModel(OUT);
        MSnew.set(NumMark,M);
        NumMark=NumMark+1;
        MarkerInfo{1,NumMark}=M.getName().toCharArray()';
        MarkerInfo{2,NumMark}=mbody;
        MarkerInfo{3,NumMark}=[Mpos.get(0);Mpos.get(1);Mpos.get(2)];
    end
end
varargout{1}=MarkerInfo;
varargout{2}=SingleMarkerInfo;
OUT.print(NOMEOUT);
OUT.disownAllComponents();
oldModello.disownAllComponents();
clear OUT oldModello
end

function [varargout]=splitBody(varargin)
if strcmpi(varargin{1},'pelvis')
    for q=2:length(varargin)
        try
            varargout{q-1}=varargin{q}/2;
        catch
            I=varargin{q};
            for w=0:I.size()-1
                I.set(w,I.get(w)/2)
            end
            varargout{q-1}=I;
        end
    end
else
    for q=2:length(varargin)
        varargout{q-1}=varargin{q};
    end
end
end
