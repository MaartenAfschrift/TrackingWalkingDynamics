function [] = MuscleFunctions_Right(S)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Path information
SubjPath    = S.OutPath;
RightMusclePath = fullfile(SubjPath,'Right');
%% Get kinematics

%% Load all the information
d=load(fullfile(SubjPath,'Right','MuscleInfo_.mat'));
e=load(fullfile(SubjPath,'Right','muscle_spanning_joint_INFO.mat'));
MuscleInfoR=d.MuscleInfo;
muscle_spanning_joint_INFOR=e.muscle_spanning_joint_INFO;
index_q_locked=[7 13];

% we need to get the number of DOFs => read from Left and Right model ?
Tracking=load(fullfile(SubjPath,'TrackingSim.mat'));

%% Get the symbolic equations
ndof=length(Tracking.output.result.setup.auxdata.NomiCoord_R);
qAll=sym('Q', [1 ndof]);
syms q x0 x1;
ctmat=1;
ctmatdiff=1;
for m=1:length(MuscleInfoR.muscle)
    disp(num2str(m));
    coeff_nr=1;
    nr_coefficients = 0;
    muscle_index=m;
    
    % get the DOF that are spanned by this muscle
    index_dof_crossing = find(muscle_spanning_joint_INFOR(muscle_index,:)==1);
    index_dof_crossing = setdiff(index_dof_crossing, index_q_locked); % Exclude the locked dofs
    nr_dof_crossing=length(index_dof_crossing);
    n_dof=nr_dof_crossing;
    order=MuscleInfoR.muscle(m).order;
    
    % n_art_mat3 function
    q=qAll(index_dof_crossing);
    if nr_dof_crossing==1
        q=[q x0 x0 x0];
    elseif nr_dof_crossing==2
        q=[q x0 x0];
    elseif nr_dof_crossing==3
        q=[q x0];
    end
    for n_q1 = 0:order
        if n_dof<2
            n_q2s = 0;
        else
            n_q2s = 0:order-n_q1;
        end
        for n_q2 = n_q2s
            if n_dof<3
                n_q3s = 0;
            else
                n_q3s = 0:order-n_q1-n_q2;
            end
            for n_q3 = n_q3s
                if n_dof<4
                    n_q4s = 0;
                else
                    n_q4s = 0:order-n_q1-n_q2-n_q3;
                end
                for n_q4 = n_q4s
                    nr_coefficients = nr_coefficients + 1;
                end
            end
        end
    end
    for n_q1 = 0:order
        if n_dof<2
            n_q2s = 0;
        else
            n_q2s = 0:order-n_q1;
        end
        for n_q2 = n_q2s
            if n_dof<3
                n_q3s = 0;
            else
                n_q3s = 0:order-n_q1-n_q2;
            end
            for n_q3 = n_q3s
                if n_dof<4
                    n_q4s = 0;
                else
                    n_q4s = 0:order-n_q1-n_q2-n_q3;
                end
                for n_q4 = n_q4s
                    Index_storeMat(ctmat,:)=[coeff_nr m];
                        mat_temp = q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4;
                        if mat_temp ==0
                            mat(ctmat) = x0;
                        elseif mat_temp == 1
                            mat(ctmat) = x1;
                        else                            
                            mat(ctmat) = mat_temp;
                        end
                        ctmat=ctmat+1;
                        
                    for dof_nr = 1:n_dof                        
                        Index_storeMatDiff(ctmatdiff,:)=[m,coeff_nr,dof_nr];
                        if dof_nr==1
                            out_temp=n_q1*q(:,1).^(n_q1-1).*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4;
                            if out_temp==0
                                diff_mat_q(ctmatdiff)=x0;
                            elseif out_temp==1
                                diff_mat_q(ctmatdiff)=x1;
                            else
                                diff_mat_q(ctmatdiff)=n_q1*q(:,1).^(n_q1-1).*q(:,2).^n_q2.*q(:,3).^n_q3.*q(:,4).^n_q4;
                            end
                            ctmatdiff=ctmatdiff+1;
                        elseif dof_nr==2
                            out_temp=q(:,1).^n_q1.*n_q2.*q(:,2).^(n_q2-1).*q(:,3).^n_q3.*q(:,4).^n_q4;
                            if out_temp==0
                                diff_mat_q(ctmatdiff)=x0;
                            elseif out_temp==1
                                diff_mat_q(ctmatdiff)=x1;
                            else
                                diff_mat_q(ctmatdiff)=q(:,1).^n_q1.*n_q2.*q(:,2).^(n_q2-1).*q(:,3).^n_q3.*q(:,4).^n_q4;
                            end
                            ctmatdiff=ctmatdiff+1;
                        elseif dof_nr==3
                            out_temp=q(:,1).^n_q1.*q(:,2).^n_q2.*n_q3.*q(:,3).^(n_q3-1).*q(:,4).^n_q4;
                            if out_temp==0
                                diff_mat_q(ctmatdiff)=x0;
                            elseif out_temp==1
                                diff_mat_q(ctmatdiff)=x1;
                            else
                                diff_mat_q(ctmatdiff)=q(:,1).^n_q1.*q(:,2).^n_q2.*n_q3.*q(:,3).^(n_q3-1).*q(:,4).^n_q4;
                            end
                            ctmatdiff=ctmatdiff+1;
                        elseif dof_nr==4
                            out_temp=q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*n_q4.*q(:,4).^(n_q4-1);
                            if out_temp==0
                                diff_mat_q(ctmatdiff)=x0;
                            elseif out_temp==1
                                diff_mat_q(ctmatdiff)=x1;
                            else
                                diff_mat_q(ctmatdiff)=q(:,1).^n_q1.*q(:,2).^n_q2.*q(:,3).^n_q3.*n_q4.*q(:,4).^(n_q4-1);
                            end
                            ctmatdiff=ctmatdiff+1;
                        end
                    end
                    coeff_nr = coeff_nr + 1;
                end
            end
        end
    end
end


%% Export the symbolic equations to matlab functions files
%  (This optimizes the equations automatically --> speed-up!)
f = matlabFunction(mat(:),'File',fullfile(RightMusclePath,'Get_Mat_Right.m'));
f3= matlabFunction(diff_mat_q(:),'File',fullfile(RightMusclePath,'Get_Mat_Q_Right.m'));
save(fullfile(RightMusclePath,'MomentArmInfoR.mat'),'f','f3','Index_storeMatDiff','Index_storeMat');

%% Store all the information in structure
SymbolicInfo=load(fullfile(RightMusclePath,'MomentArmInfoR.mat'));
MuscleInfoR.SymbolicInfo=SymbolicInfo;
for m=1:length(MuscleInfoR.muscle)
    muscle_index=m;    
    % get the DOF that are spanned by this muscle
    index_dof_crossing = find(muscle_spanning_joint_INFOR(muscle_index,:)==1);
    index_dof_crossing = setdiff(index_dof_crossing, index_q_locked); % Exclude the locked dofs
    nr_dof_crossing=length(index_dof_crossing);    
    MuscleInfoR.muscle(m).index_dof_crossing = index_dof_crossing;
    MuscleInfoR.muscle(m).nr_dof_crossing = nr_dof_crossing;
end

save(fullfile(RightMusclePath,'MuscleInfoR.mat'),'MuscleInfoR');




end

