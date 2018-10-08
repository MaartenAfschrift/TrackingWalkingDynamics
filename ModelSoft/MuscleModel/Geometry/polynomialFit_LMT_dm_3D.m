function polynomialFit_LMT_dm_3D(model_file_name, root_folder, dummy_motion_fileName, max_order, threshold, do_muscle_analysis, do_polynomial_fitting, if_construct_model)


if_reuse_dummy_motion = 1;
check_implementation = 0;

model_sel=model_file_name;
prefix = '';

%% Running the muscle analysis for a dummy motion spanning the range motion for the degrees of motion
% get the opensim model
if do_muscle_analysis
    if if_reuse_dummy_motion
        dummy_motion = importdata(dummy_motion_fileName);
        Angles = dummy_motion.data(:,2:end);
        time = dummy_motion.data(:,1);
    end
    
    
    % do a muscle analysis
    % Old implementation:
%     output_path=fullfile(root_folder, 'MuscleAnalysis');mkdir(output_path);
%     OpenSim_Muscle_Analysis(fullfile(pwd, dummy_motion_fileName),model_sel,output_path,[time(1) time(end)])

    % New implementation:
    output_path=fullfile(root_folder, 'MuscleAnalysis');mkdir(output_path);

    MA_setup_fileName = which('settings_Muscle_analysis.xml');
    MA_setup = xml_read(MA_setup_fileName);
    MA_setup.AnalyzeTool.model_file = model_file_name; %fullfile(pwd, model_file_name);
    MA_setup.AnalyzeTool.results_directory = output_path;
    MA_setup.AnalyzeTool.initial_time = time(1);
    MA_setup.AnalyzeTool.final_time = time(end);
    MA_setup.AnalyzeTool.ATTRIBUTE.name = 'dummy_motion';
    MA_setup.AnalyzeTool.coordinates_file = fullfile(dummy_motion_fileName);
       
    xml_write(MA_setup_fileName, MA_setup, 'OpenSimDocument');
    
    dos(['C:\OpenSim33\bin\analyze -S ', MA_setup_fileName]);
    
    % get the results from the analysis
    if ~isfield(dummy_motion,'colheaders'); dummy_motion.colheaders=strsplit(dummy_motion.textdata{end}); end;
    dof_names = strtok(dummy_motion.colheaders(2:end));
    dof_nr_index = 1;
    dof_nr_index_vector = [];
    for dof_nr = 1:length(dof_names)
        if exist([output_path, '/dummy_motion_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto'], 'file')
            dM=importdata(fullfile(output_path,['dummy_motion_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
            MuscleData.dM(:,:,dof_nr_index) = dM.data(:,2:end);
            dof_nr_index_vector = [dof_nr_index_vector, dof_nr];
            dof_nr_index = dof_nr_index + 1;
        else
            disp('problem');
        end
    end
    Length=importdata(fullfile(output_path,'dummy_motion_MuscleAnalysis_Length.sto'));
    
    MuscleData.Length=Length.data(:,2:end);
    MuscleData.Angles=Angles(:,dof_nr_index_vector);
    MuscleData.dof_names = dof_names(dof_nr_index_vector);
    MuscleData.muscle_names = Length.colheaders(2:end);
    
    save([root_folder '\MuscleData_', prefix, '.mat'], 'MuscleData');
end
if do_polynomial_fitting    
    %% Construct the polynomials for the moment arms and muscle length
    load([root_folder  '\MuscleData_', prefix, '.mat'])
    
    muscle_sel=[];
    
    nr_muscles = length(MuscleData.muscle_names);
    
    for m_nr = 1:nr_muscles
%         if strcmp(MuscleData.muscle_names{m_nr}(end-1:end), '_l')
            muscle_sel = [muscle_sel m_nr];
%         end
    end
    
    muscle_spanning_joint_INFO = squeeze(sum(MuscleData.dM, 1));
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO<=0.0001 & muscle_spanning_joint_INFO>=-0.0001) = 0;
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO~=0) = 1;
    
    save([root_folder '\muscle_spanning_joint_INFO.mat'],'muscle_spanning_joint_INFO');
    
    q_all = MuscleData.Angles/180*pi;
    index_q_locked = find(mean(q_all)==q_all(1,:));
    
%     max_order = 6; %3; %6;
%     threshold = 0.003; % 3mm
    nr_samples = length(q_all(:,1));
    
    LMT_all_error = zeros(length(muscle_sel), 1);
    DM_all_error = zeros(length(muscle_sel), length(q_all(1,:)));
    order_all = zeros(length(muscle_sel), 1);
    
    for m_nr=1:length(muscle_sel);
        muscle_index = muscle_sel(m_nr);
        
        index_dof_crossing = find(muscle_spanning_joint_INFO(muscle_index,:)==1);
        index_dof_crossing = setdiff(index_dof_crossing, index_q_locked); % Exclude the locked dofs
        nr_dof_crossing = length(index_dof_crossing);
        
        LMT = MuscleData.Length(:,muscle_index);
        dM = zeros(nr_samples, nr_dof_crossing);
        dM_recon = dM;
        for dof_nr = 1:nr_dof_crossing
            dM(:,dof_nr) = MuscleData.dM(:,muscle_index,index_dof_crossing(dof_nr));
        end
        
        criterion_full_filled = 0;
        order = 2;
        while criterion_full_filled==0
            try
            [mat,diff_mat_q] = n_art_mat_3(q_all(:,index_dof_crossing), order);
            catch
                disp('error');
            end
            nr_coeffs = length(mat(1,:));
            
            %         LMT = MuscleData.Length(:,muscle_index);
            %         dM = zeros(nr_samples, nr_dof_crossing);
            diff_mat_q_all = zeros(nr_samples*nr_dof_crossing, nr_coeffs);
            for dof_nr = 1:nr_dof_crossing
                diff_mat_q_all(nr_samples*(dof_nr-1)+1:nr_samples*dof_nr,:) = -squeeze(diff_mat_q(:,:,dof_nr));
            end
            
            coeff=[mat ; diff_mat_q_all]\[LMT; dM(:)];
            dM_recon = zeros(nr_samples, nr_dof_crossing);
            for dof_nr = 1:nr_dof_crossing
                dM_recon(:,dof_nr) = (-squeeze(diff_mat_q(:,:,dof_nr)))*coeff;
            end
            LMT_recon=mat*coeff;
            
            LMT_error_rms = sqrt(mean((LMT - LMT_recon).^2));
            dm_error_rms = sqrt(mean((dM - dM_recon).^2));
            
            criterion_full_filled = LMT_error_rms<=threshold & max(dm_error_rms)<=threshold;
            if order==max_order
                criterion_full_filled = 1;
            end
            if criterion_full_filled==0;
                order = order+1;
            end
        end
        
        MuscleInfo.muscle(m_nr).DOF = MuscleData.dof_names(index_dof_crossing);
        MuscleInfo.muscle(m_nr).m_name = MuscleData.muscle_names{muscle_index};
        MuscleInfo.muscle(m_nr).coeff = coeff;
        MuscleInfo.muscle(m_nr).order = order;
        MuscleInfo.muscle(m_nr).LMT_error_rms = LMT_error_rms;
        MuscleInfo.muscle(m_nr).dm_error_rms = dm_error_rms;
        
        LMT_all_error(m_nr) = LMT_error_rms;
        DM_all_error(m_nr, index_dof_crossing) = dm_error_rms;
        order_all(m_nr) = order;
        
%         if strcmp(MuscleInfo.muscle(m_nr).m_name, 'tib_ant_l') || strcmp(MuscleInfo.muscle(m_nr).m_name, 'tib_post_l'); %m_nr==1023 %10
%             dof_nr_index = 1;
%             figure;
%             hold on;
%             plot(q_all(:,index_dof_crossing(dof_nr_index)), dM(:,dof_nr_index), '*b')
%             plot(q_all(:,index_dof_crossing(dof_nr_index)), dM_recon(:,dof_nr_index), '*r')
%             
%             figure;
%             hold on;
%             plot(q_all(:,index_dof_crossing(dof_nr_index)), LMT, '*b')
%             plot(q_all(:,index_dof_crossing(dof_nr_index)), LMT_recon, '*r')
%         end       
        
        for dof_nr_index = 1:length(index_dof_crossing)
            if index_dof_crossing(dof_nr_index)==3
%                 dof_nr_index = 3;
                figure;
                hold on;
                plot(q_all(:,index_dof_crossing(dof_nr_index)), dM(:,dof_nr_index), '*b')
                plot(q_all(:,index_dof_crossing(dof_nr_index)), dM_recon(:,dof_nr_index), '*r')
                suptitle(['Moment arm of muscle ' MuscleData.muscle_names(m_nr)])
            end
        end
    end
    
    figure('Position', [610 428 1307 525]);
%     subplot(2,1,1)
    hold on;
    plot(LMT_all_error)
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold, threshold], 'r', 'linewidth', 2)
    suptitle('RMS error on the approximated muscle-tendon length')
    ylabel('RMS error (m)')
%     set(gca, 'XTickLabel', [])
%     xticklabel_rotate(1:length(muscle_sel),90,MuscleData.muscle_names(muscle_sel))
    
    figure('Position', [610 428 1307 525]);
%     subplot(2,1,2)
    hold on;
    plot(max(DM_all_error, [], 2))
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold, threshold], 'r', 'linewidth', 2)
    suptitle('maximal RMS error on the approximated muscle moment arm')
    ylabel('RMS error (m)')
%     xticklabel_rotate(1:length(muscle_sel),90,MuscleData.muscle_names(muscle_sel))
%     title('maximal RMS error on the approximated muscle moment arm')

    figure('Position', [610 428 1307 525]);
    hold on;
    plot(order_all)
    ylim([0 max_order+1])
    xlimits = get(gca, 'XLim');
    plot(xlimits, [max_order, max_order], 'r', 'linewidth', 2)
    suptitle('Order of the polynomial approximation')
    ylabel('Order')
%     xticklabel_rotate(1:length(muscle_sel),90,MuscleData.muscle_names(muscle_sel))
    
    save([root_folder '\MuscleInfo_', prefix, '.mat'], 'MuscleInfo');
end

if if_construct_model
    load(['MuscleInfo_', prefix, '.mat'])
        
    modelXML = xml_readOSIM(model_sel);
    
    if isfield(modelXML.Model.ForceSet.objects, 'Thelen2003Muscle')
        field_name = 'Thelen2003Muscle';
    elseif isfield(modelXML.Model.ForceSet.objects, 'Millard2012EquilibriumMuscle')
        field_name = 'Millard2012EquilibriumMuscle';
    end
    
    for m_nr = 1:length(modelXML.Model.ForceSet.objects.(field_name))
        model_muscle_name = modelXML.Model.ForceSet.objects.(field_name)(m_nr).ATTRIBUTE.name;
        for i=1:length(MuscleInfo.muscle)
            m_name = MuscleInfo.muscle(i).m_name;
            if strcmp(model_muscle_name(1:end-1), m_name(1:end-1))
                side_string = model_muscle_name(end);
                if strcmp(side_string, 'r') 
                    side_string2 = 'l';
                else side_string2 = 'r';
                end
                dof_names = MuscleInfo.muscle(i).DOF;
                dof_name_string = [];
                for k=1:length(dof_names)
                    if strcmp(dof_names{k}(end), side_string2)
                        dof_names{k}(end) = side_string;
                    end
                    dof_name_string = [dof_name_string, ' ', dof_names{k}];
                end
                if isempty(strfind(dof_names{k}, 'lumbar'))==0 && strcmp(side_string, 'r') % Only do it for the right leg
                    [mat_sign,~] = n_art_mat_3([1 -1 -1], MuscleInfo.muscle(i).order); % The sign of the coefficient for lumbar bending and lumbar rotation needs to be reversed
                    coeff = mat_sign'.*MuscleInfo.muscle(i).coeff;
                else
                    coeff = MuscleInfo.muscle(i).coeff;
                end
                modelXML.Model.ForceSet.objects.(field_name)(m_nr).GeometryPath.coefficients_polynomial=coeff;
                modelXML.Model.ForceSet.objects.(field_name)(m_nr).GeometryPath.polynomial_DOF1='3D';
                modelXML.Model.ForceSet.objects.(field_name)(m_nr).GeometryPath.polynomial_DOFs=dof_name_string;
                modelXML.Model.ForceSet.objects.(field_name)(m_nr).GeometryPath.polynomial_order=MuscleInfo.muscle(i).order;
                break;
            end
        end
    end
    
    output_model_file_name = [model_file_name(1:end-5), '_poly.osim'];
    xml_writeOSIM(output_model_file_name, modelXML, 'OpenSimDocument')
end

% if check_implementation
%     % do a muscle analysis
%     dummy_motion = importdata('3D_results/dummy_motion.mot');
%     time = dummy_motion.data(:,1);
%     output_path=fullfile(pwd,'3D_results', 'MuscleAnalysis_poly');mkdir(output_path);
% %     OpenSim_Muscle_Analysis(fullfile(pwd, '3D_results','dummy_motion.mot'),[model_file_name(1:end-5), '_poly.osim'],output_path,[time(1) time(end)])
% 
% %     data_folder = 'C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\MuscleSmoothing\MuscleLengthPolynomials\3D_results';
% %     data_KS = importdata('C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\Data_overground_3D\KS\KS_motion_07.mot');
% %     dof_names = data_KS.colheaders(2:end);
%     
% %     left_leg_DOFs = find(~cellfun('isempty', strfind(dof_names,'_l')));
% %     pelvis_DOFs = find(~cellfun('isempty', strfind(dof_names,'pelvis')));
% %     lumbar_DOFs = find(~cellfun('isempty', strfind(dof_names,'lumbar')));
% %     all_DOFS = [left_leg_DOFs, lumbar_DOFs];
% %     all_DOFS = unique(all_DOFS);
% %     all_DOFS = setdiff(all_DOFS, pelvis_DOFs); % exclude the pelvis DOFS
% %     dof_names = dof_names(all_DOFS);
%     
% %     for dof_nr = 1:length(dof_names)
% %         dM=importdata(fullfile(data_folder,['MuscleAnalysis\dummy_motion_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
% %         MuscleData_normal.dM(:,:,dof_nr) = dM.data(:,2:end);
% %     end
% %     Length=importdata(fullfile(data_folder,'MuscleAnalysis\dummy_motion_MuscleAnalysis_Length.sto'));
% %     MuscleData_normal.LMT = Length.data(:,2:end);
% %     
% %     for dof_nr = 1:length(dof_names)
% %         dM=importdata(fullfile(data_folder,['MuscleAnalysis_poly\dummy_motion_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
% %         MuscleData_poly.dM(:,:,dof_nr) = dM.data(:,2:end);
% %     end
% %     Length=importdata(fullfile(data_folder,'MuscleAnalysis_poly\dummy_motion_MuscleAnalysis_Length.sto'));
% %     MuscleData_poly.LMT = Length.data(:,2:end);
%     
% %     data_folder = 'C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\Data_overground_3D\TestMuscleAnalysisPoly\';
% %     data_KS = importdata('C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\Data_overground_3D\KS\KS_motion_07.mot');
% %     dof_names = data_KS.colheaders(2:end);
% 
%     data_folder = 'C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\DataAmputee\MuscleAnalysis\';
%     data_KS = importdata('C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\DataAmputee\KS_contact\KS_contact_lockedStump_HC\KS_motion_contact_HC.mot');
%     dof_names = data_KS.colheaders(2:end);
%     
%     for dof_nr = 1:length(dof_names)
%         dM=importdata(fullfile(data_folder,['Results_normal\p2-scaled_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
%         MuscleData_normal.dM(:,:,dof_nr) = dM.data(:,2:end);
%     end
%     Length=importdata(fullfile(data_folder,'Results_normal\p2-scaled_MuscleAnalysis_Length.sto'));
%     MuscleData_normal.LMT = Length.data(:,2:end);
%     
%     for dof_nr = 1:length(dof_names)
%         dM=importdata(fullfile(data_folder,['Results_poly\p2-scaled_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
%         MuscleData_poly.dM(:,:,dof_nr) = dM.data(:,2:end);
%     end
%     Length=importdata(fullfile(data_folder,'Results_poly\p2-scaled_MuscleAnalysis_Length.sto'));
%     MuscleData_poly.LMT = Length.data(:,2:end);
% 
%     rms_LMT = sqrt(mean((MuscleData_poly.LMT - MuscleData_normal.LMT).^2));
%     rms_dM = squeeze(sqrt(mean((MuscleData_poly.dM-MuscleData_normal.dM).^2)));
%     
%     figure;
%     plot(rms_LMT)
%     suptitle('RMS error on the approximated muscle-tendon length')
%     ylabel('RMS error (m)')
%     xticklabel_rotate(1:length(Length.data(1,2:end)),90,Length.colheaders(2:end))
%     
%     figure;
%     plot(max(rms_dM, [], 2))
%     suptitle('maximal RMS error on the approximated muscle moment arm')
%     ylabel('RMS error (m)')
%     xticklabel_rotate(1:length(Length.data(1,2:end)),90,Length.colheaders(2:end))
%     
%     right_muscles = find(~cellfun('isempty', strfind(Length.colheaders(2:end),'_r')));
%     left_muscles = setdiff(1:nr_muscles, right_muscles);
%     
%     figure;
%     plot(rms_LMT(right_muscles))
%     suptitle('Right')
%     
%     figure;
%     plot(rms_LMT(left_muscles))
%     suptitle('Left')
%     
% end

% setenv('PATH', old_PATH);

end
