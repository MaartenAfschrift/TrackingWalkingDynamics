function [] = GetMusclePolynomials_subtalar(S)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%% Path information
import org.opensim.modeling.*;
ModelPathL=fullfile(S.OutPath,'LeftSideModel.osim');        % this should point to the left and right side model
ModelPathR=fullfile(S.OutPath,'RightSideModel.osim');        % this should point to the left and right side model
DummyMotionPath=fullfile(S.InstallFolder,'Data','FitMuscleFunctions');

PathLeft=fullfile(S.OutPath,'Left');
PathRight=fullfile(S.OutPath,'Right');

if ~isdir(PathLeft);    mkdir(PathLeft);    end;
if ~isdir(PathRight);   mkdir(PathRight);   end;

%% Polynomail fit
max_order = 4;  threshold = 0.002;
if_construct_model=0;

% run for Left model
polynomialFit_LMT_dm_3D(ModelPathL, PathLeft,fullfile(DummyMotionPath,'dummy_motion_LeftSubtalar.mot'), max_order, threshold,...
    S.do_muscle_analysis, S.do_polynomial_fitting, if_construct_model);
% run for Right model
polynomialFit_LMT_dm_3D(ModelPathR, PathRight,fullfile(DummyMotionPath,'dummy_motion_RightSubtalar.mot'), max_order, threshold,...
    S.do_muscle_analysis, S.do_polynomial_fitting, if_construct_model);


end

