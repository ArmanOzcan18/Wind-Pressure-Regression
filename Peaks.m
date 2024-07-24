% clear
% clc
% close all
% % % % % % load %time history
% load('D:\FIU\Others\Diana\FIU\Research\Project1\WOW_Data\Diana\204n_InnovationGenerator\Scanivalve\20231027_Case2_RoofWithPV\TranferFuctionCorrected_MATLAB_Data\TF_204n_Case2_RoofwithPV_044_5perc_000deg_Test01.mat')
% V= 22.21; % velocity in m/s  use it for Throttle ratio 44%    (files with 044)
% CP= pDataCorrected./(.5*airDensitySI(76,87,1010.3)*V^2)*6894.76;
% for i = 1: size(CP,2)
% [~,  Cp_rmax, ~,  Cp_rmin] = blue4pressure_modified(CP(:,i), 100, 0.78, 500); %generate the peaks based on EV distribution
% CpMax(1,i) = Cp_rmax;
% CpMin(1,i) = Cp_rmin;
% end
% Cp_mean=mean(CP); %Mean Cp
% Cp_rms=rms(CP); %RMS Cp
% Cp_std=std(CP); %STD Cp


% function [CpMin, Cp_mean, Cp_std] = Peaks(filename)
%      clc
%      close all
%      load(filename)
%      for i = 1:size(CP, 2)
%          [~,~,~, Cp_rmin] = blue4pressure_modified(CP(:,i), 100, 0.78, 500);
%          CpMin(1,i) = Cp_rmin;
%          Cp_mean(1,i) = mean(CP(:,i));
%          Cp_std(1,i) = std(CP(:,i));
%      end
% end

function [CpMin, Cp_mean, Cp_std] = Peaks(filename)
    load(filename)
    num_epoch = 100; % Number of samples
    leng=length(CP)/num_epoch;
    for i = 1:size(CP, 2)
        [~,~,~, Cp_rmin] = blue4pressure_modified(CP(:,i), 100, 0.78, 500);
        CpMin(1,i) = Cp_rmin;
        for k = 1:num_epoch
            startIndex = round((k - 1) * leng) + 1;
            endIndex = round(k * leng);
            mins(k) = min(CP(startIndex:endIndex,i));
        end
        Cp_mean(1,i) = mean(mins);
        Cp_std(1,i) = std(mins);
    end