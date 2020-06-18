function [] = getkivalues_02(projdir, subject)
% function to calculate Ki values for one subject

% This script is divided into several parts:
% 1 - loading striatal and cerebelum masks, using the masks to get time
% activity curves (TACs) for each striatal and cerebellum regions from each
% frame of the PET scan
% 2 - I'm doing the analysis on the region-wise mean TAC instead of
% voxelwise, as this is what Howes et al have done. Doing it voxelwise is
% also an option, I've set up both analysis here using the WhichMode
% switch
% 3 - get the Blood Input Function (BIF) by using the cerebellum as a
% reference region. This uses the cumtrapz function.
% 4 -the first frame of the cumulative BIF( CBIF) output is masked because the value is 0
% 5 - Area Under Curve (AUC) is divided by BIF; so I think this is to
% estimate concentration of tracer in the reference region (aka cerebellum)
% 6 - Clearance rate in striatal region is estimated by dividing the TACS
% in each region by the cerebellum BIF

Mode = {'voxelwise','regionwise'};
WhichMode = Mode{2};

% create duration of frames; there are 26 frames
dur = [.5 1 1 1 1 2 2 2 3 3 3 5 * ones(1,15)];

frames = [1:26];
frames = reshape(frames,26,1);

framesdir = [projdir,'datadir/PET_files/',subject,'/t1-space/'];

cd(framesdir)

frames_files=dir(framesdir);
frames_files(1:2) = [];


roidir = [projdir,'datadir/derivatives/fmriprep_resample2mni_2mm/fmriprep/',subject,'/rois/'];

% load cerebellum mask
cer_mask = niftiread([roidir,'Cerebellum-MNIflirt-maxprob-thr50-2mm_t1.nii']);
cer_mask_index = find(cer_mask==1);

% load striatal mask in subject's T1 space
% associative striatum
str_assoc_mask = niftiread([roidir,'Lk32_space-T1.nii']);
str_assoc_mask_index = find(str_assoc_mask==1);
%str_assoc = reshape(str_assoc,[],length(str_assoc));

% ventral striatum
str_vent_mask = niftiread([roidir,'Lk33_space-T1.nii']);
str_vent_mask_index = find(str_vent_mask==1);
% str_vent = reshape(str_assoc,[],length(str_assoc));

% sensorimotor striatum
str_sens_mask = niftiread([roidir,'Lk31_space-T1.nii']);
str_sens_mask_index = find(str_sens_mask == 1);
% str_sens = reshape(str_sens,[],length(str_sens));

% whole striatum
str_whole_mask = niftiread([roidir,'Lk3_space-T1.nii']);
str_whole_mask_index = find(str_whole_mask == 1);
% str_sens = reshape(str_sens,[],length(str_sens));

tacs_str_assoc = [];
tacs_str_vent = [];
tacs_str_sens = [];
tacs_cer = [];

tacs_str_whole = [];

for i = 1:26
    frame = niftiread([framesdir,frames_files(i).name]);
    tacs_str_assoc(:,i) = frame(str_assoc_mask_index);
    tacs_str_vent(:,i) = frame(str_vent_mask_index);
    tacs_str_sens(:,i) = frame(str_sens_mask_index);
    tacs_cer(:,i) = frame(cer_mask_index);
    tacs_str_whole(:,i) = frame(str_whole_mask_index);
end

% Get mean tacs for each striatal region for regionwise analysis
str_assoc_mean = mean(tacs_str_assoc);
str_sens_mean = mean(tacs_str_sens);
str_vent_mean = mean(tacs_str_vent);

str_whole_mean = mean(tacs_str_whole);

% Get cerebellum blood input function
% bif{ii} = mean(t_tac(xroi{1},:)); What was in the script
bif = mean(tacs_cer);
cbif = cumtrapz(cumsum(dur),bif);
% 
n_early_frames_to_mask = -1; 
% 
if n_early_frames_to_mask == -1
    n_early_frames_to_mask = find(cbif == 0,1,'last');
end

k = n_early_frames_to_mask + 1:numel(bif);
% 
% 
% % Divide AUC by BIF
Bn = [cbif(k)./bif(k); ones(1,numel(bif(k)))]; 
% 

switch WhichMode
    case 'voxelwise'
    % Divide TACS by BIF: voxelwise
        % Associative striatum
        Rn_assoc = bsxfun(@rdivide, double(tacs_str_assoc(:,k)),bif(k));
        param_assoc = Rn_assoc/Bn;
        % Sensorimotor striatum
        Rn_sens = bsxfun(@rdivide, double(tacs_str_sens(:,k)),bif(k));
        param_sens = Rn_sens/Bn;
        % Ventral striatum
        Rn_vent = bsxfun(@rdivide, double(tacs_str_vent(:,k)),bif(k));
        param_vent = Rn_vent/Bn;
        % Whole striatum
        Rn_whole = bsxfun(@rdivide, double(tacs_str_whole(:,k)),bif(k));
        param_whole = Rn_whole/Bn;
    case 'regionwise'
    % Divide TACS by BIF: division using the mean tac of the region instead of voxel-wise
        % get ki value for associative striatum
        Rn_str_assoc = bsxfun(@rdivide, double(str_assoc_mean(:,k)),bif(k));
        param_str_assoc = Rn_str_assoc/Bn;
        % get ki value for sensorimotor striatum
        Rn_str_sens = bsxfun(@rdivide, double(str_sens_mean(:,k)),bif(k));
        param_str_sens = Rn_str_sens/Bn;
        % get ki value for ventral striatum
        Rn_str_vent = bsxfun(@rdivide, double(str_vent_mean(:,k)),bif(k));
        param_str_vent = Rn_str_vent/Bn;
        % get ki value for whole striatum
        Rn_str_whole = bsxfun(@rdivide, double(str_whole_mean(:,k)),bif(k));
        param_str_whole = Rn_str_whole/Bn;
end       
% the output for the mode of analysis above is a two column matrois, first
% column is clearance rate (Kicer) and the second is the volume of
% distribution

% Plot clearance rate (Ki) for assoc str
figure(1)
subplot(2,2,1)
scatter(Bn(1,:),Rn_str_assoc)
xlabel('auc/bif')
ylabel('tacs/bif')
title([subject,': assoc str patlak plot; ki value ',sprintf('%.6f',param_str_assoc(1,1))])

subplot(2,2,2)
scatter(Bn(1,:),Rn_str_vent)
xlabel('auc/bif')
ylabel('tacs/bif')
title([subject,': ventral str patlak plot; ki value ',sprintf('%.6f',param_str_vent(1,1))])

subplot(2,2,3)
scatter(Bn(1,:),Rn_str_sens)
xlabel('auc/bif')
ylabel('tacs/bif')
title([subject,': sens str patlak plot; ki value ',sprintf('%.6f',param_str_sens(1,1))])

subplot(2,2,4)
scatter(Bn(1,:),Rn_str_whole)
xlabel('auc/bif')
ylabel('tacs/bif')
title([subject,': whole str patlak plot; ki value ',sprintf('%.6f',param_str_whole(1,1))])


% Plot TACs
figure(2)
plot(str_assoc_mean,'Marker','o','Color','m')
hold on
plot(str_sens_mean,'Marker','o','Color','c')
plot(str_vent_mean,'Marker','o','Color','g')
plot(bif,'Marker','o','Color','r');

hold off
legend('assoc str','sens str','vent str','cerebellum')

title([subject,': mean TACs'])
xlabel('frames')
ylabel('MBq')

cd([projdir,'datadir/PET_files/',subject,'/'])
mkdir('patlak')
cd('patlak')
saveas(figure(1),'str_ki.fig')
saveas(figure(2),'mean_tacs.fig') 
 
regions = {'assoc','ventral','sens','whole_str','cereb'};
regions_tacs = {str_assoc_mean, str_vent_mean, str_sens_mean, str_whole_mean, bif};
regions_ki = {param_str_assoc, param_str_vent, param_str_sens, param_str_whole};
regions_ki_names = {'assoc','ventral','sens','whole_str'};

save([subject,'_tacs_stri_patlak.mat'],'regions','regions_tacs','regions_ki','regions_ki_names')

close all
end