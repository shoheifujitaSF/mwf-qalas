% Main script for MWF mapping using MWF-QALAS

%--------------------------------------------------------------------------
%% Set up
%--------------------------------------------------------------------------
clear; close all; clc;

load_dict   = 1;
save_nii    = 0;

addpath(genpath('./utils/'));


turbo_factor    = 128;
esp             = 5.8 * 1e-3;
alpha_deg       = 4;
num_acqs        = 6;
fprintf('done\n')

%% generating dictionaries for conventional and subspace mapping

TR                      = 5400e-3;
num_reps                = 5;
echo2use                = 3;
gap_between_readouts    = 900e-3;
time2relax_at_the_end   = 0;


% Entries for each compartment
t1_entries  = [150, 1000, 1200, 4500]; 
t2_entries  = [15, 60, 80, 500]; 


inv_eff     = 0.5:0.05:1.0; % 0.8;
inv_eff_sub = 0.5:0.05:1.0; % 0.8;
b1_val      = 0.65:0.05:1.35; % 1;
b1_val_sub  = 0.65:0.05:1.35; % 1;

T1_entries  = t1_entries; % PV
T1_entries  = T1_entries(:);

T2_entries  = t2_entries; % PV
T2_entries  = T2_entries(:);

t1t2_lut    = cat(2, T1_entries, T2_entries);

% remove cases where T2>T1
idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) < t1t2_lut(t,2)
        idx = idx+1;
    end
end

t1t2_lut_prune = zeross([length(t1t2_lut) - idx, 2]);

idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) >= t1t2_lut(t,2)
        idx = idx+1;
        t1t2_lut_prune(idx,:) = t1t2_lut(t,:);
    end 
end

fprintf('dictionary size: %d \n', length(t1t2_lut_prune)*length(b1_val)*length(inv_eff));

E = turbo_factor * num_acqs;

signal_conv_fit     = zeross([length(t1t2_lut_prune), num_acqs, length(b1_val), length(inv_eff)]);
signal_sub          = zeross([E, length(t1t2_lut_prune), length(b1_val_sub), length(inv_eff_sub)]);
signal_sub_fit      = zeross([E, length(t1t2_lut_prune), length(b1_val), length(inv_eff)]);

length_b1_val       = length(b1_val);
length_b1_val_sub   = length(b1_val_sub);
length_inv_eff      = length(inv_eff);
length_inv_eff_sub  = length(inv_eff_sub);


% for fitting
cnt         = 0;
iniTime1    = clock;
for b1 = 1:length_b1_val
    for ie = 1:length_inv_eff
        cnt             = cnt + 1;
        [Mz_, Mxy_]     = sim_qalas_pd_b1_eff_T2_mwf(TR, alpha_deg, esp, ...
                            turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, ...
                            num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, ...
                            b1_val(b1), inv_eff(ie));
        
        temp_conv       = abs(Mxy_(:,:,end).');
        
        for n = 1:size(temp_conv,1)
            temp_conv(n,:)      = temp_conv(n,:) / sum(abs(temp_conv(n,:)).^2)^0.5;
        end
        
        signal_conv_fit(:,:,b1,ie)  = temp_conv;

    end
end
delete(gcp('nocreate'))
fprintf('total elapsed time: %.1f sec\n\n',etime(clock,iniTime1));

figure; plot(Mz_(:,:,end),'LineWidth',2); title('MWF-QALAS Mz');
figure; plot(Mxy_(:,:,end),'LineWidth',2); title('MWF-QALAS Mxy');


length_dict = length(t1t2_lut_prune);
dict_cnv    = zeross([length_dict * length(inv_eff), num_acqs, length(b1_val)]);
for t = 1:length(inv_eff)
    dict_cnv(1 + (t-1)*length_dict : t*length_dict, :, :) = signal_conv_fit(:,:,:,t);
end
figure; plot(dict_cnv(:,:)','LineWidth',2); title('QALAS dictionary for conventional fitting');

fprintf('done\n')


%% PV dictionary-based

dict_cnv_    = permute(dict_cnv,[2,1]);
dict_pv      = inv(dict_cnv_.' * dict_cnv_) * dict_cnv_.';


PV_dict_1   = 0.02:0.02:0.3; % myelin water
PV_dict_2   = 0.02:0.02:0.5; % intra/extracellular water
PV_dict_3   = 0.02:0.02:0.5; % intra/extracellular water
PV_dict_4   = 0.05:0.05:1.0; % free water

PV_dict     = zeros(4,length(PV_dict_1)*length(PV_dict_2)*length(PV_dict_3)*length(PV_dict_4));

cnt = 0;

for pp = 1:length(PV_dict_1)
    for qq = 1:length(PV_dict_2)
        for rr = 1:length(PV_dict_3)
            for ss = 1:length(PV_dict_4)
                cnt = cnt + 1;
                PV_dict(:,cnt) = [PV_dict_1(pp),PV_dict_2(qq),PV_dict_3(rr),PV_dict_4(ss)]./ ...
                    (PV_dict_1(pp)+PV_dict_2(qq)+PV_dict_3(rr)+PV_dict_4(ss));
            end
        end
    end
end

PV_dict = unique(PV_dict','rows')';
dict_pv_full = dict_cnv_ * PV_dict;

for cc = 1:size(dict_pv_full,2)
    dict_pv_full(:,cc) = dict_pv_full(:,cc) / sum(abs(dict_pv_full(:,cc)).^2)^0.5;
end

% 6 source voumes of mwf qalas
img_path = '/mwf_qalas_hbcd_1mm_R5_nomoco';

save_path = [img_path,'_mwf.mat'];

load([img_path,'.mat'], 'img_mwf_qalas')

img_contrasts = permute(img_mwf_qalas, [1 2 4 3]);


PV1_map = zeros([size(img_contrasts,1),size(img_contrasts,2),size(img_contrasts,4)]);
PV2_map = zeros([size(img_contrasts,1),size(img_contrasts,2),size(img_contrasts,4)]);
PV3_map = zeros([size(img_contrasts,1),size(img_contrasts,2),size(img_contrasts,4)]);
PV4_map = zeros([size(img_contrasts,1),size(img_contrasts,2),size(img_contrasts,4)]);

for ss = 1:size(img_contrasts,4)
    tic
    fprintf('slice: %d/%d... ',ss,size(img_contrasts,4));
    img_contrasts_ = img_contrasts(:,:,:,ss);
    
    [PV1_map(:,:,ss),PV2_map(:,:,ss),PV3_map(:,:,ss),PV4_map(:,:,ss)] = ...
        dict_fit_qalas_sub_pv_recon(permute(img_contrasts_,[1,2,4,3]),permute(dict_pv_full,[2,1]),permute(PV_dict,[2,1]));
    toc
end
fprintf('done\n')

imagesc3d2( PV1_map(:,:,:), s(PV1_map)/2+[28 5 -10], 1, [0,180,180], [0,0.5]), setGcf(.5)

imagesc3d2( PV2_map(:,:,:), s(PV1_map)/2+[28 5 -10], 2, [0,180,180], [0,1]), setGcf(.5)

imagesc3d2( PV3_map(:,:,:), s(PV1_map)/2+[28 5 -10], 3, [0,180,180], [0,1]), setGcf(.5)

imagesc3d2( PV4_map(:,:,:), s(PV1_map)/2+[28 5 -10], 4, [0,180,180], [0,1]), setGcf(.5)


figure(10);
tiledlayout(1,4, 'TileSpacing', 'compact');
nexttile;imshow(sq(PV1_map(end/2,:,:)),[0,0.3]);title('Short T2 component');
nexttile;imshow(sq(PV2_map(end/2,:,:)),[0,1]);title('Intermediate T2 component');
nexttile;imshow(sq(PV3_map(end/2,:,:)),[0,1]);title('Intermediate T2 component');
nexttile;imshow(sq(PV4_map(end/2,:,:)),[0,1]);title('Long T2 component');


figure(21);
tiledlayout(1,6, 'TileSpacing', 'compact'); 
for i=1:6 
    nexttile;
imshow(flip(abs(sq(img_contrasts(:,:,i,end/2))),1),[0,4e-4]);title(['Contrast ',num2str(i)])

end

% save(save_path, 'PV1_map','PV2_map','PV3_map','img_contrasts')