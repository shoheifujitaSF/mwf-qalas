%--------------------------------------------------------------------------
%% load qalas acquisition:  
%--------------------------------------------------------------------------

addpath utils

% standard acquisition:
dicom_path = [pwd, '/qalas_dicoms/'];

dicom_name = dir(dicom_path);
dicom_name = dicom_name(end).name;

img = dicomread_dir(dicom_path);

N = s(img) ./ [1,1,5]

img = single(reshape(img, [N, 5]));

for t = 1:5
    imagesc3d2( img(:,:,:,t), s(img)/2, t, [-0,-0,0], [0,500]), 
end


%--------------------------------------------------------------------------
%% ROI mask
%--------------------------------------------------------------------------

% uncomment if on linux: precompiled mex file should work for BET
%msk = BET(rsos(img,4), N, [1,1,1]*1.2, 0.5);

% load precomputed bet mask:
load msk


t = 5;
imagesc3d2( img(:,:,:,t), s(img)/2, t, [-0,-0,0], [0,350]), 
imagesc3d2( img(:,:,:,t) .* msk, s(img)/2, 1+t, [-0,-0,0], [0,350]), 
imagesc3d2( msk, s(img)/2, 2+t, [-0,-0,0], [0,1]), 


%--------------------------------------------------------------------------
%% b1 map: interpolate
%--------------------------------------------------------------------------

b1_path = [pwd, '/tfl_b1map/9/'];
img_b1_path =  [pwd, '/tfl_b1map/8/'];

% uncomment if imresize3 exists (from matlab2017a onwards)
%img_b1 = dicomread_dir(b1_path);
%img_b1 = imresize3(slice_uninterleave( img_b1 ), N) / 800;
%im_b1 = dicomread_dir(img_b1_path);
%im_b1 = imresize3(slice_uninterleave( im_b1 ), N) / 800;

% load precomputed b1 map and tfl image:
load img_b1
load im_b1


imagesc3d2( img_b1 .* msk, s(img)/2, 1, [-0,-0,0], [0.5,1.5]), colormap jet
imagesc3d2( im_b1, s(img)/2, 2, [-0,-0,0], [0,3]), 
imagesc3d2( img(:,:,:,5), s(img)/2, 3, [-0,-0,0], [0,350]), 


%--------------------------------------------------------------------------
%% b1 map: polynomial fitting
%--------------------------------------------------------------------------

order_3d = 4;
I_fitted4 = polyfit3D_NthOrder(img_b1, msk, order_3d);


imagesc3d2( img_b1 .* msk, s(img)/2, 11, [-0,-0,0], [0.75,1.35]), colormap jet, colorbar, setGcf(.5)
imagesc3d2( I_fitted4 .* msk, N/2, 12, [-0,-0,0], [0.75,1.35]), colormap jet, colorbar, setGcf(.5)

disp([min(img_b1(msk==1)), max(img_b1(msk==1))])
disp([min(I_fitted4(msk==1)), max(I_fitted4(msk==1))])


%--------------------------------------------------------------------------
%% dicom header
%--------------------------------------------------------------------------

di = dicominfo([dicom_path, dicom_name]);

esp = di.RepetitionTime * 1e-3     % echo spacing in sec

turbo_factor = di.EchoTrainLength  % ETL


%--------------------------------------------------------------------------
%% create look up table 
%--------------------------------------------------------------------------

t1_entries = [5:2:1000, 1002:5:3000, 3100:100:5000];
t2_entries = [1:1:700, 710:20:1000, 1100:100:2500];


T1_entries = repmat(t1_entries.', [1,length(t2_entries)]).';
T1_entries = T1_entries(:);
  

T2_entries = repmat(t2_entries.', [1,length(t1_entries)]);
T2_entries = T2_entries(:);


t1t2_lut = cat(2, T1_entries, T2_entries);


% remove cases where T2>T1
idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) < t1t2_lut(t,2)
        idx = idx+1;
    end
end

t1t2_lut_prune = zeross( [length(t1t2_lut) - idx, 2] );

idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) >= t1t2_lut(t,2)
        idx = idx+1;
        t1t2_lut_prune(idx,:) = t1t2_lut(t,:);
    end
end

disp(['dictionary entries: ', num2str(length(t1t2_lut_prune))])


TR = 4500e-3;       % excluding dead time at the end of the sequence 
alpha_deg = 4;

num_reps = 5;       % number of repetitions to simulate to reach steady state
echo2use = 1;

gap_between_readouts = 900e-3;
time2relax_at_the_end = 0;  % actual TR on the console = TR + time2relax_at_the_end


% simulate Mz and Mxy for T1 T2 values in the look up table
%b1_val = 1;     % use this to ignore b1 effects
b1_val = linspace( min(I_fitted4(msk==1)), max(I_fitted4(msk==1)), 20 )

% inversion efficiency
inv_eff = 0.96;
% inv_eff = 1;

signal = zeross([length(t1t2_lut_prune), 5, length(b1_val)]);

for b1 = 1:length(b1_val)
    [Mz, Mxy] = sim_qalas_pd_b1_eff(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b1), inv_eff);

    signal(:,:,b1) = abs( Mxy(:,:,end).' );     % transverse magnitude signal

    % normalize every entry of dict.signal
    for n = 1:size(signal,1)
        signal(n,:,b1) = signal(n,:,b1) / sum(abs(signal(n,:,b1)).^2)^0.5;
    end

    figure(b1), subplot(211), plot(Mz(:,1:1e3:end,end), 'linewidth', 3), title('Mz', 'fontsize', 48), ylim([-1,1])
    figure(b1), subplot(212), plot(abs(Mxy(:,1:1e3:end,end)), 'linewidth', 3), title('Mxy', 'fontsize', 48), ylim([-.05,.1])
end


%--------------------------------------------------------------------------
%% create masks for each b1 value
%--------------------------------------------------------------------------

use_poly_fit = 1;

if use_poly_fit
    img_b1 = I_fitted4 .* msk;
end

if length(b1_val) == 1
    % do not use b1 correction
    msk_b1 = msk;
else

    msk_b1 = zeross([N,length(b1_val)]);
    
    for t = 1:length(b1_val)
        if t > 1
            msk_b1(:,:,:,t) = (img_b1 <= b1_val(t)) .* (img_b1 > b1_val(t-1));
        else
            % first bin: contains all pixels less than b1_val(1)
            msk_b1(:,:,:,t) = img_b1 <= b1_val(t);
        end
    
        if t == length(b1_val)
            % last bin: contains all pixels greater than b1_val(end-1)
            msk_b1(:,:,:,t) = img_b1 > b1_val(t-1);
        end
        
        msk_b1(:,:,:,t) = msk_b1(:,:,:,t) .* msk;
    
        imagesc3d2( img(:,:,:,5) .* sum(msk_b1(:,:,:,1:t),4), N/2, 1, [-0,-0,0], [0,500], [], num2str(b1_val(t))), 
    end
end

imagesc3d2( img(:,:,:,5) .* sum(msk_b1,4), N/2, 2, [-0,-0,0], [0,500]), 
 


%--------------------------------------------------------------------------
%% dictionary fit -> in each slice, bin voxels based on b1 value
%--------------------------------------------------------------------------

% can use parfor when b1 map is used, uncomment following 4 lines:
delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object
total_cores = c.NumWorkers; 
parpool( min(ceil(total_cores*.5), N(3)) ) 

T1_map = zeross(N);
T2_map = zeross(N);
PD_map = zeross(N);


parfor slc_select = 1:N(3)    % uncomment to use parfor
%for slc_select = 1:N(3)
    disp(num2str(slc_select))
    
    for b = 1:length(b1_val)
        msk_slc = msk_b1(:,:,slc_select,b);
        num_vox = sum(msk_slc(:)~=0);

        if num_vox > 0
            img_slc = zeross([5, num_vox]);
        
            for t = 1:5
                temp = sq(img(:,:,slc_select,t));
                img_slc(t,:) = temp(msk_slc~=0);
            end
        
            tic   
                res = signal(:,:,b) * img_slc; % dot product    
        
                [~, max_idx] = max(abs(res), [], 1);
        
                res_map = t1t2_lut_prune(max_idx,:);
        
                [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), inv_eff);
            toc
        
            t1_map = zeross(N(1:2));
            t1_map(msk_slc==1) = res_map(:,1);
        
            t2_map = zeross(N(1:2));
            t2_map(msk_slc==1) = res_map(:,2);
        
            Mxy_sim_use = abs(Mxy_sim(:,:,end));
        
            scl = zeross([num_vox,1]);
        
            for idx = 1:size(Mxy_sim_use,2)
                scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
            end
    
            pd_map = zeross(N(1:2));
            pd_map(msk_slc~=0) = scl;     
    
            T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
            T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
            PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
        end
    end

end


delete(gcp('nocreate'))   % uncomment when using parfor


%--------------------------------------------------------------------------
%% save
%--------------------------------------------------------------------------

imagesc3d2( T1_map, N/2, 1, [-0,-0,0], [0,2500]), colormap jet
imagesc3d2( T2_map, N/2, 2, [-0,-0,0], [0,250]), colormap jet
imagesc3d2( PD_map, N/2, 3, [-0,-0,0], [0,11000]), 


save_path = [pwd, '/results/'];

mkdir(save_path)

save([save_path, 'T1_map.mat'], 'T1_map')
save([save_path, 'T2_map.mat'], 'T2_map')
save([save_path, 'PD_map.mat'], 'PD_map')


