%--------------------------------------------------------------------------
%% load qalas acquisition:  
%--------------------------------------------------------------------------
set(0,'DefaultFigureWindowStyle','docked')

addpath utils/

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

imagesc3d2( rsos(img,4), s(img)/2, 1+t, [-0,-0,0], [0,500]), 


%--------------------------------------------------------------------------
% ROI mask
%--------------------------------------------------------------------------

msk = zeross(N);

for slc_select = 1:N(3)
    msk(:,:,slc_select) = imerode(imfill(rsos(img(:,:,slc_select,:),4) > 100, 'holes'), strel('disk', 3));
end

t = 5;
imagesc3d2( img(:,:,:,t), s(img)/2, t, [-0,-0,0], [0,350]), 
imagesc3d2( img(:,:,:,t) .* msk, s(img)/2, 1+t, [-0,-0,0], [0,350]), 
imagesc3d2( msk, s(img)/2, 2+t, [-0,-0,0], [0,1]), 



%--------------------------------------------------------------------------
%% b1 map: interpolate
%--------------------------------------------------------------------------

b1_path = [pwd, '/TFL_B1MAP_SAGITAL_0005/'];


img_b1_load = dicomread_dir(b1_path);

% requires imresize3 (available in later matlab versions):
%img_b1_load = imresize3(slice_uninterleave( img_b1_load ), N) / 800;

% directly load .mat file instead:
load img_b1_load

img_b1 = img_b1_load;


imagesc3d2( img_b1 .* msk, s(img)/2, 1, [-0,-0,0], [0.5,1.5]), colormap jet, colorbar, setGcf(.5)
imagesc3d2( img(:,:,:,5), s(img)/2, 2, [-0,-0,0], [0,350]), 



%--------------------------------------------------------------------------
%% threshold high and low B1 values: use raw b1 map without polyfit
%--------------------------------------------------------------------------

thre_high = 1.35;
thre_low = 0.65;        

temp = img_b1_load .* msk;

temp(temp > thre_high) = thre_high;
temp(temp < thre_low) = thre_low;

temp = temp .* msk;
img_b1 = temp .* msk;    


rmse(temp .* msk, img_b1_load .* msk)

imagesc3d2( temp .* msk, N/2, 12, [-0,-0,0], [0.75,1.35]), colormap jet, colorbar, setGcf(.5)

disp([min(temp(msk==1)), max(temp(msk==1))])



%--------------------------------------------------------------------------
%% create masks for each b1 value
%--------------------------------------------------------------------------


num_b1_bins = 20;

b1_val = linspace( min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins )


sum_msk = sum(msk(:));

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
            msk_b1(:,:,:,t) = msk.*(img_b1 <= b1_val(t));
        end
        
        percent_bin = sum(sum(sum(msk_b1(:,:,:,t),1),2),3) / sum_msk;

    
        if t == length(b1_val)
            % last bin: contains all pixels greater than b1_val(end-1)
            msk_b1(:,:,:,t) = img_b1 > b1_val(t-1);
        end
        
        msk_b1(:,:,:,t) = msk_b1(:,:,:,t) .* msk;
    
        imagesc3d2( img(:,:,:,5) .* msk_b1(:,:,:,t), N/2, 1, [-0,-0,0], [0,500], [], [num2str(b1_val(t)), '  voxels in this bin: ', num2str(100*percent_bin), '%']), pause(.1)
    end
end

imagesc3d2( img(:,:,:,5) .* sum(msk_b1,4), N/2, 1, [-0,-0,0], [0,500]), 


%--------------------------------------------------------------------------
%% dicom header
%--------------------------------------------------------------------------

di = dicominfo([dicom_path, dicom_name]);

esp = di.RepetitionTime * 1e-3     % echo spacing in sec

turbo_factor = di.EchoTrainLength  % ETL


%--------------------------------------------------------------------------
%% create look up table 
%--------------------------------------------------------------------------


t1_entries = [5:5:3000, 3100:100:5000];
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

% inversion efficiency
inv_eff = 0.7:0.05:.9
% inv_eff = 1;

signal = zeross([length(t1t2_lut_prune), 5, length(b1_val), length(inv_eff)]);

delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

total_cores = c.NumWorkers; 
parpool( min(ceil(total_cores*.5), length(inv_eff)) ) 


length_b1_val = length(b1_val);

parfor ie = 1:length(inv_eff)
    for b1 = 1:length_b1_val                    
        [Mz, Mxy] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b1), inv_eff(ie));

        temp = abs( Mxy(:,:,end).' );     % transverse magnitude signal

        % normalize every entry of dict.signal
        for n = 1:size(temp,1)
            temp(n,:) = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
        end
        
        signal(:,:,b1,ie) = temp;
    end    
end

delete(gcp('nocreate'))



%--------------------------------------------------------------------------
%%  concat dictionary across inv_Eff dimension
%--------------------------------------------------------------------------

length_dict = length(t1t2_lut_prune);

dict = zeross([length_dict * length(inv_eff), 5, length(b1_val)]);

for t = 1:length(inv_eff)
    dict(1 + (t-1)*length_dict : t*length_dict, :, :) = signal(:,:,:,t);
end



%--------------------------------------------------------------------------
%% dictionary fit -> in each slice, bin voxels based on b1 value
%--------------------------------------------------------------------------

estimate_pd_map = 0;     % set to 1 to estiamte PD map (makes it slower)

T1_map = zeross(N);
T2_map = zeross(N);

PD_map = zeross(N);     % proton density-> won't be estimated if estimate_pd_map is set to 0
IE_map = zeross(N);     % inversion efficiency


if length(b1_val) > 1
    % use parfor
    
    delete(gcp('nocreate'))
    c = parcluster('local');    % build the 'local' cluster object
 
    total_cores = c.NumWorkers; % 48 cores for marduk
    parpool( min(ceil(total_cores*.5), N(3)) ) 

    parfor slc_select = 1:N(3)
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
                    res = dict(:,:,b) * img_slc; % dot product    
    
                    % find T1, T2 values    
                    [~, max_idx] = max(abs(res), [], 1);

                    max_idx_t1t2 = mod(max_idx, length_dict);
                    max_idx_t1t2(max_idx_t1t2==0) = length_dict;

                    res_map = t1t2_lut_prune(max_idx_t1t2,:);

                    % find inv_eff                         
                    max_idx_ie = 1 + (max_idx - max_idx_t1t2) / length_dict;

                    ie_to_use = inv_eff(max_idx_ie);

                    if estimate_pd_map
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                    end
                toc

                t1_map = zeross(N(1:2));
                t1_map(msk_slc==1) = res_map(:,1);

                t2_map = zeross(N(1:2));
                t2_map(msk_slc==1) = res_map(:,2);

                ie_map = zeross(N(1:2));
                ie_map(msk_slc==1) = ie_to_use;

                if estimate_pd_map
                    Mxy_sim_use = abs(Mxy_sim(:,:,end));
    
                    scl = zeross([num_vox,1]);
    
                    for idx = 1:size(Mxy_sim_use,2)
                        scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                    end
    
                    pd_map = zeross(N(1:2));
                    pd_map(msk_slc~=0) = scl;     
                    PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
                end

                T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
                T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
                IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;

            end
        end
    end
    
    delete(gcp('nocreate'))
    
else
    for slc_select = 1:N(3)
        % don't use parfor with b1 map is not used
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
                    res = dict(:,:,b) * img_slc; % dot product    

                    % find T1, T2 values    
                    [~, max_idx] = max(abs(res), [], 1);

                    max_idx_t1t2 = mod(max_idx, length_dict);
                    max_idx_t1t2(max_idx_t1t2==0) = length_dict;

                    res_map = t1t2_lut_prune(max_idx_t1t2,:);

                    % find inv_eff                         
                    max_idx_ie = 1 + (max_idx - max_idx_t1t2) / length_dict;

                    ie_to_use = inv_eff(max_idx_ie);
    
                    if estimate_pd_map
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                    end
                toc

                t1_map = zeross(N(1:2));
                t1_map(msk_slc==1) = res_map(:,1);

                t2_map = zeross(N(1:2));
                t2_map(msk_slc==1) = res_map(:,2);

                ie_map = zeross(N(1:2));
                ie_map(msk_slc==1) = ie_to_use;

                if estimate_pd_map
                    Mxy_sim_use = abs(Mxy_sim(:,:,end));
    
                    scl = zeross([num_vox,1]);
    
                    for idx = 1:size(Mxy_sim_use,2)
                        scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                    end
    
                    pd_map = zeross(N(1:2));
                    pd_map(msk_slc~=0) = scl;     
                    PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
                end

                T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
                T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
                IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;
            end
        end
    end
end




%--------------------------------------------------------------------------
%% show results
%--------------------------------------------------------------------------

imagesc3d2( T1_map, N/2, 1, [-0,-0,0], [0,2000]), colormap jet
imagesc3d2( T2_map, N/2, 2, [-0,-0,0], [0,200]), colormap jet
imagesc3d2( PD_map, N/2, 3, [-0,-0,0], [0,23000]), 
imagesc3d2( IE_map, N/2, 3, [-0,-0,0], [0,1]), 


