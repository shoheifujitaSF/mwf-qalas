%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------

clear; clc;
disp('Loading data');

addpath(genpath('./utils/'));

% example data available at https://zenodo.org/records/15546728
load('kspace_unsorted_mwf_qalas_hbcd_1mm_R5.mat');


% For Siemens data
% data_file_path = '/mwf_qalas_hbcd_1mm_R5.dat';
% [p,n,e] = fileparts(data_file_path);
% basic_file_path = fullfile(p,n);
% twix_obj = mapVBVD2(data_file_path);
% data_unsorted = twix_obj{end}.image.unsorted();

[adc_len,ncoil,readouts] = size(data_unsorted);

% Read params from seq file
% pulseq_file_path = [p, '/',regexprep(n, 'meas_MID\d+_FID\d+_', ''), '.seq'];
pulseq_file_path            =       'REPOSITORY/sequences/mwf_qalas_hbcd/mwf_qalas_hbcd_1mm_R5.seq';

seq = mr.Sequence();
seq.read(pulseq_file_path);

N = seq.getDefinition('Matrix');
nTR = seq.getDefinition('nTR');
nETL = seq.getDefinition('nETL');

os_factor = seq.getDefinition('os_factor');

traj_y = seq.getDefinition('traj_y');
traj_z = seq.getDefinition('traj_z');
step_size = 6 * nETL;


%--------------------------------------------------------------------------
%% patref scan
%--------------------------------------------------------------------------
disp('Ref scan');

Ny_acs=32;
Nz_acs=32;


temp = 1:N(2);
iY_acs_indices = temp(1+end/2-Ny_acs/2:end/2+Ny_acs/2);

temp = 1:N(3);
iZ_acs_indices = temp(1+end/2-Nz_acs/2:end/2+Nz_acs/2);

ref = zeros(N(1),N(2),N(3),ncoil);

data_unsorted_ref = data_unsorted(:,:,1:Ny_acs*Nz_acs);

index=1;
for iZ = iZ_acs_indices
    for iY = iY_acs_indices

            ref(:,iY,iZ,:) = data_unsorted_ref(:,:,index);    
            index=index+1;
    end
end


img_ref = ifft3call(ref);
imagesc3d2(rsos(img_ref,4), s(img_ref)/2, 1, [0,180,180], [-0,1e-3]), setGcf(.5)



data_unsorted_img = data_unsorted(:,:,Ny_acs*Nz_acs+1:end);
kspace = zeros([N(1)*os_factor N(2) N(3) 6 ncoil]);
for TR=1:nTR
    for contrast=1:6
        index_start = (contrast-1)*nETL + (TR-1)*step_size +1;

        for index=(index_start):(index_start+nETL-1)
            ky = traj_y(index);
            kz = traj_z(index);
            kspace(:,ky,kz,contrast,:) = data_unsorted_img(:,:,index);     
        end
    end
end



%--------------------------------------------------------------------------
% coil compression
%--------------------------------------------------------------------------

num_chan = 16;  % num channels to compress to

[ref_svd, cmp_mtx] = svd_compress3d(ref, num_chan, 1);

rmse(rsos(ref_svd,4), rsos(ref,4))

N = size(kspace(:,:,:,1,1));

kspace = permute(kspace, [1,2,3,5,4]);

num_eco = size(kspace,5);

kspace_svd = zeross([N,num_chan,num_eco]);

% apply compression to imaging data
for t = 1:size(kspace,5)
    kspace_svd(:,:,:,:,t) = svd_apply3d(kspace(:,:,:,:,t), cmp_mtx);
end

rmse(rsos(kspace_svd,4), rsos(kspace,4))



patref_pad = ref_svd;

img_patref_pad = ifft3c(patref_pad);

imagesc3d2( rsos(img_patref_pad,4), s(img_patref_pad)/2, 10, [0,180,180], [0,2e-4],1,'img_patref_pad'), setGcf(.5)
imagesc3d2( rsos(ifft3call(kspace_svd(:,:,:,:,end)),4), s(img_patref_pad)/2, 11, [0,180,180], [0,2e-4],1,'kspace_svd ifft'), setGcf(.5)


%--------------------------------------------------------------------------
%% calculate sens map using ESPIRiT: parfor
%--------------------------------------------------------------------------
disp('Calculating sens map');

num_acs = min(Ny_acs,Nz_acs);
kernel_size = [6,6];
eigen_thresh = 0.7;

receive = zeross(size(kspace_svd(:,:,:,:,1)));


delete(gcp('nocreate'))
c = parcluster('local');    

total_cores = c.NumWorkers;  
parpool(ceil(total_cores/2))


tic
parfor slc_select = 1:s(img_patref_pad,1)     

    [maps, weights] = ecalib_soft( fft2c( sq(img_patref_pad(slc_select,:,:,:)) ), num_acs, kernel_size, eigen_thresh );

    receive(slc_select,:,:,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
end 
toc

delete(gcp('nocreate'))

% remove the phase of first channel to eliminate slice to slice phase jumps
rec = abs(receive) .* exp(1i * angle( receive .* repmat(conj(receive(:,:,:,1)), [1,1,1,s(receive,4)]) ));

% save([basic_file_path, '_rec_svd_', num2str(num_chan), 'ch.mat'], 'rec', '-v7.3')


%--------------------------------------------------------------------------
%% load coil sensitivities
%--------------------------------------------------------------------------

% tic
%     load([data_path, 'rec_svd_', num2str(num_chan), 'ch_test.mat'])
%     load([data_path, 'cmp_mtx_N.mat'])
% toc


%--------------------------------------------------------------------------
%% lsqr 3d-Sense recon each time-binned under-sampled k space data
%--------------------------------------------------------------------------
disp('Reconstructing low-res image from each TR');

% whosbetter -b
clear kspace kspace_svd patref_pad receive

window=1; %nTR; %num of TR to use for each time frame. Must be devisor of 48 in current implementation
nTR_abbr = nTR/window;        

kspace_c345_combined = single(zeros([N(1)*os_factor N(2) N(3) num_chan 1 nTR_abbr])); %changed to single, and use only 3 contrasts


for time_frame = 1:nTR_abbr
    kspace_tmp = single(zeros([N(1)*os_factor N(2) N(3) ncoil]));
    kspace_tmp_svd = single(zeros([N(1)*os_factor N(2) N(3) num_chan]));
    for TR=((time_frame-1)*window +1):(time_frame*window)
        for contrast=4:6
            index_start = (contrast-1)*nETL + (TR-1)*step_size +1;

            for index=(index_start):(index_start+nETL-1)
                ky = traj_y(index);
                kz = traj_z(index);
                if all(kspace_tmp(:,ky,kz,:) == 0)
                    kspace_tmp(:,ky,kz,:) = data_unsorted_img(:,:,index); %v2
                    %break;
                end
            end
        end
    end

    kspace_tmp_svd = single(svd_apply3d(kspace_tmp(:,:,:,:), cmp_mtx));

    kspace_c345_combined(:,:,:,:,1,time_frame) = kspace_tmp_svd;
end


% Quick check with ifft
% for time_frame=1%nTR_abbr
%     imagesc3d2( rsos(ifft3call(kspace_c345_combined(:,:,:,:,1,time_frame)),4), s(img_patref_pad)/2, time_frame, 180+[0,0,0], [0,1e-3],1,['time frame:',num2str(time_frame)]), setGcf(.5)
% end


msk = single(zeros(N(1),N(2),N(3),nTR_abbr));
msk(:,1+end/2-16:end/2+16,1+end/2-16:end/2+16,:) = ones(N(1),32,32,nTR_abbr); %32x32 msk
msk = logical(reshape(msk, [N(1) N(2) N(3) 1 1 nTR_abbr]));
msk = repmat(msk,[1,1,1,num_chan,1]);

kspace2use = kspace_c345_combined.*msk;
clear kspace_c345_combined %release memory

%lsqr_iter = 20;             % max num of iterations
lsqr_iter = 50;             % max num of iterations

%lsqr_tol = 1e-2;            % tolerance to terminate iterations
lsqr_tol = 1e-3;            % tolerance to terminate iterations

m3d = kspace2use~=0;

img_sense = zeros([size(kspace2use,1), size(kspace2use,2), size(kspace2use,3), size(kspace2use,5)]);


param = [];
param.lambda = 1e-2;        % L2 regularization amount
param.lambda = 1e-3;        % L2 regularization amount

param.N = size(m3d(:,:,:,1,1));
param.num_chan = size(m3d,4);

param.sens = rec;

tic
for tr = 1:nTR_abbr
    disp(['Tr: ', num2str(tr)])

    param.m3d = m3d(:,:,:,:,tr);         

    kspace_coils = kspace2use(:,:,:,:,tr);

    res = lsqr(@apply_sense_tikc3d, cat(1, kspace_coils(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  

    Res = reshape(res, param.N);

    img_sense(:,:,:,tr) = Res;


end
toc
clear kspace2use
% save([data_path, 'img_sense_lambda1Eneg3_iter20_test.mat'], 'img_sense', '-v7.3')
% save([data_path, 'img_sense_R_lambda1Eneg3_iter20_test.mat'], 'img_sense_R', '-v7.3')

% for tr = 45%:nTR_abbr
%     imagesc3d2( img_sense(:,:,:,tr), s(img_sense)/2, tr, [180,180,180], [0,2e-4],1,['time frame:',num2str(tr)]), setGcf(.5)
%     %imagesc3d2( angle(img_sense(:,:,:,tr)), s(img_sense)/2, tr+10, [180,180,180], [-pi,pi]), setGcf(.5)
% 
%     % imagesc3d2( img_sense_R(:,:,:,t), s(img_sense)/2, t+100, [180,180,180], [0,2e-4]), setGcf(.5)
% end



%--------------------------------------------------------------------------
%% estimate motion
%--------------------------------------------------------------------------
disp('Estimating motion');

param = [];

[optimizer, metric] = imregconfig('monomodal');
transformType = 'rigid';
interpType = 'cubic';

for tr=1:nTR_abbr
%    img_N = img_sense(:,:,:,1);
    img_N = mean(abs(img_sense(:,:,:,1:2)),4); %v2
    img_N = img_N / prctile(abs(img_N(:)),75);  % signal normalization assumed
     
    img_R = img_sense(:,:,:,tr);
    img_R = img_R / prctile(abs(img_R(:)),75); % signal normalization assumed

    if tr == 1 %identity
        tform_N2N = imregtform(abs(img_N), abs(img_N), transformType, optimizer, metric);
        param.tform{tr} = tform_N2N;
        param.tformInv{tr} = invert(param.tform{tr});

    else
        tform_N2R = imregtform(abs(img_N), abs(img_R), transformType, optimizer, metric);
        param.tform{tr} = tform_N2R;
        param.tformInv{tr} = invert(param.tform{tr});

        %img_R2N = imwarp(img_R, param.tformInv{tr}, interpType, 'OutputView', imref3d(size(img_N)), 'FillValues', 0);
        %imagesc3d2( img_R2N, s(img_sense)/2, 20+tr, [180,180,180], [0,2],1,['Current to initial:',num2str(tr)]), setGcf(.5)
    end
end

tform_matrix = param.tform; %save this for cs
tform_matrixInv = param.tformInv; %save this for cs

%--------------------------------------------------------------------------
%% lsqr 3d-Sense recon [without CS] -> with moco
%--------------------------------------------------------------------------

img_sense_mc = zeros([N, num_eco]);


for t = 4%1:5
    disp(['lsqr 3d-Sense recon with moco contrast: ',num2str(t)]);
    tic
    % sort kspace to be used for reconstruction
    kspace2use = single(zeros([N(1)*os_factor N(2) N(3) num_chan 1 nTR_abbr])); % only 1 contrast at a time

    for time_frame = 1:nTR_abbr
        kspace_tmp = single(zeros([N(1)*os_factor N(2) N(3) 1 ncoil]));
        kspace_tmp_svd = single(zeros([N(1)*os_factor N(2) N(3) num_chan 1]));
        for TR=((time_frame-1)*window +1):(time_frame*window)
            for contrast=t
                index_start = (contrast-1)*nETL + (TR-1)*step_size +1;

                % for index=(index_start):(index_start+nETL-1)
                for index=(index_start+8):(index_start+nETL-1)
                    ky = traj_y(index);
                    kz = traj_z(index);
                    kspace_tmp(:,ky,kz,1,:) = data_unsorted_img(:,:,index);
                end

            end
        end
        kspace_tmp = permute(kspace_tmp, [1,2,3,5,4]);

        kspace_tmp_svd(:,:,:,:,1) = single(svd_apply3d(kspace_tmp(:,:,:,:,1), cmp_mtx));

        kspace2use(:,:,:,:,:,time_frame) = kspace_tmp_svd;
    end


    idx = 1;

    m3d_tmp = sq(kspace2use(:,:,:,1,:,:));
    param.m3d = reshape(m3d_tmp,[N(1), N(2), N(3),1,nTR_abbr])~= 0;  %dim N(1) N(2) N(3) 1 nTR

    param.nTR = size(param.m3d,5);

    param.num_chan = num_chan;

    kspace_coils = zeros([sum(param.m3d(:)) * param.num_chan, 1]);

    for tr = 1:param.nTR
        % sampling mask for the current TR
        msk = repmat(param.m3d(:,:,:,:,tr), [1,1,1,param.num_chan,1]);

        num_points_current_TR = sum(msk(:));

        temp = kspace2use(:,:,:,:,1,tr) .* msk;

        % concatenate all k-space data
        kspace_coils(idx:idx+num_points_current_TR-1) = temp( msk == 1 );

        idx = idx + num_points_current_TR;
    end


    lsqr_iter = 20;         % max num of iterations
    %lsqr_iter = 1;         % max num of iterations

    lsqr_tol = 1e-3;        % tolerance to terminate iterations
    %lsqr_tol = 1e-1;        % tolerance to terminate iterations




    param.lambda = 1e-3;        % L2 regularization amount

    param.N = [N(1) N(2) N(3)];%size(m3d(:,:,:,1,1));
    % param.num_chan = size(m3d,4);

    param.sens = rec;

    param.interpType = 'cubic';

    tic
    res = lsqr(@apply_sense_tikc3d_TR_v0, cat(1, kspace_coils(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);
    toc

    Res = reshape(res, param.N);

    img_sense_mc(:,:,:,t) = Res;



    % save([data_path, 'img_mcsense_lambda1Eneg3_iter50.mat'], 'img_sense_mc', '-v7.3')
    % %

    for t = 1:5
        imagesc3d2( img_sense_mc(:,:,:,t), s(img_sense_mc)/2, 50+t, [180,180,180], [0,2e-4],1,'img sense moco lambdaEneg2 iter20'), setGcf(.5)
        %    imagesc3d2( angle(img_sense_mc(:,:,:,t)), s(img_sense_mc)/2, t+60, [180,180,180], [-pi,pi]), setGcf(.5)
    end
    toc
end

%--------------------------------------------------------------------------
%% 3d-Sense recon with wavelet and TV regularization --> with moco
%--------------------------------------------------------------------------

% check if wavelet toolbox exists -> use only TV if not present
v = ver;
has_wavelet = any(strcmp(cellstr(char(v.Name)), 'Wavelet Toolbox'));

%mtx_size = size(m3d(:,:,:,1,1));
mtx_size = [N(1), N(2), N(3)];
size_powerof2 = nextpow2( mtx_size );

Nwav = 2.^size_powerof2;        % zero padded mtx size for wavelet transform to the next power of 2

param = init;

% 3d-wavelet parameters:
param.wav_scale = min(size_powerof2);
param.wav_type = 'db1';
param.wav_mode = 'per';


param.TV = TV3D;

param.TVWeight = 1e-6;                                  % 3D TV penalty
param.WavWeight = has_wavelet * param.TVWeight;         % 3D Wavelet penalty: set to 0 if wavelet toolbox not present. use same lambda as TV


param.num_chan = num_chan;
param.wav_pad = (Nwav - mtx_size) / 2;


param.Receive = rec;
param.Ct = conj(param.Receive);


%param.Itnlim = 5;           % num max inner iters
param.Itnlim = 10;           % num max inner iters
%num_outer_iters = 2;        % num outer iters
num_outer_iters = 5;        % num outer iters

param.tol = 1e-2;           % tolerance to terminate inner loop
param.pNorm = 1;            % use L1 constraint

param.nTR = nTR_abbr;
param.tform = tform_matrix;
param.tformInv = tform_matrixInv;

param.interpType = 'cubic';
param.N = mtx_size;

img_cs_mc = zeros(N(1), N(2), N(3), num_eco);


for t = 4%1:5
    tic
    disp(['3d-Sense recon with wavelet and TV regularization with moco contrast: ',num2str(t)]);

    % sort kspace to be used for reconstruction
    kspace2use = single(zeros([N(1)*os_factor N(2) N(3) num_chan 1 nTR_abbr])); % only 1 contrast

    for time_frame = 1:nTR_abbr
        kspace_tmp = single(zeros([N(1)*os_factor N(2) N(3) 1 ncoil]));
        kspace_tmp_svd = single(zeros([N(1)*os_factor N(2) N(3) num_chan 1]));
        for TR=((time_frame-1)*window +1):(time_frame*window)
            for contrast=t
                index_start = (contrast-1)*nETL + (TR-1)*step_size +1;

                % for index=(index_start):(index_start+nETL-1)
                for index=(index_start+7):(index_start+nETL-1)
                    ky = traj_y(index);
                    kz = traj_z(index);
                    %kspace_tmp(:,ky,kz,contrast-2,:) = data_unsorted(:,:,index);
                    kspace_tmp(:,ky,kz,1,:) = data_unsorted_img(:,:,index);
                end
            end
        end
        kspace_tmp = permute(kspace_tmp, [1,2,3,5,4]);

        kspace_tmp_svd(:,:,:,:,1) = single(svd_apply3d(kspace_tmp(:,:,:,:,1), cmp_mtx));

        kspace2use(:,:,:,:,:,time_frame) = kspace_tmp_svd;
    end


    m3d_tmp = sq(kspace2use(:,:,:,1,:,:));
    param.M3d = reshape(m3d_tmp,[N(1), N(2), N(3),1,nTR_abbr])~= 0;  %dim N(1) N(2) N(3) 1 nTR

    idx = 1;

    kspace_coils = zeros([sum(param.M3d(:)) * param.num_chan, 1]);

    for tr = 1:param.nTR
        % sampling mask for the current TR
        msk = repmat(param.M3d(:,:,:,:,tr), [1,1,1,param.num_chan,1]);

        num_points_current_TR = sum(msk(:));

        temp = kspace2use(:,:,:,:,1,tr) .* msk;

        % concatenate all k-space data
        kspace_coils(idx:idx+num_points_current_TR-1) = temp( msk == 1 );

        idx = idx + num_points_current_TR;
    end

    param.data = kspace_coils;

    res = zeros(mtx_size);

    for n = 1:num_outer_iters
        res = fnlCg_ics_sense_TR_v0(res, param);

        imagesc3d2(res, mtx_size/2, t+20, [0,0,0], [0,1e-4], 0, ['CS recon iter: ', num2str(n)]), setGcf(.5)
    end

    img_cs_mc(:,:,:,t) = res;
    toc
end


%save([data_path, 'img_mccs_lambda1Eneg6_iter5x2.mat'], 'img_cs_mc', '-v7.3')


for t = 1:num_eco
    imagesc3d2( img_cs_mc(:,:,:,t), s(img_sense_mc)/2, 20+t, [180,180,180], [0,2e-4],1,'img sense wavelet+TV moco iter5'), setGcf(.5)
    %imagesc3d2( angle(img_cs_mc(:,:,:,t)), s(img_sense_mc)/2, 30+t, [180,180,180], [-pi,pi]), setGcf(.5)
end

save_file_name = [basic_file_path, '_moco.mat'];
save(save_file_name, 'img_sense_mc', 'img_cs_mc' ,'img_sense', 'tform_matrix', '-v7.3')