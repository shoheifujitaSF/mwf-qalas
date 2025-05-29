%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------
clear; clc;

disp('Loading data');

addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'));
addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/pulseq/pulseq-develop_1.4.0'))
addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code/TOOLBOXES/SENSE_LSQR_Toolbox'))
addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/TOOLBOXES/CS_Wave_Toolbox'))
addpath '/autofs/cluster/berkin/berkin/Matlab_Code_New/NEATR/NEATR_Wip_936_Toolbox/neatr_wip_936_sense'

data_file_path = '/autofs/space/marduk_002/users/shohei/2024_11_03_BCH_pediatric_invivo/raw/meas_MID00023_FID38437_mwf_qalas_hbcd_1mm_R5.dat';


[p,n,e] = fileparts(data_file_path);
basic_file_path = fullfile(p,n);

twix_obj = mapVBVD2(data_file_path);
data_unsorted = twix_obj{end}.image.unsorted();
[adc_len,ncoil,readouts] = size(data_unsorted);

% Read params from seq file
pulseq_file_path = [p, '/',regexprep(n, 'meas_MID\d+_FID\d+_', ''), '.seq'];

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
%%  display zero filled images
%--------------------------------------------------------------------------

kspace = permute(kspace, [1 2 3 5 4]);
img = ifft3call(kspace);

img_rsos = sq(rsos(img, 4));

for t = 1:s(img_rsos, 4)
    imagesc3d2( img_rsos(:,:,:,t), s(img_rsos)/2, t, [0,0,0], [0,3e-4]), setGcf(.5)
end


%--------------------------------------------------------------------------
%% patref scan
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% coil compression
%--------------------------------------------------------------------------

num_chan = 20;  % num channels to compress to

[ref_svd, cmp_mtx] = svd_compress3d(ref, num_chan, 1);

rmse(rsos(ref_svd,4), rsos(ref,4))

N = size(kspace(:,:,:,1,1));
num_eco = size(kspace,5);

kspace_svd = zeross([N,num_chan,num_eco]);

for t = 1:size(kspace,5)
    kspace_svd(:,:,:,:,t) = svd_apply3d(kspace(:,:,:,:,t), cmp_mtx);
end

rmse(rsos(kspace_svd,4), rsos(kspace,4))


%--------------------------------------------------------------------------
%% interpolate patref by zero padding to the high res matrix size
%--------------------------------------------------------------------------


size_data = size(kspace_svd(:,:,:,1,1));
size_patref = size(ref_svd(:,:,:,1,1));

patref_pad = padarray( ref_svd, [size_data-size_patref, 0, 0, 0]/2 );

img_patref_pad = ifft3c(patref_pad);

imagesc3d2( rsos(img_patref_pad,4), s(img_patref_pad)/2, 10, [0,0,0], [0,2e-4])


%--------------------------------------------------------------------------
%% calculate sens map using ESPIRiT: parfor
%--------------------------------------------------------------------------

num_acs = min(size_patref);
kernel_size = [6,6];
eigen_thresh = 0.7;

receive = zeross(size(kspace_svd(:,:,:,:,1)));


delete(gcp('nocreate'))
c = parcluster('local');    

total_cores = c.NumWorkers;  
%parpool(ceil(total_cores/2))
parpool(ceil(total_cores/8))

    
tic
parfor slc_select = 1:s(img_patref_pad,1)     
    disp(num2str(slc_select))
    
    [maps, weights] = ecalib_soft( fft2c( sq(img_patref_pad(slc_select,:,:,:)) ), num_acs, kernel_size, eigen_thresh );

    receive(slc_select,:,:,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
end 
toc

delete(gcp('nocreate'))
 
% save([data_path, 'receive_svd_', num2str(num_chan), 'ch.mat'], 'receive', '-v7.3')


%--------------------------------------------------------------------------
%% lsqr 3d-Sense recon [without CS]
%--------------------------------------------------------------------------

lsqr_iter = 100;        % max num of iterations
lsqr_tol = 1e-3;        % tolerance to terminate iterations

m3d = kspace_svd~=0;
img_sense = zeros([size(kspace_svd,1), size(kspace_svd,2), size(kspace_svd,3), size(kspace_svd,5)]);
        

param = [];
param.lambda = 1e-3;        % L2 regularization amount

param.N = size(m3d(:,:,:,1,1));
param.num_chan = size(m3d,4);

rec = abs(receive) .* exp(1i * angle( receive .* repmat(conj(receive(:,:,:,1)), [1,1,1,param.num_chan]) ));
param.sens = rec;

tic
for t = 1:size(img,5)
    disp(['eco: ', num2str(t)])

    param.m3d = m3d(:,:,:,:,t);         
    kspace_coils = kspace_svd(:,:,:,:,t);
    
    res = lsqr(@apply_sense_tikc3d, cat(1, kspace_coils(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  

    Res = reshape(res, param.N);
   
    img_sense(:,:,:,t) = Res;
end
toc


%--------------------------------------------------------------------------
%% save as nifti
%--------------------------------------------------------------------------


for t = 1:size(img_sense,4)
    imagesc3d2( img_sense(:,:,:,t), s(img_sense)/2, t, [-0,180,180], [0,2e-4]), setGcf(.5)
    %imagesc3d2( angle(img_sense(:,:,:,t)), s(img_sense)/2, t+10, [-0,-0,0], [-pi,pi]), setGcf(.5)
end
 
%--------------------------------------------------------------------------
%% 3d-Sense recon with wavelet and TV regularization
%--------------------------------------------------------------------------

kspace_svd = kspace_svd ./ norm(kspace_svd(:));

% check if wavelet toolbox exists -> use only TV if not present
v = ver;
has_wavelet = any(strcmp(cellstr(char(v.Name)), 'Wavelet Toolbox'));
% has_wavelet = 0;
 
mtx_size = size(m3d(:,:,:,1,1));
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


param.num_chan = size(kspace_svd,4);
param.wav_pad = (Nwav - mtx_size) / 2;

% use channel 1 as phase reference to remove phase discontinuities across slices
rec = abs(receive) .* exp(1i * angle( receive .* repmat(conj(receive(:,:,:,1)), [1,1,1,param.num_chan]) ));

param.Receive = rec;
param.Ct = conj(param.Receive);


param.Itnlim = 10;          % num max inner iters
num_outer_iters = 5;        % num outer iters

param.tol = 1e-2;           % tolerance to terminate inner loop
param.pNorm = 1;            % use L1 constraint

img_cs = zeros([size(kspace_svd,1), size(kspace_svd,2), size(kspace_svd,3), size(kspace_svd,5)]);


tic
for t = 1:size(kspace,5)
    disp(['contrast: ', num2str(t)])
    
    param.data = kspace_svd(:,:,:,:,t);
    param.M3d = m3d(:,:,:,:,t);         

    res = zeros(mtx_size);

    for n = 1:num_outer_iters
        res = fnlCg_ics_sense(res, param);

        imagesc3d2(res, mtx_size/2, t+20, [0,0,0], [0,1e-4], 0, ['CS recon iter: ', num2str(n)]), setGcf(.5)
    end
    
    img_cs(:,:,:,t) = res;
end
toc
 

for t = 1:size(img_sense,4)
    % imagesc3d2( img_sense(:,:,:,t), s(img_sense)/2, t, [-0,-0,0], [0,2e-4]), setGcf(.5)
%     imagesc3d2( angle(img_sense(:,:,:,t)), s(img_sense)/2, t+10, [-0,-0,0], [-pi,pi]), setGcf(.5)

    imagesc3d2( img_cs(:,:,:,t), s(img_sense)/2+[0 0 0 0], t+100, [-0,-0,0], [0,3e-4]), setGcf(.5)
%     imagesc3d2( angle(img_cs(:,:,:,t)), s(img_sense)/2, t+110, [-0,-0,0], [-pi,pi]), setGcf(.5)
end


tiledlayout(1,size(img_sense,4))
for t = 1:size(img_sense,4)
    nexttile;
    % imshow(flip(abs(sq(img_cs(:,:,end/2-30,t))),1),[0,3e-4]);
    imshow(abs(sq(img_cs(end/2+20,:,:,t))),[0,3e-4]);
end


% save_file_name = [basic_file_path, '_nomoco.mat'];
% save(save_file_name, 'img_sense', 'img_cs' ,'receive','-v7.3')