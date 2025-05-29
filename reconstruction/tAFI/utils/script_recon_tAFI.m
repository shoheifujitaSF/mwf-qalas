%--------------------------------------------------------------------------
%% Recon Pulseq tAFI
%--------------------------------------------------------------------------

addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'));
addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code/TOOLBOXES/SENSE_LSQR_Toolbox'))
close all; clear; clc;

flag.save_nii               =       0;
flag.save_kspace            =       1;
msk_thres                   =       0;

%--------------------------------------------------------------------------
%% Data sorting and simple fft
%--------------------------------------------------------------------------

% data_file_path = '/autofs/space/marduk_002/users/shohei/2024_10_03_mwf_invivo/afi/meas_MID00336_FID75744_tAFI_nTR203_TR14.dat';
data_file_path = '/autofs/space/marduk_002/users/shohei/2024_10_03_mwf_invivo/afi/meas_MID00335_FID75743_tAFI_nTR203.dat';
traj            = readmatrix('/autofs/space/marduk_002/users/shohei/2024_10_03_mwf_invivo/afi/fov_240x208_msize_60x52_tl_4_ncal_0_spiral_4_acc_1.8x1.8_nTR203_nph1.txt');

[p,n,e]                     =       fileparts(data_file_path);
basic_file_path             =       fullfile(p,n);

twix_obj                    =       mapVBVD2(data_file_path);
data_unsorted               =       twix_obj{end}.image.unsorted();
[adc_len,ncoil,readouts]    =       size(data_unsorted);

% Read params from seq file
pulseq_file_path            =       [p, '/',regexprep(n, 'meas_MID\d+_FID\d+_', ''), '.seq'];
seq                         =       mr.Sequence();
seq.read(pulseq_file_path);
N                           =       seq.getDefinition('Matrix');
os_factor                   =       seq.getDefinition('os_factor');
etl                         =       4;


% Sort kspace data 
kspace                      =       zeros([N(1)*os_factor N(2) N(3) 2 ncoil]);


for idx = 1:etl:floor(size(traj,1)/etl)*etl
    index = (idx-1)*2;

    for i = 1:etl
    kspace(:,traj(idx+i-1, 2)  ,traj(idx+i-1, 1),  1,:) = data_unsorted(:,:,index+i);
    kspace(:,traj(idx+i-1, 2)  ,traj(idx+i-1, 1),  2,:) = data_unsorted(:,:,index+etl+i);
    end
end

S1      =       rsos(ifft3call(kspace(:,:,:,1,:)),5);
S2      =       rsos(ifft3call(kspace(:,:,:,2,:)),5);

%--------------------------------------------------------------------------
%% Fitting for B1 based on Yarnykh et al. MRM 2007
%--------------------------------------------------------------------------

r       =       S2 ./ S1;
n       =       5; % TR2/TR1 in the acquisition
msk     =       single(S1(:,:,:,1) > 0); %1e-3 for 4mm3291

theta   =       acos(min((r*n - 1)./(n-r),1)); % actual flip angle in radian

b1      =       rad2deg(real(theta))/60.*msk; % nominal FA = 60 degrees

imagesc3d2(S1 , s(S1)/2+[0 0 0], 1, 180+[0,0,0], [0,3e-3],1,'S1'), setGcf(.5)
imagesc3d2(S2 , s(S2)/2+[0 0 0], 2, 180+[0,0,0], [0,3e-3],1,'S2'), setGcf(.5)
imagesc3d2(b1 , s(b1)/2+[0 0 0], 3, 180+[0,0,0], [0.5,1.3],1,'b1'), setGcf(.5), colormap jet

%--------------------------------------------------------------------------
%% patref scan
%--------------------------------------------------------------------------
disp('Ref scan');

ref = sq(kspace(:,:,:,1,:)); % use TR1 for coil sens

num_acs = 24;



ref = ref(:,end/2-num_acs/2+1:end/2+num_acs/2,end/2-num_acs/2+1:end/2+num_acs/2,:);

img_ref = ifft3call(ref);
imagesc3d2(rsos(img_ref,4), s(img_ref)/2, 1, [0,0,0], [-0,5e-3]), setGcf(.5)

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


%--------------------------------------------------------------------------
% interpolate patref by zero padding to the high res matrix size
%--------------------------------------------------------------------------

size_data = size(kspace_svd(:,:,:,1,1));
size_patref = size(ref_svd(:,:,:,1,1));

patref_pad = padarray( ref_svd, [size_data-size_patref, 0, 0, 0]/2 );

img_patref_pad = ifft3c(patref_pad);

imagesc3d2( rsos(img_patref_pad,4), s(img_patref_pad)/2, 10, 180+[0,0,0], [0,2e-3],1,'img patref pad'), setGcf(.5)
imagesc3d2( rsos(ifft3call(kspace_svd(:,:,:,:,end)),4), s(img_patref_pad)/2, 11, 180+[0,0,0], [0,2e-3],1,'kspace svd ifft'), setGcf(.5)


%--------------------------------------------------------------------------
%% calculate sens map using ESPIRiT: parfor
%--------------------------------------------------------------------------
disp('Calculating sens map');

num_acs = min(size(kspace(:,:,:,1,1)));
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


mosaic(abs(sq(rec(end/2,:,:,:))), 2, 4, 2, '(simulated) coil sensitivities', [0,.5]), set(gcf, 'color', 'k')

%--------------------------------------------------------------------------
%% lsqr 3d-Sense recon [without CS]
%--------------------------------------------------------------------------

lsqr_iter = 30;        % max num of iterations
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
for t = 1:size(kspace_svd,5)
    disp(['eco: ', num2str(t)])

    param.m3d = m3d(:,:,:,:,t);         
    kspace_coils = kspace_svd(:,:,:,:,t);
    
    res = lsqr(@apply_sense_tikc3d, cat(1, kspace_coils(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  

    Res = reshape(res, param.N);
   
    img_sense(:,:,:,t) = Res;
end
toc

% imagesc3d2( img_sense(:,:,:,1), s(img_sense)/2, 11, 180+[0,0,0], [0,1e-2],1,'S1 SENSE'), setGcf(.5)
% imagesc3d2( img_sense(:,:,:,2), s(img_sense)/2, 12, 180+[0,0,0], [0,1e-2],1,'S2 SENSE'), setGcf(.5)


%--------------------------------------------------------------------------
%% Fitting for B1 based on Yarnykh et al. MRM 2007
%--------------------------------------------------------------------------

S1 = abs(img_sense(:,:,:,1));
S2 = abs(img_sense(:,:,:,2));

r       =       S2 ./ S1;
n       =       5; % TR2/TR1 in the acquisition
msk     =       single(S1(:,:,:,1) > msk_thres);

theta   =       acos(min((r*n - 1)./(n-r),1)); % actual flip angle in radian

b1      =       rad2deg(real(theta))/60.*msk; % nominal FA = 60 degrees

imagesc3d2(S1 , s(S1)/2+[0 0 0], 1, 180+[0,0,0], [0,2e-3],1,'S1 SENSE'), setGcf(.5)
imagesc3d2(S2 , s(S2)/2+[0 0 0], 2, 180+[0,0,0], [0,2e-3],1,'S2 SENSE'), setGcf(.5)
imagesc3d2(b1 , s(b1)/2+[5 0 0], 3, 180+[0,0,0], [0.5,1.3],1,'Pulseq AFI b1 SENSE'), setGcf(.5), colormap jet

