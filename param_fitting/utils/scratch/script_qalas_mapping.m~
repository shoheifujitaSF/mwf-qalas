%--------------------------------------------------------------------------
%% qalas bloch sim T1, T2, PD and inversion efficiency estimation
%% NOTE THAT THE SCRIPT REQUIRES FSL TO BE IN THE PATH 
%--------------------------------------------------------------------------

addpath utils

nii_folder = 'nii_input/';

% dicoms from 3d-Qalas c2p have been converted to nifti and json using dcm2niix:
qalas_nifti_fname = [nii_folder, 'dicoms_1p15mm_QALAS3d_1p15mm_R2_20220823120716_10.nii'];
qalas_json_fname = [nii_folder, 'dicoms_1p15mm_QALAS3d_1p15mm_R2_20220823120716_10.json'];


% dicoms from AFI B1 mapping sequence have been converted to nifti using dcm2niix:
% first nifti in afi series for registration onto qalas volume:
afi_img_nifti_fname = [nii_folder, 'dicoms_1p15mm_AFI3D_RS_Spoilx35_R3_20220823120716_18.nii'];

% third nifti in afi series is the B1 map:
afi_b1_nifti_fname = [nii_folder, 'dicoms_1p15mm_AFI3D_RS_Spoilx35_R3_20220823120716_20.nii'];


output_path = 'nii_output/';    % where to save the output maps

num_b1_bins = 25;               % number of discrete bins to use to quantize the b1 map (in the range [0.65, 1.35] at 3T)
inv_eff = 0.75:0.025:1;         % inversion efficiency range to use for the Bloch sim dictionary

estimate_pd = 1;                % set to 0 to avoid estimating proton density -> estimation becomes much faster
use_parfor = 1;                 % use parallel processing with parfor
num_workers = 4;                % number of parfor workers: recommended total_ram / 64 gb -> e.g. use 4 workers for 256 gb ram

use_bet = 1;                    % set to 0 if this is a phantom scan, where FSL BET will not work


tic
    fn_qalas_blochsim_map_v3( qalas_nifti_fname, qalas_json_fname, afi_img_nifti_fname, afi_b1_nifti_fname, output_path, num_b1_bins, inv_eff, estimate_pd, use_parfor, num_workers, use_bet );
toc
 


%--------------------------------------------------------------------------
%% load results
%--------------------------------------------------------------------------


t1 = load_nifti([output_path, 'T1_map.nii']);
t2 = load_nifti([output_path, 'T2_map.nii']);
ie = load_nifti([output_path, 'IE_map.nii']);
pd = load_nifti([output_path, 'PD_map.nii']);

qalas_rsos = load_nifti([output_path, 'qalas_rsos.nii.gz']);
afi = load_nifti([output_path, 'afi2qalas.nii.gz']);
b1 = load_nifti([output_path, 'b1_afi_flirt.nii.gz']);


imagesc3d2( t1.vol, s(t1.vol)/2, 1, [180,-0,90], [0,2500]),
imagesc3d2( t2.vol, s(t1.vol)/2, 2, [180,-0,90], [0,200]),
imagesc3d2( ie.vol, s(t1.vol)/2, 3, [180,-0,90], [0,1]),
imagesc3d2( pd.vol, s(t1.vol)/2, 4, [180,-0,90], [0,10000]),

imagesc3d2( qalas_rsos.vol, s(t1.vol)/2, 5, [180,-0,90], [0,400]), setGcf(.5)
imagesc3d2( afi.vol, s(t1.vol)/2, 6, [180,-0,90], [0,1000]), setGcf(.5)
imagesc3d2( b1.vol/60, s(t1.vol)/2, 7, [180,-0,90], [0.,1.35]), setGcf(.5), colormap jet


