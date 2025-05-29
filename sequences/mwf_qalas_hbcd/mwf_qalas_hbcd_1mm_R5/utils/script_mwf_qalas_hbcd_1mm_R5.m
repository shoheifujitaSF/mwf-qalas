%--------------------------------------------------------------------------
%% Define high-level parameters
%--------------------------------------------------------------------------

addpath(genpath('/autofs/cluster/berkin/berkin/Matlab_Code_New/qalas_7T/pulseq-master/'));

gMax = 28;
sMax = 150;

% Set system parameters
sys             =   mr.opts('MaxGrad', gMax, 'GradUnit', 'mT/m', ...
                            'MaxSlew', sMax, 'SlewUnit', 'T/m/s', ...
                            'rfRingdownTime', 20e-6, ...
                            'rfDeadTime', 100e-6, ...
                            'adcDeadTime', 10e-6, ...
                            'B0',2.89);

%--------------------------------------------------------------------------
% Imaging sequence
%--------------------------------------------------------------------------


path2use = '/autofs/space/marduk_002/users/shohei/sharing/mwf_qalas_ismrm2025_hbcd/';

% Load ky-kz trajectory data %[kz_position ky_position index]

traj = readmatrix([path2use, 'SeqSmplPattern_CSGRE_228x176_R5_idx_nph6_v2_dc6_v3.csv']); 

nETL = 128;
nTR = size(traj,1)/nETL/6;

% Create sequence object and define FOV, resolution, and other high-level parameters
seq = mr.Sequence(sys);         % Create a new sequence object
fov=[228 228 176]*1e-3;         % Define FOV and resolution
N = [228 228 176];

Nx = N(1);
Ny = N(2);
Nz = N(3);

flip_angle = 4;                 % degrees
rfSpoilingInc=117;              % RF spoiling increment
dwell = 14e-6;                  % 320
Tread = dwell*Nx;               % bandwidth = 1 / Tread Hz/pixel
os_factor = 1;                  % oversampling
deltak = 1 ./ fov;



%--------------------------------------------------------------------------
% Prepare calibration sequence blocks
%--------------------------------------------------------------------------

% Currently hard-coded
Ny_acs = 32;    %8mm
Nz_acs = 32;    %8mm
Ndummy_acs = 50;

% Calculate timing
TE_acs = 6e-3;
TR_acs = 13e-3;

flip_angle_acs = 15;        % degrees

gx_acs = mr.makeTrapezoid('x',sys,'Amplitude',Nx*deltak(1)/Tread,'FlatTime',ceil(Tread/sys.gradRasterTime)*sys.gradRasterTime ,'system',sys);    % readout gradient
gxPre_acs = mr.makeTrapezoid('x',sys,'Area',-gx_acs.area/2);        % Gx prewinder
gxSpoil_acs = mr.makeTrapezoid('x',sys,'Area',gx_acs.area);         % Gx spoiler

rf_acs = mr.makeBlockPulse(flip_angle_acs*pi/180, sys, 'Duration', 0.1e-3);
adc_acs = mr.makeAdc(Nx * os_factor,'Duration',Tread,'Delay',gx_acs.riseTime,'system',sys);

delayTE_acs = ceil( (TE_acs - mr.calcDuration(rf_acs) + mr.calcRfCenter(rf_acs) + rf_acs.delay - mr.calcDuration(gxPre_acs)  ...
    - mr.calcDuration(gx_acs)/2)/seq.gradRasterTime)*seq.gradRasterTime;

if delayTE_acs < 0
    disp('error: acs delay TE is negative')    
else
    disp(['acs delay TE: ', num2str(1e3*delayTE_acs), ' ms'])
end


delayTR_acs = ceil((TR_acs - mr.calcDuration(rf_acs) - mr.calcDuration(gxPre_acs) ...
    - mr.calcDuration(gx_acs) - mr.calcDuration(gxSpoil_acs) - delayTE_acs(1))/seq.gradRasterTime)*seq.gradRasterTime;

if delayTR_acs < 0
    disp('error: acs delay TR is negative')    
else
    disp(['acs delay TR: ', num2str(delayTR_acs*1e3), ' ms'])
end

areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);


%--------------------------------------------------------------------------
% ACS data: dummies
%--------------------------------------------------------------------------

gpe1Pre_acs = mr.makeTrapezoid('y','Area',areaY(floor(Ny/2)),'Duration',mr.calcDuration(gxPre_acs));
gpe1Reph_acs = mr.makeTrapezoid('y','Area',-areaY(floor(Ny/2)),'Duration',mr.calcDuration(gxSpoil_acs));

% Drive magnetization to steady state
for iY = 1:Ndummy_acs

    rf_acs.phaseOffset = mod(117*(iY^2+iY+2)*pi/180,2*pi);  % add RF phase spoiling with 117 degrees
    seq.addBlock(rf_acs);

    % Gradients    
    seq.addBlock(gxPre_acs,gpe1Pre_acs);                  % add Gx pre-winder, go to desired ky
    seq.addBlock(mr.makeDelay(delayTE_acs));    % add delay needed before the start of readout
    
    seq.addBlock(gx_acs);                           % add readout Gx
    
    seq.addBlock(gpe1Reph_acs,gxSpoil_acs);               % add Gx spoiler, and go back to DC in ky
    seq.addBlock(mr.makeDelay(delayTR_acs))     % add delay to the end of TR
end


%--------------------------------------------------------------------------
% ACS data: Loop over phase encodes and define sequence blocks
%--------------------------------------------------------------------------

mask_acs = zeros([Ny,Nz]);

temp = 1:Ny;
iY_acs_indices = temp(1+end/2-Ny_acs/2:end/2+Ny_acs/2);

temp = 1:Nz;
iZ_acs_indices = temp(1+end/2-Nz_acs/2:end/2+Nz_acs/2);

% Make trapezoids for inner loop to save computation
for iY = iY_acs_indices
    gpe1Pre_acs(iY) = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',mr.calcDuration(gxPre_acs));
    gpe1Reph_acs(iY) = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',mr.calcDuration(gxSpoil_acs));
end

for iZ = iZ_acs_indices
    % Gz blips to go desired kz, and to come back to DC in kz
    gpe2Pre_acs = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',mr.calcDuration(gxPre_acs));        
    gpe2Reph_acs = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',mr.calcDuration(gxSpoil_acs));

    for iY = iY_acs_indices
        disp(['kZ: ', num2str(iZ), ' kY: ', num2str(iY)])
        mask_acs(iY,iZ) = 1;

        % RF spoiling
        rf_acs.phaseOffset = mod(117*(iY^2+iY+2)*pi/180,2*pi);
        adc_acs.phaseOffset = rf_acs.phaseOffset;

        % Excitation
        seq.addBlock(rf_acs);

        % Encoding
        seq.addBlock(gxPre_acs,gpe1Pre_acs(iY),gpe2Pre_acs);        % Gz, Gy blips, Gx pre-winder
        seq.addBlock(mr.makeDelay(delayTE_acs));    % delay until readout

        seq.addBlock(gx_acs,adc_acs);                       % readout
        
        seq.addBlock(gpe1Reph_acs(iY),gpe2Reph_acs,gxSpoil_acs);% -Gz, -Gy blips, Gx spoiler
        seq.addBlock(mr.makeDelay(delayTR_acs))     % wait until end of TR
    end
end


%--------------------------------------------------------------------------
%% Prepare QALAS sequence blocks
%--------------------------------------------------------------------------

deltak = 1 ./ fov;
gx = mr.makeTrapezoid('x',sys,'Amplitude',Nx*deltak(1)/Tread,'FlatTime',ceil(Tread/sys.gradRasterTime)*sys.gradRasterTime ,'system',sys);    % readout gradient
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2);        % Gx prewinder
gxSpoil = mr.makeTrapezoid('x',sys,'Area',gx.area);         % Gx spoiler

% Make trapezoids for inner loop to save computation

gyPre = mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxPre), 'system',sys);
gzPre = mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxPre), 'system',sys);

gyReph = mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxSpoil), 'system',sys);
gzReph = mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxSpoil), 'system',sys);

stepsY=((traj(:,2)-1)-Ny/2)/Ny*2;
stepsZ=((traj(:,1)-1)-Nz/2)/Nz*2;

% Create non-selective pulse
rf = mr.makeBlockPulse(flip_angle*pi/180, sys, 'Duration', 0.1e-3);

% Spoilers after t2prep and IR prep
gslSp_t2prep = mr.makeTrapezoid('z','Amplitude',-42.58*8*1e3,'Risetime',0.84e-3,'Duration',0.84e-3+8e-3+0.84e-3,'system',sys); %Amplitude in Hz/m
gslSp_IRprep = mr.makeTrapezoid('z','Amplitude',-42.58*8*1e3,'Risetime',1e-3,'Duration',1e-3+8e-3+1e-3,'system',sys); %Amplitude in Hz/m

adc = mr.makeAdc(Nx * os_factor,'Duration',Tread,'Delay',gx.riseTime,'system',sys);

% Prep pulses
% T2 prep and IR prep pulse are imported from external txt files
% txt file is in mag and phase, while mr.makeArbitraryRf assumes real and imag
t2prep = readmatrix([path2use, 'csv_MocoQalas/T2prep.txt']);
Re = t2prep(:,1) .* cos(t2prep(:,2));
Im = t2prep(:,1) .* sin(t2prep(:,2));
%t2prep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 360*pi/180, 'system',sys, 'dwell', 1e-6);
t2prep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 380.4*pi/180, 'system',sys, 'dwell', 1e-6);

text = readmatrix([path2use, 'csv_MocoQalas/rf90.txt']);
Re = text(:,1) .* cos(text(:,2));
Im = text(:,1) .* sin(text(:,2));
rf90 = mr.makeArbitraryRf((Re+Im*1i).', pi/2, 'system',sys, 'dwell', 1e-6);
rf90_180PhaseOffset = mr.makeArbitraryRf((Re+Im*1i).', pi/2, 'system',sys, 'dwell', 1e-6, 'PhaseOffset',-180*pi/180);

IRprep = readmatrix([path2use, 'csv_MocoQalas/invpulse.txt']);
Re = IRprep(:,1) .* cos(IRprep(:,2));
Im = IRprep(:,1) .* sin(IRprep(:,2));
IRprep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 750*pi/180, 'system',sys, 'dwell', 1e-5);
%IRprep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 1500*pi/180, 'system',sys, 'dwell', 1e-5);

%-------------------------------------------------------------------------
% Sequence timings to match IDEA implementation
%--------------------------------------------------------------------------

esp = 5.8e-3;
gap_between_readouts = 900e-3;

% delay_1_t2prep  =   11e-3 + 80e-6;      % 11.08 ms
% delay_2_t2prep  =   25e-3;              % 25 ms
% delay_3_t2prep  =   14e-3 - 80e-6;      % 13.92 ms


% TE20
delay_pre_t2prep_TE20 = 80e-3;
delay_1_t2prep_TE20  =   11e-3 + 80e-6 - 10e-3;      % 11.08 ms
delay_2_t2prep_TE20  =   5e-3;              % 25 ms
delay_3_t2prep_TE20  =   14e-3 - 80e-6 - 10e-3;      % 13.92 ms

% TE80
delay_pre_t2prep_TE80 = 20e-3;
delay_1_t2prep_TE80  =   11e-3 + 80e-6 - 2.5e-3;      % 11.08 ms
delay_2_t2prep_TE80  =   20e-3;              % 25 ms
delay_3_t2prep_TE80  =   14e-3 - 80e-6 - 2.5e-3;      % 13.92 ms


delay_IRprep    =   100e-3 - mr.calcDuration(IRprep_pulse)/2;         % gap between end of inversion and start of readout#2 
delay_TE        =   0;
delay_TRinner   =   esp - (mr.calcDuration(rf) + delay_TE + mr.calcDuration(gxPre)+mr.calcDuration(gx)+mr.calcDuration(gxSpoil));         
delay_TRouter   =   gap_between_readouts - esp*nETL;
delT_M3_M4      =   gap_between_readouts - esp*nETL - mr.calcDuration(IRprep_pulse) - delay_IRprep;     % between end of readout#1 and start of inversion
delT_M3_M4      =   delT_M3_M4 - 0.22e-3;   % this 0.22 ms is needed to match IDEA
delT_M13_2end   =   53.5e-3;


%--------------------------------------------------------------------------
%% QALAS data: dummies
%--------------------------------------------------------------------------

nDummies = 1;
useAdc = 0;

for iZ = 1:nDummies
    disp(['dummy ', num2str(iZ), '/', num2str(nDummies)])

    rf_phase=0;
    rf_inc=0;

    % T2 prep pulse TE 20
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE20));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE20));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    contrast = 1;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc,...
        rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph, useAdc);

    % T2 prep pulse TE 80
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE80));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE80));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    % Contrast 2
    contrast = 2;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc,...
        rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph, useAdc);
    seq.addBlock(mr.makeDelay(delT_M3_M4));

    % IR prep
    seq.addBlock(IRprep_pulse);
    seq.addBlock(gslSp_IRprep,mr.makeDelay(delay_IRprep));

    % Contrast 3
    contrast = 3;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc,...
        rfSpoilingInc, rf_phase, rf_inc,stepsZ,stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph, useAdc);
    seq.addBlock(mr.makeDelay(delay_TRouter))

    % Contrast 4
    contrast = 4;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, ...
        rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil,delay_TRinner, gyPre,gyReph,gzPre,gzReph, useAdc);
    seq.addBlock(mr.makeDelay(delay_TRouter));

    % Contrast 5
    contrast = 5;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, ...
        rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph, useAdc);
    seq.addBlock(mr.makeDelay(delay_TRouter));

    % Contrast 6
    contrast = 6;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, ...
        rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph, useAdc);
    seq.addBlock(mr.makeDelay(delT_M13_2end));

end


%--------------------------------------------------------------------------
% QALAS data
%--------------------------------------------------------------------------

for iZ = 1:nTR
    disp([num2str(iZ), '/', num2str(nTR)])

    rf_phase=0;
    rf_inc=0;

    % T2 prep pulse TE 20
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE20));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE20));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    contrast = 1;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph);


    % T2 prep pulse TE 80
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE80));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE80));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    % Contrast 2
    contrast = 2;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph);

    seq.addBlock(mr.makeDelay(delT_M3_M4));

    % IR prep
    seq.addBlock(IRprep_pulse);
    seq.addBlock(gslSp_IRprep,mr.makeDelay(delay_IRprep));

    % Contrast 3
    contrast = 3;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph);

    seq.addBlock(mr.makeDelay(delay_TRouter))

    % Contrast 4
    contrast = 4;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph);

    seq.addBlock(mr.makeDelay(delay_TRouter));

    % Contrast 5
    contrast = 5;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph);

    seq.addBlock(mr.makeDelay(delay_TRouter));

    % Contrast 6
    contrast = 6;
    [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ,  stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph);

    seq.addBlock(mr.makeDelay(delT_M13_2end));

end
%--------------------------------------------------------------------------
%% Check timing and write sequence
%--------------------------------------------------------------------------

% check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end


%--------------------------------------------------------------------------
%% Check PNS and forbidden frequency
%--------------------------------------------------------------------------

addpath /autofs/cluster/berkin/berkin/Matlab_Code_New/pulseq/jon_3d_epi/ISMRM-Virtual-Meeting--November-15-17-2023-main/tutorials/day2_SMS-EPI/safe_pns_prediction-master

[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('/autofs/cluster/berkin/fujita/matlab_codes/MP_GPA_K2309_2250V_951A_AS82_XA30A_mod.asc'); % prisma xa30

if (pns_ok)
     fprintf('PNS check passed successfully\n');
else
     fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end

% gradSpectrum

%--------------------------------------------------------------------------
%% Plot and write seq file
%--------------------------------------------------------------------------

% Set definitions for recon
seq.setDefinition('FOV', fov);
seq.setDefinition('Matrix', N);
seq.setDefinition('nETL', nETL);
seq.setDefinition('nTR', nTR);
seq.setDefinition('traj_y', traj(:,2)); %ky-kz sampling pattern saved for recon 
seq.setDefinition('traj_z', traj(:,1));
seq.setDefinition('os_factor', os_factor);

acs_duration = (Ndummy_acs + Ny_acs * Nz_acs) * TR_acs;


% plot
% seq.plot('TimeRange',[0 acs_duration],'timeDisp','ms');
% 
% seq.plot('TimeRange',[acs_duration acs_duration+4.5],'timeDisp','ms');
% 
% seq.plot('TimeRange',[acs_duration+4.5 acs_duration+9],'timeDisp','ms');
% 
% seq.plot('TimeRange',[acs_duration+13.5 acs_duration+18],'timeDisp','ms');

% Write to pulseq file
filename = strrep(mfilename, 'script_', '');
% filename = '/autofs/space/marduk_002/users/shohei/2024_10_XX_bay4_mesoscale/mwf_qalas_0p72_R7';

seq.write([filename,'.seq']);


%--------------------------------------------------------------------------
%% scan time
%--------------------------------------------------------------------------


disp(['acs duration ', num2str(acs_duration), ' sec'])

scan_duration = (nDummies + nTR) * 5.4;

disp(['scan duration ', num2str(scan_duration), ' sec'])


disp(['total duration ', num2str((scan_duration+acs_duration)/60), ' min'])
