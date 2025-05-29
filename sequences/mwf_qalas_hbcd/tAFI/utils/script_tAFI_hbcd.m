%--------------------------------------------------------------------------
%% Pulseq AFI
%--------------------------------------------------------------------------

close all; clear; clc;
addpath(genpath('/autofs/cluster/berkin/fujita/Pulseq'))

checkPNS        =   0;
checkFreq       =   0;
writeSeq        =   0;

%--------------------------------------------------------------------------
%% Define high-level parameters
%--------------------------------------------------------------------------

sys             =   mr.opts('MaxGrad', 64, 'GradUnit', 'mT/m', ...
                            'MaxSlew', 160, 'SlewUnit', 'T/m/s', ...
                            'rfRingdownTime', 20e-6, ...
                            'rfDeadTime', 100e-6, ...
                            'adcDeadTime', 10e-6, ...
                            'B0',2.89);

seq             =   mr.Sequence(sys);            % Create a new sequence object

traj            =   readmatrix('fov_228x176_msize_56x44_tl_4_ncal_0_spiral_4_acc_1.5x1.5_nTR237_nph1.txt'); %size(traj) = [kz_position, ky_position, index]

fov             =   [228 228 176]*1e-3;         % Define FOV and resolution
N               =   [56 56 44]; % 2x2x2

Nx              =   N(1);
Ny              =   N(2);
Nz              =   N(3);

flip_angle_exc      =   60;                % degrees
rfSpoilingInc   =   39;              % RF spoiling increment
dwell           =   54e-6;
Tread           =   dwell*Nx;
os_factor       =   1;                  % readout oversampling amount
Ndummy          =   30;

TE              =   2.5e-3;
TR1             =   19e-3;
TR2             =   95e-3;
etl             =   4;
nTR             =   floor(size(traj,1)/etl);

%--------------------------------------------------------------------------
%% Prepare sequence blocks
%--------------------------------------------------------------------------

% Create non-selective pulse
rf              =   mr.makeBlockPulse(flip_angle_exc*pi/180, sys, 'Duration', 0.2e-3);


% Define other gradients and ADC events
deltak          =   1 ./ fov;

gx              =   mr.makeTrapezoid('x',sys,'Amplitude',Nx*deltak(1)/Tread,'FlatTime',ceil(Tread/sys.gradRasterTime)*sys.gradRasterTime ,'system',sys);    % readout gradient
gxPre           =   mr.makeTrapezoid('x','Area',-gx.area/2, 'system',sys);        % Gx prewinder
gxSpoil         =   mr.makeTrapezoid('x','Area',gx.area, 'system',sys);         % Gx spoiler
gxflyback         =   mr.makeTrapezoid('x','Area',-gx.area,'system',sys);         % Gx flyback

adc             =   mr.makeAdc(Nx * os_factor,'Duration',Tread,'Delay',gx.riseTime,'system',sys);

% 1 mT/m is 42.58 kHz/m
% 1 mT/m is 10 G/m
% momentum of gslSp2_x is 20% of that of gslSp2_z
sp_momentum = 300;
gslSp1_x        =   mr.makeTrapezoid('x','Area',42.58*sp_momentum*2,'system',sys); %Amplitude in Hz/m
gslSp2_z        =   mr.makeTrapezoid('z','Area',42.58*5*sp_momentum*5/6,'system',sys); %Amplitude in Hz/m
gslSp2_x        =   mr.makeTrapezoid('x','Area',42.58*5*sp_momentum*1/6,'Duration',mr.calcDuration(gslSp2_z) ,'system',sys); %Amplitude in Hz/m


gyPre           =       mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxPre), 'system',sys);
gzPre           =       mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxPre), 'system',sys);

gyReph          =       mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxflyback), 'system',sys);
gzReph          =       mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxflyback), 'system',sys);

gySpoil          =       mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gslSp1_x), 'system',sys);
gzSpoil          =       mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gslSp1_x), 'system',sys);


stepsY          =       ((traj(:,2)-1)-Ny/2)/Ny*2;
stepsZ          =       ((traj(:,1)-1)-Nz/2)/Nz*2;



delayTE         =   ceil( (TE - mr.calcDuration(rf) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre) ...
                    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR1        =  TR1 - mr.calcDuration(rf) - mr.calcDuration(gxPre) - delayTE - mr.calcDuration(gx)*etl- mr.calcDuration(gxflyback)*(etl-1);
delayTR2        =  TR2 - mr.calcDuration(rf) - mr.calcDuration(gxPre) - delayTE - mr.calcDuration(gx)*etl- mr.calcDuration(gxflyback)*(etl-1);


if delayTE < 0
    disp('error: delay TE is negative')
else
    disp(['delay TE: ', num2str(1e3*delayTE), ' ms'])
end

%--------------------------------------------------------------------------
%% Build sequence
%--------------------------------------------------------------------------
% 
% % Drive magnetization to steady state
% for iY = 1:Ndummy
% 
%     % Readout 1
%     seq.addBlock(rf);
%     seq.addBlock(gxPre); 
%     seq.addBlock(mr.makeDelay(delayTE));  
%     seq.addBlock(gx);
%     seq.addBlock(gslSp1_x, mr.makeDelay(delayTR1));
% 
%     % Readout 2
%     seq.addBlock(rf);
%     seq.addBlock(gxPre); 
%     seq.addBlock(mr.makeDelay(delayTE)); 
%     seq.addBlock(gx);
%     seq.addBlock(gslSp2_x, mr.makeDelay(delayTR2));
% end
% 



% actual imaging
counter=1;
for index = 1:etl:floor(size(traj,1)/etl)*etl


    % Readout 1
    tmp = 2*counter -1;
    rf.phaseOffset = mod(rfSpoilingInc*(tmp^2+tmp+2)*pi/180,2*pi);
    adc.phaseOffset = rf.phaseOffset;

    for i=1:etl
        if i==1
            seq.addBlock(rf);
            seq.addBlock(mr.makeDelay(delayTE));            % delay until readout

            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre,-stepsY(index)), ...
                mr.scaleGrad(gzPre,-stepsZ(index)));    % Gz, Gy blips, Gx pre-winder

            seq.addBlock(gx, adc);

        else
            gyPre_blip   =   mr.scaleGrad(gyReph,-stepsY(index+i-1) - -stepsY(index+i-2));
            gzPre_blip   =   mr.scaleGrad(gzReph,-stepsZ(index+i-1) - -stepsZ(index+i-2));

            seq.addBlock(gxflyback, gyPre_blip, gzPre_blip);
            seq.addBlock(gx, adc);
        end
    end

    seq.addBlock(gslSp1_x, mr.scaleGrad(gySpoil,stepsY(index+etl-1)), mr.scaleGrad(gzSpoil,stepsZ(index+etl-1)), mr.makeDelay(delayTR1));            % delay until TR1

    % Readout 2
    tmp = 2*counter;
    rf.phaseOffset = mod(rfSpoilingInc*(tmp^2+tmp+2)*pi/180,2*pi);
    adc.phaseOffset = rf.phaseOffset;


    for i=1:etl
        if i==1
            seq.addBlock(rf);
            seq.addBlock(mr.makeDelay(delayTE));            % delay until readout

            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre,-stepsY(index)), ...
                mr.scaleGrad(gzPre,-stepsZ(index)));    % Gz, Gy blips, Gx pre-winder

            seq.addBlock(gx, adc);

        else
            gyPre_blip   =   mr.scaleGrad(gyReph,-stepsY(index+i-1) - -stepsY(index+i-2));
            gzPre_blip   =   mr.scaleGrad(gzReph,-stepsZ(index+i-1) - -stepsZ(index+i-2));

            seq.addBlock(gxflyback, gyPre_blip, gzPre_blip);
            seq.addBlock(gx, adc);
        end
    end

    seq.addBlock(gslSp2_x, mr.scaleGrad(gySpoil,stepsY(index+etl-1)), gslSp2_z, mr.makeDelay(delayTR2));     % wait until TR2

    counter = counter+1;
end


%--------------------------------------------------------------------------
%% Check timing, PNS, forbidden freq.
%--------------------------------------------------------------------------

% Check whether the timing of the sequence is correct
[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% PNS check
if (checkPNS)
    [pns_ok, pns_n, pns_c, tpns] = seq.calcPNS('/autofs/cluster/berkin/fujita/matlab_codes/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
    %[pns_ok, pns_n, pns_c, tpns] = seq.calcPNS_C2('/autofs/cluster/berkin/fujita/matlab_codes/MP_GradSys_P034_c_CX600.asc'); % Connectom.X

    if (pns_ok)
        fprintf('PNS check passed successfully\n');
    else
        fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
    end
end

% Mechanical resonance check
if (checkFreq)
    gradSpectrum
end

% Plot
seq.plot('TimeRange',[2 3],'timeDisp','ms');

% Write to pulseq file with definitions for recon
seq.setDefinition('FOV', fov);
seq.setDefinition('Matrix', N);
seq.setDefinition('os_factor', os_factor);

filename = strrep(mfilename, 'script_', '');
seq.write([filename,'.seq']);
