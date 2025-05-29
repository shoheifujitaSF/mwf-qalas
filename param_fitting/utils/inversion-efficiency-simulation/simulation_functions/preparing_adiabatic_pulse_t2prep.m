function event_blocks = preparing_adiabatic_pulse_t2prep(dt,adiabatic_raw,pulse_time_ha,pulse_time_ad,ha_flip,ad_flip)

Gz = @(t) 0;  %preparation module never turns on gradients during RF.

gamma = 42.577478518e6; %[Hz/T]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Scaling Adiabatic and Hard Pulses
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-adiabatic pulse
time_grid_ad = (dt:dt:pulse_time_ad).';

pulse_ad       = adiabatic_raw(:,1) .* exp(1i * adiabatic_raw(:,2));
pulse_scale_ad = (2 * pi * ad_flip / 360)/sum(pulse_ad * gamma * dt);
pulse_real_ad  = real(pulse_scale_ad*pulse_ad);
pulse_imag_ad  = imag(pulse_scale_ad*pulse_ad);

%-real and imaginary of hard pulses
time_grid_ha = (dt:dt:pulse_time_ha).';

pulse_ha       = ones(length(time_grid_ha),1);
pulse_scale_ha = (2 * pi * ha_flip / 360)/sum(pulse_ha * gamma * dt);
pulse_real_ha  = real(pulse_scale_ha*pulse_ha);
pulse_imag_ha  = imag(pulse_scale_ha*pulse_ha);

Bx_ad = @(t) interp1(time_grid_ad,pulse_real_ad,t);
By_ad = @(t) interp1(time_grid_ad,pulse_imag_ad,t);

Bx_ha = @(t) interp1(time_grid_ha,pulse_real_ha,t);
By_ha = @(t) interp1(time_grid_ha,pulse_imag_ha,t);

noBx = @(t) 0;
noBy = @(t) 0;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Hard-coding information about each event block.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
event_blocks = {};

%-(1) 90 (phase 0)
event_blocks{1}.time_grid = time_grid_ha;
event_blocks{1}.Bx        = Bx_ha;
event_blocks{1}.By        = By_ha;
event_blocks{1}.Gz        = Gz;
event_blocks{1}.pulseflag = 1;

%-(2) Relaxation between 90 (phase 0) and first adiabatic pulse
event_blocks{2}.time_grid = (dt:dt:(11179e-6 - (100e-6 + pulse_time_ha))).';
event_blocks{2}.Bx        = noBx;
event_blocks{2}.By        = noBy;
event_blocks{2}.Gz        = Gz;
event_blocks{2}.pulseflag = 0;

%-(3) First Adiabatic Pulse
event_blocks{3}.time_grid = time_grid_ad;
event_blocks{3}.Bx        = Bx_ad;
event_blocks{3}.By        = By_ad;
event_blocks{3}.Gz        = Gz;
event_blocks{3}.pulseflag = 1;

%-(4) Relaxation between first and second adiabatic pulses
event_blocks{4}.time_grid = (dt:dt:(36179e-6 - (11179e-6 + pulse_time_ad))).';
event_blocks{4}.Bx        = noBx;
event_blocks{4}.By        = noBy;
event_blocks{4}.Gz        = Gz;
event_blocks{4}.pulseflag = 0;

%-(5) Second Adiabatic Pulse
event_blocks{5}.time_grid = time_grid_ad;
event_blocks{5}.Bx        = Bx_ad;
event_blocks{5}.By        = By_ad;
event_blocks{5}.Gz        = Gz;
event_blocks{5}.pulseflag = 1;

%-(6) Relaxation between second and third adiabatic pulses
event_blocks{6}.time_grid = (dt:dt:(61179e-6 - (36179e-6 + pulse_time_ad))).';
event_blocks{6}.Bx        = noBx;
event_blocks{6}.By        = noBy;
event_blocks{6}.Gz        = Gz;
event_blocks{6}.pulseflag = 0;

%-(7) Third Adiabatic pulse
event_blocks{7}.time_grid = time_grid_ad;
event_blocks{7}.Bx        = Bx_ad;
event_blocks{7}.By        = By_ad;
event_blocks{7}.Gz        = Gz;
event_blocks{7}.pulseflag = 1;

%-(8) Relaxation between third and fourth adiabatic pulses
event_blocks{8}.time_grid = (dt:dt:(86179e-6 - (61179e-6 + pulse_time_ad))).';
event_blocks{8}.Bx        = noBx;
event_blocks{8}.By        = noBy;
event_blocks{8}.Gz        = Gz;
event_blocks{8}.pulseflag = 0;

%-(9) Fourth Adiabatic pulse
event_blocks{9}.time_grid = time_grid_ad;
event_blocks{9}.Bx        = Bx_ad;
event_blocks{9}.By        = By_ad;
event_blocks{9}.Gz        = Gz;
event_blocks{9}.pulseflag = 1;

%-(10) Relaxation between fourth adiabitic pulse and hard 90 (-180 phase)
event_blocks{10}.time_grid = (dt:dt:(100100e-6 - (86179e-6 + pulse_time_ad))).';
event_blocks{10}.Bx        = noBx;
event_blocks{10}.By        = noBy;
event_blocks{10}.Gz        = Gz;
event_blocks{10}.pulseflag = 0;

%-(11) 90 (phase -180)
event_blocks{11}.time_grid = time_grid_ha;
event_blocks{11}.Bx        = @(t) -Bx_ha(t);
event_blocks{11}.By        = @(t) -By_ha(t);
event_blocks{11}.Gz        = Gz;
event_blocks{11}.pulseflag = 1;
end