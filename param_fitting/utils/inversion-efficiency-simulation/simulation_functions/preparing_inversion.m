function event_blocks = preparing_inversion(dt,inversion_raw,pulse_time,flip)

Gz = @(t) 0;  %preparation module never turns on gradients during RF.

gamma = 42.577478518e6; %[Hz/T]

%~~~~~~~~~~~~~~~~~~~~~~~~
%Scaling Inversion Pulse
%~~~~~~~~~~~~~~~~~~~~~~~~
%-adiabatic pulse
time_grid = (dt:dt:pulse_time).';

pulse_scale = (2 * pi * flip / 360)/sum(inversion_raw * gamma * dt);
pulse_real  = real(pulse_scale*inversion_raw);
pulse_imag  = imag(pulse_scale*inversion_raw);

Bx = @(t) interp1(time_grid,pulse_real,t);
By = @(t) interp1(time_grid,pulse_imag,t);

event_blocks = {};

%-(1) inversion
event_blocks{1}.time_grid = time_grid;
event_blocks{1}.Bx        = Bx;
event_blocks{1}.By        = By;
event_blocks{1}.Gz        = Gz;
event_blocks{1}.pulseflag = 1;

end