function [rf_phase, rf_inc] = addAcq_v0_contrast6(contrast, seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, ...
    stepsY, gxPre, gx, gxSpoil, delay_TRinner, gyPre,gyReph,gzPre,gzReph, useAdc)

if nargin < 20
    useAdc = 1;
end

for iY = 1:nETL

    
    % Calculate index for ky-kz look up table
    %index = iY + (iZ-1)*nETL;
    index = (nETL*6)*(iZ-1)+(contrast-1)*nETL+iY;

    % RF spoiling
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);       %increment RF phase

    % Excitation
    seq.addBlock(rf);

    % Delay included to attain the desired echo time
    %seq.addBlock(mr.makeDelay(delay_TE));

    % Encoding
    seq.addBlock(gxPre, ...
        mr.scaleGrad(gyPre,stepsY(index)), ...
        mr.scaleGrad(gzPre,stepsZ(index)));    % Gz, Gy blips, Gx pre-winder

    if useAdc
        seq.addBlock(gx, adc);                      % Gx readout
    else
        seq.addBlock(gx);                           % Gx readout
    end

    seq.addBlock(gxSpoil, ...
        mr.scaleGrad(gyReph,-stepsY(index)), ...
        mr.scaleGrad(gzReph,-stepsZ(index)));    % -Gz, -Gy blips, Gx spoiler

    seq.addBlock(mr.makeDelay(delay_TRinner));  % wait until desired echo spacing
end
end


