function [ res, tflag ] = apply_sense_tikc3d_TR_v0( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp')
    % Transposed SENSE operator:
  
    Res = zeros(params.N);
    
    idx = 1;
    
    for t = 1:params.nTR
        % extract k-space data for the current TR
        msk = repmat(params.m3d(:,:,:,:,t), [1,1,1,params.num_chan,1]);

        num_points_current_TR = sum(msk(:));
     
        temp = zeros([params.N, params.num_chan]);
        
        temp(msk == 1) = in(idx:idx+num_points_current_TR-1);
        
        
        img_coils = ifft3call( temp .* msk );
        res = sum(img_coils .* conj(params.sens), 4);

        Res = Res + imwarp(res, params.tformInv{t}, params.interpType, 'OutputView', imref3d(params.N), 'FillValues', 0);
        
        idx = idx + num_points_current_TR;
    end
         
    b = in(idx:end);    

    res = Res(:) + sqrt(params.lambda) * b;
    
else
    fprintf('+')
    % Forward SENSE operator:
  
    img = reshape(in, params.N);
    
    kspace_coils = zeros([sum(params.m3d(:)) * params.num_chan, 1]);

    idx = 1;
    
    for t = 1:params.nTR
        img_reg = imwarp(img, params.tform{t}, params.interpType, 'OutputView', imref3d(params.N), 'FillValues', 0);

        img_coils = repmat(img_reg, [1,1,1,params.num_chan]) .* params.sens;
        
        % sampling mask for the current TR
        msk = repmat(params.m3d(:,:,:,:,t), [1,1,1,params.num_chan,1]);

        num_points_current_TR = sum(msk(:));
        
        temp = fft3call(img_coils) .* msk;
        
        % concatenate all k-space data 
        kspace_coils(idx:idx+num_points_current_TR-1) = temp( msk == 1 );
  
        idx = idx + num_points_current_TR;
    end
    
    res = cat(1, kspace_coils(:), sqrt(params.lambda)*in);
    
end

end

    