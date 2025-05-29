function [ img ] = slice_uninterleave( img_int )
%SLICE_UNINTERLEAVE_ODD Summary of this function goes here
%   Detailed explanation goes here


img = zeros(size(img_int), 'single');

if mod(size(img,3),2)

    img(:,:,1:2:end,:,:) = img_int(:,:,1:(end+1)/2,:,:);
    img(:,:,2:2:end,:,:) = img_int(:,:,1+(end+1)/2:end,:,:);
 
else

    img(:,:,2:2:end,:,:) = img_int(:,:,1:end/2,:,:);
    img(:,:,1:2:end,:,:) = img_int(:,:,1+end/2:end,:,:);
    
end

end

