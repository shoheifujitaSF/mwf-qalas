function [ res ] = svd_apply3d( in, comp_mtx )
% svd coil compression for 3d data
% assumes that coil axis is the 4th dimension

mtx_size = size(in(:,:,:,1));
 
if size(in,3)==1
    mtx_size = [mtx_size, 1];
end

temp = reshape(in, prod(mtx_size), []);

res = reshape(temp * comp_mtx, [mtx_size, size(comp_mtx,2)]);

end