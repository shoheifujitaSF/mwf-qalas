function [PV1_map,PV2_map,PV3_map,PV4_map] = ...
    dict_fit_qalas_sub_pv_recon_4compartments(imgs,dict,pv)

img         = abs(imgs);
N           = [size(img,1),size(img,2),size(img,3)];

% length_dict = length(dict);
length_dict = length(pv);

msk         = zeross(N);

for slc_select = 1:N(3)
    thresh = rsos(img(:,:,slc_select,:),4);
    thresh = mean(thresh(:)).*0.5;
    msk(:,:,slc_select) = imerode(imfill(rsos(img(:,:,slc_select,:),4) > thresh, 'holes'), strel('disk', 3));
end

% w/o msk
msk         = oness(N); 


% b1
b1          = ones(size(msk));


thre_high   = 1.35;
thre_low    = 0.65;

% NO b1 correction
% sigma = 200;
% X = 1:N(2);
% Y = [1:N(1)]';
% G = exp(-1/(sigma^2)*((Y-N(1)/2).^2 + (X-N(2)/2).^2));
% G = (G-min(G(:)))/(max(G(:))-min(G(:)))*(thre_high-thre_low)+ thre_low;
% b1 = repmat(G,[1,1,N(3)]);
%

temp        = b1 .* msk;
temp(temp > thre_high)  = thre_high;
temp(temp < thre_low)   = thre_low;
img_b1      = temp .* msk;

num_b1_bins = size(dict,3);
b1_val      = linspace(min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins);
% b1_val      = 1; % no b1 correction

if length(b1_val) == 1
    msk_b1 = msk;
else
    msk_b1 = zeross([N,length(b1_val)]);
    
    for t = 1:length(b1_val)
        if t > 1
            msk_b1(:,:,:,t)   = (img_b1 <= b1_val(t)) .* (img_b1 > b1_val(t-1));
        else
            msk_b1(:,:,:,t)   = msk.*(img_b1 <= b1_val(t));
        end
        
        if t == length(b1_val)
            msk_b1(:,:,:,t)   = img_b1 > b1_val(t-1);
        end
        
        msk_b1(:,:,:,t)       = msk_b1(:,:,:,t) .* msk;
    end
end

% img_b1  = ones(size(msk));
% inv_eff = 1;
% b1_val  = 1;
% msk_b1  = msk;

PV1_map = zeross(N);
PV2_map = zeross(N);
PV3_map = zeross(N);
PV4_map = zeross(N);

for slc_select = 1:N(3)
    % tic
    % fprintf('slice: %d/%d... ',slc_select,N(3));
    
    for b = 1:length(b1_val)
        msk_slc = msk_b1(:,:,slc_select,b);
        num_vox = sum(msk_slc(:)~=0);
        
        if num_vox > 0
            img_slc = zeross([size(dict,2), num_vox]);
            
            for t = 1:size(dict,2)
                temp            = squeeze(img(:,:,slc_select,t));
                img_slc(t,:)    = temp(msk_slc~=0);
            end
            
            res             = dict(:,:,b) * img_slc;
            [~, max_idx]    = max(abs(res), [], 1);
            
            max_idx_pv      = mod(max_idx, length_dict);
            max_idx_pv(max_idx_pv==0) = length_dict;
            
            res_map         = pv(max_idx_pv,:);
            
            pv1_map = zeross(N(1:2));
            pv1_map(msk_slc==1) = res_map(:,1);
            
            pv2_map = zeross(N(1:2));
            pv2_map(msk_slc==1) = res_map(:,2);

            pv3_map = zeross(N(1:2));
            pv3_map(msk_slc==1) = res_map(:,3);

            pv4_map = zeross(N(1:2));
            pv4_map(msk_slc==1) = res_map(:,4);
            
            PV1_map(:,:,slc_select) = PV1_map(:,:,slc_select) + pv1_map;
            PV2_map(:,:,slc_select) = PV2_map(:,:,slc_select) + pv2_map;
            PV3_map(:,:,slc_select) = PV3_map(:,:,slc_select) + pv3_map;
            PV4_map(:,:,slc_select) = PV4_map(:,:,slc_select) + pv4_map;
        end
    end
    % toc
end

end