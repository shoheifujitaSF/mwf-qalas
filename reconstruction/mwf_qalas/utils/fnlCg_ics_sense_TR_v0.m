function [x, con, wtv] = fnlCg_ics_sense_TR_v0(x0,params)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F* W' *x - y||^2 + lambda1*|x|_1 + lambda2*TV(W'*x) 
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% (c) Michael Lustig 2007
%
% 
% This function has been modified to perform CS-Wave reconstruction
% 
%-------------------------------------------------------------------------

x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
alpha = params.lineSearchAlpha; 
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;

% compute g0  = grad(Phi(x))
g0 = wGradient(x,params);
dx = -g0;

% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTx, FTdx, Gx, Gdx, Wx, Wdx] = preobjective(x, dx, params);
	f0 = objective(FTx, FTdx, Gx, Gdx, Wx, Wdx, 0, params);
	t = t0;
    
    [f1, con, wtv] = objective(FTx, FTdx, Gx, Gdx, Wx, Wdx, t, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 && (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, con, wtv] = objective(FTx, FTdx, Gx, Gdx, Wx, Wdx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
    end
    
    x_update = norm(t*dx(:)) / norm(x(:));
    disp(['Iter: ', num2str(k), '  Update: ', num2str(100*x_update), ' %', '  Obj: ', num2str(f1), '   Con: ', num2str(con), '   Reg: ', num2str(wtv)])
    
    if (x_update < params.tol) || (k == params.Itnlim)
        break
    end
    
	x = (x + t*dx);
 
    %conjugate gradient calculation    
	g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;	
end

end
 

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [FTx, FTdx, Gx, Gdx, Wx, Wdx] = preobjective(x, dx, params)
% precalculates transforms to make line search cheap


idx = 1;

FTx = zeros([sum(params.M3d(:)) * params.num_chan, 1]);
FTdx = zeros([sum(params.M3d(:)) * params.num_chan, 1]);

for t = 1:params.nTR
    img_reg_x = imwarp(x, params.tform{t}, params.interpType, 'OutputView', imref3d(size(x)), 'FillValues', 0);
    img_coils_x = repmat(img_reg_x, [1,1,1,params.num_chan]) .* params.Receive;

    
    img_reg_dx = imwarp(dx, params.tform{t}, params.interpType, 'OutputView', imref3d(size(x)), 'FillValues', 0);
    img_coils_dx = repmat(img_reg_dx, [1,1,1,params.num_chan]) .* params.Receive;

    
    % sampling mask for the current TR
    msk = repmat(params.M3d(:,:,:,:,t), [1,1,1,params.num_chan,1]);

    num_points_current_TR = sum(msk(:));

    
    temp = fft3call(img_coils_x) .* msk;
    FTx(idx:idx+num_points_current_TR-1) = temp( msk == 1 );


    temp = fft3call(img_coils_dx) .* msk;
    FTdx(idx:idx+num_points_current_TR-1) = temp( msk == 1 );

    idx = idx + num_points_current_TR;
end

    
if params.TVWeight
    Gx = params.TV*x;
    Gdx = params.TV*dx;
else
    Gx = 0;
    Gdx = 0;
end


if params.WavWeight
    Wx = wavedec3( padarray(x, params.wav_pad), params.wav_scale, params.wav_type, 'mode', params.wav_mode );
    Wdx = wavedec3( padarray(dx, params.wav_pad), params.wav_scale, params.wav_type, 'mode', params.wav_mode );
else
    Wx = 0;
    Wdx = 0;
end


end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [res, obj, WTV] = objective(Fx, Fdx, Gx, Gdx, Wx, Wdx, t, params)
%calculates the objective function

p = params.pNorm;

obj = Fx + t*Fdx - params.data;
obj = obj(:)'*obj(:);


TV = 0;
Wav = 0;


if params.TVWeight
    w = Gx(:) + t*Gdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
end


if params.WavWeight    
    for w = 1:length(Wx.dec)    
        wav = Wx.dec{w} + t * Wdx.dec{w};   
        Wav = Wav + sum( sum( sum( (wav.*conj(wav)+params.l1Smooth).^(p/2) ) ) );
    end
end

TV = sum(TV) * params.TVWeight;
 
Wav = Wav * params.WavWeight;

WTV = Wav + TV;

res = obj + WTV;

end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grad = wGradient(x,params)
gradTV = 0;
gradW = 0;

gradObj = gOBJ(x,params);
 
if params.TVWeight
    gradTV = gTV(x,params);
end

if params.WavWeight
    gradW = gWav(x,params);
end

grad = gradObj + params.TVWeight.*gradTV + params.WavWeight.*gradW;

end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

idx = 1;

FTx = zeros([sum(params.M3d(:)) * params.num_chan, 1]);

for t = 1:params.nTR
    img_reg_x = imwarp(x, params.tform{t}, params.interpType, 'OutputView', imref3d(size(x)), 'FillValues', 0);
    img_coils_x = repmat(img_reg_x, [1,1,1,params.num_chan]) .* params.Receive;

    % sampling mask for the current TR
    msk = repmat(params.M3d(:,:,:,:,t), [1,1,1,params.num_chan,1]);

    num_points_current_TR = sum(msk(:));
    
    temp = fft3call(img_coils_x) .* msk;
    FTx(idx:idx+num_points_current_TR-1) = temp( msk == 1 );

    idx = idx + num_points_current_TR;
end


kspace_coils = 2*(FTx - params.data);    

gradObj = zeros(params.N);
    
idx = 1;
    
for t = 1:params.nTR
    % extract k-space data for the current TR
    msk = repmat(params.M3d(:,:,:,:,t), [1,1,1,params.num_chan,1]);

    num_points_current_TR = sum(msk(:));
     
    temp = zeros([params.N, params.num_chan]);
        
    temp(msk == 1) = kspace_coils(idx:idx+num_points_current_TR-1);
        
        
    img_coils = ifft3call( temp .* msk );
    res = sum(img_coils .* params.Ct, 4);

    gradObj = gradObj + imwarp(res, params.tformInv{t}, params.interpType, 'OutputView', imref3d(params.N), 'FillValues', 0);
        
    idx = idx + num_points_current_TR;
end
         

end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grad = gTV(x,params)
% compute gradient of TV operator
p = params.pNorm;

Dx = params.TV*(x);

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);

grad = (params.TV'*G);

end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grad = gWav(x,params)
% compute gradient of wavelet operator
p = params.pNorm;

Dx = wavedec3( padarray(x, params.wav_pad), params.wav_scale, params.wav_type, 'mode', params.wav_mode );

dx = Dx;

for w = 1:length(Dx.dec)    
    dx.dec{w} = p*Dx.dec{w} .* (Dx.dec{w}.*conj(Dx.dec{w}) + params.l1Smooth).^(p/2-1);
end

grad = croparray(waverec3(dx), params.wav_pad);

end
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


