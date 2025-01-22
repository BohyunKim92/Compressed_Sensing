function [SNR]=runReferenceAlgo(cs_level, NUM_OF_ITERATIONS, f_VD, f_R, im0,im)
percentage_of_k_space = cs_level/NUM_OF_ITERATIONS;
epsilon1=0.1;
SNR = 0;

%% generate fully sampled k-spaced of follow-up (Y_full) and reference (Y_0) images
size_im=size(im,1);
dft_mat=dftmtx(size_im)/sqrt(size_im);
Y_0=dft_mat*im0*dft_mat';
Y_full=dft_mat*im*dft_mat';
W=Wavelet;
X0=W*im0;

inverse_wav_baseline_representation=(1./(1+abs(X0)));% used to update W1 later

%% Initialize pdf parameters
f_VD_SUM = sum(sum(f_VD));
f_R_SUM = sum(sum(f_R));


%% initializing parameters
W1=ones(size(X0));
W2=zeros(size(X0));
Y = zeros(size_im);
S = zeros(size_im); % sample matrix


%% make sure gamma isn't making the pdf 0 matrix
% Note that when we only use f_VD as a pdf (i.e. f_VD is ward or eldar pdf
% but f_ND is a zero matrix), setting gamma =0 will make the pdf 0 matrix.
% This can cause trouble in sampling so we prevent this

if f_VD_SUM == 0
    gamma=1;
else 
    gamma = 0;
end

%% Sample for first iteration

for j=1:NUM_OF_ITERATIONS
    
    
    %% Sampling in fourier k-space
    [S,Sj]=generate_samples_without_rep(f_R,f_VD,gamma,percentage_of_k_space,S); % updating the sample matrix
    Y=Y_full.*S; % updating Y with samples that has been taken
    
    
    %% use FISTA based solver to calculate optimal X
    h = waitbar(0,'Please wait...','Name',['Minimization solver, ',num2str(percentage_of_k_space*100),'% of data']);
    X=FISTA_based_solver(Y,dft_mat,W,X0,W1,W2,h);
    close(h);
    
    
    %% Compute weighting coeficients for next iteration
    reconstructed_image=abs((W'*X));
    diff_in_wav=abs(W*(reconstructed_image-im0));
    map_diff=diff_in_wav./(diff_in_wav+1);
    W1(map_diff>epsilon1)=1;
    W1(map_diff<=epsilon1)=inverse_wav_baseline_representation(map_diff<=epsilon1);  
    W2=1./(1+abs(reconstructed_image-im0));
    gamma=mean(W2(:));  
    
    
    %% Recover Image and calculate SNR
    im_adaptive=abs(W'*X);
    SNR = 20*log(norm(flipud(im'), 'fro')/norm(flipud(im_adaptive') - flipud(im'), 'fro'));
end
end