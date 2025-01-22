%
% RunParallelSimulation
%
% Modified by: Bohyun Kim (hyunbrianna@gmail.com)
% 
% Department of Mathematics
% University of California Irvine -> Now University of Utah
% 
% Last revision: July 2025
%
% This code uses 32x32 image to simulate compressive sensing with
% mixture of variable density pdf and adaptive sampling pdf using
% the software ``Reference-Based MRI'' and POCS L1 minimization. 


clear; clc; close all;


tic;
%% Prompting users to choose data. You may choose portion of brain image or use phantom image
fprintf('Run case simulations.....\n');
whichdata = 1;


%% Define parameters for simulation
% Here, we consider two data set to test. For the case 1, we do simulation
% using 32x32 phantom image. For the case 2, we do simulation using 512x512
% brain image provided in Reference-Based MRI


switch whichdata


% For case1, we use 32x32 phantom image do cases simulation.
case 1 
load('phantom32x32.mat')
% define brain images
% We define what im0 and im are in the following section.
% im0 = reference image and im = target image. Note that we only use im0 to recover im when 
% we are using Eldar algorithm
im0 = X0;
im = X;

% Load various pdfs you may use.
% Load variable density pdfs
load('32X32reference_pdf.mat'); 
reference_fVD = pdf./sum(sum(pdf)); % load Reference-Based MRI variable density pdf in 32x32 form 

% creating a new fvd
W_fVD = zeros(size(im,1));
C = 1;
p = 0.3;

for j=1:32
    for k=1:32
        rj = j-16;
        rk = k-16;
        myk = 1/((rj^2+rk^2)^p);
        W_fVD(j,k) = min(C,myk);
    end
end
W_fVD = ifftshift(W_fVD); % load k-space variable density pdf in 32x32 form
W_fVD = W_fVD./sum(sum(W_fVD));

% Load adaptive sampling pdfs
load('DMNW_pdf.mat'); 
DMNW = adaptive_pdf; % load Dr. Deanna Needell's adaptive sampling pdf
DMNW = DMNW./sum(sum(DMNW));

size_im = size(im,1);
dft_mat=dftmtx(size_im)/sqrt(size_im);
Y_0=dft_mat*im0*dft_mat';
E_fND = abs(Y_0)./sum(sum(abs(Y_0))); % load Eldar's adaptive pdf in the newer direction. 



% For case2, we use 512x512 brain image from ``Reference Based MRI''  do cases simulation.
case 2
load('Brain_image.mat')
reference_fVD = pdf_vardens./sum(sum(pdf_vardens)); % variable density pdf in 512 form 

% creating a ward pdf
W_fVD = zeros(size(im,1));
C = 1;
p = 1;

for j=1:512
    for k=1:512
        rj = j-256;
        rk = k-256;
        myk = 1/((rj^2+rk^2)^p);
        W_fVD(j,k) = min(C,myk);
    end
end
W_fVD = ifftshift(W_fVD); % load Ward's variable density pdf in 32x32 form
W_fVD = W_fVD./sum(sum(W_fVD));

size_im = size(im,1);
dft_mat=dftmtx(size_im)/sqrt(size_im);
Y_0=dft_mat*im0*dft_mat';
E_fND = abs(Y_0)./sum(sum(abs(Y_0))); % construct adaptive pdf. 
DMNW = zeros(size_im);
end



%% initializing parameters

TOTAL_ITERATIONS = -1;
CS_LEVEL = 0.01;

while((TOTAL_ITERATIONS <1 || TOTAL_ITERATIONS > 30) && (CS_LEVEL < 0.03 || CS_LEVEL > 0.2))
    TOTAL_ITERATIONS = str2num(input(['How many iterations do you want to test? (TOTAL_ITERATIONS)'...
        '\nPlease type integers between 1 ~ 30 \n '],'s')); 
    CS_LEVEL = str2num(input(['What is the CS level of your data? (CS_LEVEL)'...
        '\nPlease type doubles in between 0.03 ~ 0.2 \n '],'s'));
    Num_of_core = str2num(input(['What is the number of core in your computer? (Num_of_core)'...
        '\nPlease type integers in between 2~4 \n '],'s'));
end
% TOTAL_ITERATIONS is how many times you are solving for L1 minimization problem 
% CS_LEVEL is the percentage of the samples you are getting
NUM_OF_ITERATIONS = 3; % This is the number of iterations in a single eldar's and eldar's newer algorithm


fprintf('\nParameter initialized! Here we start simulations.\n');
% From here we start case simulation. On the very top, we have defined pdf
% that we are using i.e. defined variable density and adaptive density.
% We also defined the algorithm that we are using. You can also check out 
% case_simulation_refence in braindump/figs to look up each cases.  

parpool('local',Num_of_core);


%% case 8. fvd = ward, fnd: DMNW. algorithm: Reference Based MRI
disp(strcat('simulating case 8 with cs level = ', num2str(CS_LEVEL)));
if sum(sum(DMNW))~=0
f_VD = W_fVD;
f_ND = DMNW;
SNR8 = zeros(TOTAL_ITERATIONS,1);

parfor i = 1: TOTAL_ITERATIONS
SNR = runReferenceAlgo(CS_LEVEL, NUM_OF_ITERATIONS, f_VD, f_ND, im0, im);
SNR8(i) = SNR;
end
MSNR8 = sum(SNR8)/TOTAL_ITERATIONS;

else
fprintf('\nDMNW is not defined since we are doing brain image simulation. Skipping this case...\n')
end


delete(gcp)
toc;