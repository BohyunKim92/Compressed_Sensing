function [already_sampled,sampling_matrix]=generate_samples_without_rep(f_R,f_VD,gamma,percentage_to_sample,already_sampled)
% This function generates a sampling matrix in k-space using a pdf where the pdf is
% a mixture of f_VD and f_R. The weight of this is determined by gamma
% value where 0<=gamma<=1.

rng('shuffle');
pdf=gamma*f_R+(1-gamma)*f_VD;
r_mat = zeros(size(f_R));


%% make sure r_mat is nonzero matrix. If r_mat is a zero matrix, we do randomizing again 
while (sum(sum(r_mat ~= zeros(size(f_R))))==0) 
rng('shuffle') % using a new random matrix for every iteration
r_mat=rand(size(f_R));
end


%% using information from already_sampled data to exclude already sampled points
temp_sampling = ones(size(f_R));
temp_sampling(already_sampled>0)=0;
r_mat =(r_mat.*temp_sampling);

pdf2=(r_mat.*pdf);
pdf3=pdf2(:); % makes into a vector 
[~, b]=sort(pdf3);
b=flipud(b);
threshold_for_sampling=pdf3(b(round(percentage_to_sample*length(b))));
sampling_matrix=zeros(size(f_R));
sampling_matrix(pdf2>=threshold_for_sampling)=1;
already_sampled = sampling_matrix+already_sampled; % updating the already sampled matrix. This is to prevent sample the same space twice.
end