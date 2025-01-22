function X=FISTA_based_solver(Y,dft_mat,W,X0,W1,W2,h)
%inputs: Y: the samples of the k-space, W: Wavelet matrix, X0: the wavelet transform of the baseline
%image
%wavlet_mat=build_wavelet_matrix(size(Y,1),4);
numIter=20;
beta=0.8;
beta1=0.8;
beta2=0.8;
L=100;
lambda_bar=0.158;
lambda_bar_s=0.158;
mu_bar=0.025;
lambda=30;
lambda2=5;
mu=0.1;
t_k=1;
t_k_m1=1;
X=W*(dft_mat'*Y*dft_mat); % equivalent to W*im, wavelet transform of the original image
X_k_m1=X;
for i=1:numIter
     waitbar(i/numIter,h,['Iteration #',num2str(i),' of ',num2str(numIter)]);
    Z=X+((t_k_m1-1)/t_k)*(X-X_k_m1);
    temp_val1=W2.*(W'*(Z-X0));
    
    big_lambda=(abs(temp_val1)-lambda*mu);
    big_lambda=big_lambda./abs(temp_val1+eps).*temp_val1.*(big_lambda>0);
    
    temp_val2=W1.*Z;
    big_lambda2=(abs(temp_val2)-lambda2*mu);
    big_lambda2=big_lambda2./abs(temp_val2+eps).*temp_val2.*(big_lambda2>0);
    U=Z-(1/L)*(W*(dft_mat'*((Y~=0).*dft_mat*(W'*X)*dft_mat'-Y)*dft_mat)+(1/mu)*((W2.*(W*(temp_val1-big_lambda)))+W1.*(temp_val2-big_lambda2)));
    
    temp_val=(abs(U)-mu/L);
    
    X_k_m1=X;
    
    X=temp_val./abs(U+eps).*U.*(temp_val>0);
    temp_Y=dft_mat*(W'*X)*dft_mat';
    temp_Y(Y~=0)=Y(Y~=0);
    X=W*(dft_mat'*temp_Y*dft_mat);
    
    t_k_m1=t_k; % current t_k
    t_k=(1+sqrt(4*t_k*t_k+1))/2; % calculating t_k+1
    lambda=max(beta1*lambda,lambda_bar_s);
    mu=max(beta*mu,mu_bar);
    lambda2=max(beta2*lambda2,lambda_bar);
end



