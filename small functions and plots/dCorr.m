function dist_correlation = dCorr(X,Y)
%DCORR calculates the distance correlation (1) as defined in (2).
%   X,Y are mXn containing m observations, which are R^n vectors. Function
%   returns the correlation between observations in X and Y.
%(1) G. J. Szekely, M. L. Rizzo, ´ et al., “Brownian distance covariance,”The annals of applied statistics, vol. 3, no. 4, pp. 1236–1265, 2009.
%(2) Generalized Multiple Correlation Coefficient as a Similarity
%Measurement between Trajectories, Julen Urain and Jan Peters, https://arxiv.org/pdf/1906.09802.pdf

[m,n]=size(X);

%calc dist matrices
a=squareform(pdist(X));
b=squareform(pdist(Y));

%double normelize
A=a-mean(a,1)-mean(a,2)+mean(a(:));
B=b-mean(b,1)-mean(b,2)+mean(b(:));

%calc varinces
varAsqr=sum(A.*A,'all')/m^2;
varBsqr=sum(B.*B,'all')/m^2;
CovSqr=sum(A.*B,'all')/m^2;

dist_correlation=sqrt(CovSqr/sqrt(varAsqr*varBsqr));
% scatter3(A(:,1),A(:,2),A(:,3),'r')
% hold on
% scatter3(B(:,1),B(:,2),B(:,3),'b')
% title('Normalized')

% figure
% scatter3(a(:,1),a(:,2),a(:,3),'r')
% hold on
% scatter3(b(:,1),b(:,2),b(:,3),'b')
% title('Unnormalized')

end

