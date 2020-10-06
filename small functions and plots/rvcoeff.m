function rvc = rvcoeff(X,Y)
%RVCOEFF calculated the RV coefficient as described in https://hal-cea.archives-ouvertes.fr/cea-00371054/file/Kherifetal_NeuroImage.pdf
%   X,Y are mXn containing m observations, which are R^n vectors

[m,n]=size(X);

%mean center data
X=X-mean(X,1);
Y=Y-mean(Y,1);

ZX=X*X';
ZY=Y*Y';

rvc=trace(ZX*ZY)/sqrt(trace(ZX^2)*trace(ZY^2));

end

