function Z=eexp(X)

Z=zeros(size(X));

Z(isnan(X))=0;
Z(~isnan(X))=exp(X(~isnan(X)));

end