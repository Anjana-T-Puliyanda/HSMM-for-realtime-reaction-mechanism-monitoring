function Z=eln(X)

[m,n]=size(X);
Z=zeros(m,n);
for i=1:m
    for j=1:n
        if X(i,j)>0
            Z(i,j)=log(X(i,j));
        elseif X(i,j)==0
            Z(i,j)=NaN;
        else
            disp('Negative error input');
        end
    end
end


end