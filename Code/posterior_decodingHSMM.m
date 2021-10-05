function [elnGamma]=posterior_decodingHSMM(Data,T,K,M,Mu,Cov,P,Pi,C,lambda)

p=size(Data,2);

elndur_p(win).window=NaN(K,T);
for st=1:K
	elndur_p(st,:)=eln(poisspdf(1:T,lambda(st)));
end
  
elnalpha(win).window=NaN(T,K); elnbeta(win).window=NaN(T,K); elnalpha_star(win).window=NaN(T,K); elnbeta_star(win).window=NaN(T,K);
elnGamma(win).window=NaN(T,K);
elnU=NaN(T,K,T);
   
k1=(2*pi)^(-p/2);

iCov=[];k2=[];
  
    
    
for i=1:T
	for l=1:K
          	ss=0;
          	for m=1:M
			d=Mu(m).mixture(l,:)-Data_loop(win).window(i,:);
			if Cov(m).mixture(l)==0
				Cov(m).mixture(l)=1e-40;
			end
			iCov=inv(Cov(m).mixture(l)); k2=k1/sqrt(det( Cov(m).mixture(l) ));
			tempr(win).window(m).mixture(i,l)=C(l,m)*k2*exp(-0.5*d*iCov*d');
			ss=ss+tempr(win).window(m).mixture(i,l);
   		end
          	elnB(win).window(i,l)=eln(ss);
      	end
end
    
%scaling B
for i=1:K
	normalizer=NaN;
        for t=1:T
        	normalizer=elnsum(normalizer,elnB(win).window(t,i));
        end
        for t=1:T
        	elnB(win).window(t,i)=elnproduct(elnB(win).window(t,i),-normalizer);
        end
end
    

%forward

elnalpha_star(win).window(1,:)=eln(Pi);
for t=1:T
	for j=1:K
        	for d=1:T
                	if t-d+1>0
                    		if d==1
                    			elnU(t,j,d)=elnB(win).window(t,j);
                    		else
                    			elnU(t,j,d)=elnproduct(elnU(t-1,j,d-1),elnB(win).window(t,j));
                    		end 
                		elnalpha(win).window(t,j)=elnsum(elnalpha(win).window(t,j), elnproduct(elnalpha_star(win).window(t-d+1,j), elnproduct(elndur_p(win).window(j,d),elnU(t,j,d)) ) );
                	end
            	end
        end

    	elnScale(t)=NaN;
   	%Scaling for alpha
   	for st=1:K
        	elnScale(win).window(t)=elnsum(elnScale(win).window(t),elnalpha(win).window(t,st));
   	end
   	for st=1:K
    		elnalpha(win).window(t,st)=elnproduct(elnalpha(win).window(t,st),-elnScale(win).window(t));
   	end
   
      	if t<=T-1
      		for j=1:K
          		jj_t=[];
          		for jj=1:K
          			jj_t=[jj_t; elnproduct(elnalpha(win).window(t,jj),eln(P(jj,j)))] ;
          		end
          		elnalpha_star(win).window(t+1,j)=NaN;
          		for jj=1:K
              			elnalpha_star(win).window(t+1,j)=elnsum(elnalpha_star(win).window(t+1,j),jj_t(jj));
          		end
          
      		end
      	end
      
end

% Backward run
for st=1:K
	elnbeta(win).window(T,st)=elnproduct(1,-elnScale(win).window(T));
end

for t=T-1:-1:0
	for j=1:K
        	elnss=NaN;
        	for d=1:T
            		if t+d<=T
                		elnss=elnsum( elnss, elnproduct( elndur_p(win).window(j,d), elnproduct(elnU(t+d,j,d), elnbeta(win).window(t+d,j))) );
            		end
        	end
        	elnbeta_star(win).window(t+1,j)=elnss;
    	end
    	if t~=T-1
		for j=1:K
			jj_t=[];
			for jj=1:K
				jj_t=[jj_t; elnproduct(elnbeta_star(win).window(t+2,jj),eln(P(j,jj)))] ;
			end
			elnbeta(win).window(t+1,j)=NaN;
			for jj=1:K
				elnbeta(win).window(t+1,j)=elnsum(elnbeta(win).window(t+1,j),jj_t(jj));
			end 


		end
    	end
end 

for time=1:T-1
	for st=1:K
    		elnbeta(win).window(time,st)=elnproduct(elnbeta(win).window(time,st),-elnScale(win).window(time));  
    	end
end

normalizer=NaN;
for time=1:T
	normalizer=NaN;

    	for st=1:K
        	elnGamma(win).window(time,st)=elnproduct(elnalpha(win).window(time,st),elnbeta(win).window(time,st));
    	end
    
    	for st=1:K
        	normalizer=elnsum(normalizer,elnGamma(win).window(time,st));
    	end
    
    	for st=1:K
        	elnGamma(win).window(time,st)=elnproduct(elnGamma(win).window(time,st),-normalizer);
    	end

        
end


    
end













    
    
