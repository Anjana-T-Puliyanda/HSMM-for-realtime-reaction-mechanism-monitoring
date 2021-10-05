function [elnGamma]=posterior_decodingHSMM(Data,T,K,M,Mu,Cov,P,Pi,C,lambda)

p=size(Data,2);

elndur_p=NaN(K,T);
for st=1:K
	elndur_p(st,:)=eln(poisspdf(1:T,lambda(st)));
end
  
elnalpha=NaN(T,K); elnbeta=NaN(T,K); elnalpha_star=NaN(T,K); elnbeta_star=NaN(T,K);
elnGamma=NaN(T,K);
elnU=NaN(T,K,T);
   
k1=(2*pi)^(-p/2);

iCov=[];k2=[];
  
    
    
for i=1:T
	for l=1:K
          	ss=0;
          	for m=1:M
			d=Mu(m).mixture(l,:)-Data(i,:);
			if Cov(m).mixture(l)==0
				Cov(m).mixture(l)=1e-40;
			end
			iCov=inv(Cov(m).mixture(l)); k2=k1/sqrt(det( Cov(m).mixture(l) ));
			tempr(m).mixture(i,l)=C(l,m)*k2*exp(-0.5*d*iCov*d');
			ss=ss+tempr(m).mixture(i,l);
   		end
          	elnB(i,l)=eln(ss);
      	end
end
    
%scaling B
for i=1:K
	normalizer=NaN;
        for t=1:T
        	normalizer=elnsum(normalizer,elnB(t,i));
        end
        for t=1:T
        	elnB(t,i)=elnproduct(elnB(t,i),-normalizer);
        end
end
    

%forward

elnalpha_star(1,:)=eln(Pi);
for t=1:T
	for j=1:K
        	for d=1:T
                	if t-d+1>0
                    		if d==1
                    			elnU(t,j,d)=elnB(t,j);
                    		else
                    			elnU(t,j,d)=elnproduct(elnU(t-1,j,d-1),elnB(t,j));
                    		end 
                		elnalpha(t,j)=elnsum(elnalpha(t,j), elnproduct(elnalpha_star(t-d+1,j), elnproduct(elndur_p(j,d),elnU(t,j,d)) ) );
                	end
            	end
        end

    	elnScale(t)=NaN;
   	%Scaling for alpha
   	for st=1:K
        	elnScale(t)=elnsum(elnScale(t),elnalpha(t,st));
   	end
   	for st=1:K
    		elnalpha(t,st)=elnproduct(elnalpha(t,st),-elnScale(t));
   	end
   
      	if t<=T-1
      		for j=1:K
          		jj_t=[];
          		for jj=1:K
          			jj_t=[jj_t; elnproduct(elnalpha(t,jj),eln(P(jj,j)))] ;
          		end
          		elnalpha_star(t+1,j)=NaN;
          		for jj=1:K
              			elnalpha_star(t+1,j)=elnsum(elnalpha_star(t+1,j),jj_t(jj));
          		end
          
      		end
      	end
      
end

% Backward run
for st=1:K
	elnbeta(T,st)=elnproduct(1,-elnScale(T));
end

for t=T-1:-1:0
	for j=1:K
        	elnss=NaN;
        	for d=1:T
            		if t+d<=T
                		elnss=elnsum( elnss, elnproduct( elndur_p(j,d), elnproduct(elnU(t+d,j,d), elnbeta(t+d,j))) );
            		end
        	end
        	elnbeta_star(t+1,j)=elnss;
    	end
    	if t~=T-1
		for j=1:K
			jj_t=[];
			for jj=1:K
				jj_t=[jj_t; elnproduct(elnbeta_star(t+2,jj),eln(P(j,jj)))] ;
			end
			elnbeta(t+1,j)=NaN;
			for jj=1:K
				elnbeta(t+1,j)=elnsum(elnbeta(t+1,j),jj_t(jj));
			end 


		end
    	end
end 

for time=1:T-1
	for st=1:K
    		elnbeta(time,st)=elnproduct(elnbeta(time,st),-elnScale(time));  
    	end
end

normalizer=NaN;
for time=1:T
	normalizer=NaN;

    	for st=1:K
        	elnGamma(time,st)=elnproduct(elnalpha(time,st),elnbeta(time,st));
    	end
    
    	for st=1:K
        	normalizer=elnsum(normalizer,elnGamma(time,st));
    	end
    
    	for st=1:K
        	elnGamma(time,st)=elnproduct(elnGamma(time,st),-normalizer);
    	end

        
end


    
end













    
    
