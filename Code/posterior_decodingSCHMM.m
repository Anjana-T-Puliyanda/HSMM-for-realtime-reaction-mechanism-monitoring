function [post_state_probab]=posterior_decodingSCHMM(Data,T,K,M,Mu,Cov,P,Pi,C,lambda)

[N,p]=size(Data);



elnalpha=NaN(T,K); elnbeta=NaN(T,K); elnalpha_star=NaN(T,K); elnbeta_star=NaN(T,K);
elnU=NaN(T,K,T);

for m=1:M
gamma(m).mixture=zeros(T,K);
tempr(m).mixture=zeros(T,K);
end

elnB=NaN(T,K);
k1=(2*pi)^(-p/2);

  
  %%%% FORWARD-BACKWARD 
  for m=1:M
      
      Gamma(m).mixture=[];
      Gammasum(m).mixture=zeros(1,K);
  end
  
  
  
  for n=1:floor(N/T)
      
      elndur_p=NaN(K,T);
      for st=1:K
        elndur_p(st,:)=eln(poisspdf(1:T,lambda(st)));
      end
  
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
elnbeta(T,:)=eln(ones(1,K));%./Scale(win).window(T);
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

    
normalizer=NaN;
for st=1:K
 elnGamma(1,st)=elnproduct(eln(Pi(st)),elnbeta_star(1,st));
 normalizer=elnsum( normalizer, elnGamma(1,st) );
end

for st=1:K
    elnGamma(1,st)=elnproduct(elnGamma(1,st),-normalizer);
end

for t=2:T
    normalizer=NaN;
    for i=1:K
     elnGamma(t,i)= eln( eexp(elnGamma(t-1,i))+ eexp(elnproduct(elnalpha_star(t,i),elnbeta_star(t,i)))- eexp(elnproduct(elnalpha(t-1,i),elnbeta(t-1,i)))  ) ;
     normalizer=elnsum( normalizer, elnGamma(t,i) );
    end

    for i=1:K
        elnGamma(t,i)=elnproduct(elnGamma(t,i),-normalizer);
    end 
    
end

%changed this part to evaluate over mixtures
    for m=1:M
        %ab=alpha(win).window.*beta(win).window;
        for i=1:K
            for t=1:T
                elngamma(m).mixture(t,i)=elnproduct(elnGamma(t,i), elnproduct(eln(tempr(m).mixture(t,i)), -elnB(t,i)) ); 
            end
        end
      gamma(m).mixture=elngamma(m).mixture;
      Gamma(m).mixture=[Gamma(m).mixture; gamma(m).mixture];
        
    end
    
    
    
    
    
  end
  
 post_state_probab=NaN(size(Gamma(1).mixture));
   
for m=1:M
    for tt=1:T
        for st=1:K
   post_state_probab(tt,st)=elnsum(post_state_probab(tt,st),Gamma(m).mixture(tt,st));
        end
    end
end

% post_state_probab=zeros(size(Gamma(1).mixture));
%    
% for m=1:M
%    post_state_probab=post_state_probab+Gamma(m).mixture;
% end
% 
% post_state_probab=rdiv(post_state_probab,rsum(post_state_probab));


end