function [P,LL,Mu,Cov,Pi,C,lambda]=HMM_BaumW(Data,K,M,Mu,Cov,P,Pi,C,lambda,cyc,tol,p)


num_windows=length(Data);
elngamma=[];gammasum=[];
Xi=[];elnScale=[];elnxi_prior=[];
tempr=[];elnB=[];elndur_p=[];
%cycle loop starts here
lik=0; LL=[];

Pcycle=[]; LLcycle=[]; Mucycle=[]; Covcycle=[]; Picycle=[]; Ccycle=[]; lambdacycle=[];

for cycle=1:cyc

%%%% FORWARD-BACKWARD 
  for win=1:num_windows

    
  T=size(Data(win).window,1);
  elndur_p(win).window=NaN(K,T);
  for st=1:K
    elndur_p(win).window(st,:)=eln(poisspdf(1:T,lambda(st)));
  end
  
  elnScale(win).window=NaN; 
  elnXsi(win).window=NaN(T-1,K*K);
  elnalpha(win).window=NaN(T,K); elnbeta(win).window=NaN(T,K); elnalpha_star(win).window=NaN(T,K); elnbeta_star(win).window=NaN(T,K);
  elnxi_prior(win).window=NaN(K,1); elnGamma(win).window=NaN(T,K);
  elnU=NaN(T,K,T);
   
  k1=(2*pi)^(-p/2);

    iCov=[];k2=[];
%     for m=1:M
%         for st=1:K
%         iCov(m).mixture(st).state=inv(Cov(m).mixture(st).state);
%         k2(m).mixture(st).state=k1/sqrt(det(Cov(m).mixture(st).state));
%         end
%     end    
    
    
    for i=1:T
      for l=1:K
          ss=0;
          for m=1:M
	d=Mu(m).mixture(l,:)-Data(win).window(i,:);
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
    
    %scale=zeros(T,1);
%     alpha(win).window(1,:)=Pi.*B(win).window(1,:);
%     Scale(win).window(1)=sum(alpha(win).window(1,:));
%     alpha(win).window(1,:)=alpha(win).window(1,:)/Scale(win).window(1);

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
%         elnScale(win).window(t)=NaN;
%    %Scaling for alpha
%    for st=1:K
%         elnScale(win).window(t)=elnsum(elnScale(win).window(t),elnalpha(win).window(t,st));
%    end
%     alpha(win).window(t,:)=alpha(win).window(t,:)/Scale(win).window(t);
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
elnbeta(win).window(T,:)=eln(ones(1,K));%./Scale(win).window(T);
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

%evaluating neta
normalizer=NaN;
for st=1:K
 elnGamma(win).window(1,st)=elnproduct(eln(Pi(st)),elnbeta_star(win).window(1,st));
 normalizer=elnsum( normalizer, elnGamma(win).window(1,st) );
end

for st=1:K
    elnGamma(win).window(1,st)=elnproduct(elnGamma(win).window(1,st),-normalizer);
end


for t=1:T-1
    tmp=zeros(K,K);
    for i=1:K
  
         for j=1:K
        tmp(i,j)=elnproduct( elnproduct(elnalpha(win).window(t,i), eln(P(i,j))), elnbeta_star(win).window(t+1,j)) ;
         end     
    end
    elnXsi(win).window(t,:)=reshape(tmp,1,K*K);
end

for t=2:T
    normalizer=NaN;
    for i=1:K
     elnGamma(win).window(t,i)= eln( eexp(elnGamma(win).window(t-1,i))+ eexp(elnproduct(elnalpha_star(win).window(t,i),elnbeta_star(win).window(t,i)))- eexp(elnproduct(elnalpha(win).window(t-1,i),elnbeta(win).window(t-1,i)))  ) ;
     normalizer=elnsum( normalizer, elnGamma(win).window(t,i) );
    end

    for i=1:K
        elnGamma(win).window(t,i)=elnproduct(elnGamma(win).window(t,i),-normalizer);
    end 
    
end
    
    %changed this part to evaluate over mixtures
    for m=1:M
        %ab=alpha(win).window.*beta(win).window;
        for i=1:K
            for t=1:T
                elngamma(win).window(m).mixture(t,i)=elnproduct(elnGamma(win).window(t,i), elnproduct(eln(tempr(win).window(m).mixture(t,i)), -elnB(win).window(t,i)) ); 
            end
        end
        %gammasum(win).window(m).mixture=sum(gamma(win).window(m).mixture);
    end
    
    %xi=zeros(T-1,K*K);
%     for i=1:T-1
%       t=P.*( alpha(win).window(i,:)' * (beta(win).window(i+1,:).*B(win).window(i+1,:)));
%       Xi(win).window(i,:)=t(:)'/sum(t(:));
%     end
    
    
    %Scale=Scale+log(scale);
%     Gamma=[Gamma; gamma];
%     Gammasum=Gammasum+gammasum;
 %   Xi=Xi+xi;
    
    
    %log the value of xi for prior in M step
    tmp=reshape(elnXsi(win).window(1,:),K,K);
    for st=1:K
        rsumval=NaN;
        for j=1:K
            rsumval=elnsum(rsumval,tmp(st,j));
        end
        elnxi_prior(win).window(st)=rsumval;
    end
    
    %evaluate the likelihood of the sequence
    for t=1:T
    for i=1:K
        for j=1:K
            for d=1:T
                if t-d>0
                    elnScale(win).window=elnsum(elnScale(win).window, elnproduct(elnproduct(elnproduct(elnproduct(elnalpha(win).window(t-d,i),eln(P(i,j))), elndur_p(win).window(j,d)),elnU(t,j,d)), elnbeta(win).window(t,j))     );
                end
            end
        end
    end
    end
    
    
    
  end

  
%%%% M STEP 
  
  % outputs
  for m=1:M
      num=NaN(K,p);denom=NaN(K,1);
  for win=1:num_windows
      T=size(Data(win).window,1);
%       nss=eexp(elngamma(win).window(m).mixture)'*Data(win).window;
      dataw=eln(Data(win).window);  
      nss=NaN(K,p); dss=NaN(K,1); 
      for st=1:K
          for pp=1:p
              for tt=1:T
                  nss(st,pp)=elnsum(nss(st,pp),elnproduct(elngamma(win).window(m).mixture(tt,st),dataw(tt,pp)) );
               end
          end
          
          for tt=1:T
              dss(st,1)=elnsum(dss(st,1),elngamma(win).window(m).mixture(tt,st) );
          end
      end
      
      
      for st=1:K
            for pp=1:p
                num(st,pp)=elnsum(num(st,pp), elnproduct(nss(st,pp),elnScale(win).window));
            end
            denom(st,1)=elnsum(denom(st,1), elnproduct(dss(st,1),elnScale(win).window) );
        end
  end
  
  for st=1:K
      for pp=1:p
          Mu(m).mixture(st,pp)=eexp(elnproduct(num(st,pp),-denom(st)));%rdiv(num,denom');
      end
  end
  end
  
 % C = bsxfun(@rdivide, A, B);

  % transition matrix 
  num=NaN(1,K*K); denom=NaN(K,1);
  for win=1:num_windows
   T=size(Data(win).window,1);
%    sxi=sxi+prod(eexp(elnScale(win).window(:))).*sum(eexp(elnXsi(win).window)); 
    tsum=NaN(1,K*K);
    for dm=1:K*K
        for t=1:T-1
            tsum(1,dm)=elnsum(tsum(1,dm),elnXsi(win).window(t,dm));
        end
    end
    
    tsum=reshape(tsum,K,K);
    rsum=NaN(K,1);
    for i=1:K
        for j=1:K
            rsum(i,1)=elnsum(rsum(i,1),tsum(i,j));
        end
    end
    
    tsum=reshape(tsum, 1, K*K);
   for dm=1:K*K
       num(1,dm)=elnsum(num(1,dm),elnproduct(tsum(1,dm),elnScale(win).window) );
   end
   
   for st=1:K
       denom(st,1)=elnsum(denom(st,1), elnproduct(rsum(st),elnScale(win).window) );
   end
    
  end
  
  sxi=NaN(K,K);
  num=reshape(num, K,K);
  for i=1:K
      for j=1:K
          sxi(i,j)=elnproduct(num(i,j),-denom(i));
      end
  end
 
  
%     for i=1:K
%     sxi(i,i) = 0;
%     end
for i=1:K
    for r=1:K
    if r==i
        sxi(i,r)=NaN;
    end
%     if r>i+1
%         sxi(i,r)=0;
%     end
%     
    end
end

rwsm=NaN(K,1);
for st=1:K
    for j=1:K
    rwsm(st,1)=elnsum(rwsm(st,1),sxi(st,j));
    end
   
end

for st=1:K
    for j=1:K
        P(st,j)=eexp(elnproduct(sxi(st,j),-rwsm(st,1)));
    end
end
% P=rdiv(sxi,rsum(sxi));
  
 % priors
  Pi=NaN(1,K);
 % Pi(1,1)=1;
%   for i=1:batches
%       for m=1:M
%         Pi=Pi+Gamma((i-1)*T+1,:);
%       end
%   end
for st=1:K
for win=1:num_windows
    Pi(1,st)=elnsum(Pi(1,st),elnproduct(elnxi_prior(win).window(st),elnScale(win).window));
end
end

normalizer=NaN;
for st=1:K
    normalizer=elnsum(normalizer,Pi(st));
end
for st=1:K
    Pi(st)=elnproduct(Pi(st),-normalizer);
end

Pi=eexp(Pi);

%   Pi=Pi/batches;
  
% covariance
   
  for m=1:M
      num=NaN(K,1); denom=NaN(K,1);
  for win=1:num_windows
      T=size(Data(win).window,1);nss=NaN(K,1); dss=NaN(K,1);
      for st=1:K
      for t=1:T
         nss(st,1)=elnsum(nss(st,1), elnproduct( elngamma(win).window(m).mixture(t,st),eln((Data(win).window(t,:)-Mu(m).mixture(st,:))*(Data(win).window(t,:)-Mu(m).mixture(st,:))') )   );
         dss(st,1)=elnsum(dss(st,1),elngamma(win).window(m).mixture(t,st) );
      end
      
      end
      
      for st=1:K
          num(st,1)=elnsum(num(st,1),elnproduct(nss(st,1),elnScale(win).window));
          denom(st,1)=elnsum(denom(st,1),elnproduct(dss(st,1),elnScale(win).window));
      end
      
  end
  
  for st=1:K
      Cov(m).mixture(st)=eexp(elnproduct(num(st),-denom(st)));
  end
  end
  
  %mixture weights of gaussian probabilities
           
 for st=1:K
     num=NaN(M,1); denom=NaN;
     for win=1:num_windows
        T=size(Data(win).window,1); 
        gssm=NaN(M,1);
        for m=1:M
            for t=1:T
                gssm(m,1)=elnsum(gssm(m,1),elngamma(win).window(m).mixture(t,st));
            end
        end
        
        dss=NaN;
        for m=1:M
            dss=elnsum(dss,gssm(m,1));
        end
        
        for m=1:M
            num(m,1)=elnsum(num(m,1),elnproduct(gssm(m,1),elnScale(win).window));
            denom=elnsum(denom, elnproduct(dss, elnScale(win).window));
        end
     end
        
     for m=1:M
         C(st,m)=elnproduct(num(m,1),-denom);
     end
 end
 
 rwsm=NaN(K,1);
for st=1:K
    for m=1:M
    rwsm(st,1)=elnsum(rwsm(st,1),C(st,m));
    end
   
end

for st=1:K
    for m=1:M
        C(st,m)=eexp(elnproduct(C(st,m),-rwsm(st,1)));
    end
end
  
 %re-estimating duration probabilities
for j=1:K
    num=NaN; denom=NaN;   
    for win=1:num_windows
    T=size(Data(win).window,1); 
    nss=NaN; dss=NaN; 
    for t0=1:T
        for t1=t0:T
            product=1;
            for tt=t0:t1
                product=product*eexp(elnB(win).window(tt,j));
            end
            product=eln(product);
            if t0==1
                xss= elnproduct(elnproduct(elnproduct(elnproduct(eln(sum(Pi.*P(:,j)')),product), elndur_p(win).window(j,t1-t0+1)),elnbeta(win).window(t1,j)),-elnScale(win).window)   ;
            else
                sumval=NaN;
                for st=1:K
                    sumval=elnsum(sumval, elnproduct(elnalpha(win).window(t0-1,st),eln(P(st,j))) );
                end
                xss= elnproduct(elnproduct(elnproduct(elnproduct(sumval,product), elndur_p(win).window(j,t1-t0+1)),elnbeta(win).window(t1,j)), -elnScale(win).window)  ;
            end
      nss=elnsum(nss, elnproduct(xss,eln(t1-t0+1))); dss=elnsum(dss, xss);      
    end
   % num=num+xss*(t1-t0+1);den=den+xss;
    end
    num=elnsum(num, elnproduct(nss,elnScale(win).window)); denom=elnsum(denom, elnproduct(dss,elnScale(win).window));
    end
    lambda(j)=eexp(elnproduct(num,-denom));
end


  
 

oldlik=lik;
  lik=NaN;
  for win=1:num_windows
      lik=elnsum(lik,elnScale(win).window);
  end
  
  if isnan(lik)
    P=Pcycle(cycle-1).cycle;LL=LLcycle(cycle-1).cycle; Mu=Mucycle(cycle-1).cycle;
    Cov=Covcycle(cycle-1).cycle;Pi=Picycle(cycle-1).cycle;C=Ccycle(cycle-1).cycle;
    lambda=lambdacycle(cycle-1).cycle;
    break;
  end
    
  
  %lik=sum(Scale);
  LL=[LL lik];
  fprintf('cycle %i log likelihood = %f ',cycle,lik);  
  
  if (cycle<=2)
    likbase=lik;
  elseif (lik<oldlik) 
    fprintf('log lik is dec');
  elseif ((lik-likbase)<(1 + tol)*(oldlik-likbase)|~isfinite(lik)) 
    fprintf('\n');
    break;
  end
  fprintf('\n');
 
  Pcycle(cycle).cycle=P;LLcycle(cycle).cycle=LL; Mucycle(cycle).cycle=Mu; 
  Covcycle(cycle).cycle=Cov; Picycle(cycle).cycle=Pi; Ccycle(cycle).cycle=C;
  lambdacycle(cycle).cycle=lambda;
end

end