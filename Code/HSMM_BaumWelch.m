function [P,LL,Mu,Cov,Pi,C,lambda]=HSMM_BaumWelch(Data,K,M,Mu,Cov,P,Pi,C,lambda,cyc,tol,p)


elngamma=[];gammasum=[];
Xi=[];elnScale=[];elnxi_prior=[];
tempr=[];elnB=[];elndur_p=[];
%cycle loop starts here
lik=0; LL=[];


for cycle=1:cyc
        
%%%% FORWARD-BACKWARD 
    
    T=size(Data,1);
    elndur_p=NaN(K,T);
    for st=1:K
        elndur_p(st,:)=eln(poisspdf(1:T,lambda(st)));
    end
  
    elnScale=NaN; 
    elnXsi=NaN(T-1,K*K);
    elnalpha=NaN(T,K); elnbeta=NaN(T,K); elnalpha_star=NaN(T,K); elnbeta_star=NaN(T,K);
    elnxi_prior=NaN(K,1); elnGamma=NaN(T,K);
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

    %Backward run
    elnbeta(T,:)=eln(ones(1,K));
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

    %evaluating neta
    normalizer=NaN;
    for st=1:K
        elnGamma(1,st)=elnproduct(eln(Pi(st)),elnbeta_star(1,st));
        normalizer=elnsum( normalizer, elnGamma(1,st) );
    end

    for st=1:K
        elnGamma(1,st)=elnproduct(elnGamma(1,st),-normalizer);
    end


    for t=1:T-1
        tmp=zeros(K,K);
        for i=1:K
            for j=1:K
                tmp(i,j)=elnproduct( elnproduct(elnalpha(t,i), eln(P(i,j))), elnbeta_star(t+1,j)) ;
            end     
        end
        elnXsi(t,:)=reshape(tmp,1,K*K);
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
        for i=1:K
            for t=1:T
                elngamma(m).mixture(t,i)=elnproduct(elnGamma(t,i), elnproduct(eln(tempr(m).mixture(t,i)), -elnB(t,i)) ); 
            end
        end
    end
    
    
    
    %log the value of xi for prior in M step
    tmp=reshape(elnXsi(1,:),K,K);
    for st=1:K
        rsumval=NaN;
        for j=1:K
            rsumval=elnsum(rsumval,tmp(st,j));
        end
        elnxi_prior(st)=rsumval;
    end
    
    %evaluate the likelihood of the sequence
    for t=1:T
        for i=1:K
            for j=1:K
                for d=1:T
                    if t-d>0
                        elnScale=elnsum(elnScale, elnproduct(elnproduct(elnproduct(elnproduct(elnalpha(t-d,i),eln(P(i,j))), elndur_p(j,d)),elnU(t,j,d)), elnbeta(t,j))  );
                    end
                end
            end
        end
    end
    
    
    
 

  
    %%%% M STEP 
  
    % outputs
    for m=1:M
        num=NaN(K,p);denom=NaN(K,1);
        %T=size(Data,1);
        dataw=eln(Data);  
        nss=NaN(K,p); dss=NaN(K,1); 
        for st=1:K
            for pp=1:p
                for tt=1:T
                    nss(st,pp)=elnsum(nss(st,pp),elnproduct(elngamma(m).mixture(tt,st),dataw(tt,pp)) );
                end
            end
          
            for tt=1:T
                dss(st,1)=elnsum(dss(st,1),elngamma(m).mixture(tt,st) );
            end
        end
      
      
        for st=1:K
            for pp=1:p
                num(st,pp)=elnsum(num(st,pp), elnproduct(nss(st,pp),elnScale));
            end
            denom(st,1)=elnsum(denom(st,1), elnproduct(dss(st,1),elnScale) );
        end
  
  
        for st=1:K
            for pp=1:p
                Mu(m).mixture(st,pp)=eexp(elnproduct(num(st,pp),-denom(st)));%rdiv(num,denom');
            end
        end
    end
  

    % transition matrix 
    num=NaN(1,K*K); denom=NaN(K,1);
    %T=size(Data,1);
    tsum=NaN(1,K*K);
    for dm=1:K*K
        for t=1:T-1
            tsum(1,dm)=elnsum(tsum(1,dm),elnXsi(t,dm));
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
        num(1,dm)=elnsum(num(1,dm),elnproduct(tsum(1,dm),elnScale) );
    end
   
    for st=1:K
        denom(st,1)=elnsum(denom(st,1), elnproduct(rsum(st),elnScale) );
    end
    
  
  
    sxi=NaN(K,K);
    num=reshape(num, K,K);
    for i=1:K
        for j=1:K
            sxi(i,j)=elnproduct(num(i,j),-denom(i));
        end
    end
 
  

    for i=1:K
        for r=1:K
            if r==i
                sxi(i,r)=NaN;
            end
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
  
    % priors
    Pi=NaN(1,K);
 
    for st=1:K
        Pi(1,st)=elnsum(Pi(1,st),elnproduct(elnxi_prior(st),elnScale));
    end

    normalizer=NaN;
    for st=1:K
        normalizer=elnsum(normalizer,Pi(st));
    end

    for st=1:K
        Pi(st)=elnproduct(Pi(st),-normalizer);
    end

    Pi=eexp(Pi);

  
    % covariance
   
    for m=1:M
        num=NaN(K,1); denom=NaN(K,1);
        nss=NaN(K,1); dss=NaN(K,1);
        for st=1:K
            for t=1:T
                nss(st,1)=elnsum(nss(st,1), elnproduct( elngamma(m).mixture(t,st),eln((Data(t,:)-Mu(m).mixture(st,:))*(Data(t,:)-Mu(m).mixture(st,:))') )   );
                dss(st,1)=elnsum(dss(st,1),elngamma(m).mixture(t,st) );
            end
      
        end
      
        for st=1:K
            num(st,1)=elnsum(num(st,1),elnproduct(nss(st,1),elnScale));
            denom(st,1)=elnsum(denom(st,1),elnproduct(dss(st,1),elnScale));
        end
      
  
        for st=1:K
            Cov(m).mixture(st)=eexp(elnproduct(num(st),-denom(st)));
        end
    end
  
    %mixture weights of gaussian probabilities
           
    for st=1:K
        num=NaN(M,1); denom=NaN;
        %T=size(Data,1); 
        gssm=NaN(M,1);
        for m=1:M
            for t=1:T
                gssm(m,1)=elnsum(gssm(m,1),elngamma(m).mixture(t,st));
            end
        end
        
        dss=NaN;
        for m=1:M
            dss=elnsum(dss,gssm(m,1));
        end
        
        for m=1:M
            num(m,1)=elnsum(num(m,1),elnproduct(gssm(m,1),elnScale));
            denom=elnsum(denom, elnproduct(dss, elnScale));
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
        %T=size(Data,1); 
        nss=NaN; dss=NaN; 
        for t0=1:T
            for t1=t0:T
                product=1;
                for tt=t0:t1
                    product=product*eexp(elnB(tt,j));
                end
                product=eln(product);
                if t0==1
                    xss= elnproduct(elnproduct(elnproduct(elnproduct(eln(sum(Pi.*P(:,j)')),product), elndur_p(j,t1-t0+1)),elnbeta(t1,j)),-elnScale);
                else
                    sumval=NaN;
                    for st=1:K
                        sumval=elnsum(sumval, elnproduct(elnalpha(t0-1,st),eln(P(st,j))) );
                    end
                    xss= elnproduct(elnproduct(elnproduct(elnproduct(sumval,product), elndur_p(j,t1-t0+1)),elnbeta(t1,j)), -elnScale)  ;
                end
                nss=elnsum(nss, elnproduct(xss,eln(t1-t0+1))); dss=elnsum(dss, xss);      
            end
        end
        num=elnsum(num, elnproduct(nss,elnScale)); denom=elnsum(denom, elnproduct(dss,elnScale));
    
        lambda(j)=eexp(elnproduct(num,-denom));
    end


  
 

    oldlik=lik;
    lik=NaN;
    lik=elnsum(lik,elnScale);

  
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
 
  
end

end
