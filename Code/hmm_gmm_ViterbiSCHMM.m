function [Jstar,Dstar]=hmm_gmm_ViterbiSCHMM(Data,T,K,M,Mu,Cov,P,Pi,C,lambda)
[N,p]=size(Data);
k1=(2*pi)^(-p/2);

state_seq=[];
for n=1:N/T
    
    elndur_p=NaN(K,T);
      for st=1:K
        elndur_p(st,:)=eln(poisspdf(1:T,lambda(st)));
      end
      
iCov=[];k2=[];  elnB=NaN(T,K); 

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
        
%     B=eexp(elnB);
%     B=rdiv(B,rsum(B));

%% write the Viterbi algorithm for decoding state sequences from here


 elndelta_init=NaN(K,T);
      for st=1:K
        elndelta_init(st,:)=eln(poisspdf(1:T,Pi(st)));
      end

 elndelta.time.state.dur=[]; phi.time.state.dur=[];
 
      for t=1:T
          for j=1:K
              for d=1:T
                 
                  elndelta(t).time(j).state(d).dur=NaN;
                  phi(t).time(j).state(d).dur=zeros(3,1);
                  if t-d==0
                      valmatrix=zeros(K,T);
                      for i=1:K
                      for h=1:T
                     %valmatrix(i,h)= elndelta_init(i,h)+eln(P(i,j))+elndur_p(j,d)+sum(elnB(t-d+1:t,j));
                     valmatrix(i,h)= elndelta_init(i,h)+eln(P(j,i))+elndur_p(j,d)+sum(elnB(t-d+1:t,j));

                      end
                      end
                      elndelta(t).time(j).state(d).dur=max(valmatrix(:));
                      %[istar, hstar]=find(valmatrix==elndelta(t).time(j).state(d).dur);
                      [istar, hstar] = find(ismember(valmatrix, max(valmatrix(:))));
                      if isempty(istar)
                      istar=0;
                      end
                      if isempty(hstar)
                      hstar=0;
                      end
                      %phi(t).time(j).state(d).dur(1)=t-d;phi(t).time(j).state(d).dur(2)=istar; phi(t).time(j).state(d).dur(3)=hstar;
                      phi(t).time(j).state(d).dur=[t-d;istar;hstar];
                  end
                  
                  if t-d>0
                      valmatrix=zeros(K,T);
                      for i=1:K
                      for h=1:T
                     valmatrix(i,h)= elndelta(t-d).time(i).state(h).dur+eln(P(j,i))+elndur_p(j,d)+sum(elnB(t-d+1:t,j));
                      end
                      end
                      elndelta(t).time(j).state(d).dur=max(valmatrix(:));
                      [istar, hstar] = find(ismember(valmatrix, max(valmatrix(:))));
                      if isempty(istar)
                      istar=0;
                      end
                      if isempty(hstar)
                      hstar=0;
                      end
                      %phi(t).time(j).state(d).dur(1)=t-d;phi(t).time(j).state(d).dur(2)=istar; phi(t).time(j).state(d).dur(3)=hstar;
                      phi(t).time(j).state(d).dur=[t-d;istar;hstar];
                  end
                      
              end
          end
      end

%Trace back part of code
nn=1; tt(nn)=T;
valmatrix=NaN(K,T);      
for j=1:K
    for d=1:T
        valmatrix(j,d)=elndelta(tt(nn)).time(j).state(d).dur;
    end
end
[jj,dd]=find(valmatrix==max(valmatrix(:)));

if isempty(jj)& isempty(dd)
[~,jj]=max(Pi);
dd=1;
end


Jstar(nn)=jj; Dstar(nn)=dd;
nn=nn+1;

while (nn>0)
    
    temp=phi(tt(nn-1)).time(Jstar(nn-1)).state(Dstar(nn-1)).dur;
    tt(nn)=temp(1);Jstar(nn)=temp(2);Dstar(nn)=temp(3);
     if Jstar(nn)==0 & Dstar(nn)==0
         Jstar(nn)=Jstar(nn-1); 
         Dstar(nn)=Dstar(nn-1);
     end
     
   if tt(nn)-Dstar(nn)+1<=1
        break;
    else
        nn=nn+1;
    end
end





% delta=zeros(T,K);phi.dur=zeros(T,K);phi.state=zeros(T,K);
% p_star=0; q_star=zeros(T,1);
% 
% %initialize
% for st=1:K
% delta(1,st)=Pi(st)*B(1,st)*dur_p(st,1);
% phi.dur(1,st)=0; phi.state(1,st)=0;
% end
% 
% %recursion and forward propagation step
% for i=1:K
%     
%     for t=2:T
%         sc=zeros(1,t-1);
%         for tt=1:t-1
%            [val,arg]=max(delta(tt,:).*P(:,i)');
%            phi.state(t,i)=arg;
%            product=1;
%            for pt=tt+1:t
%                product=product*B(pt,i);
%            end
%            sc(tt)=val*dur_p(i,t-tt)*product;
%         end
%         [delta(t,i),phi.dur(t,i)]=max(sc);
%         
%     end
% end
%         
%             
% %termination
% [p_star,q_star(T)]=max(delta(T,:));
% tau_star=T-phi.dur(T,q_star(T));
% for tau=1:tau_star
%     q_star(T-tau)=q_star(T);
% end
% %backtracking recursively
% for t=T-tau_star:-1:1
%     q_star(t)=phi.state(t+tau_star,q_star(t+tau_star));
%     tau_star=t-phi.dur(t,q_star(t));
%     for tt=1:tau_star-1
%         q_star(t-tt)=q_star(t);
%     end
%     t=t-tau_star;
% end
% state_seq=[state_seq;q_star];


% viterbi end 
end







end