%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code loads the data and train a Hidden Semi-Markov model on the data
% This is followed by using the trained model for mode idnetification by using the Viterbi algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%First let us read the data

filename = 'Datasets/Decreasing_temperature_data.csv';
DATA = csvread(filename);

temperature= DATA(:,1); % the first column includes the temperatures
residence_time= DATA(:,2); % the second column includes the residence time across which spectra are sampled
Data_spectra= DATA(:,3:end);

%Define parameters that control the HSMM model structure

K=4; % K is the number of states
M=2; % M is the number of mixture components used to model emission probabilitites using a Gaussian Mixture Model
cyc=10; % cyc is the number of Expectation Maximization iterations
tol=1e-3; % convergence tolerance of the EM algorithm

%%Parameters learned via the EM algorithm are initialized

Pi=[]; % Pi is the initial state distribution
Pi=rand(1,K);
Pi=Pi/sum(Pi);

P=[];  % P is the transition probability 
P=rand(K);
for i=1:K
  for r=1:K
    if r==i
      P(i,r)=0;
    end

  end
end
P=rdiv(P,rsum(P));

Mu=[]; Cov=[]; C=[]; % Mu is the mean, Cov is the covariance, C is the mixing weights of Gaussian components used in modeling the emission probability distribution
%k means initialization for the emission probability represented by GMM
d=size(Data_spectra,2);
idx=kmeans(Data_spectra,M);
for i=1:M
    data=Data_spectra(find(idx ==i),:);
    Mu(i).mixture=randn(K,d)*sqrtm(diag(diag(cov(data))))+ones(K,1)*mean(data);
    Cov(i).mixture=rand(K,1);
end

C=rand(K,M);
for i=1:K
   C(i,:)=C(i,:)./sum(C(i,:));
end


lambda=[]; % parameter characterizing the state duration distribution
segment_length=round(size(Data_spectra,1)/K); % segment the observation sequence into equal parts by the hidden states
lambda=sort(randsample(segment_length:5:2*segment_length,K)); % randomly initialize the average state duration parameter by sampling



%%Book-keeping of inferences using the estimated HSMM parameters from EM

LL=[]; % LL is the log-likelihood recorded across all EM iterations
Gamma=[]; % posterior state probabailities along the observation sequence
Viterbi_seq=[]; % optimal sequence of states from the Viterbi algorithm


[P,LL,Mu,Cov,Pi,C,lambda]=HSMM_BaumWelch(Data_spectra,K,M,Mu,Cov,P,Pi,C,lambda,cyc,tol,d);

%%Use the HSMM parameters in Viterbi state decoding

disp('Viterbi state decoding')

T=size(Data_spectra,1);
[Jstar, Dstar]=hmm_gmm_ViterbiHSMM(Data_spectra,T,K,M,Mu,Cov,P,Pi,C,lambda);
for ii=1:length(Jstar)
  Viterbi_seq=[repmat([Jstar(ii)],Dstar(ii),1)];
end

%%Use the HSMM parameters for posterior decoding

disp('Posterior state decoding')

Gamma=posterior_decodingHSMM(Data_spectra,T,K,M,Mu,Cov,P,Pi,C,lambda);



