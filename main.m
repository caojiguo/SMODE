% Simulations in Semiparametric mixed-effects ODE models with time-varying coefficients.
% Note: this program is built on Matlab (at least 2014 version) due to
% using some truncated distributions.

% load the fda package developed by James O. Ramsay
addpath('fdaM')

odeoptions = odeset('RelTol',1e-7,'AbsTol',1e-7);

%% Generate simulated gene expression data.

ncohort = 50;
tEnd    = 8;
ngrid   = 201;
nobs    = 201;
N = ncohort*nobs;

sigma_epsilon = 0.3;
tobs  = linspace(0,tEnd,nobs);
tfine   = linspace(0,tEnd,201);

input1.alpha = 1.2;
input1.beta = 3.5;
input1.delta = 1.0;

[tfine,X] = ode45(@Sim1Ode, tfine, 1,odeoptions,input1);

figure(1)
plot(tfine,X);


%% set up basis functions for eta(t).

% B-splines for eta(t).

knots_eta = linspace(0,tEnd,18);
norder   = 4;
nbasis_eta   = length(knots_eta) + norder - 2;
basis_eta   = create_bspline_basis([0,tEnd],nbasis_eta,norder,knots_eta);
 
% basis functions evaluated at tobs
Phimat_eta   = eval_basis(tobs, basis_eta);
K            = size(Phimat_eta,2);

% first derivative of basis functions evaluated at tobs
D1Phimat_eta = eval_basis(tobs, basis_eta, 1);
% second derivative of basis functions evaluated at tobs
D2Phimat_eta = eval_basis(tobs, basis_eta, 2);

AA = zeros((K-3),(K-1));
for i=1:(K-3)
    AA(i,i) = 1;
    AA(i,i+1) = -2;
    AA(i,i+2) = 1;
end

D_eta = AA'*AA;


%% Calculate penalty term using Simpson's rule
h = 0.01;
quadts = (0:h:tEnd)';
nquad = length(quadts);
quadwts = ones(nquad,1);
quadwts(2:2:(nquad-1)) = 4;
quadwts(3:2:(nquad-2)) = 2;
quadwts = quadwts.*(h/3);

D0Qbasismat_eta = eval_basis(quadts, basis_eta);
D1Qbasismat_eta = eval_basis(quadts, basis_eta, 1);
D2Qbasismat_eta = eval_basis(quadts, basis_eta, 2);
R1mat_eta = D1Qbasismat_eta'*(D1Qbasismat_eta.*(quadwts*ones(1,nbasis_eta)));
R2mat_eta = D2Qbasismat_eta'*(D2Qbasismat_eta.*(quadwts*ones(1,nbasis_eta)));

%% Generating the simulation data.

Nreps = 20; 

Sim1_n50_TrueParArray = cell(1,Nreps);
Sim1_n50_muSoluFArray = cell(1,Nreps);
Sim1_n50_mumatArray = cell(1,Nreps);
Sim1_n50_YobsmatArray = cell(1,Nreps);


Simmumat = zeros(nobs,ncohort); 
truepar = [1.2;3.5;1.0]; 
X0 = 1.0; 

truetheta = log(truepar); 
SimCorr_True = [1.0000   0.45  -0.25;
                0.45    1.0000    -0.45;
                -0.25   -0.45   1.0000];
SimSigma_True = 0.04*SimCorr_True;
%AD = chol(SimSigma_True,'lower'); 

thetaTru = log([1.2;3.5;1.0]);
    
for tt=1:Nreps
   display(num2str(tt));
 
    truealpha =  truepar(1);
    truebeta  =  truepar(2);
    truedelta =  truepar(3);
    
    input2.alpha = truealpha;
    input2.beta = truebeta;
    input2.delta = truedelta;
    
    [tobs,fit_pop] = ode45(@Sim1Ode, tobs, X0,odeoptions,input2);
    
%      figure(1)
%      plot(tobs,fit_pop);
        
    thetaImat=zeros(3,ncohort); 
 
    X0I = X0+0.1*trnd(3,ncohort,1);
    
    SimthetaImat = zeros(3,ncohort);
    U = gamrnd(3/2,2/3,ncohort,1);

    for i=1:ncohort
        SimthetaImat(:,i) = thetaTru+1/sqrt(U(i))*mvnrnd(zeros(3,1),SimSigma_True)';
    end
    
    for i=1:ncohort
   %     display(num2str(i))
        
        mytheta  = exp(SimthetaImat(:,i));
        input2.alpha = mytheta(1);
        input2.beta = mytheta(2);
        input2.delta = mytheta(3); 
        [tbs,myfit] = ode45(@Sim1Ode, tobs, X0I(i),odeoptions,input2);
        Simmumat(:,i)  = myfit;
    end
        
%     figure(11)
%     plot(tobs,Simmumat,'o-');
%     
%     figure(3)
%     hist(reshape(log(Simmumat),nobs*ncohort,1));
%     
%     figure(4)
%     qqplot(reshape(log(Simmumat),nobs*ncohort,1)) 
%     
%     figure(12)
%     for i=1:20
%         subplot(4,5,i);
%         plot(tobs,Simmumat(:,i));
%     end
    
    thetaITruM = SimthetaImat;
    
    SimthetaImat_True = zeros(4,ncohort);
    SimthetaImat_True(1,:) = X0I;
    SimthetaImat_True(2:4,:) = SimthetaImat;
    Sim1_n50_TrueParArray{tt} = SimthetaImat_True;
    Sim1_n50_mumatArray{tt}=Simmumat;
 
    Yobsmat = Simmumat+0.3*trnd(3,nobs,ncohort);
    
%     figure(5)
%     hist(reshape(Yobsmat,nobs*ncohort,1));
%     figure(6)
%     qqplot(reshape(Yobsmat,nobs*ncohort,1))
%     
%     figure(7)
%     for i=1:20
%         subplot(4,5,i)
%         plot(tobs,Simmumat(:,i),'k-',tobs,exp(Yobsmat(:,i)),'o')
%     end 

    Sim1_n50_YobsmatArray{tt} = Yobsmat; 
end 

figure(5)
hist(reshape(Yobsmat,nobs*ncohort,1));
figure(6)
qqplot(reshape(Yobsmat,nobs*ncohort,1)) 
figure(34) 
hist(reshape(Yobsmat,ncohort*nobs,1)); 

% save Sim1_n50_YobsmatArray Sim1_n50_YobsmatArray;
% save Sim1_n50_mumatArray Sim1_n50_mumatArray;
% save Sim1_n50_TrueParArray Sim1_n50_TrueParArray;

  
%% Starting values of mu_g(t). 

tfine = linspace(0,tEnd,201);
trueeta = exp(-(tobs-2).^2/1.44)+exp(-(tobs-5).^2)+0.5*exp(-(tobs-6).^2/0.64)-0.062; 
trueeta_fine = exp(-(tfine-2).^2/1.44)+exp(-(tfine-5).^2)+0.5*exp(-(tfine-6).^2/0.64)-0.062; 

basis_eta   = create_bspline_basis([0,tEnd],nbasis_eta,norder,knots_eta);
Phimat_eta_fine = eval_basis(tfine, basis_eta);

zetaInt    = (Phimat_eta'*Phimat_eta+0.001*R2mat_eta)\(Phimat_eta'*trueeta); 
zetaInt(1) = 0;
oldeta = D0Qbasismat_eta*zetaInt; 

figure(4)
plot(tfine,trueeta_fine,'k-',quadts,oldeta,'r--') 


%%
options1 = optimset( ...
    'Algorithm', 'levenberg-marquardt', ...
    'LargeScale', 'off', ...
    'Display','iter', ...
    'DerivativeCheck','off', ...
    'MaxIter', 1000, ...
    'Jacobian','off',...
    'MaxFunEvals', 100000,...
    'TolFun', 1e-10, ...
    'TolX',   1e-10);
options2 = optimset( ...
    'Algorithm', 'levenberg-marquardt', ...
    'LargeScale', 'off', ...
    'Display','off', ...
    'DerivativeCheck','off', ...
    'MaxIter', 1000, ...
    'Jacobian','off',...
    'MaxFunEvals', 100000,...
    'TolFun', 1e-10, ...
    'TolX',   1e-10); 

%% Metroplis-Hastings within Gibbs Samplepler, fixed lambda=10^4.
% Setting up. 
%%

Sim1_n50_thetaSam_TTArray=cell(Nreps,1);
Sim1_n50_zetaSam_TTArray=cell(Nreps,1);

Sim1_n50_thetaSam_NNArray=cell(Nreps,1);
Sim1_n50_zetaSam_NNArray=cell(Nreps,1);


for tt=1:Nreps
    display(num2str(tt));

    Yobsmat=Sim1_n50_YobsmatArray{tt};
    thetaITru=Sim1_n50_TrueParArray{tt};
    Simmumat = Sim1_n50_mumatArray{tt};
    SimthetaImat = thetaITru(2:4,:);
    
    BurnIn = 1e4;
    Final_size = 1e4; 

    p  = 4; 
    df  = p+1;

    Psi01 = 100*eye(p);
    S01   = inv(Psi01); 
    eta0 = [0;log(1.2);log(3.5);log(1.0)];
    Omega = 0.1*eye(p); 

    InvOmega = inv(Omega); 

    thetaInt = [0;log(1.2);log(3.5);log(1.0)]+normrnd(0,1,p,1); 
    
    thetaIInt = thetaITru+normrnd(0,1,p,ncohort);
    thetaIInt(1,:) = thetaITru(1,:);
    
    InvSigma1Int = eye(p); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% T/T model.                                                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tauESam_TT  = zeros(Final_size,1); 
    zetaSam_TT       = zeros(K,Final_size);
    thetaSam_TT      = zeros(p,Final_size);
    thetaISam_TT     = zeros(p,ncohort,Final_size);
    InvSigmaSam_TT  = zeros(p,p,Final_size);
    muSam_TT = zeros(nobs,ncohort,Final_size);
    lambdaSam_TT = zeros(Final_size,1); 
    USam_TT = zeros(ncohort,Final_size);
    WSam_TT = zeros(ncohort,Final_size);
    kappaSam_TT = zeros(Final_size,1);
    nuSam_TT = zeros(Final_size,1); 
    
    oldmumat  = Simmumat+normrnd(0,0.5,nobs,ncohort);

    oldtheta  = thetaInt;
    oldzeta   = zetaInt+normrnd(0,0.001,length(zetaInt),1); 
    oldzeta(1) = zetaInt(1);
    oldtauE  = 1;
    oldlambda_eta = 0.01; 

    oldInvSigma  = InvSigma1Int; 
    
    oldthetaImat = thetaIInt;
 
    oldkappa = 3;
    oldnu    = 3; 
    oldUmat = ones(ncohort,1);
    oldWmat = ones(ncohort,1);
   
    
    %  proposal; 
    sig_thetaI = cell(1,ncohort);
    for i=1:ncohort
         sig_thetaI{i}=  0.2*diag([0.20;0.20; 0.20; 0.20]);
    end
    sig_zeta  = 0.02*abs(zetaInt);
    sig_kappa = 0.2;
    sig_nu = 0.2; 
    
    countthetaI = zeros(ncohort,1);
    countzeta   = 0;
    countkappa = 0;
    countnu = 0; 

    input1.ncohort = ncohort;
    input1.oldthetaImat = oldthetaImat;
    input1.oldtheta     = oldtheta;
    input1.oldmumat    = oldmumat;
    input1.oldzeta      = oldzeta;
    input1.oldtauE      = oldtauE; 
    input1.oldUmat      = oldUmat;
    input1.oldWmat      = oldWmat;
    input1.oldkappa     = oldkappa;
    input1.oldnu        = oldnu; 

    input1.sig_thetaI = sig_thetaI;
    input1.sig_zeta = sig_zeta;
    input1.sig_kappa = sig_kappa;
    input1.sig_nu = sig_nu; 

    input1.p    = p;
    input1.K    = K;
    input1.tobs = tobs';
    input1.Yobsmat = Yobsmat;
    input1.quadwts = quadwts;
    input1.quadts  = quadts; 
    input1.eta0    = eta0;
    input1.S01     = S01;
    input1.df      = df;
    input1.shape_tauE       = 1;
    input1.scale_tauE       = 10; 
    input1.shape_lambda = 1;
    input1.scale_lambda = 100;
    input1.D0Qbasismat_eta = D0Qbasismat_eta;
    input1.oldInvSigma = oldInvSigma;
    input1.InvOmega    = InvOmega;
    input1.basis_eta   = basis_eta;
    input1.D_eta        = D_eta; 
    input1.oldlambda_eta = oldlambda_eta; 
    input1.countthetaI = countthetaI;
    input1.countzeta   = countzeta;
    input1.countkappa  = countkappa;
    input1.countnu     = countnu; 
    
    input1.odeoptions = odeoptions;

    out1 = Sim1MCMC_funTT(input1);

    oldmu_k = zeros(p,ncohort);
    oldSigma_k = sig_thetaI;
    newSigma_k = cell(1,ncohort);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Burn-in 
    for iter=1:BurnIn
 %       display(num2str(iter))
 
        input1 = out1;
        out1 = Sim1MCMC_funTT(input1);

        oldthetaImat = out1.oldthetaImat;

% Adaptive proposal.

       newmu_k = oldmu_k+1/iter*(oldthetaImat-oldmu_k);
       for i = 1:ncohort
           newSigma_k{i}=oldSigma_k{i}+1/iter*((oldthetaImat(:,i)-oldmu_k(:,i))*(oldthetaImat(:,i)-oldmu_k(:,i))'-oldSigma_k{i})+0.000001*eye(p); 
        end   
        oldmu_k = newmu_k;
        oldSigma_k=newSigma_k;

        for i=1:ncohort
              out1.sig_thetaI{i}=0.05*oldSigma_k{i};	
        end

    end
    
    countthetaI  = out1.countthetaI;
    countzeta    = out1.countzeta;
    countkappa   = out1.countkappa;
    countnu      = out1.countnu;

    countthetaI/BurnIn;
    countzeta/BurnIn;
    countkappa/BurnIn;
    countnu/BurnIn; 
    
    %% Final Samples. 
    
    countthetaI = zeros(ncohort,1);
    countzeta   = 0;
    countkappa = 0;
    countnu = 0; 
    out1.countthetaI = countthetaI;
    out1.countzeta = countzeta;
    out1.countkappa = countkappa;
    out1.countnu = countnu;        
 
    for iter=1:Final_size
  %      display(num2str(iter)) 
         input1 = out1;
        out1 = Sim1MCMC_funTT(input1);

        oldthetaImat = out1.oldthetaImat;

% Adaptive proposal.

       newmu_k = oldmu_k+1/(BurnIn+iter)*(oldthetaImat-oldmu_k);
       for i = 1:ncohort
           newSigma_k{i}=oldSigma_k{i}+1/(BurnIn+iter)*((oldthetaImat(:,i)-oldmu_k(:,i))*(oldthetaImat(:,i)-oldmu_k(:,i))'-oldSigma_k{i})+0.000001*eye(p); 
        end   
        oldmu_k = newmu_k;
        oldSigma_k=newSigma_k;

        for i=1:ncohort
              out1.sig_thetaI{i}=0.05*oldSigma_k{i};	
        end
        
        oldthetaImat = out1.oldthetaImat;
        oldtheta     = out1.oldtheta;
        oldmumat     = out1.oldmumat;
        oldInvSigma  = out1.oldInvSigma;
        oldtauE      = out1.oldtauE;
        oldlambda_eta= out1.oldlambda_eta; 
        oldUmat      = out1.oldUmat;
        oldWmat      = out1.oldWmat;
        oldkappa     = out1.oldkappa;
        oldnu        = out1.oldnu; 
        oldzeta      = out1.oldzeta;
        
        countthetaI  = out1.countthetaI;
        countzeta    = out1.countzeta;
        countkappa   = out1.countkappa;
        countnu      = out1.countnu; 
        
        thetaISam_TT(:,:,iter) = oldthetaImat;
        zetaSam_TT(:,iter)     = oldzeta; 
        thetaSam_TT(:,iter)    = oldtheta;
        InvSigmaSam_TT(:,:,iter) = oldInvSigma;
        tauESam_TT(iter)       = oldtauE;
        lambdaSam_TT(iter)     = oldlambda_eta; 
        USam_TT(:,iter)  = oldUmat;
        WSam_TT(:,iter)  = oldWmat;
        kappaSam_TT(iter) = oldkappa;
        nuSam_TT(iter) = oldnu; 

        muSam_TT(:,:,iter)    = oldmumat; 
    end 
    countthetaI/Final_size
    countzeta/Final_size
    countkappa/Final_size
    countnu/Final_size 

    display(num2str('End of Student-t/Student-t model')); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Normal/Normal model.                                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tauESam_NN  = zeros(Final_size,1);  
    zetaSam_NN       = zeros(K,Final_size);
    thetaSam_NN      = zeros(p,Final_size);
    thetaISam_NN     = zeros(p,ncohort,Final_size);
    InvSigmaSam_NN  = zeros(p,p,Final_size);
    muSam_NN = zeros(nobs,ncohort,Final_size);
    lambdaSam_NN = zeros(Final_size,1); 

    countthetaI = zeros(ncohort,1);
    countzeta   = 0;

    input2.ncohort = ncohort;
    input2.oldthetaImat = oldthetaImat;
    input2.oldtheta     = oldtheta;
    input2.oldmumat    = oldmumat;
    input2.oldzeta      = oldzeta;
    input2.oldtauE      = oldtauE; 

    input2.sig_thetaI = sig_thetaI;
    input2.sig_zeta = sig_zeta;

    input2.p    = p;
    input2.K    = K;
    input2.tobs = tobs';
    input2.Yobsmat = Yobsmat;
    input2.quadwts = quadwts;
    input2.quadts  = quadts; 
    input2.eta0    = eta0;
    input2.S01     = S01;
    input2.df      = df;
    input2.shape_tauE  = 1;
    input2.scale_tauE  = 10; 
    input2.shape_lambda = 1;
    input2.scale_lambda = 100;
    input2.D0Qbasismat_eta = D0Qbasismat_eta;
    input2.oldInvSigma = oldInvSigma;
    input2.InvOmega    = InvOmega;
    input2.basis_eta   = basis_eta;
    input2.D_eta        = D_eta; 
    input2.oldlambda_eta = oldlambda_eta; 
    input2.countthetaI = countthetaI;
    input2.countzeta   = countzeta;

    input2.odeoptions = odeoptions;

    out2 = Sim1MCMC_funNN(input2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Burn-in
    
    for iter=1:BurnIn
%        display(num2str(iter))

        input2 = out2;
        out2 = Sim1MCMC_funNN(input2);

        oldthetaImat = out2.oldthetaImat;

% Adaptive proposal.

       newmu_k = oldmu_k+1/iter*(oldthetaImat-oldmu_k);
       for i = 1:ncohort
           newSigma_k{i}=oldSigma_k{i}+1/iter*((oldthetaImat(:,i)-oldmu_k(:,i))*(oldthetaImat(:,i)-oldmu_k(:,i))'-oldSigma_k{i})+0.000001*eye(p); 
        end   
        oldmu_k = newmu_k;
        oldSigma_k=newSigma_k;

        for i=1:ncohort
              out2.sig_thetaI{i}=0.1*oldSigma_k{i};	
        end	

    end

    countthetaI  = out2.countthetaI;
    countzeta    = out2.countzeta;

    countthetaI/BurnIn;
    countzeta/BurnIn;

    %% Final Samples. 

    countthetaI = zeros(ncohort,1);
    countzeta   = 0;

    out2.countthetaI = countthetaI;
    out2.countzeta = countzeta;

    for iter=1:Final_size
  %      display(num2str(iter)) 
        input2 = out2; 
        out2 = Sim1MCMC_funNN(input2);

% Adaptive proposal.

       newmu_k = oldmu_k+1/(BurnIn+iter)*(oldthetaImat-oldmu_k);
       for i = 1:ncohort
           newSigma_k{i}=oldSigma_k{i}+1/(BurnIn+iter)*((oldthetaImat(:,i)-oldmu_k(:,i))*(oldthetaImat(:,i)-oldmu_k(:,i))'-oldSigma_k{i})+0.000001*eye(p); 
        end   
        oldmu_k = newmu_k;
        oldSigma_k=newSigma_k;

        for i=1:ncohort
              out2.sig_thetaI{i}=0.1*oldSigma_k{i};	
        end	  

        oldthetaImat = out2.oldthetaImat;
        oldtheta     = out2.oldtheta;
        oldmumat     = out2.oldmumat;
        oldInvSigma  = out2.oldInvSigma;
        oldtauE      = out2.oldtauE;
        oldlambda_eta= out2.oldlambda_eta; 
        oldzeta      = out2.oldzeta;
        countthetaI  = out2.countthetaI;
        countzeta    = out2.countzeta;

        thetaISam_NN(:,:,iter) = oldthetaImat;
        zetaSam_NN(:,iter)     = oldzeta; 
        thetaSam_NN(:,iter)    = oldtheta;
        InvSigmaSam_NN(:,:,iter) = oldInvSigma;
        tauESam_NN(iter)       = oldtauE;
        lambdaSam_NN(iter)     = oldlambda_eta; 

        muSam_NN(:,:,iter)    = oldmumat; 
    end 
    countthetaI/Final_size
    countzeta/Final_size

    display(num2str('End of Normal/Normal model')); 
  
%% 


    Sim1_n50_thetaSam_TTArray{tt}=thetaSam_TT;
    Sim1_n50_zetaSam_TTArray{tt}=zetaSam_TT;
    Sim1_n50_thetaSam_NNArray{tt}=thetaSam_NN;
    Sim1_n50_zetaSam_NNArray{tt}=zetaSam_NN;
    
    

    save Sim1_n50_thetaSam_NNArray Sim1_n50_thetaSam_NNArray
    save Sim1_n50_thetaSam_TTArray Sim1_n50_thetaSam_TTArray
    save Sim1_n50_zetaSam_NNArray Sim1_n50_zetaSam_NNArray
    save Sim1_n50_zetaSam_TTArray Sim1_n50_zetaSam_TTArray
    
%     save Sim1_n50_thetaISampleArray Sim1_n50_thetaISampleArray
%     save Sim1_n50_thetaISam_TTArray Sim1_n50_thetaISam_TTArray
%     save Sim1_n50_tauESampleArray Sim1_n50_tauESampleArray
%     save Sim1_n50_tauESam_TTArray Sim1_n50_tauESam_TTArray
%     save Sim1_n50_InvSigmaSampleArray Sim1_n50_InvSigmaSampleArray
%     save Sim1_n50_InvSigmaSam_TTArray Sim1_n50_InvSigmaSam_TTArray
    
end

