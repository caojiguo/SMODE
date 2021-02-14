function out1 = Sim1MCMC_funNN(input1)
% Note: this program is built on Matlab (at least 2014 version) due to
% using some truncated distributions.

        ncohort      = input1.ncohort;
        oldthetaImat = input1.oldthetaImat;
        oldtheta     = input1.oldtheta;
        oldmumat     = input1.oldmumat;
        oldzeta      = input1.oldzeta;
        oldtauE      = input1.oldtauE;
 
        oldInvSigma = input1.oldInvSigma;
        oldlambda_eta = input1.oldlambda_eta;  

        sig_thetaI    = input1.sig_thetaI;
        sig_zeta    = input1.sig_zeta;

        p       = input1.p;
        K       = input1.K;
        tobs    = input1.tobs;
        nobs    = length(tobs);

        Yobsmat = input1.Yobsmat;
        
        basis_eta = input1.basis_eta;

        eta0    = input1.eta0;
        S01     = input1.S01;
        df      = input1.df;
        shape_tauE = input1.shape_tauE;
        scale_tauE = input1.scale_tauE;   
        shape_lambda = input1.shape_lambda;
        scale_lambda = input1.scale_lambda;          
        

        InvOmega       = input1.InvOmega;
        D_eta           = input1.D_eta;
        
       
        countthetaI = input1.countthetaI;
        countzeta   = input1.countzeta;
       
        odeoptions  = input1.odeoptions;
 
        %% Update theta_i
        
    for i=1:ncohort
       
        oldthetaI   = oldthetaImat(:,i);   
        %newthetaI   = oldthetaI+normrnd(0,sig_thetaI,p,1);
        newthetaI   = mvnrnd(oldthetaI, sig_thetaI{i})';


        Yobs = Yobsmat(:,i);

        old.y0      = oldthetaI(1);
        myfitstr1.alpha = exp(oldthetaI(2));
        myfitstr1.beta = exp(oldthetaI(3));
        myfitstr1.delta = exp(oldthetaI(4));
         
        myfitstr1.basis_eta = basis_eta;
        myfitstr1.zeta = oldzeta;       
 
        [tbs,oldfit] = ode45(@Sim1WorkingOde, tobs, old.y0,odeoptions,myfitstr1);
        oldmumat(:,i) = oldfit;

%       oldfit = oldmumat(:,i); 
        
        olddiff  = (Yobs-oldfit);  
        oldSSE   = -0.5*oldtauE*sum(olddiff.^2);
        oldprior = -0.5*(oldthetaI-oldtheta(1:p))'*oldInvSigma*(oldthetaI-oldtheta(1:p));
        
        den1     = oldSSE + oldprior;

        new.y0      = newthetaI(1);
        myfitstr1.alpha = exp(newthetaI(2));
        myfitstr1.beta = exp(newthetaI(3));
        myfitstr1.delta = exp(newthetaI(4));
        
        [tbs,newfit] = ode45(@Sim1WorkingOde, tobs, new.y0,odeoptions,myfitstr1);
        
        newdiff  = (Yobs-newfit); 
        newSSE   = -0.5*oldtauE*sum(newdiff.^2);
        newprior = -0.5*(newthetaI-oldtheta(1:p))'*oldInvSigma*(newthetaI-oldtheta(1:p));
        
        num1     = newSSE + newprior;
       
        if log(rand(1))<=(num1-den1)
            oldthetaImat(:,i) = newthetaI;
            countthetaI(i)    = countthetaI(i)+1;
            oldmumat(:,i)   = newfit;
        else
            oldthetaImat(:,i) = oldthetaI;
            oldmumat(:,i)   = oldfit;
        end
    end

    %% Update zeta.
    
    oldSSE = 0; newSSE = 0;
    
    newzeta = oldzeta+normrnd(0,sig_zeta,K,1);
    newzeta(1) = 0;

    newfitmat = oldmumat;

    for i=1:ncohort

        oldthetaI = oldthetaImat(:,i);
        Yobs = Yobsmat(:,i);
        
        old.y0      = oldthetaI(1);
        myfitstr1.alpha = exp(oldthetaI(2));
        myfitstr1.beta = exp(oldthetaI(3));
        myfitstr1.delta = exp(oldthetaI(4));
         
        myfitstr1.basis_eta = basis_eta;
        myfitstr1.zeta = oldzeta;       
 
%         [tbs,oldfit] = ode45(@Sim1WorkingOde, tobs, old.y0,odeoptions,myfitstr1);
%          oldmumat(:,i) = oldfit;      
       oldfit = oldmumat(:,i);

        myfitstr1.zeta = newzeta;
        [tbs,newfit] = ode45(@Sim1WorkingOde, tobs, old.y0,odeoptions,myfitstr1);
        newfitmat(:,i) = newfit;
        
        oldSSE  = oldSSE-0.5*oldtauE*sum((Yobs-oldfit).^2);
        newSSE  = newSSE-0.5*oldtauE*sum((Yobs-newfit).^2);

    end

    oldprior  = -0.5*oldlambda_eta*(oldzeta(2:K)'*D_eta*oldzeta(2:K));
    newprior  = -0.5*oldlambda_eta*(newzeta(2:K)'*D_eta*newzeta(2:K));

 
    num1    = newSSE+newprior;
    den1    = oldSSE+oldprior;
    
    if log(rand(1))<=(num1-den1)
        oldzeta    = newzeta;
        countzeta  = countzeta+1;
        oldmumat = newfitmat;
    end

    
    %% Sample lambda

    oldprior  = oldzeta(2:K)'*D_eta*oldzeta(2:K);
    oldlambda_eta = gamrnd(shape_lambda+(K-2-1)/2, 1/(0.5*oldprior+1/scale_lambda));
    
    %% Sample theta;

    MuThetaI      = sum(oldthetaImat,2);
    W             = ncohort*oldInvSigma+InvOmega;
    postW         = inv(W);
    temp1         = oldInvSigma*MuThetaI+InvOmega*reshape(eta0,p,1);
    postmu        = W\temp1;
 
    oldtheta      = (mvnrnd(postmu, postW))';

    %% Sample inverse of SigmaThetaI;

    diffthetaI   = oldthetaImat-oldtheta*ones(1,ncohort);
    newOmega     = diffthetaI*diffthetaI'+S01;
    oldSigma = iwishrnd(newOmega, ncohort+df);
    oldInvSigma = inv(oldSigma);
    

    %% Sample tauE

    diffmumat = Yobsmat-oldmumat;
    SSE      = sum(sum(diffmumat.^2));
    oldtauE = gamrnd(0.5*ncohort*nobs+shape_tauE, 1/(0.5*SSE+1/scale_tauE));
     
    
    %% Output.
    
    out1 = input1;
    
    out1.oldtheta     = oldtheta;
    out1.oldthetaImat = oldthetaImat;
    out1.oldmumat     = oldmumat;
    out1.oldzeta      = oldzeta;
    out1.oldtauE      = oldtauE;
    out1.oldlambda    = oldlambda_eta;
    
    out1.oldInvSigma = oldInvSigma;
    
    out1.countthetaI  = countthetaI;
    out1.countzeta    = countzeta; 
    
     

 end

