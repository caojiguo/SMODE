function out1 = Sim1MCMC_funTT(input1)
% Note: this program is built on Matlab (at least 2014 version) due to
% using some truncated distributions.

        ncohort      = input1.ncohort;
        oldthetaImat = input1.oldthetaImat;
        oldtheta     = input1.oldtheta;
        oldmumat     = input1.oldmumat;
        oldzeta      = input1.oldzeta;
        oldUmat      = input1.oldUmat;
        oldWmat      = input1.oldWmat;  
        oldtauE      = input1.oldtauE;
 
        oldInvSigma = input1.oldInvSigma;
        oldlambda_eta  = input1.oldlambda_eta;  
        oldkappa     = input1.oldkappa;
        oldnu        = input1.oldnu;  

        sig_thetaI    = input1.sig_thetaI;
        sig_zeta    = input1.sig_zeta;
        sig_kappa    = input1.sig_kappa;
        sig_nu    = input1.sig_nu;
        
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
        countkappa  = input1.countkappa;
        countnu     = input1.countnu;
        
        odeoptions  = input1.odeoptions;
        

        %% Update theta_i
        
    for i=1:ncohort

        oldthetaI   = oldthetaImat(:,i);   
        %newthetaI   = oldthetaI+normrnd(0,sig_thetaI,p,1);
        newthetaI   = mvnrnd(oldthetaI, sig_thetaI{i})';

        oldU = oldUmat(i);
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
        oldSSE   = -0.5*oldU*oldtauE*sum(olddiff.^2);
        oldprior = -0.5*oldWmat(i)*(oldthetaI-oldtheta(1:p))'*oldInvSigma*(oldthetaI-oldtheta(1:p));
        den1     = oldSSE + oldprior;

        new.y0      = newthetaI(1);
        myfitstr1.alpha = exp(newthetaI(2));
        myfitstr1.beta = exp(newthetaI(3));
        myfitstr1.delta = exp(newthetaI(4));
        
        
        [tbs,newfit] = ode45(@Sim1WorkingOde, tobs, new.y0,odeoptions,myfitstr1);
        
        newdiff  = (Yobs-newfit);  
        newSSE   = -0.5*oldU*oldtauE*sum(newdiff.^2);
        newprior = -0.5*oldWmat(i)*(newthetaI-oldtheta(1:p))'*oldInvSigma*(newthetaI-oldtheta(1:p));
        num1     = newSSE + newprior;
        
%         display(num2str(num1))
%         
%        display(num2str(den1))
%         
%        display(num2str(num1-den1))
        
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
        oldU = oldUmat(i);
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
        
        oldSSE  = oldSSE-0.5*oldU*oldtauE*sum((Yobs-oldfit).^2);
        newSSE  = newSSE-0.5*oldU*oldtauE*sum((Yobs-newfit).^2);

    end

    oldprior  = -0.5*oldlambda_eta*(oldzeta(2:K)'*D_eta*oldzeta(2:K));
    newprior  = -0.5*oldlambda_eta*(newzeta(2:K)'*D_eta*newzeta(2:K));

 
    num1    = newSSE+newprior;
    den1    = oldSSE+oldprior;
    
%        display(num2str(num1))
%         
%        display(num2str(den1))
%         
%        display(num2str(num1-den1))
    
    if log(rand(1))<=(num1-den1)
        oldzeta    = newzeta;
        countzeta  = countzeta+1;
        oldmumat = newfitmat;
    end

    %% Sample ui
    
    olddiff = Yobsmat-oldmumat;
    resids2 = diag(olddiff'*olddiff);
    oldUmat = gamrnd(oldkappa/2+nobs/2,2./(oldkappa+oldtauE.*resids2),ncohort,1);
    
    %% Sample wi
    temp1 = oldthetaImat-oldtheta*ones(1,ncohort);
    temp2 = diag(temp1'*oldInvSigma*temp1);
    
    oldWmat = gamrnd(oldnu/2+p/2,2./(oldnu+temp2),ncohort,1);
 
    
    %% Sample kappa
    
    pd = makedist('normal',log(oldkappa),sig_kappa);
    newnormpdf = truncate(pd,log(2.5),inf);
    newkappa = exp(random(newnormpdf,1)); 
    
    num1 = ncohort*(0.5*newkappa*log(0.5*newkappa)-gammaln(0.5*newkappa))+(0.5*newkappa-1)*sum(log(oldUmat))-0.5*newkappa*sum(oldUmat)-1/3*newkappa;
    den1 = ncohort*(0.5*oldkappa*log(0.5*oldkappa)-gammaln(0.5*oldkappa))+(0.5*oldkappa-1)*sum(log(oldUmat))-0.5*oldkappa*sum(oldUmat)-1/3*oldkappa;
    if log(rand(1))<=(num1-den1)
        oldkappa = newkappa;
        countkappa = countkappa+1;
    end

   %% Sample nu
   
    pd = makedist('normal',log(oldnu),sig_nu);
    newnormpdf = truncate(pd,log(2.5),inf);
    newnu = exp(random(newnormpdf,1)); 

    num1 = ncohort*(0.5*newnu*log(0.5*newnu)-gammaln(0.5*newnu))+(0.5*newnu-1)*sum(log(oldWmat))-0.5*newnu*sum(oldWmat)-1/3*newnu;
    den1 = ncohort*(0.5*oldnu*log(0.5*oldnu)-gammaln(0.5*oldnu))+(0.5*oldnu-1)*sum(log(oldWmat))-0.5*oldnu*sum(oldWmat)-1/3*oldnu;
    if log(rand(1))<=(num1-den1)
        oldnu = newnu;
        countnu = countnu+1;
    end

       
    %% Sample lambda_eta

    oldprior  = oldzeta(2:K)'*D_eta*oldzeta(2:K);
    oldlambda_eta = gamrnd(shape_lambda+(K-2-1)/2, 1/(0.5*oldprior+1/scale_lambda));
    
    %% Sample theta;

    MuThetaI      = oldthetaImat*oldWmat;
    W             = sum(oldWmat)*oldInvSigma+InvOmega;
    postW         = inv(W);
    temp1         = oldInvSigma*MuThetaI+InvOmega*reshape(eta0,p,1);
    postmu        = W\temp1;
    oldtheta      = (mvnrnd(postmu, postW))';

    %% Sample inverse of SigmaThetaI;

    diffthetaI   = oldthetaImat-oldtheta*ones(1,ncohort);
    newOmega     = diffthetaI*diag(oldWmat)*diffthetaI'+S01;
    oldSigma = iwishrnd(newOmega, ncohort+df);
    oldInvSigma = inv(oldSigma);

    %% Sample tauE

    diffmumat = Yobsmat-oldmumat;
    SSE      = sum(oldUmat'.*sum(diffmumat.^2));

    oldtauE = gamrnd(0.5*ncohort*nobs+shape_tauE, 1/(0.5*SSE+1/scale_tauE));
    
    %% Output.
    
    out1 = input1;
    
    out1.oldtheta     = oldtheta;
    out1.oldthetaImat = oldthetaImat;
    out1.oldmumat     = oldmumat;
    out1.oldzeta      = oldzeta;
    out1.oldUmat      = oldUmat;
    out1.oldtauE      = oldtauE;
    out1.oldkappa     = oldkappa;
    out1.oldlambda    = oldlambda_eta;
    out1.oldWmat      = oldWmat;
    out1.oldnu        = oldnu;
    
    out1.oldInvSigma = oldInvSigma;
    
    out1.countthetaI  = countthetaI;
    out1.countzeta    = countzeta; 
    out1.countkappa   = countkappa; 
    out1.countnu      = countnu; 
  
        

 end

