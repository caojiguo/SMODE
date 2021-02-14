function dy=Sim1WorkingOde(t,y,input)

     basis_eta = input.basis_eta;
     zeta      = input.zeta;
     
     alpha   = input.alpha;
     beta   = input.beta;
     delta   = input.delta;
     
     A   = eval_basis(t, basis_eta);
     eta = A*zeta;
     
     dy  = -delta*y+alpha+beta*eta;
     
end