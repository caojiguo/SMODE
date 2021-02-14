function dy=Sim1Ode(t,y,input)

     alpha   = input.alpha;
     beta   = input.beta;
     delta   = input.delta;

     eta = exp(-(t-2).^2/1.44)+exp(-(t-5).^2)+0.5*exp(-(t-6).^2/0.64)-0.062;
     
     dy  = -delta*y+alpha+beta*eta;
end