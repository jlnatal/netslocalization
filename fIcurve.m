function y=fIcurve(x)

    y=18*(log(1+log(1+exp(.5*x-8)))).^1.5;
    %y=18*(log(1+log(1+exp(.5*x-8*ones(1,length(x)))))).^1.5;
    %y=(1.715/100)*18*(log(1+log(1+exp(.5*x-8)))).^1.5; new?
    %y=.02*(length(x)*y)/norm(y,1);
    %y=3*tanh(.05*(x-30)+1.5);
    
    % No longer needed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%a_new=0.1+0.1*((rand-0.5)/5);
    %a_new=normrnd(.02,.0075);
    %y=(a_new)*(length(x)*y)/norm(y,1);
    %%y=(0.2*N*y)/norm(y,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end