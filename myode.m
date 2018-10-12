function drdt = myode(t,r,J,N,input_IDs)
%I = interp1(It,I,t); % Interpolate the data set (gt,g) at time t

activa=fIcurve(J*r+input_funct(t,N,input_IDs));
%drdt = -r+activa+.02-sum(activa)/N; % Evaluate ODE at time t
%drdt = -r+fIcurve(J*r+input_funct(t,N,input_IDs)); % fullwording
%.02*N*(max(activa))/sum(activa)
drdt = -r+.02*N*(activa)/sum(activa);