function drdt = myode_NOinput(t,r,J,N)
%I = interp1(It,I,t); % Interpolate the data set (gt,g) at time t

drdt = -r+fIcurve(J*r); % Evaluate ODE at time t