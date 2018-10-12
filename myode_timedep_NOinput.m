function drdt = myode_timedep_NOinput(t,r,J,N,noisestream_smooth)

activa=fIcurve(J*r);
activa=activa*noisestream_smooth(t);
%drdt = -r+activa; % Evaluate ODE at time t
%drdt = ceil(-r+activa); % Evaluate ODE at time t