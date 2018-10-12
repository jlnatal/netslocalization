function drdt = myode_timedep(t,r,J,N,input_IDs,noisestream_smooth)

activa=fIcurve(J*r+input_funct(t,N,input_IDs));
activa=activa*noisestream_smooth(t);
drdt = -r+activa; % Evaluate ODE at time t