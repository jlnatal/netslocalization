function I=input_funct(t,N,input_IDs)
offset=5;
strength=100;
I=zeros(N,1);
%I=repmat(I,[50 1]);
%I=[I; zeros(N-50,1)];
I(input_IDs)=strength*(1-heaviside(t-offset));