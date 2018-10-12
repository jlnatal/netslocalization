function [J,positions,distances,N,dale,cutoff]=configure_topology(N,cutoff)

MU=-0.702;                  
V=0.8752;
SIGMA=sqrt(V);
  
       M = log(1)+exp(MU + SIGMA^2/2);                         %mean
       V = exp(2*MU + SIGMA^2) * (exp(SIGMA^2) - 1);    %var
       
positions=5*rand(N,2);  % 5mm x 5mm !


%%% UNIFORM TEST
%[X,Y] = meshgrid(linspace(0,5,32),linspace(0,5,32));
%positions=[X(:) Y(:)];

%distances=squareform(pdist(positions));
%scatter(positions(:,1),positions(:,2))

Xs=squareform(pdist(positions(:,1)));
XsADJUSTED=5-squareform(pdist(positions(:,1)));
Ys=squareform(pdist(positions(:,2)));
YsADJUSTED=5-squareform(pdist(positions(:,2)));
distances=sqrt((min(Xs,XsADJUSTED)).^2+(min(Ys,YsADJUSTED)).^2);



%cutoff=5/sqrt(4900)+.005;
%cutoff=5/sqrt(4900)-.005;
%cutoff=0.3;

%Connection matrix
%J=4900*lognrnd(MU,SIGMA,N)/N;
J=lognrnd(MU,SIGMA,N);

dale = 2*(rand(1,N)<1)-1;
%dale = 2*(rand(1,N)<0.8)-1;
%signs=2*(rand(N)<0.5)-1;
signs=repmat(dale,N,1);


excitatory_max=0.3;
inhibitory_max=0.9;

%J=J.*signs.*((distances>0.3 & distances<inhibitory_max).*(signs<0)+(distances<=excitatory_max).*(signs>0));
J=J.*signs.*(distances<cutoff);


J=J-diag(diag(J));
%J=J*10;

%J(J>0)=1;
%J(J<0)=-1;






%%% NO...
%J=J./sum(J(1,:)~=0);


%J=rand(N);

%figure;
%gplot(J>0,positions,'b')
%hold on
%gplot(J<0,positions,'r')
