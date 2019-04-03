%CGpipeline
%close all; clear; clc; CGpipeline

 
%Load neural data
load bint_fishmovie32_100
data=reshape(bint(297,:,:),[160,953]);
data=data*2-1;
mean(sum(data>0,2)/953);
mean(mean(data,2));

 
%Pre-allocate space for saving each GC step
CGevents=3;            %Number of iterations
C{1,CGevents+1}=[];
C{1}=data;
Rec=data;
%R{CGevents}=[];

%Break down each neuron's trajectory as member in a cell array
for g=1:size(C{1},1)
c{g}=C{1}(g,:);
end
D{1}=c;
%masterlist{1}=c;

%"Raster plot"-style image of original data matrix
figure,imagesc(C{1})
colormap(flipud(gray))



%Main loop

for k=1:CGevents
k

 
%Find most highly correlated pair
[mosthigh,hival]=find_most_highly(data);
MostHi(k,:)=mosthigh;

%Calculate stats for most highly correlated pair
v1=data(mosthigh(1),:); v2=data(mosthigh(2),:);
means=[mean(data(mosthigh(1),:)) mean(data(mosthigh(2),:))];
pairwise=mean(data(mosthigh(1),:).*data(mosthigh(2),:));
covarnce=pairwise-mean(data(mosthigh(1),:))*mean(data(mosthigh(2),:));

%Break the loop if <s_1*s_2> is negative or zero
if pairwise<=0
    fprintf('Non-positive pairwise correlation')
    k
    break   %Violates assumptions; stop here
end


%Analytical solution equations
tanhJ0sqd=means(1)*means(2)/pairwise;
tanhJ1sqd=means(1)*pairwise/means(2);
tanhJ2sqd=means(2)*pairwise/means(1);

 
%Break the loop if any of above is >1
if sum([tanhJ0sqd tanhJ1sqd tanhJ2sqd]>1)
    fprintf('Error: Component > 1')
    k
    break
end

 %Choose values of couplings (as interdependent, from each other)
J0=atanh(-tanhJ0sqd^.5);       %choose minus root as default
%J1=atanh(tanhJ1sqd^.5);
J1=atanh(means(1)/tanh(J0));
%J2=atanh(tanhJ2sqd^.5);
J2=atanh(pairwise/tanh(J1));
%mean(v1)^2+mean(v2)^2

%Correlation for current pair -- check
corry=(pairwise-mean(v1)*mean(v2))/sqrt(var(v1)*var(v2));
[J0 J1 J2]

%In case any of above didn't catch it...
if sum((imag([J0 J1 J2]))>0)
    fprintf('Imaginary couplings')
    k
    break

end

%Run script to propagate one step forward and backward
%Calculate sigma (hidden) from v1,v2; then in reverse
logistic_through
Meansig(k)=meansig;
data(mosthigh,:)=[];
data=[data; samples];
newsigs{k}=samples;
%masterlist{k+1}=[masterlist{k}(setdiff(1:length(masterlist{k}),MostHi(:))) newsigs];


%Not in use currently; reconstruction of v1, v2 statistics
[insertLoHi,indexes]=sort(mosthigh);
    if (sum(insertLoHi)~=1) & (sum(insertLoHi)~=size(data,1))
        rec=[data(1:insertLoHi(1)-1,:); vis(indexes(1),:); data(insertLoHi(1):insertLoHi(2)-1,:); vis(indexes(2),:); data(insertLoHi(2):end-1,:)];
    else
        fprintf('Error in reconstruction')
        k
    end

 

%Save couplings each iteration
couplings(k,:)=[J0 J1 J2];
covthrutime(k)=covarnce;

 

%Compare probability distribution ((1,1),(1,-1),(-1,1), & (-1,-1))
%s_1 and s_2 with reconstruction for a given CG iteration
%fprintf('Individual reconstructions: sample stats')
JPD_compar(:,:,k)=[[sum((v1+v2)==2) sum((v1==1) & (v2==-1)) sum((v1==-1) & (v2==1)) sum((v1+v2)==-2)]/953; [sum((vis(1,:)+vis(2,:)==2)) sum((vis(1,:)==1) & (vis(2,:)==-1)) sum((vis(1,:)==-1) & (vis(2,:)==1)) sum((vis(1,:)+vis(2,:))==-2)]/953];
%figure
%bar(JPD_compar(:,:,k)')
%abs(JPD_compar(1,:,k)-JPD_compar(2,:,k))
%errorbar(mean(JPD_compar(:,:,k),3),mean(abs(JPD_compar(1,:,:)-JPD_compar(2,:,:)),3))

%errorbar(mean(JPD_compar(1,:,k),3),mean(abs(JPD_compar(1,:,:)-JPD_compar(2,:,:)),3))
%hold on
%errorbar(mean(JPD_compar(2,:,k),3),mean(abs(JPD_compar(1,:,:)-JPD_compar(2,:,:)),3))
 

%figure
%bar(1:4,[sum((v1+v2)==2) sum((v1==1) & (v2==-1)) sum((v1==-1) & (v2==1)) sum((v1+v2)==-2)]/953,0.5,'FaceColor',[0.2 0.2 0.5]); hold on
%bar(1:4,[sum((vis(1,:)+vis(2,:)==2)) sum((vis(1,:)==1) & (vis(2,:)==-1)) sum((vis(1,:)==-1) & (vis(2,:)==1)) sum((vis(1,:)+vis(2,:))==-2)]/953,.25,'FaceColor',[0 0.7 0.7])



%Save the coarse-grained data matrix

    C{1,k+1}=data;
    
    %size(D{k})
    included{k}=setdiff(1:length(D{k}),MostHi(k,:));
    D{k+1}=[D{k}(included{k}) newsigs{k}];  %add new neurons to array
    %setdiff(1:length(D{k}),MostHi(k,:))
    

end

    



   %Building the Reconstructions for each iteration backwards
R=D;
    for k=CGevents:-1:1
        R{k}(included{k})=R{k+1}(1:end-1);
        activ=[cell2mat(R{k+1}(end)'); cell2mat(R{k+1}(end)')].*couplings(k,2:3)';
        recons=(exp((activ))./(2.*cosh(activ))>rand(2,num_samps))*2-1;
        R{k}(MostHi(k,:))={recons(1,:),recons(2,:)};
        %imagesc(cell2mat(R{k}')==cell2mat(D{k}'))
        %find(~ismember(cell2mat(D{1}'),cell2mat(R{1}'),'rows'))
    end
    
    
        figure
        imagesc(cell2mat(R{1}')==C{1})
        title('Entries that differ between data and reconstruction')
        colormap('gray')
    

    
    
        
    %[covC_raw,covC]=cov_matrix(C{1});
    [covDat_raw,covDat]=cov_matrix(cell2mat(D{1}'));
    [covRec_raw,covRec]=cov_matrix(cell2mat(R{1}'));
    [covCG_raw,covCG]=cov_matrix(cell2mat(D{end}'));
    
    %figure,imagesc(covC); title('Original Data Matrix')
    figure,imagesc(covDat); title('Original Data'); colormap(flipud(gray))
    figure,imagesc(covRec); title('Reconstructed Data'); colormap(flipud(gray))
    figure,imagesc(covCG); title('Coarse-Grained Data'); colormap(flipud(gray))
        %xlim([0 160]); ylim([0 160]);
    
    %figure,imagesc(cell2mat(R{3}')==C{3})
    %figure,imagesc(cov_matrix(cell2mat(R{3}'))==cov_matrix(C{3}))

    
    
    
    
    howmanytop=10;  %Top __ covariance entries (in data/recons)
    
    [covDat_unique indD]=sort(unique(covDat(~isnan(covDat))),'descend');
    covDat_top=covDat_unique(1:howmanytop);
    [rD,cD]=find(ismember(triu(covDat),covDat_top(1:howmanytop)));
    
    [covRec_unique]=sort(unique(covRec(~isnan(covRec))),'descend');
    covRec_top=covRec_unique(1:howmanytop);
    [rR,cR]=find(ismember(triu(covRec),covRec_unique(1:howmanytop)));

    
    %Lists of top entries and locations
    
    for topindsD=1:howmanytop
        topsD(topindsD)=covDat(rD(topindsD),cD(topindsD));
    end
        [topsD,indsD]=sort(topsD','descend');
        rD=rD(indsD); cD=cD(indsD);
    
    for topindsR=1:howmanytop
        topsR(topindsR)=covDat(rR(topindsR),cR(topindsR));
    end
        [topsR,indsR]=sort(topsR','descend');
        rR=rR(indsR); cR=cR(indsR);
    
        
    %Side-by-side plot of top cov's for original data vs. recons
    figure
    subplot(1,2,1)
    imagesc(covDat>=covDat_top(end)); colormap(flipud(gray));
    subplot(1,2,2)
    imagesc(covRec>=covRec_top(end)); colormap(flipud(gray));
    figure,imagesc(covRec>=covRec_top(end) == covDat>=covDat_top(end)); colormap(flipud(gray));
       

    %Joint probability distributions w/ average subtractive errors
bar(mean(JPD_compar(1,:,k),3),0.5,'FaceColor',[0.2 0.2 0.5]); hold on; bar(mean(JPD_compar(2,:,k),3),.25,'FaceColor',[0 0.7 0.7])
errorbar(mean(JPD_compar(1,:,k),3),mean(abs(JPD_compar(1,:,:)-JPD_compar(2,:,:)),3),'LineStyle','none','MarkerSize', 20)



%Reconstructed (different) elements vs data matrix
    %for k=CGevents+1:-1:1
    for k=1:CGevents+1
        figure
        imagesc(cell2mat(R{k}')~=C{k})
    end
        
    
%Reconstructed elements of covariance matrix
    %for k=CGevents+1:-1:1
    for k=1:CGevents+1
        figure
        imagesc(cov_matrix(cell2mat(R{k}'))~=cov_matrix(cell2mat(D{k}')))
    end
        