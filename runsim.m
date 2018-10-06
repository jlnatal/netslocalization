N=2^12;
restr=0.3;

animateon=0;

%[a,center_of_excitation,cents,culprits,cutoff,dale,distances,i,input_IDs,J,N,positions,r,radii,t,target,FINAL_STATE,Firing_at_rad,howfar,original,asmany,cee,cutstring,filenaming,endtime_in_mins]=get_data(2^10,0.3);
%%%


cutoff=0.3;
noiseon=0
if noiseon
    %noiselevel=0.001;
    noiselevel=0.1;
else
    noiselevel=0;
end


number_trials=250;
[X,Y] = meshgrid(0:.1:4.9,0:.1:4.9);
allgridpts=[X(:) Y(:)];
gridpositions=randi(length(allgridpts),number_trials,1);
% At number_trials=2500, hits most of the 2500 gridpoints
% At number_trials=15000, hits consistently in the 2490's

asmany=length(gridpositions)
%asmany=input('How many trials? (Default 250)  ');
%if isempty(asmany)
%    asmany=250;
%end

secondz=200;
tspan=[0:0.1:secondz];
Trajects=zeros(length(tspan),N,asmany);

%%% Moved to bottom:
%%%%
%cee=clock;
%if cutoff<1
%    cutstring=['0pt' num2str(cutoff*10)];
%elseif (cutoff>1 & cutoff<2)
%    cutstring=['1pt' num2str(10*rem(cutoff,1))];
%end
%filenaming=['Data' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5)) num2str(floor(cee(6))) '__N' num2str(N)  '__cutoff' num2str(cutstring) '___' num2str(asmany) 'runs']
%save(filenaming)
%%%

%Gain for dynamical threshold function
g=0.5;

%Sparsity (fraction of active neurons in each memory)
%a=0.02;


%Generate E/I topology with cutoff distance 0.3 on 5x5 pexperne
[J,positions,distances,N,dale,cutoff]=configure_topology(N,cutoff);



target=dsearchn(positions,[2.5 2.5]);
cents=rangesearch(positions,[positions(target,1) positions(target,2)],restr);
culprits=cents{1,1};
input_IDs=culprits;
%figure,scatter(positions(culprits,1),positions(culprits,2))
%ylim([0 5])
%xlim([0 5])


%a=length(input_IDs)/N;
%a=0.2;
a=0.02;


vi=zeros(1,N);
%shoots=randperm(N,floor(a*N));
%vi(shoots)=1/(a*N);
for i=1:size(input_IDs,2)
%vi(input_IDs(i))=memory(1,input_IDs(i));
%vi(input_IDs(i))=N*a/length(input_IDs);
%vi(input_IDs(i))=1/size(input_IDs,2);
vi(input_IDs(i))=1+vi(i);
end
%vi=memory(1,:);
%vi=vi/norm(vi,1);
%vi=vi*(a/length(input_IDs));


%gain=g*ones(N,1);
%gain([input_IDs])=1.5*g;

%%
%[r,t,center_of_excitation]=dyn_contin(J,N,gain,vi,positions,target,a,restr);
%stable_state=r(end,:);
%%


%FINAL_STATE=zeros(length(input_IDs),N);
%FINAL_STATE=zeros(N);


howfar=0;
%howfar=[0:.01:1];
%howfar=0:0.01:.5;
%howfar=0:0.025:1;
%howfar=0:0.5:1.2;
%asmany=10;
corr_at_range=zeros(length(howfar),asmany);
corr_at_range_END=zeros(length(howfar),asmany);

tic


%original=zeros(401,N);

%Initialize system with all rates = a, then let it run without input for t=1000
r0=a*ones(1,N);
if noiseon
    noisestream=1+cumsum(noiselevel*randn(1000,1));
    noisestream_smooth=fit((1:1000)',noisestream,'linearinterp');
    [t,r] = ode45(@(t,r) myode_timedep_NOinput(t,r,J,N,noisestream_smooth), tspan,r0);
else
    [t,r] = ode45(@(t,r) myode_NOinput(t,r,J,N), [0:0.1:1000],r0);
end



%Run the system from the ``equilibration" of random initialization
exper=0

    %ah=(lo-1)/stepper;
   
    r0=r(end,:);
    %r0=a*ones(1,N);
    %r0=stable_state;
    
    
    
    reach=0;
    
    %target=dsearchn(positions,rand(1,2)*2.5+1.25);
    %target=tohit(lo)
    %target=1
    %cents=rangesearch(positions,[positions(target,1) positions(target,2)],.3);
%    cents=rangesearch(positions,[1.75+howfar 2.5],.3);
%    culprits=cents{1,1};
%    input_IDs=culprits;
    %target=input_IDs(lo)

    ang=2*pi*rand;
    x_coord=reach*cos(ang);
    y_coord=reach*sin(ang);

    
    cents=rangesearch(positions,[2.5+x_coord 2.5+y_coord],.3);
    culprits=cents{1,1};
    input_IDs=culprits;
    
    %BLOCK
    %hold on
    %figure(103),scatter(positions(input_IDs,1),positions(input_IDs,2),50,'k')
    %pause(0.5)
    %%hold off
    
   
    
    

%tspan=[0 100];


%Running just once after ``equilibration" as described above
if noiseon
    noisestream=1+cumsum(noiselevel*randn(tspan(end),1));
    noisestream_smooth=fit((1:tspan(end))',noisestream,'linearinterp');
    [t,r] = ode45(@(t,r) myode_timedep(t,r,J,N,input_IDs,noisestream_smooth), tspan,r0);
else
    [t,r] = ode45(@(t,r) myode(t,r,J,N,input_IDs), tspan,r0);
end
%figure,plot(t,r)

    
        for la=1:10:length(tspan)


            if animateon==1
                animate
            end

            if la==length(tspan)
                cutthresh=0.1;
                %viscircles([positions(target,1) positions(dsearchn(positions,[2.5 2.5]),2)],restr,'Color','k');
                center_of_excitation=[sum((r(end,:)>cutthresh).*positions(:,1)')/sum((r(end,:)>cutthresh)) sum((r(end,:)>cutthresh).*positions(:,2)')/sum((r(end,:)>cutthresh))];
                %drawnow
                %xlim([0 5])
                %ylim([0 5])


            elseif la==1
                %figure(101),hold on;
                %scatter(positions(input_IDs,1),positions(input_IDs,2),100,'k','filled')
                %pause(.5)
                %hold off
                %set(gca,'FontSize',18)


            end


        end
    

        original=r;
        %fprintf('r becoming original')
        %r([1:5])

        top=find(r(end,:)==max(r(end,:)));

        radii_at_which_to_measure=0.1:0.1:7;
        firing_at_rad=zeros(1,length(radii_at_which_to_measure));

        for ay=1:length(radii_at_which_to_measure)
            within_ring=rangesearch(positions,positions(top,:),radii_at_which_to_measure(ay));

            in_ring=within_ring{1,1};
            firing_at_rad(ay)=mean(r(end,in_ring));
        end

        %disty=[[sum((r(end,:)>a).*positions(:,1)')/sum((r(end,:)>a)) sum((r(end,:)>a).*positions(:,2)')/sum((r(end,:)>a))]] - [2.5+x_coord 2.5+y_coord];

        
    %corr_at_range(1,exper)=sum(mean(original.*r))/N;
    %corr_at_range_END(1,exper)=mean(original(end,:).*r(end,:));
    
    
    
    



for lo=1:length(howfar)
    
    howfar(lo)
    
    reach=howfar(lo);
    
    toc
  
%BLOCK    
%figure(103),scatter(positions(:,1),positions(:,2),50,zeros(1,N),'filled')























final_states=zeros(asmany,N);



stim_position=zeros(asmany,2);
%Vector_dists=


%centerz=zeros(asmany,2);




    for exper=1:asmany
        [fprintf('exper') num2str(exper)]
        

    %ah=(lo-1)/stepper;
   
    
    %Run the system from the previous attractor it settled into
    r0=r(end,:);
    %r0=a*ones(end);
    %r0=stable_state;


    %target=dsearchn(positions,rand(1,2)*2.5+1.25);
    %target=tohit(lo)
    %target=1
    %cents=rangesearch(positions,[positions(target,1) positions(target,2)],.3);
%    cents=rangesearch(positions,[1.75+howfar 2.5],.3);
%    culprits=cents{1,1};
%    input_IDs=culprits;
    %target=input_IDs(lo)

    cents=rangesearch(positions,[allgridpts(gridpositions(exper),1) allgridpts(gridpositions(exper),2)],.3);
    culprits=cents{1,1};
    input_IDs=culprits;
    
    %BLOCK
    %hold on
    %figure(103),scatter(positions(input_IDs,1),positions(input_IDs,2),50,'k')
    %pause(0.5)
    %%hold off
    
    
    
    
    
    
    
    %fprintf('r checking still same')
    %r([1:5])
    

%tspan=[0 100];

if noiseon
    noisestream=1+cumsum(noiselevel*randn(tspan(end),1));
    noisestream_smooth=fit((1:tspan(end))',noisestream,'linearinterp');
    [t,r] = ode45(@(t,r) myode_timedep(t,r,J,N,input_IDs,noisestream_smooth), tspan,r0);
else
    [t,r] = ode45(@(t,r) myode(t,r,J,N,input_IDs), tspan,r0);
end
%figure,plot(t,r)
%fprintf('r should be changed now')
%r([1:5])
    
        for la=1:10:length(tspan)


            if animateon==1
                animate
            end

            if la==length(tspan)
                cutthresh=0.1;
                %viscircles([positions(target,1) positions(dsearchn(positions,[2.5 2.5]),2)],restr,'Color','k');
                center_of_excitation=[sum((r(end,:)>cutthresh).*positions(:,1)')/sum((r(end,:)>cutthresh)) sum((r(end,:)>cutthresh).*positions(:,2)')/sum((r(end,:)>cutthresh))];
                %drawnow
                %xlim([0 5])
                %ylim([0 5])


            elseif la==1
                %figure(101),hold on;
                %scatter(positions(input_IDs,1),positions(input_IDs,2),100,'k','filled')
                %pause(.5)
                %hold off
                %set(gca,'FontSize',18)


            end


        end

 
    
    
   
 
    
    
    
    
    
    
    %FINAL_STATE(lo,:)=r(end,:);
    
    %if lo==1
        %if exper==1
        %%original=r(end,:);
        %original=r;
        %end
   % end
    
    
        %BLOCK
        %figure(102)

        
        top=find(r(end,:)==max(r(end,:)));

        radii_at_which_to_measure=0.1:0.1:7;
        firing_at_rad=zeros(1,length(radii_at_which_to_measure));

        for ay=1:length(radii_at_which_to_measure)
            within_ring=rangesearch(positions,positions(top,:),radii_at_which_to_measure(ay));
            in_ring=within_ring{1,1};
            firing_at_rad(ay)=mean(r(end,in_ring));
        end

        %scatter(radii_at_which_to_measure,firing_at_rad)
        %hold on
        
        
        
        
        %Firing_at_rad(lo,:)=firing_at_rad;

        %runname=['distval' num2str(lo) 'expno' num2str(exper)]
        %eval([runname '= r;']);
        %R(:,:,lo,exper)=r;
%       %save(filenaming,runname,'-append');
        Trajects(:,:,exper)=r;
        
        %figure(104)
        %errorbar(0.1:0.1:7,mean(Firing_at_rad),var(Firing_at_rad))
      
        %Old???
        %disty=[[sum((r(end,:)>a).*positions(:,1)')/sum((r(end,:)>a)) sum((r(end,:)>a).*positions(:,2)')/sum((r(end,:)>a))]] - [2.5+x_coord 2.5+y_coord];
        

        
    corr_at_range(lo,exper)=sum(mean(original.*r))/N;
    corr_at_range_END(lo,exper)=mean(original(end,:).*r(end,:));
    %corr_at_range_END(lo,exper)=mean(original(end,:)>a.*r(end,:)>a);
    %[original(end,1:5); r(end,1:5)]
    
    final_states(exper,:)=r(end,:);
    %centerz(exper,:)=[sum((final_states(exper,:)>a).*positions(:,1)')/sum((final_states(exper,:)>a),2) sum((final_states(exper,:)>a).*positions(:,2)')/sum((final_states(exper,:)>a),2)];
    COMperiodic_circlemethod; centaz(exper,:)=centa;
    
    
    end 


radii=[0:0.1:sqrt(2)*5];
%ranges




end

    curve=mean(corr_at_range,2);
    curvy=mean(corr_at_range_END,2);
    var_temporal=var(corr_at_range,0,2);
    var_stable=var(corr_at_range_END,0,2);

    %corr_at_range_END
    %sum(sum(r==original))
    
    raw_deviations=[centaz(:,1)-allgridpts(gridpositions,1) centaz(:,2)-allgridpts(gridpositions,2)];
    distances_traveled=[min(abs(raw_deviations(:,1)),5-abs(raw_deviations(:,1))) min(abs(raw_deviations(:,2)),5-abs(raw_deviations(:,2)))];
    
        
endtime_in_mins=toc/60
%save(filenaming)

cee=clock;
if cutoff<1
    cutstring=['0pt' num2str(cutoff*10)];
elseif (cutoff>1 & cutoff<2)
    cutstring=['1pt' num2str(10*rem(cutoff,1))];
end
filenaming=['Data' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5)) num2str(floor(cee(6))) '__N' num2str(N)  '__cutoff' num2str(cutstring) '___' num2str(asmany) 'runs' '___' 'noiselevel' num2str(noiselevel*1000) 'thousandths']
%save(filenaming)









%b=combnk(1:4,2)
e=rand(asmany,N)<(mean(sum(final_states>1,2))/N);

D = pdist(final_states>1,'hamming');
%eucd = pdist(centerz);
eucd = pdist(centaz);
E=pdist(e,'hamming');

overlap=(N-2*D)/N;


figure(101),scatter(eucd,E,[],[0.5 0.5 0.5])
hold on
scatter(eucd,D,[],'b')
%title('Geometrically distant stable states represent distinct firing patterns')
title('Distances between stable bump pairs')
xlabel('Euclidean distance')
ylabel('Hamming distance')
set(gca,'FontSize',18)
hold off

figure(102),scatter(eucd,(N-2*E)/N,[],[0.5 0.5 0.5])
hold on
scatter(eucd,(N-2*D)/N,[],'b')
%title('Geometrically distant stable states represent distinct firing patterns')
title('Overlap for pairs of stable firing patterns')
xlabel('Euclidean distance')
ylabel('Hamming distance')
set(gca,'FontSize',18)
hold off


%figure,scatter(positions(:,1),positions(:,2),50,sum(final_states>1),'filled')
%colormap(flipud(gray))
%colorbar
%hold on
%%%scatter(centerz(:,1),centerz(:,2))
%[Au,~,ic] = unique(centerz, 'rows');

centaz_abbreviated=round(centaz,1);

[Au,~,ic] = unique(centaz_abbreviated, 'rows');
tally = accumarray(ic, 1);
Result = [Au tally];



[unique_gridpositions,~,label_gridpositions] = unique(gridpositions, 'rows');
counts_gridpositions = accumarray(label_gridpositions, 1);
[unique_attractors,~,label_attractors] = unique(centaz_abbreviated, 'rows');
counts_attractors = accumarray(label_attractors, 1);
countz=accumarray([label_attractors label_gridpositions],1);
JPD=countz/sum(sum(countz));


for k=1:size(JPD,1)
    PX_givenY(k,:)=JPD(k,:)./sum(JPD(k,:));
end

for l=1:size(JPD,2)
    PY_givenX(:,l)=JPD(:,l)./sum(JPD(:,l));
end


P_X=sum(JPD,1);
H_X=sum(-P_X.*log2(P_X));

P_Y=sum(JPD,2);
H_Y=sum(-P_Y.*log2(P_Y));


intermed_X_givenY=-log2(PX_givenY);
        intermed_X_givenY(isnan(intermed_X_givenY))=0;
        intermed_X_givenY(isinf(intermed_X_givenY))=0;
intermed_X_givenY=PX_givenY.*intermed_X_givenY;

H_X_givenY=sum(sum(intermed_X_givenY,2).*P_Y);


intermed_Y_givenX=-log2(PY_givenX);
        intermed_Y_givenX(isnan(intermed_Y_givenX))=0;
        intermed_Y_givenX(isinf(intermed_Y_givenX))=0;
intermed_Y_givenX=PY_givenX.*intermed_Y_givenX;

H_Y_givenX=sum(sum(intermed_Y_givenX,1).*P_X);


MI=H_X-H_X_givenY
%MI=H_Y-H_Y_givenX






save(filenaming)










%%%FIgures


%scatter(Result(:,1),Result(:,2),50,Result(:,3),'filled')
%colormap(jet)
%hold off

%Hiow Many identically 0?
%How many below a certain val? Sort; plot by sort, need unique tho
%hist(D,100)
%hist(eucd) seems a hump? NOT uniform..?


figure(103)
%% Create two axes
ax1 = axes;
scatter(positions(:,1),positions(:,2),50,sum(final_states>1),'filled')
view(2)
ax2 = axes;
scatter(Result(:,1),Result(:,2),50,Result(:,3),'filled')
%% Link them together
linkaxes([ax1,ax2])
%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax1,flipud(gray))
colormap(ax2,flipud('autumn'))
%% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);


figure(104)
%% Create two axes
ax1 = axes;
scatter(positions(:,1),positions(:,2),50,sum(final_states>1)>0,'filled')
view(2)
ax2 = axes;
scatter(Result(:,1),Result(:,2),50,Result(:,3),'filled')
%% Link them together
linkaxes([ax1,ax2])
%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax1,flipud(gray))
colormap(ax2,flipud('autumn'))
%% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
hold on
scatter(allgridpts(gridpositions,1),allgridpts(gridpositions,2))
hold off

%% Save figures etc.
saveas(figure(101),['Dists_hamm_vs_eucd__' num2str(N) '_' num2str(asmany) 'trials' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5))],'epsc')
saveas(figure(102),['Overlaps_hamm_vs_eucd__' num2str(N) '_' num2str(asmany) 'trials' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5))],'epsc')
figure(103)
title('Stable firing patterns: do they span the space?')
%Estimating the number of attractors
%xlabel('Horizontal position [mm]')
%ylabel('Horizontal position [mm]')
set(gca,'FontSize',18)
saveas(figure(103),['Space_of_attractors__' num2str(N) '_' num2str(asmany) 'trials' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5))],'epsc')





%dists_old=squareform(pdist(positions));
%J_new=J.*dists_old>0.3;
%figure,gplot(J_new,positions);

%figure,errorbar(howfar,curvy,sqrt(var_stable),'LineWidth',2)
%title('Stable pattern correlation vs. stimulation distance')
%xlabel('Distance from trial stimulation center to original stimulation center')
%ylabel('C = \langle r_original(t)\cdot r_trial(t) \rangle_neurons')
%set(gca,'FontSize',18)



for i=101:103
h(i-100)=figure(i)
end
%savefig(h,[filenaming '_figs'])









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scatter(positions(:,1),positions(:,2),degrees-min(degrees)+1,'filled')
degrees_with_attractors
eva=evalclusters(Result(:,1:2),'kmeans','silhouette','Klist',[1:40]);
eva.OptimalK
%overlay_clusters



[idx,C]=kmeans(Result(:,1:2),ceil(2^MI));
figure(105)
%% Create two axes
ax1 = axes;
scatter(positions(:,1),positions(:,2),50,sum(final_states>1)>0,'filled')
view(2)
ax2 = axes;
scatter(Result(:,1),Result(:,2),50,Result(:,3),'filled')
%% Link them together
linkaxes([ax1,ax2])
%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax1,flipud(gray))
colormap(ax2,flipud('autumn'))
%% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
hold on
scatter(allgridpts(gridpositions,1),allgridpts(gridpositions,2))
scatter(C(:,1),C(:,2),30,'filled','b')
hold off



mean_dists_traveled=sqrt(sum(distances_traveled.^2,2));
mean_dists_traveled=min(mean_dists_traveled,5-mean_dists_traveled);
%mean(mean_dists_traveled)==0.3047
