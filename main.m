
%%Set up system
    
    tic
    fprintf('Setting up system... \n')
    
    N=2^12;         %Number of neurons
    cutoff=0.3;     %Cutoff for synapses
    restr=0.3;      %Size of input patch
    a=0.02 ;        %Sparsity (fraction of active neurons in each memory)



%%Animation and stochasticity switches

    animateon=0;
    
    if animateon
        fprintf('(Animations on) \n')
        
    else
        fprintf('(Animations off) \n')
        
    end

    noiseon=0;
    if noiseon
        %noiselevel=0.001;
        noiselevel=0.1;
        fprintf(['(Noise level ' num2str(noiselevel) ') \n'])
    else
        noiselevel=0;
        fprintf('(Noiseless) \n')
    end


%%Grid-style simulations
    
    %Construct grid
    [X,Y] = meshgrid(0:.1:4.9,0:.1:4.9);
    allgridpts=[X(:) Y(:)];
    
    %How many points to sample from grid?
    number_trials=5;
    %number_trials=input('How many trials? (Default 250)  ');
    %if isempty(number_trials)
    %    number_trials=250;
    %end
    
    %Sample number_trials grid points for number_trials experimental simulations
    gridpositions=randi(length(allgridpts),number_trials,1);
    % At number_trials=2500, hits most of the 2500 gridpoints
    % At number_trials=15000, hits consistently in the 2490's
    %For spacing 0.1 at least...

    fprintf(['(' num2str(number_trials) ' trials) \n'])
    
    
%%Set up simulations

%Duration of experiments
    secondz=200;
    tspan=[0:0.1:secondz];
    devtimes=[1:25:751];         %Times at which to measure deviations of COE's (bump size?)
    fprintf(['(Trial duration ' num2str(secondz) ' seconds) \n'])

%Large array initialization
    Trajects=zeros(length(tspan),N,number_trials);

COE_temporal=zeros(number_trials,2,length(devtimes));

%Generate E/I topology with cutoff distance 0.3 on 5x5 plane
    [J,positions,distances,N,dale,cutoff]=configure_topology(N,cutoff);
    fprintf('     Connectivity generated \n')



%%Prepare system for experimental trials


%Initialize system with all rates = a, then let it run without input for a time

    fprintf('Initializing firing rates... \n')
    r0=a*ones(1,N);
    if noiseon
        noisestream=1+cumsum(noiselevel*randn(1000,1));
        noisestream_smooth=fit((1:1000)',noisestream,'linearinterp');
        [t,r] = ode45(@(t,r) myode_timedep_NOinput(t,r,J,N,noisestream_smooth), tspan,r0);
    else
        [t,r] = ode45(@(t,r) myode_NOinput(t,r,J,N), [0:0.1:1000],r0);
    end


%Run system once from the end of above dynamics for uniformly initialized conditions
% i.e. No external input

    fprintf('Running to equilibration (no external input)... \n')
    r0=r(end,:);
    reach=0;

    ang=2*pi*rand;
    x_coord=reach*cos(ang);
    y_coord=reach*sin(ang);

    cents=rangesearch(positions,[2.5+x_coord 2.5+y_coord],.3);
    culprits=cents{1,1};
    input_IDs=culprits;
    

    if noiseon
        noisestream=1+cumsum(noiselevel*randn(tspan(end),1));
        noisestream_smooth=fit((1:tspan(end))',noisestream,'linearinterp');
        [t,r] = ode45(@(t,r) myode_timedep(t,r,J,N,input_IDs,noisestream_smooth), tspan,r0);
    else
        [t,r] = ode45(@(t,r) myode(t,r,J,N,input_IDs), tspan,r0);
    end
    %figure,plot(t,r)
    
    %SHOULD HAVE A CHECK IF IT FAILS ON THIS FIRST GO

    
    for la=1:10:length(tspan)

            if animateon==1
                animate
            end

            if la==length(tspan)
                cutthresh=0.1;
                %viscircles([positions(target,1) positions(dsearchn(positions,[2.5 2.5]),2)],restr,'Color','k');
                %find_COE
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

        top=(find(r(end,:)==max(r(end,:))));

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
    
    
    
    
howfar=0;


for lo=1:length(howfar)
    
    howfar(lo)
    
    reach=howfar(lo);
    
    toc
  
%BLOCK    
%figure(103),scatter(positions(:,1),positions(:,2),50,zeros(1,N),'filled')



final_states=zeros(number_trials,N);



stim_position=zeros(number_trials,2);
%Vector_dists=


%centerz=zeros(number_trials,2);




%% Run experimental simulations

    for exper=1:number_trials
        [fprintf('Experiment no. ') num2str(exper)]
        

    %Run the system directly from the previous attractor it settled into
    r0=r(end,:);

    cents=rangesearch(positions,[allgridpts(gridpositions(exper),1) allgridpts(gridpositions(exper),2)],.3);
    culprits=cents{1,1};
    input_IDs=culprits;


    if noiseon
        noisestream=1+cumsum(noiselevel*randn(tspan(end),1));
        noisestream_smooth=fit((1:tspan(end))',noisestream,'linearinterp');
        [t,r] = ode45(@(t,r) myode_timedep(t,r,J,N,input_IDs,noisestream_smooth), tspan,r0);
    else
        [t,r] = ode45(@(t,r) myode(t,r,J,N,input_IDs), tspan,r0);
    end

    
        for la=1:10:length(tspan)

            if animateon==1
                animate
            end

            if la==length(tspan)
                cutthresh=0.1;
                %viscircles([positions(target,1) positions(dsearchn(positions,[2.5 2.5]),2)],restr,'Color','k');
%                center_of_excitation=[sum((r(end,:)>cutthresh).*positions(:,1)')/sum((r(end,:)>cutthresh)) sum((r(end,:)>cutthresh).*positions(:,2)')/sum((r(end,:)>cutthresh))];
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
    
        
        top=(find(r(end,:)==max(r(end,:))));

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
        %Trajects(:,:,exper)=r;
        
        %figure(104)
        %errorbar(0.1:0.1:7,mean(Firing_at_rad),var(Firing_at_rad))
      
        %Old???
        %disty=[[sum((r(end,:)>a).*positions(:,1)')/sum((r(end,:)>a)) sum((r(end,:)>a).*positions(:,2)')/sum((r(end,:)>a))]] - [2.5+x_coord 2.5+y_coord];
        

        
    corr_at_range(lo,exper)=sum(mean(original.*r))/N;
    corr_at_range_END(lo,exper)=mean(original(end,:).*r(end,:));
    %corr_at_range_END(lo,exper)=mean(original(end,:)>a.*r(end,:)>a);
    %[original(end,1:5); r(end,1:5)]
    
    final_states(exper,:)=r(end,:);
    centa=find_COE(a,positions,final_states(exper,:));
    centaz(exper,:)=centa;
    
    
    %devtimes=[0:50:750];  ALREADY DEFINED ABOVE
    
    for devmeas=1:length(devtimes)
        current_state=r(devtimes(devmeas),:);
        COE=find_COE(a,positions,current_state);
        COE_temporal(exper,:,devmeas)=COE;    
    end 
    
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
    
    deviations_raw=[centaz(:,1)-allgridpts(gridpositions,1) centaz(:,2)-allgridpts(gridpositions,2)];
    %distances_traveled=[min(abs(deviations_raw(:,1)),abs(5-abs(deviations_raw(:,1)))) min(abs(deviations_raw(:,2)),abs(5-abs(deviations_raw(:,2))))];
    distances_traveled=convert_distances(deviations_raw);
        
endtime_in_mins=toc/60
%save(filenaming)

cee=clock;
if cutoff<1
    cutstring=['0pt' num2str(cutoff*10)];
elseif (cutoff>1 & cutoff<2)
    cutstring=['1pt' num2str(10*rem(cutoff,1))];
end
filenaming=['Data' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5)) num2str(floor(cee(6))) '__N' num2str(N)  '__cutoff' num2str(cutstring) '___' num2str(number_trials) 'runs' '___' 'noiselevel' num2str(noiselevel*1000) 'thousandths']
%save(filenaming)







%% Processessing

%b=combnk(1:4,2)
e=rand(number_trials,N)<(mean(sum(final_states>1,2))/N);

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
saveas(figure(101),['Dists_hamm_vs_eucd__' num2str(N) '_' num2str(number_trials) 'trials' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5))],'epsc')
saveas(figure(102),['Overlaps_hamm_vs_eucd__' num2str(N) '_' num2str(number_trials) 'trials' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5))],'epsc')
figure(103)
title('Stable firing patterns: do they span the space?')
%Estimating the number of attractors
%xlabel('Horizontal position [mm]')
%ylabel('Horizontal position [mm]')
set(gca,'FontSize',18)
saveas(figure(103),['Space_of_attractors__' num2str(N) '_' num2str(number_trials) 'trials' num2str(cee(1)) num2str(cee(2)) num2str(cee(3)) num2str(cee(4)) num2str(cee(5))],'epsc')





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





mean_dists_traveled=sqrt(sum(distances_traveled.^2,2));
mean_dists_traveled=min(mean_dists_traveled,5-mean_dists_traveled);
mindists



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


%Moved above
%mean_dists_traveled=sqrt(sum(distances_traveled.^2,2));
%mean_dists_traveled=min(mean_dists_traveled,5-mean_dists_traveled);
%%mean(mean_dists_traveled)==0.3047
