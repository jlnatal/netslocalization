
%load Data201810218231__N4096__cutoff0pt3___250runs___noiselevel0thousandths
    

    %deviations_raw=[centaz(:,1)-allgridpts(gridpositions,1) centaz(:,2)-allgridpts(gridpositions,2)];
    %distances_traveled=[min(abs(deviations_raw(:,1)),abs(5-abs(deviations_raw(:,1)))) min(abs(deviations_raw(:,2)),abs(5-abs(deviations_raw(:,2))))];
    

Distances=distances;
Distances(Distances==0)=NaN;

%figure,hist(min(Distances),100)
%figure,hist(min(Distances))
figure,hist(min(Distances),13)
[mean(min(Distances)) std(min(Distances))]
