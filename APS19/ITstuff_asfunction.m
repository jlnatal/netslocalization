
centaz_abbreviated=round(centaz,2);

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
