%Counts the number of cells with various relaxation times
%Used to generate minor subplot of Figure 3 C

paramtable=readtable('Data/finalparameters.csv');

tau_thresh=0:.125:2.5;
cell_counts=zeros(size(tau_thresh));
for itr=1:length(tau_thresh)
    unit_list_indexes=find(table2array(paramtable(:,4))>=tau_thresh(itr));
    unit_list=table2array(paramtable(unit_list_indexes,1));
    cell_counts(itr)=length(unit_list);
end
figure(4)
plot(tau_thresh,cell_counts,'.k-','MarkerSize',26)
xlabel('Relaxation Time Threshold (s)')
ylabel('Number of Cells')
