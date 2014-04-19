A = importdata('dat/square.1.output');

figure
for i = 1:3:length(A(:,1))
    trisurf([i,i+1,i+2],A(:,1),A(:,2),A(:,3));
    hold on
end