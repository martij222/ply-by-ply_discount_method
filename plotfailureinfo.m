function [] = plotfailureinfo(failure_loads,failure_plies)
%plotfailureinfo Creates a visualization of failure loads and plies
%   failure_loads is an array of failure loads. failure_plies is a cell array of the corresponding failed plies.
n = length(failure_loads);

figure; grid on; hold on;
for i=1:n
    X = failure_plies{i};
    X = categorical( sort(X) );

    Y = failure_loads(i);
    Y = repelem(Y, length(X));
    b=bar(X,Y);
    xtips = b(1).XEndPoints; ytips = b(1).YEndPoints;
    labels = string( round( (b(1).YData)/(1e6),2 ) );
    text(xtips,ytips,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
end
hold off;
title("Failure Loads and Plies");
xlabel("Ply Number"); ylabel("F");

end

