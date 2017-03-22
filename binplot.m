% Function binplot, Created by Alexis, 08/12/2016

% plot y = f(x) by bins of x. x is separated in n quantile, for each
% quantile, the function plot the median of y for each quantile of x, and
% quartiles as shade. y is ploted at the location x = median of the
% quantile. 

% To do: add optional input, such as dashed line and color, etc... with
% default option

function [s,p] = binplot(x,y,n,c,sized) % s is shade (quartiles) and p is plot (median)
xbinedges = quantile(x,0:1/n:1); %* n bins of driver class, by quantile (n bins with same amount of data)
Quartilesy_xbins = nan(n,3); Median_xbins = nan(n,1);
for i = 1:n 
    xbins = find(x >= xbinedges(i) & x <= xbinedges(i+1)); %* iteration change x class (x bins)
    Quartilesy_xbins(i,[1 2 3]) = quantile(y(xbins),[0.25, 0.5, 0.75]);
    Median_xbins(i) = median(x(xbins));
end
X = [Median_xbins' fliplr(Median_xbins')];
y1 = Quartilesy_xbins(:,1)'; y2 = Quartilesy_xbins(:,3)';
Y = [y1 fliplr(y2)];
s = fill(X,Y,c,'FaceAlpha',0.2,'EdgeAlpha',0.0); hold on;  % !! x1,x2,y1 and y2 must be dim 1-n and 1-m (not n-1 and m-1)
p = plot(Median_xbins,Quartilesy_xbins(:,2),'LineWidth',3,'Marker','o','MarkerEdgeColor', ...
    'none','MarkerFaceColor',c,'Color',c,'MarkerSize',sized);
end