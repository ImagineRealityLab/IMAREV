function plotCI(x,y,kind,linecolor,fillcolor,type)
% function plotCI(x,y,kind,linecolor,fillcolor,type)
%
% x is a p x n data matrix with n measurements of p features. y is a 1 x p
% vector such as time. 'kind' can be 'CI' or 'SEM', default is 'CI';
% 'linecolor' is the color of the average line and'fillcolor' is the color
% of the shaded area. If 'type' is 'over' you can plot multiple shaded
% areas over each other.
%
% plots the 95% confidence interval as shaded area.

if nargin < 3 || isempty(kind)
    kind = 'CI';
end

if nargin < 4 || isempty(linecolor)
    linecolor = 'black';
end

if nargin < 5 || isempty(fillcolor)
    fillcolor = [0.5 0.5 0.5];
end

if nargin < 6 || isempty(type)
    type = 'one';
end

% remove NaNs
nanIdx = isnan(x(:,1));
x      = x(~nanIdx,:);
y      = y(~nanIdx);

n = size(x,2);

SEM = std(x,[],2)/sqrt(n); 
ts = tinv([0.025  0.975],n-1); 

mX = mean(x,2);

if strcmp(kind,'CI')
    z1 = mX+(SEM*ts(2));
    z2 = mX-(SEM*abs(ts(1)));
elseif strcmp(kind,'SEM')
    z1 = mX+SEM;
    z2 = mX-SEM;
end

if strcmp(type,'over')
h = fill([y';flipud(y')],[z1;flipud(z2)],fillcolor,'LineStyle','none');
set(h,'facealpha',.5)
else
% plotting main structure
plot(y,mX,y,z1,y,z2);

% add SEM shaded area
lims = ylim;
a1 = area(y,z1,lims(1));
hold on;
set(a1,'LineStyle','none');     set(a1,'FaceColor',fillcolor);
alpha 0.5
a2 = area(y,z2,lims(1));
set(a2,'LineStyle','none');     set(a2,'FaceColor',[1 1 1]);

hold on; plot(y,mX,'Color',linecolor); hold off;
end