%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will attempt to find the location of corners on a discrete
% curve represented by xy points in a text file.
%
% args:
%
%   file:   data file handle
%   plotit: bool value (true = plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corner_points = find_corners_special(x,y,Nx,z)

%%%%%%%%%
% Setup %
%%%%%%%%%

% Compute m, d_lim and bd from the number of data points 
% (Jacob Staples' algorithm)
m = round(sqrt(Nx)/2);
d_lim = round(sqrt(m));
bd = m / d_lim;


x_av = mean(x);                     % average of x values
y_av = mean(y);                     % average of y values

tol = 2*round(Nx/100);              % Allowed closeness of corners

s = 201;                            % Steps for calculating frequency of Kt

%%%%%%%%%%%%%%%%%%%%
% Run computations %
%%%%%%%%%%%%%%%%%%%%

% Kt is a measure of curvature relative to neighboring points.  To find it,
% we perform the following steps:
%
% 1) we find the distance between points and the average value for the 
%    points (the center of the shape) and subtract off the average offset
%    distance = sigma((x-x_avg)^2 + (y-y_avg)^2) - mu(deviation)
Kt = ((x'-x_av).^2+(y'-y_av).^2);
Kt = Kt - mean(Kt);

% Pad the values of Kt to allow modular wrapping about the closed curve
Kt = [Kt(Nx+1-3*m:Nx),Kt,Kt(1:3*m)];
N = length(Kt);

% Define left-slope and right-slope
ml = zeros(1,N);
mr = zeros(1,N);

% 2) now we compute the "average derivative" of Kt, which is done in the
%    following manner:
%
% Compute the average of m slopes to left and right of center point (i),
% then set Kt equal to the absolute value of the difference in LR slopes
for i = (m+1):(N-m)
    holderl = 0;
    holderr = 0;
    for j = 1:m
        %holderx += m(i,j)
        holderl = holderl+((Kt(i)-Kt(i-j))/(j));
        holderr = holderr+((Kt(i+j)-Kt(i))/(j));
    end
    ml(i)= holderl/m;
    mr(i)= holderr/m;   
end

% set the curvature variable equal to the difference in average slope from
% the m points to the left and the m points to the right
Kt = ml-mr;

% Smooth the data (thrice) using averaging to eliminate horrendous spikes
for j = 1:3
    for i = 1:N-1
        Kt(i)=(Kt(i)+Kt(i+1))/2;
    end
    Kt(N) = (Kt(1)+Kt(N))/2;
end

% Determine the candidate corner points
cand = find(Kt>bd*mean(abs(Kt)));

% Refine the candidates to find the proper corners
corners = [];
lef = 1;
rt = 1;

for i = 2:length(cand)
    holder = 0;
    
    % if the corner candidate isnt too close to another candidate that has
    % already been chosen
    if (cand(i)-cand(i-1)<tol)
        rt = i;
         
        if (i == length(cand))
            % find the largest candidate corner in this region
            holder = find(Kt == max(Kt(cand(lef):cand(rt))))-3*m;
        end
        
    % otherwise, begin a new region here
    else
        holder = find(Kt == max(Kt(cand(lef):cand(rt))))-3*m;
        lef = i;
        rt = i;
    end
    
    % now that we have the curve split up into regions with each region
    % containing a corner somewhere, we add the best candidate corner such
    % that d_limit is exceeded by a normal line to the line between points
    % P1,P3 that passes through P2
    for j = 1:length(holder)
        if ((holder(j)>0)&&(holder(j)<(Nx+1)))
            
            %possible corner found...test d_lim to be sure
            if(isempty(find(holder(j)==corners,1)))
                %use a smaller radius for the secondary check
                r = round(m/3);
                
                x1 = x(rem(holder(j)-r+Nx-1,Nx)+1);
                y1 = y(rem(holder(j)-r+Nx-1,Nx)+1);
                xa = x(holder(j));
                ya = y(holder(j));
                x2 = x(rem(holder(j)+r+Nx-1,Nx)+1);
                y2 = y(rem(holder(j)+r+Nx-1,Nx)+1);
                
                dd = new(x1,y1,xa,ya,x2,y2);
                
                %use secondary check to verify we have a corner
                if (dd > 1.25)
                    corners=[corners,holder(j)];
                    %new(x2,y2,xa,ya,x1,y1,1)
                    %find_dr(x2,y2,xa,ya,x1,y1,1);
                end
            end
        end
    end
end

%corners = [18 288 490 731]

% % Remove the wrap-around fringe - only view one cycle of closed loop
Kt = Kt(3*m+1:(N-(3*m)));
% 
% % Normalize the data
Kt = Kt/max(abs(Kt));
t=1:length(Kt);
% errrr = [100,200,475,712,937,1155,1610,1697,1939,2014,2245,2363]
%subplot(2,1,1)
%plot(t,Kt,'-')%,t(errrr),Kt(errrr),'or');
%title('Curvature Kt')
 
%Test section to find frequency of Kt
Kf = linspace(-1,1,s);
Kf2 = linspace(0,0,s);
for i = 1:Nx
    for j = 1:s-1
        if Kf(1,j) <= Kt(1,i) && Kt(1,i) < Kf(1,j+1)
            Kf2(1,j)= Kf2(1,j)+1;
        end
    end
end

%subplot(2,1,2)
%h2 = plot(Kf,Kf2);
%title('Frequency Kf2 of Kt')
% saveas(h2,sprintf('piece%02d_plot.png',z))
% close all

mn = mean(Kt);
st = std(Kt);

%10-16-10
 bends = [];
if length(corners)>0
    if (corners(1)-m) > 0
    psh = Kt(1:corners(1)-m);
        if (length(find(psh>mn+st|psh<mn-st))>1)
        holder = psh_func(psh,m);
        bends=[bends,holder];
        end
    end
    for i = 1:(length(corners) - 1)
    psh = Kt(corners(i)+m:corners(i+1)-m);
        if (length(find(psh>mn+st|psh<mn-st))>1)
        holder = psh_func(psh,m);
        bends=[bends,holder+corners(i)+m];
        end
    end
    if (corners(length(corners))+m) < Nx
    psh = Kt(corners(length(corners))+m:Nx);
        if (length(find(psh>mn+st|psh<mn-st))>1)    
        holder = psh_func(psh,m);
        bends=[bends,holder+corners(length(corners))+m];
        end
    end
end
 
%%%%%%%%%%
% Output %
%%%%%%%%%%

%d = zeros(N,1);
%for j = 1:(length(corners)-1)
    %x1 = x(corners(j));
    %x2 = x(corners(j+1));
    %y1 = y(corners(j));
    %y2 = y(corners(j+1));
    %m = (y2-y1)/(x2-x1);
    %ctr = 1;
    %for k = corners(j):corners(j+1)
        %d(ctr) = abs(-m*x(k)+y(k)+(m*x1-y1))/sqrt(m^2+1);
        %ctr = ctr+1;
    %end
%end
%x1 = x(corners(4));
%x2 = x(corners(1));
%y1 = y(corners(4));
%y2 = y(corners(1));
%m = (y2-y1)/(x2-x1);
%ctr = 1;
%for k = corners(4):Nx
    %d(ctr) = abs(-m*x(k)+y(k)+(m*x1-y1))/sqrt(m^2+1);
    %ctr = ctr+1;
%end
%for k = 1:corners(1)
    %d(ctr) = abs(-m*x(k)+y(k)+(m*x1-y1))/sqrt(m^2+1);
    %ctr = ctr+1;
%end

holder = zeros(2,max(length(corners),length(bends)));

holder(1,1:length(corners)) = corners;
holder(2,1:length(bends)) = bends;
corner_points = holder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d_ratio = new(x1,y1,xa,ya,x2,y2)

zeta = atan((y1 + y2 - 2*ya)/(x1 + x2 - 2*xa));

m = (sqrt((x2 - xa).^2 + (y2 - ya).^2) + sqrt((xa - x1).^2 + (ya - y1).^2))/2;

x_disp = m*cos(zeta);
y_disp = m*sin(zeta);

p_x1 = xa - x_disp;
p_y1 = ya - y_disp;
p_x2 = xa + x_disp;
p_y2 = ya + y_disp;

d1 = sqrt((p_x1 - x1).^2 + (p_y1 - y1).^2);
d2 = sqrt((p_x1 - x2).^2 + (p_y1 - y2).^2);
d3 = sqrt((p_x2 - x1).^2 + (p_y2 - y1).^2);
d4 = sqrt((p_x2 - x2).^2 + (p_y2 - y2).^2);

d_ratio = ((d1 + d2)/(d3 + d4));

if(d_ratio < 1)
    d_ratio = 1/d_ratio;
end

%d_ratio


