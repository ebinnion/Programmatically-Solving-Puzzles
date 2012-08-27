function Bez_Spline(filename,z,plotit)

data_hold = load(filename);
x = data_hold(:,1);
y = data_hold(:,2);
Nx = length(x);
% step_size = max(round(sqrt(Nx/100)),2);

holder = find_corners_special(x,y,Nx,z);
corners = holder(1,holder(1,:)~=0)+2;
bends = holder(2,(holder(2,:)~=0)&(holder(2,:)<Nx));
    

oops = find(corners>Nx);
if length(oops)
    corners(oops) = corners(oops)-Nx;
end

figure
errrr = [100,200,475,712,937,1155,1610,1697,1939,2014,2245,2363];
plot(x,y,'.')%,x(errrr),y(errrr),'or')

% For data sheet 12/02/10
% Nx
% width=max(x)-min(x)
% heigth=max(y)-min(y)
% no_corners=length(corners)

plot(x,y,'.k',x(corners),y(corners),'oc',x(bends),y(bends),'*r');
axis([min(x)-5,max(x)+5,min(y)-5,max(y)+5])
title(filename)

%Start Eric Code
linList = [];
bendVar = 1;
cornerVar = 1;
i = 1;
listIndex = 1;
%Create a single list of all bends and corners
while i < Nx
   if bendVar <= length(bends)
        if i ==  bends(bendVar) 
            linList(listIndex,1) = bends(bendVar);
            linList(listIndex,2) = 1;
            bendVar = bendVar + 1;
            listIndex = listIndex + 1;
        end
   end
   if cornerVar <= length(corners)
        if i == corners(cornerVar)
            linList(listIndex,1) = corners(cornerVar);
            linList(listIndex,2) = -1;
            cornerVar = cornerVar + 1;
            listIndex = listIndex + 1;
        end
   end
   i = i+1;
end
%Let's print the list now ;)
linList;

%Now let's split into cubics
i = 2; 
cubicsList = [];
cubicsList(1,1) = linList(1,1);
bendCount = 0;
colVar = 2;
cubicsIndex = 1;

while i < length(linList)
    %This gives us our middle bends
    if linList(i,2) > 0 & bendCount < 2
        bendCount = bendCount + 1;
        cubicsList(cubicsIndex, colVar) = linList(i,1);
        colVar = colVar + 1;
    %This case is for having 3 bends and finding mid point
    elseif linList(i,2) > 0 & bendCount >= 2
        cubicsList(cubicsIndex, 4) = (linList(i,1) + cubicsList(cubicsIndex,3)) / 2;
        %Reset colVar to keep track of middle bends
        colVar = 2;
        %Move to next row of cubicsList
        cubicsIndex = cubicsIndex + 1;
        %Ending point is starting point for next cubic
        cubicsList(cubicsIndex,1) = cubicsList(cubicsIndex-1 ,4);
        bendCount = 0;
        %Decrement i by 1 to account for going past the end point
        i = i - 1;
    elseif linList(i,2) < 0 
        cubicsList(cubicsIndex, 4) = linList(i,1);
        %Reset colVar to keep track of middle bends
        colVar = 2;
        %Move to next row of cubicsList
        cubicsIndex = cubicsIndex + 1;
        %Ending point is starting point for next cubic
        cubicsList(cubicsIndex,1) = cubicsList(cubicsIndex-1 ,4);
        bendCount = 0;
    end
    i = i + 1;
end
%The last point is also the first point
cubicsList(cubicsIndex,4) = cubicsList(1,1);
%Print the list of cubics
cubicsList = round(cubicsList)

%Initialize some variables to use for correcting interpolating curve


%This first while loop will do curve matching on the first half of the
%cubic
u = linspace(0,1,100);
cubicsIndex = 1;
correctionMatrix = [];
figure
while cubicsIndex <= length(cubicsList)
    correct1 = 0;
    correct2 = 0;
    flag = 1;
    tempDiff1 = 0;
    %This while loop is for the first part of the cubic
    while flag == 1;
        x0 = x(cubicsList(cubicsIndex,1));
        x1 = x(cubicsList(cubicsIndex,2)+correct1);
        x2 = x(cubicsList(cubicsIndex,3)+correct2);
        x3 = x(cubicsList(cubicsIndex,4));
        y0 = y(cubicsList(cubicsIndex,1));
        y1 = y(cubicsList(cubicsIndex,2)+correct1);
        y2 = y(cubicsList(cubicsIndex,3)+correct2);
        y3 = y(cubicsList(cubicsIndex,4));

        %This is a matrix for the interpolating curve
        xPoly = ((-4.5*u.^3+9*u.^2-5.5*u+1)*x0+(13.5*u.^3-22.5*u.^2+9*u)*x1+(-13.5*u.^3+18*u.^2-4.5*u)*x2+(4.5*u.^3-4.5*u.^2+u)*x3) - x0;
        yPoly = ((-4.5*u.^3+9*u.^2-5.5*u+1)*y0+(13.5*u.^3-22.5*u.^2+9*u)*y1+(-13.5*u.^3+18*u.^2-4.5*u)*y2+(4.5*u.^3-4.5*u.^2+u)*y3) - y0;

        %Original data bounded by cubic
        xTemp = x(cubicsList(cubicsIndex,1):cubicsList(cubicsIndex,4)) - x0;
        yTemp = y(cubicsList(cubicsIndex,1):cubicsList(cubicsIndex,4)) - y0;

        xx = x1 - x0;
        yy = y1 - y0;
        Theta = -atan(yy/xx);
        A = [cos(Theta), -sin(Theta); sin(Theta), cos(Theta)];

        %This is to rotate the original discrete data
        rot = A*[xTemp';yTemp'];
        xrot=rot(1,:);
        yrot=rot(2,:);

        %This is to rotate the Interpolating curve
        rotPoly = A*[xPoly;yPoly];
        xRotPoly = rotPoly(1,:);
        yRotPoly = rotPoly(2,:);

        rSum1 = 0;
        rPolySum1 = 0;

        %These next two loops will sum up and average the area under/above a
        %curve
        for i = 1:(cubicsList(cubicsIndex,2)-cubicsList(cubicsIndex,1))
           rSum1 = rSum1 + yrot(i);
        end
        rSum1 = rSum1 / (cubicsList(cubicsIndex,2) - cubicsList(cubicsIndex,1))

        i = 1;
        while xRotPoly(i) <= xrot(cubicsList(cubicsIndex,2)-cubicsList(cubicsIndex,1))
            rPolySum1 = rPolySum1 + yRotPoly(i);
            i = i +1;
        end
        rPolySum1 = rPolySum1 / (i-1)
        
        tempDiff1 = abs(rSum1 - rPolySum1)
        %We use absolute value > 'X' to set our threshold
        if rSum1 > rPolySum1 && tempDiff1 > 1 && rSum1 < 0
            correct1 = correct1 + 1
        elseif rSum1 < rPolySum1 && tempDiff1 > 1 && rSum1 < 0
            correct1 = correct1 - 1
        elseif rSum1 > rPolySum1 && tempDiff1 > 1 && rSum1 > 0
            correct1 = correct1 - 1
        elseif rSum1 < rPolySum1 && tempDiff1 > 1 && rSum1 > 0
            correct1 = correct1 + 1
        else
            flag = 0;
        end

        plot(xrot,yrot,'*', xRotPoly, yRotPoly,'.')
        pause
    end

    flag = 1;
    tempDiff2 = 0;
    %This while loop is for the second part of each cubic
    while flag == 1
        x0 = x(cubicsList(cubicsIndex,1));
        x1 = x(cubicsList(cubicsIndex,2)+correct1);
        x2 = x(cubicsList(cubicsIndex,3)+correct2);
        x3 = x(cubicsList(cubicsIndex,4));
        y0 = y(cubicsList(cubicsIndex,1));
        y1 = y(cubicsList(cubicsIndex,2)+correct1);
        y2 = y(cubicsList(cubicsIndex,3)+correct2);
        y3 = y(cubicsList(cubicsIndex,4));

        %This is a matrix for the interpolating curve
        xPoly = ((-4.5*u.^3+9*u.^2-5.5*u+1)*x0+(13.5*u.^3-22.5*u.^2+9*u)*x1+(-13.5*u.^3+18*u.^2-4.5*u)*x2+(4.5*u.^3-4.5*u.^2+u)*x3) - x2;
        yPoly = ((-4.5*u.^3+9*u.^2-5.5*u+1)*y0+(13.5*u.^3-22.5*u.^2+9*u)*y1+(-13.5*u.^3+18*u.^2-4.5*u)*y2+(4.5*u.^3-4.5*u.^2+u)*y3) - y2;

        %Original data bounded by cubic
        xTemp = x(cubicsList(cubicsIndex,1):cubicsList(cubicsIndex,4)) - x2;
        yTemp = y(cubicsList(cubicsIndex,1):cubicsList(cubicsIndex,4)) - y2;

        xx = x3 - x2;
        yy = y3 - y2;
        Theta = -atan2(yy,xx);
        A = [cos(Theta), -sin(Theta); sin(Theta), cos(Theta)];

        %This is to rotate the original discrete data
        rot = A*[xTemp';yTemp'];
        xrot=rot(1,:);
        yrot=rot(2,:);

        %This is to rotate the Interpolating curve
        rotPoly = A*[xPoly;yPoly];
        xRotPoly = rotPoly(1,:);
        yRotPoly = rotPoly(2,:);

        rSum2 = 0;
        rPolySum2 = 0;

        %These next two loops will sum up and average the area under/above a
        %curve
        for i = cubicsList(cubicsIndex,3)-cubicsList(cubicsIndex,1):(cubicsList(cubicsIndex,4)-cubicsList(cubicsIndex,1))
           rSum2 = rSum2 + yrot(i);
        end
        rSum2 = rSum2 / (cubicsList(cubicsIndex,4) - cubicsList(cubicsIndex,3))
        
        i = 100;
        while xRotPoly(i) >= 0
            rPolySum2 = rPolySum2 + yRotPoly(i);
            i = i - 1;
        end
        rPolySum2 = rPolySum2 / (i-1)
        
        %Let's break the loop if the previous change makes the difference
        %greater
        abs(rSum2 -rPolySum2)
        if tempDiff2 < abs(rSum2 -rPolySum2) && tempDiff2 ~= 0
            break;
        end
        
        tempDiff2 = abs(rSum2 -rPolySum2)
        %We use absolute value > 'X' to set our threshold
        if rSum2 > rPolySum2 && tempDiff2 > 1 && rSum2 < 0
            correct2 = correct2 + 1
        elseif rSum2 < rPolySum2 && tempDiff2 > 1 && rSum2 < 0
            correct2 = correct2 - 1
        elseif rSum2 > rPolySum2 && tempDiff2 > 1 && rSum2 > 0
            correct2 = correct2 - 1
        elseif rSum2 < rPolySum2 && tempDiff2 > 1 && rSum2 > 0
            correct2 = correct2 + 1
        else
            flag = 0;
        end

        plot(xrot,yrot,'*', xRotPoly, yRotPoly,'.')
        pause
    end
    correctionMatrix(cubicsIndex,1) = correct1
    correctionMatrix(cubicsIndex,2) = correct2
    cubicsIndex = cubicsIndex + 1
end



%End Eric code

% saveas(1,sprintf('piece%02d_figure.png',z))
% close all;

% for i = 1:length(corners)-1
%     xs = x(corners(i):corners(i+1));
%     ys = y(corners(i):corners(i+1));
%     holdit = Find_end_splines(xs,ys,step_size);
%     figure
%     plot(xs,ys,'.',xs(holdit),ys(holdit),'or')
%     axis([min(x)-5,max(x)+5,min(y)-5,max(y)+5])
% end
% xs = [x(corners(length(corners)):Nx);x(1:corners(1))];
% ys = [y(corners(length(corners)):Nx);y(1:corners(1))];
% holdit = Find_end_splines(xs,ys,step_size);
% figure
% plot(xs,ys,'.',xs(holdit),ys(holdit),'or')
% axis([min(x)-5,max(x)+5,min(y)-5,max(y)+5])
% 
% if (plotit==0)
%     close all
% end