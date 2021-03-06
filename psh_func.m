function bends = psh_func(Kt,m)

tol = round(m/4);     % how close bends can be to another
min_wid = round(m/3); % how large the region has to be to qualify as a bend
mn = mean(Kt);
std_prop = .8;        % threshold proportion
psh = (Kt-mn)/max(abs(Kt-mn));
xes = 1:length(psh);
both = find((psh<(mean(psh)-std_prop*std(psh)))|(psh>(mean(psh)+std_prop*std(psh))));

% Refine the candidates to find the proper bendies :)
gumby = [];
lef = 1;
rt = 1;
isup = psh(both)./abs(psh(both));

for i = 2:length(both)
    holder = [];
    
    % if the bend candidate isn't too close to another candidate that has
    % already been chosen
    if ((both(i)-both(i-1)<tol) & (isup(i)==isup(i-1)))
        rt = i;
         
        if (i == length(both))
            % find the largest candidate bends in this region
            if ((rt - lef)>=min_wid)
                if (isup(lef) == 1)
                    holder = find(psh == max(psh(both(lef):both(rt))));
                else
                    holder = find(psh == min(psh(both(lef):both(rt))));
                end
            end
        end
        
    % otherwise, begin a new region here
    else
        if ((rt - lef)>=min_wid)
            if (isup(lef) == 1)
                holder = find(psh == max(psh(both(lef):both(rt))));
            else
                holder = find(psh == min(psh(both(lef):both(rt))));
            end
        end
        lef = i;
        rt = i;
    end
    gumby=[gumby,holder];
    
end

%upper = find((psh>(mean(psh)+std(psh))));
%lower = find((psh<(mean(psh)-std(psh))));
%figure;
%plot(xes,psh,'-k', xes, mean(psh),'-g',xes,(mean(psh)+std_prop*std(psh)),'-c',...
    %xes,(mean(psh)-std_prop*std(psh)),'-c',xes(both),psh(both),'ob',xes(gumby),psh(gumby),'r*');
%title(['Corner Graph ',num2str(psh(1))]);

%flag = (psh(both>1))-(psh(both)<0);
bends = gumby;

% should probably use if loop instead of the case function?

% switch
% case1 top1, bottom2 within 30 pixels
%   look for top3, bottom4, bottom5 to the right
%   call top1, bottom4, bottom5 bends
%   if not top1, bottom2 qualify as bends
% case2 bottom1, bottom2, top3
%   look for bottom4, top5 withing 70 pixels to the right of bottom2
%   call bottom1, bottom2, top5 bends
%   if not bottom1, bottom2, top3 qualify as bends
% otherwise,
%   call bend

% if no cases are satisfied, then look for top&top or bottom&bottom and
% call those a bend
