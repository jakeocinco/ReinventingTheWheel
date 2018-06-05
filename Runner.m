clc;close;clear;
Rebirth = [6 4 4 4 4 4 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
res = 0;
for i = 1:100
    len = randi(15, 1) + 2;
    for k = 1:len
        ls(k) = 1.5;
    end
    shapes(i) = Polygon(zeros(1,len),zeros(1,len));
    shapes(i).plotPolygonD
    pause(.01);
    res(i) = shapes(i).Residual;
    
end


newListofShapes = shapes;
k = 1;
max = 20;
lst = order(newListofShapes);
old = lst(99);

while newListofShapes(1).Residual > .001 && k < max
    
    lst = order(newListofShapes);
    l = sort(killIndex(50, 100));
    alive = kill(lst,l);
    
    
    storeR = newListofShapes(randi(100));
    
    x = 50 / max;
    
    new = regenerate(alive,Rebirth, 50 - floor(x * k));
   
    
    
    newListofShapes = order(new);
    store(k) = newListofShapes(1);
%     storeR(k) = newListofShapes();
    
    
    k = k + 1;
end

disp('Got here');
% newListofShapes(1).movePrint;
%n(95).movePrint;


%%%
% 1st   - 6 Offspring
% 2-5   - 4 Offspring
% 6-11  - 3 Offspring
% 12-24 - 2 Offsrping
% 24-50 - 1 Offspring
%   5   - random
%%%

function a = randomPolygon()
    len = randi(100, 1) + 2;
    a = Polygon(zeros(1,len),zeros(1,len));
end
function new = regenerate(alive,rebirth, var)
    new(1) = alive(1);
    count = 1;
    for i = 1:length(alive)
        for k = 1:rebirth(i)
            newAng = Polygon.varyAngles(var, alive(i).Angles);
            newLen = Polygon.varyLengths(var, alive(i).Lengths);
            fprintf('About to do shape %i, gen %i', i, k);
            new(count) = Polygon(newAng, newLen);
            count = count + 1;
        end
        
    end
    for k = i+1:100
        new(k) = randomPolygon();
    end
end
function killedInd = killIndex(n,range)
    killedInd(1) = 0;
    i = 1;
    while i <= n
        
        v1 = randi(range - 1) + 1;
        v2 = randi(range - 1) + 1;
        while ismember(v1,killedInd) == 1
            v1 = randi(range - 1) + 1;
        end
        while ismember(v2,killedInd) == 1
            v2 = randi(range - 1) + 1;
        end
        
        if v1 > v2
            killedInd(i) = v1;
        else
            killedInd(i) = v2;
        end
        i = i +1;
    end
end
function alive = kill(lst, ind)
    alive(1) = lst(1);
    count = 1;
    for i = 1:length(lst)
        if ismember(i,ind) == 0
            alive(count) = lst(i);
            count = count + 1;
        end
    end
    
end
function ordered = order(list)
    ordered(1) = list(1);
    cont(1) = 0;
    for i = 1:length(list)
        location = 1;
        min = 100;
        for r = 1:length(list)
        
            if list(r).Residual < min && ismember(r, cont) == 0
                min = list(r).Residual;
                location = r;
            end
        end
        
    ordered(i) = list(location);
    cont(i) = location;
    
    end
end

