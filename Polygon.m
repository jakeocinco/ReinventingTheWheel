classdef Polygon
    
    properties
        Sides % number of Sides
        Angles % Value for angles
        XPos % X position of variables
        YPos % Y position of variables
        Scale % Max length of sides
        Type % 1-3 -- 1 all the same; 2 slighly different; 3 mahor difference.
        Lengths % Array of sides
        CoM % Center of Mass, CoM(1) - x, CoM(2) - y
        Residual
        
    end
    methods
        %Plots polygon
        function plotPolygon(obj)
            
            plot(obj.XPos,obj.YPos);
            axis([-1.375*3 1.375*3 -3 3]);
            
            
            %axis([-1.375*obj.Scale 1.375*obj.Scale -obj.Scale obj.Scale]);
        end
        %Plots polygon with CoM
        function plotPolygonD(obj)
            
            plot(obj.XPos,obj.YPos, obj.CoM(1),obj.CoM(2),'d');
            axis([-3 3 -3 3]);
            
            %axis([-1.375*obj.Scale 1.375*obj.Scale -obj.Scale obj.Scale]); pause(0.03);
        end
        
        function movePrint(obj)
            
            clc; disp('Generating 100 species...')
            
            sides = length(obj.Angles);
            
            ChartTitle = 'Title';
            
            CoMX = obj.CoM(1);
            CoMY = obj.CoM(2);
            PosX = obj.XPos;
            PosY = obj.YPos;
            
            [Ymin,L] = min(PosY(1,1:sides));
            Xmin = PosX(1,L); % "Xmin" is the X-value of the coordinate which has Ymin.
            plot(PosX,PosY,'k',CoMX,CoMY,'rx');
            axis([-1.588*4 1.588*4 -1.25*4 1.25*4]); title(ChartTitle); pause(0.03);
            
            clc; disp('Species generation complete.')
            
            % Extra preallocated variables (to be used for rotation).
            maxsides = sides;
            
            Xdisp = zeros(1,maxsides);
            Ydisp = zeros(1,maxsides);
            radius = zeros(1,maxsides);
            theta = zeros(1,maxsides);
            CoMXdisp = 0;
            CoMYdisp = 0;
            CoMradius = zeros(1,1);
            CoMtheta = zeros(1,1);
            
            
            CoMXdisp = CoMX - Xmin;
            CoMYdisp = CoMY - Ymin;
            CoMradius = sqrt(CoMXdisp^2 + CoMYdisp^2);
            if CoMX > Xmin && CoMY > Ymin % Quadrant I
                CoMtheta = atand(CoMYdisp/CoMXdisp);
            elseif CoMX < Xmin && CoMY > Ymin % Quadrant II
                CoMtheta = 180 + atand(CoMYdisp/CoMXdisp);
            elseif CoMX < Xmin && CoMY < Ymin % Quadrant III
                CoMtheta = 180 + atand(CoMYdisp/CoMXdisp);
            elseif CoMX > Xmin && CoMY < Ymin % Quadrant IV
                CoMtheta = 360 + atand(CoMYdisp/CoMXdisp);
            elseif CoMX == Xmin && CoMY > Ymin
                CoMtheta = 90;
            elseif CoMX == Xmin && CoMY < Ymin
                CoMtheta = 270;
            elseif CoMX > Xmin && CoMY == Ymin
                CoMtheta = 0;
            elseif CoMX < Xmin && CoMY == Ymin
                CoMtheta = 180;
            elseif CoMX == Xmin && CoMY == Ymin
                CoMtheta = 0;
            end
            
            for n = 1:sides
                
                Xdisp(n) = PosX(n) - Xmin;
                Ydisp(n) = PosY(n) - Ymin;
                radius(n) = sqrt(Xdisp(n)^2 + Ydisp(n)^2);
                if PosX(n) > Xmin && PosY(n) > Ymin % Quadrant I
                    theta(n) = atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) < Xmin && PosY(n) > Ymin % Quadrant II
                    theta(n) = 180 + atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) < Xmin && PosY(n) < Ymin % Quadrant III
                    theta(n) = 180 + atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) > Xmin && PosY(n) < Ymin % Quadrant IV
                    theta(n) = 360 + atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) == Xmin && PosY(n) > Ymin
                    theta(n) = 90;
                elseif PosX(n) == Xmin && PosY(n) < Ymin
                    theta(n) = 270;
                elseif PosX(n) > Xmin && PosY(n) == Ymin
                    theta(n) = 0;
                elseif PosX(n) < Xmin && PosY(n) == Ymin
                    theta(n) = 180;
                elseif PosX(n) == Xmin && PosY(n) == Ymin
                    theta(n) = 0;
                end
            end
            
            
            % Display Options
            Terminate = 1;
            
            k = 1;
            if k == 0
                
            elseif k >= 1 || k <= 100
                RotModelX = zeros(360,maxsides);%May need to be maxsides + 1
                RotModelY = zeros(360,maxsides);%May need to be maxsides + 1
                RotModelCoMX = zeros(360,1);
                RotModelCoMY = zeros(360,1);
                XdispNew = Xdisp;
                YdispNew = Ydisp;
                radiusNew = radius;
                thetaNew = theta;
                CoMXdispNew = CoMXdisp;
                CoMYdispNew = CoMYdisp;
                CoMradiusNew = CoMradius;
                CoMthetaNew = CoMtheta;
                XminNew = Xmin;
                YminNew = Ymin;
                
                for k2 = 1:1:360
                    if k2 ~= 1
                        if min(RotModelY(k2-1,:)) < YminNew(k)
                            [~,L] = min(RotModelY(k2-1,1:sides));
                            RotModelX(k2,L) = radiusNew(k,L)*cosd(thetaNew(k,L) - k2);
                            RotModelY(k2,L) = radiusNew(k,L)*sind(thetaNew(k,L) - k2);
                            thetaoffset = 1-atand(abs((RotModelY(k2,L))/(RotModelX(k2,L))));
                            YminNew(k) = min(RotModelY(k2-1,1:sides));
                            XminNew(k) = RotModelX(k2-1,L);
                            
                            CoMXdispNew(k) = RotModelCoMX(k2-1) - XminNew(k);
                            CoMYdispNew(k) = RotModelCoMY(k2-1) - YminNew(k);
                            CoMradiusNew(k) = sqrt(CoMXdispNew(k)^2 + CoMYdispNew(k)^2);
                            if RotModelCoMX(k2-1) > XminNew(k) && RotModelCoMY(k2-1) > YminNew(k) % Quadrant I
                                CoMthetaNew(k) = atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) < XminNew(k) && RotModelCoMY(k2-1) > YminNew(k) % Quadrant II
                                CoMthetaNew(k) = 180 + atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) < XminNew(k) && RotModelCoMY(k2-1) < YminNew(k) % Quadrant III
                                CoMthetaNew(k) = 180 + atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) > XminNew(k) && RotModelCoMY(k2-1) < YminNew(k) % Quadrant IV
                                CoMthetaNew(k) = 360 + atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) == XminNew(k) && RotModelCoMY(k2-1) > YminNew(k)
                                CoMthetaNew(k) = 90 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) == XminNew(k) && RotModelCoMY(k2-1) < YminNew(k)
                                CoMthetaNew(k) = 270 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) > XminNew(k) && RotModelCoMY(k2-1) == YminNew(k)
                                CoMthetaNew(k) = 0 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) < XminNew(k) && RotModelCoMY(k2-1) == YminNew(k)
                                CoMthetaNew(k) = 180 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) == XminNew(k) && RotModelCoMY(k2-1) == YminNew(k)
                                CoMthetaNew(k) = 0 + k2 + sum(thetaoffset);
                            end
                            for n = 1:sides
                                XdispNew(k,n) = RotModelX(k2-1,n) - XminNew(k);
                                YdispNew(k,n) = RotModelY(k2-1,n) - YminNew(k);
                                radiusNew(k,n) = sqrt(XdispNew(k,n)^2 + YdispNew(k,n)^2);
                                if RotModelX(k2-1,n) > XminNew(k) && RotModelY(k2-1,n) > YminNew(k) % Quadrant I
                                    thetaNew(k,n) = atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) < XminNew(k) && RotModelY(k2-1,n) > YminNew(k) % Quadrant II
                                    thetaNew(k,n) = 180 + atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) < XminNew(k) && RotModelY(k2-1,n) < YminNew(k) % Quadrant III
                                    thetaNew(k,n) = 180 + atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) > XminNew(k) && RotModelY(k2-1,n) < YminNew(k) % Quadrant IV
                                    thetaNew(k,n) = 360 + atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) == XminNew(k) && RotModelY(k2-1,n) > YminNew(k)
                                    thetaNew(k,n) = 90 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) == XminNew(k) && RotModelY(k2-1,n) < YminNew(k)
                                    thetaNew(k,n) = 270 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) > XminNew(k) && RotModelY(k2-1,n) == YminNew(k)
                                    thetaNew(k,n) = 0 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) < XminNew(k) && RotModelY(k2-1,n) == YminNew(k)
                                    thetaNew(k,n) = 180 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) == XminNew(k) && RotModelY(k2-1,n) == YminNew(k)
                                    thetaNew(k,n) = 0 + k2 + sum(thetaoffset);
                                end
                            end
                        end
                    end
                    RotModelCoMX(k2,1) = CoMradiusNew(k)*cosd(CoMthetaNew(k) - k2) + XminNew(k);
                    RotModelCoMY(k2,1) = CoMradiusNew(k)*sind(CoMthetaNew(k) - k2) + YminNew(k);
                    for n = 1:sides + 1
                        if n == sides + 1
                            RotModelX(k2,n) = RotModelX(k2,1);
                            RotModelY(k2,n) = RotModelY(k2,1);
                        else
                            RotModelX(k2,n) = radiusNew(k,n)*cosd(thetaNew(k,n) - k2) + XminNew(k);
                            RotModelY(k2,n) = radiusNew(k,n)*sind(thetaNew(k,n) - k2) + YminNew(k);
                        end
                    end
                    %just sides??
                    tempy = RotModelY(k2,1:sides+1) + ( k2/360);
                    
                        [tempComx(k2),tempComy(k2)] = Polygon.centerOfMass(RotModelX(k2,1:sides),tempy);
                    
                    
                    
                    plot(RotModelX(k2,1:sides+1), tempy,'k',tempComx(k2),tempComy(k2),'rx');
%                     plot(RotModelX(k2,1:sides+1),RotModelY(k2,1:sides+1),'k',RotModelCoMX(k2),RotModelCoMY(k2),'rx');
                    axis([-2.5*4 + max(RotModelX(k2,:)) 0.676*4 + max(RotModelX(k2,:)) -0*4 + Ymin(k) 2.5*4 + Ymin(k)]);
                    ChartTitle = 'Title';
                    title(ChartTitle); xlabel('Distance traveled'); pause(0.001);
                end
                
                t = CoMX(k):(max(RotModelCoMX)+CoMX(k))/1000:max(RotModelCoMX);
                IdealPath = 1:1:length(t);
                IdealPath = IdealPath / length(t);
%                 IdealPath = zeros(1,length(t));
                disp(IdealPath);
                for t2 = 1:length(t)
%                     IdealPath(t2) = CoMY(k);
                end
                
%                 plot(RotModelCoMX,RotModelCoMY,'r',t,IdealPath,'b--')
                plot(PosX(k,1:sides+1),PosY(k,1:sides+1),'k',RotModelCoMX,RotModelCoMY,'r',t,IdealPath,'b--');
                axis([-1.588*4 + min(PosX(k,:)) 1.588*4 + max(RotModelX(k2,:)) -0.75*4 + Ymin(k) 1.75*4 + Ymin(k)]);
                title(ChartTitle); legend('Shape','Movement of CoM','Ideal Path'); pause(0.6);
                %clc; clearvars RotModel* t t2 thetaoffset; clearvars -REGEXP New$; pause(0.01); disp('Viewing complete! View another species or type "0" if done.')
            else
                error('Invalid species selected. Must be an integer from 1 to 100.');
            end
            
            
            
            
            
            
            
            
            
            
            
            
        end
        
        function a = childPolygon(obj)
            if rand() < .95
                if randi() > .5
                    a.Sides = obj.Sides + randi(2);
                else
                    a.Sides = obj.Sides - randi(2);
                end
            else
                a.Sides = obj.Sides;
            end
            
            
        end
    end
    methods (Static)
        %Called to set up a 'polygon'
        %Requires s (number of sides)
        %   and maxScale(Maximum size of length)
        function obj = Polygon(ang, leng)
            if length(ang) >= 3
                s = length(ang);
                
                
                if sum(leng) > .5
                    obj.Lengths = leng;
                else
                    obj.Lengths = Polygon.randomLengths(s, randi(3));
                end
                
                obj.Sides = s;
                if sum(ang) < 10
                    obj.Angles = Polygon.randomAngles(s);
                else
                    obj.Angles = ang;
                end
                
                [obj.XPos, obj.YPos] = Polygon.randomPolygon(s,obj.Angles,obj.Lengths);
                [obj.CoM(1), obj.CoM(2)] = Polygon.centerOfMass(obj.XPos, obj.YPos);
                obj.Scale = sum(obj.Lengths) / length(obj.Lengths);
                
                obj.Residual = driveNoShow(obj);  %(Polygon.move(obj.XPos, obj.YPos, obj.CoM(1), obj.CoM(2), s, obj.Scale);
                
                 %disp(r);
                
            end
            
        end
        
        %Vary Angles, enter old angles and variance between 1 and 100
        %   1 being close to no change
        function a = varyAngles(var, angles)
            N = length(angles);
            for i = 1:length(angles)
                if mod(randi(2) ,2) == 0
                    a(i) = angles(i) + randi(var);
                else
                    r = randi(var);
                    if angles(i) - r > 10
                        a(i) = angles(i) - r;
                    else
                        a(i) = angles(i);
                    end
                end
            end
            
            total = sum(a);
            while total > 360
                ran = randi(N,1);
                if a(ran) - 1 > 0
                    a(ran) = a(ran) - 1;
                end
                total = sum(a);
            end
            while total < 360
                ran = randi(N,1);
                if a(ran) + 1 < 360
                    a(ran) = a(ran) + 1;
                end
                total = sum(a);
            end
            
            if length(a) ~= length(angles) 
                a = Polygon.varyAngles(var, angles);
            end
            
        end
        
        %Vary lengths, enter old lengths and variance between 1 and 100
        %   1 being close to no change
        function l = varyLengths(var, lengths)
            N = length(lengths);
            for i = 1:length(lengths)
                if mod(randi(2) ,2) == 0
                    r = lengths(i) + (var/ 100) * rand();
                    if r < 3 
                        l(i) = lengths(i) + (var/ 100) * rand();
                    else
                        l(i) = r;
                    end
                else
                    r = lengths(i) - (var/ 100) * rand();
                    if r > .5
                        l(i) = r;
                    else
                        l(i) = rand() / 2 + .3;
                    end
                    
                end
            end
            
            if length(l) ~= length(lengths) 
                l = Polygon.varyAngles(var, lengths);
            end
            
        end
        
        %%%%%%%%
        %%%%%%%%    NONE OF THESE SHOULD BE
        %%%%%%%%    CALLED FROM THE OUTSIDE
        %%%%%%%%          OF THE CLASS.
        %%%%%%%%      ONLY USED IN SET-UP
        %%%%%%%%
        %Returns N angles, adding up to 360
        function angles = randomAngles(N)
            %N = input('Number of angles: ');
            angles = 1:N;
            total = 0;
            for i = 1:N
                angles(i) = randi([50,150],1);
                total = total + angles(i);
            end
            x = 360 /total;
            threesixty = 0;
            for i = 1:N
                angles(i) = angles(i) * x;
                if angles(i) < 1
                    angles(i) = 1;
                end
            end
            
            for i = 1:N
                angles(i) = floor(angles(i));
                threesixty = threesixty + angles(i);
            end
            while threesixty > 360
                rand = randi(N,1);
                if angles(rand) + 1 > 0
                    angles(rand) = angles(rand) - 1;
                end
                threesixty = sum(angles);
            end
            while threesixty < 360
                rand = randi(N,1);
                if angles(rand) + 1 < 360
                    angles(rand) = angles(rand) + 1;
                end
                threesixty = sum(angles);
            end
        end
        
        function lengths = randomLengths(N, type)
            if type == 1 %%All the same 
                len = (rand() * 1.5) + 1;
                for i = 1:N
                    lengths(i) = len;
                end
            elseif type == 2 %slighlty random
                v = (rand() * 1.5) + 1;
                for i = 1:N
                    plusMinus = randi(2);
                    if plusMinus == 1
                        lengths(i) = v + (rand() /2);
                    else
                        lengths(i) = v - (rand() /2);
                    end
                end 
            else
                v = (rand() * 1.5) + 1;
                for i = 1:N
                    plusMinus = randi(2);
                    if plusMinus == 1
                        lengths(i) = v + (rand());
                    else
                        lengths(i) = v - (rand());
                    end
                end 
            end
           
        end
        %Returns X and Y points of Poylgon with N sides, an N sized
        %array of angles, and integer Scale
        function [posX,posY] = randomPolygon(sides,angles,scale)
            
            posX = zeros(1,sides);
            posY = zeros(1,sides);
            %angles = randomAngles(sides);
            startAng = randi(360,1);
            for n = 1:length(angles) + 1
                if n == length(angles) + 1
                    posX(n) = posX(1);
                    posY(n) = posY(1);
                else
                    ang = startAng;
                    for l = 1:n
                        ang = ang + angles(l);
                    end
                    
                    posX(n) = cosd(ang) * scale(n);
                    posY(n) = sind(ang) * scale(n);
                    %fprintf("Angle: %i X: %.03f  Y: %.03f \n",ang,cosd(ang),sind(ang));
                end
            end
        end
        %Returns CoM of shape
        function [x, y] = centerOfMassPart(x1Pos, x2Pos, y1Pos, y2Pos)
            
            p1 = [x1Pos/2, y1Pos/2];
            p2 = [x2Pos/2, y2Pos/2];
            m1 = (y2Pos - p1(2))/(x2Pos - p1(1));
            m2 = (y1Pos - p2(2))/(x1Pos - p2(1));
            
            b1 = p1(2) - (m1 * p1(1));
            b2 = p2(2) - (m2 * p2(1));
            
            x = -(b1 - b2)/(m1-m2);
            y = m1*x + b1;
            %[linesX, linesY] = drawLines([x1Pos x2Pos], [y1Pos y2Pos]);
            %plot([x1Pos x2Pos], [y1Pos y2Pos],0,0,'rx',linesX,linesY, x, y, 'd'); %[-1 p1(1) 1], [(m1 * -1 + b1) p1(2) (m1 * 1 + b1)],'-',[-1 p2(1) 1], [(m2 * -1 + b2) p2(2) (m2 * 1 + b2)],'-');
        end
        %Used in finding of Center of Mass
        function [x, y] = centerOfMass(PosX,PosY)
            N = length(PosX);
            xP = zeros(1,N);
            yP = zeros(1,N);
            areas = zeros(1,N);
            for i = 1:length(PosX)
                if i == length(PosX)
                    [xP(i), yP(i)] = Polygon.centerOfMassPart(PosX(i),PosX(1),PosY(i),PosY(1));
                else
                    [xP(i), yP(i)] = Polygon.centerOfMassPart(PosX(i),PosX(i+1),PosY(i),PosY(i+1));
                end
            end
            
            for i = 1:length(areas)
                if i == length(PosX)
                    [rX, rY] = Polygon.findRecipricol(PosX(i),PosX(1),PosY(i),PosY(1));
                    areas(i) = Polygon.findArea(rX,rY,PosX(i),PosX(1),PosY(i),PosY(1));
                else
                    [rX, rY] = Polygon.findRecipricol(PosX(i),PosX(i+1),PosY(i),PosY(i+1));
                    areas(i) = Polygon.findArea(rX,rY,PosX(i),PosX(i+1),PosY(i),PosY(i+1));
                end
            end
            
            
            totalX = 0;
            totalY = 0;
            divide = 0;
            for i = 1:length(areas)-1
                totalX = totalX + (areas(i) * xP(i));
                totalY = totalY + (areas(i) * yP(i));
                divide = divide + areas(i);
            end
            
            x = totalX / divide;
            y = totalY / divide;
            
        end
        function [x, y] = findRecipricol(x1Pos, x2Pos, y1Pos, y2Pos)
            
            
            mY = y1Pos - y2Pos;
            mX = x1Pos - x2Pos;
            
            slope = mY/mX;
            iSlope = -mX/mY;
            b = y1Pos - (slope * x1Pos);
            x = -(b / (slope - (iSlope)));
            y = iSlope * (x);
        end
        function area = findArea(recipX,recipY,x1,x2,y1,y2)
            h = sqrt(recipX^2 + recipY^2);
            l = sqrt((x1 - x2)^2 + (y1 - y2)^2);
            area = l * h * (1/2);
        end
        function residual = move(PosX, PosY, CoMX, CoMY, sides, scale)
            
            clc; disp('Generating 100 species...')
            
            ChartTitle = 'Title';
            
            [Ymin,L] = min(PosY(1,1:sides));
            Xmin = PosX(1,L); % "Xmin" is the X-value of the coordinate which has Ymin.
            
            clc; disp('Species generation complete.')
            
            % Extra preallocated variables (to be used for rotation).
            maxsides = sides;
            
            Xdisp = zeros(1,maxsides);
            Ydisp = zeros(1,maxsides);
            radius = zeros(1,maxsides);
            theta = zeros(1,maxsides);
            CoMXdisp = 0;
            CoMYdisp = 0;
            CoMradius = zeros(1,1);
            CoMtheta = zeros(1,1);
            
            
            CoMXdisp = CoMX - Xmin;
            CoMYdisp = CoMY - Ymin;
            CoMradius = sqrt(CoMXdisp^2 + CoMYdisp^2);
            if CoMX > Xmin && CoMY > Ymin % Quadrant I
                CoMtheta = atand(CoMYdisp/CoMXdisp);
            elseif CoMX < Xmin && CoMY > Ymin % Quadrant II
                CoMtheta = 180 + atand(CoMYdisp/CoMXdisp);
            elseif CoMX < Xmin && CoMY < Ymin % Quadrant III
                CoMtheta = 180 + atand(CoMYdisp/CoMXdisp);
            elseif CoMX > Xmin && CoMY < Ymin % Quadrant IV
                CoMtheta = 360 + atand(CoMYdisp/CoMXdisp);
            elseif CoMX == Xmin && CoMY > Ymin
                CoMtheta = 90;
            elseif CoMX == Xmin && CoMY < Ymin
                CoMtheta = 270;
            elseif CoMX > Xmin && CoMY == Ymin
                CoMtheta = 0;
            elseif CoMX < Xmin && CoMY == Ymin
                CoMtheta = 180;
            elseif CoMX == Xmin && CoMY == Ymin
                CoMtheta = 0;
            end
            
            for n = 1:sides
                Xdisp(n) = PosX(n) - Xmin;
                Ydisp(n) = PosY(n) - Ymin;
                radius(n) = sqrt(Xdisp(n)^2 + Ydisp(n)^2);
                if PosX(n) > Xmin && PosY(n) > Ymin % Quadrant I
                    theta(n) = atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) < Xmin && PosY(n) > Ymin % Quadrant II
                    theta(n) = 180 + atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) < Xmin && PosY(n) < Ymin % Quadrant III
                    theta(n) = 180 + atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) > Xmin && PosY(n) < Ymin % Quadrant IV
                    theta(n) = 360 + atand(Ydisp(n)/Xdisp(n));
                elseif PosX(n) == Xmin && PosY(n) > Ymin
                    theta(n) = 90;
                elseif PosX(n) == Xmin && PosY(n) < Ymin
                    theta(n) = 270;
                elseif PosX(n) > Xmin && PosY(n) == Ymin
                    theta(n) = 0;
                elseif PosX(n) < Xmin && PosY(n) == Ymin
                    theta(n) = 180;
                elseif PosX(n) == Xmin && PosY(n) == Ymin
                    theta(n) = 0;
                end
            end
            
            
            % Display Options
            Terminate = 1;
            
            k = 1;
            if k == 0
                
            elseif k >= 1 || k <= 100
                RotModelX = zeros(360,maxsides+1);
                RotModelY = zeros(360,maxsides+1);
                RotModelCoMX = zeros(360,1);
                RotModelCoMY = zeros(360,1);
                XdispNew = Xdisp;
                YdispNew = Ydisp;
                radiusNew = radius;
                thetaNew = theta;
                CoMXdispNew = CoMXdisp;
                CoMYdispNew = CoMYdisp;
                CoMradiusNew = CoMradius;
                CoMthetaNew = CoMtheta;
                XminNew = Xmin;
                YminNew = Ymin;
                
                for k2 = 1:1:360
                    if k2 ~= 1
                        if min(RotModelY(k2-1,:)) < YminNew(k)
                            [~,L] = min(RotModelY(k2-1,1:sides(k)));
                            RotModelX(k2,L) = radiusNew(k,L)*cosd(thetaNew(k,L) - k2);
                            RotModelY(k2,L) = radiusNew(k,L)*sind(thetaNew(k,L) - k2);
                            thetaoffset = 1-atand(abs((RotModelY(k2,L))/(RotModelX(k2,L))));
                            YminNew(k) = min(RotModelY(k2-1,1:sides(k)));
                            XminNew(k) = RotModelX(k2-1,L);
                            
                            CoMXdispNew(k) = RotModelCoMX(k2-1) - XminNew(k);
                            CoMYdispNew(k) = RotModelCoMY(k2-1) - YminNew(k);
                            CoMradiusNew(k) = sqrt(CoMXdispNew(k)^2 + CoMYdispNew(k)^2);
                            if RotModelCoMX(k2-1) > XminNew(k) && RotModelCoMY(k2-1) > YminNew(k) % Quadrant I
                                CoMthetaNew(k) = atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) < XminNew(k) && RotModelCoMY(k2-1) > YminNew(k) % Quadrant II
                                CoMthetaNew(k) = 180 + atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) < XminNew(k) && RotModelCoMY(k2-1) < YminNew(k) % Quadrant III
                                CoMthetaNew(k) = 180 + atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) > XminNew(k) && RotModelCoMY(k2-1) < YminNew(k) % Quadrant IV
                                CoMthetaNew(k) = 360 + atand(CoMYdispNew(k)/CoMXdispNew(k)) + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) == XminNew(k) && RotModelCoMY(k2-1) > YminNew(k)
                                CoMthetaNew(k) = 90 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) == XminNew(k) && RotModelCoMY(k2-1) < YminNew(k)
                                CoMthetaNew(k) = 270 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) > XminNew(k) && RotModelCoMY(k2-1) == YminNew(k)
                                CoMthetaNew(k) = 0 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) < XminNew(k) && RotModelCoMY(k2-1) == YminNew(k)
                                CoMthetaNew(k) = 180 + k2 + sum(thetaoffset);
                            elseif RotModelCoMX(k2-1) == XminNew(k) && RotModelCoMY(k2-1) == YminNew(k)
                                CoMthetaNew(k) = 0 + k2 + sum(thetaoffset);
                            end
                            for n = 1:sides(k)
                                XdispNew(k,n) = RotModelX(k2-1,n) - XminNew(k);
                                YdispNew(k,n) = RotModelY(k2-1,n) - YminNew(k);
                                radiusNew(k,n) = sqrt(XdispNew(k,n)^2 + YdispNew(k,n)^2);
                                if RotModelX(k2-1,n) > XminNew(k) && RotModelY(k2-1,n) > YminNew(k) % Quadrant I
                                    thetaNew(k,n) = atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) < XminNew(k) && RotModelY(k2-1,n) > YminNew(k) % Quadrant II
                                    thetaNew(k,n) = 180 + atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) < XminNew(k) && RotModelY(k2-1,n) < YminNew(k) % Quadrant III
                                    thetaNew(k,n) = 180 + atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) > XminNew(k) && RotModelY(k2-1,n) < YminNew(k) % Quadrant IV
                                    thetaNew(k,n) = 360 + atand(YdispNew(k,n)/XdispNew(k,n)) + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) == XminNew(k) && RotModelY(k2-1,n) > YminNew(k)
                                    thetaNew(k,n) = 90 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) == XminNew(k) && RotModelY(k2-1,n) < YminNew(k)
                                    thetaNew(k,n) = 270 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) > XminNew(k) && RotModelY(k2-1,n) == YminNew(k)
                                    thetaNew(k,n) = 0 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) < XminNew(k) && RotModelY(k2-1,n) == YminNew(k)
                                    thetaNew(k,n) = 180 + k2 + sum(thetaoffset);
                                elseif RotModelX(k2-1,n) == XminNew(k) && RotModelY(k2-1,n) == YminNew(k)
                                    thetaNew(k,n) = 0 + k2 + sum(thetaoffset);
                                end
                            end
                        end
                    end
                    RotModelCoMX(k2,1) = CoMradiusNew(k)*cosd(CoMthetaNew(k) - k2) + XminNew(k);
                    RotModelCoMY(k2,1) = CoMradiusNew(k)*sind(CoMthetaNew(k) - k2) + YminNew(k);
                    for n = 1:sides(k) + 1
                        if n == sides(k) + 1
                            RotModelX(k2,n) = RotModelX(k2,1);
                            RotModelY(k2,n) = RotModelY(k2,1);
                        else
                            RotModelX(k2,n) = radiusNew(k,n)*cosd(thetaNew(k,n) - k2) + XminNew(k);
                            RotModelY(k2,n) = radiusNew(k,n)*sind(thetaNew(k,n) - k2) + YminNew(k);
                        end
                    end
                end
                
%                 t = CoMX(k):(max(RotModelCoMX)+CoMX(k))/1000:max(RotModelCoMX);
                IdealPath = zeros(1,length(RotModelCoMY));
                for t2 = 1:length(RotModelCoMY)
                    IdealPath(t2) = CoMY(k);
                end
                
                %clc; clearvars RotModel* t t2 thetaoffset; clearvars -REGEXP New$; pause(0.01); disp('Viewing complete! View another species or type "0" if done.')
            
                diff(1) = 0;
                for i = 1:length(RotModelCoMY)
                        
                        ind = mod(i,length(IdealPath));
%                         fprintf('i: %i, Ind: %i', i ,ind);
                        diff(i) = abs(RotModelCoMY(i) - IdealPath(i));
                end
                
                residual = (sum(diff)/length(RotModelCoMY)) * (1 / scale);
            
            else
                error('Invalid species selected. Must be an integer from 1 to 100.');
            end
            
      
        end
    end
end
