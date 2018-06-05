function residual = driveNoShow(poly)

  %Polygon Generation Script
%  This code does the following:
%    - Asks for user-based input  regarding number of sides, scale factor,
%    and inclusion of randomization of length.
%    - Creates and displays generation of 100 species per iteration of
%    code.
%    - Displays rotational movement of each user-specified species about
%    its lowermost point, and models displacement of the center of mass.

% Obligatory clean-up.
clc; %clearvars

% Graph title functions.
intmessage = 'Title';%GraphTitleFunction;


% Track inputs
[PathType,x,Track,DerTrack] = PolygonTrackFunction;
for tra = 1:length(Track)
    Track(tra + length(Track)) = Track(tra);
end
maxsides = poly.Sides;
% Track = Track + Track;

% Preallocation of variables (for shape generation).
PosX = zeros(100,maxsides+1);
PosY = zeros(100,maxsides+1);
CoMX = zeros(100,1);
CoMY = zeros(100,1);
Ymin = zeros(100,1);
Xmin = zeros(100,1);

clc; disp('Generating 100 species...')
for k = 1:1
    sides(k) = poly.Sides;
    for n = 1:sides(k)+1
        if n == sides(k) + 1
            PosX(k,n) = PosX(k,1);
            PosY(k,n) = PosY(k,1);
        else
            PosX(k,n) = 3 + poly.XPos(n); %3+scale(k)*cosd(360*(n/sides(k)))*(randi(randlengthtol(k)+2,1)/(randlengthtol(k)+2));
            PosY(k,n) = 5 + poly.YPos(n); %5+scale(k)*sind(360*(n/sides(k)))*(randi(randlengthtol(k)+2,1)/(randlengthtol(k)+2));
        end
    end
%     ChartTitle = intmessage(k,:);
    CoMX(k,1) = mean(PosX(k,1:n-1)); % Rough approximation.
    CoMY(k,1) = mean(PosY(k,1:n-1)); % Rough approximation.
    [Ymin(k),S] = min(PosY(k,1:n));
    Xmin(k) = PosX(k,S); % "Xmin" is the X-value of the coordinate which has Ymin.
%     plot(PosX(k,1:sides(k)+1),PosY(k,1:sides(k)+1),'k',CoMX(k),CoMY(k),'rx');
%     axis([-1.588*4 1.588*4 -1.25*4 1.25*4]); title('Not important'); pause(0.03);
    
end
clc; disp('Species generation complete.')

% Extra preallocated variables (to be used for rotation).
Xdisp = zeros(100,maxsides);
Ydisp = zeros(100,maxsides);
radius = zeros(100,maxsides);
theta = zeros(100,maxsides);
CoMXdisp = zeros(100,1);
CoMYdisp = zeros(100,1);
CoMradius = zeros(100,1);
CoMtheta = zeros(100,1);


for k = 1:1
    CoMXdisp(k) = CoMX(k) - Xmin(k);
    CoMYdisp(k) = CoMY(k) - Ymin(k);
    CoMradius(k) = sqrt(CoMXdisp(k)^2 + CoMYdisp(k)^2);
    if CoMX(k) > Xmin(k) && CoMY(k) > Ymin(k) % Quadrant I
        CoMtheta(k) = atand(CoMYdisp(k)/CoMXdisp(k));
    elseif CoMX(k) < Xmin(k) && CoMY(k) > Ymin(k) % Quadrant II
        CoMtheta(k) = 180 + atand(CoMYdisp(k)/CoMXdisp(k));
    elseif CoMX(k) < Xmin(k) && CoMY(k) < Ymin(k) % Quadrant III
        CoMtheta(k) = 180 + atand(CoMYdisp(k)/CoMXdisp(k));
    elseif CoMX(k) > Xmin(k) && CoMY(k) < Ymin(k) % Quadrant IV
        CoMtheta(k) = 360 + atand(CoMYdisp(k)/CoMXdisp(k));
    elseif CoMX(k) == Xmin(k) && CoMY(k) > Ymin(k)
        CoMtheta(k) = 90;
    elseif CoMX(k) == Xmin(k) && CoMY(k) < Ymin(k)
        CoMtheta(k) = 270;
    elseif CoMX(k) > Xmin(k) && CoMY(k) == Ymin(k)
        CoMtheta(k) = 0;
    elseif CoMX(k) < Xmin(k) && CoMY(k) == Ymin(k)
        CoMtheta(k) = 180;
    elseif CoMX(k) == Xmin(k) && CoMY(k) == Ymin(k)
        CoMtheta(k) = 0;
    end
    
    for n = 1:sides(k)
        Xdisp(k,n) = PosX(k,n) - Xmin(k);
        Ydisp(k,n) = PosY(k,n) - Ymin(k);
        radius(k,n) = sqrt(Xdisp(k,n)^2 + Ydisp(k,n)^2);
        if PosX(k,n) > Xmin(k) && PosY(k,n) > Ymin(k) % Quadrant I
            theta(k,n) = atand(Ydisp(k,n)/Xdisp(k,n));
        elseif PosX(k,n) < Xmin(k) && PosY(k,n) > Ymin(k) % Quadrant II
            theta(k,n) = 180 + atand(Ydisp(k,n)/Xdisp(k,n));
        elseif PosX(k,n) < Xmin(k) && PosY(k,n) < Ymin(k) % Quadrant III
            theta(k,n) = 180 + atand(Ydisp(k,n)/Xdisp(k,n));
        elseif PosX(k,n) > Xmin(k) && PosY(k,n) < Ymin(k) % Quadrant IV
            theta(k,n) = 360 + atand(Ydisp(k,n)/Xdisp(k,n));
        elseif PosX(k,n) == Xmin(k) && PosY(k,n) > Ymin(k)
            theta(k,n) = 90;
        elseif PosX(k,n) == Xmin(k) && PosY(k,n) < Ymin(k)
            theta(k,n) = 270;
        elseif PosX(k,n) > Xmin(k) && PosY(k,n) == Ymin(k)
            theta(k,n) = 0;
        elseif PosX(k,n) < Xmin(k) && PosY(k,n) == Ymin(k)
            theta(k,n) = 180;
        elseif PosX(k,n) == Xmin(k) && PosY(k,n) == Ymin(k)
            theta(k,n) = 0;
        end
    end
end

% Display Options
L = [S,0];
% L(1) = New point of contact.
% L(2) = Prior point of contact.

    k = 1;
    if k == 0
        clc; Terminate = 2; 
    elseif (k >= 1 && k <= 100) && rem(k,1) == 0
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
        thetaoffset = 0;
        
        HorizontalOffset = Xmin(k);
        VerticalOffset = Ymin(k)*ones(1,length(Track));
        Track = Track + Ymin(k);
       
        
        for k2 = 1:1:360
            
            if k2 ~= 1
                ys(1) = 0;
                A(1) = -999;
                for i = 1:sides(k)
                    
                    CheckPos = uint32(1000*round(RotModelX(k2-1,i),3));%str2double(sprintf('%d',1000*round(RotModelX(k2-1,i),3))); % Gets rid of sci. notation.
                    
                    if CheckPos > length(Track) || CheckPos < 1
                        break;
                    end
                    Track(uint32(CheckPos));
                    ys(i) = round(RotModelY(k2-1,i),3) - Track(CheckPos);
                    
                    A(i) = abs(ys(i)) < 0.1/poly.Scale;
                end
                if A(1) == -999
                    break;
                end
                if k2 == 2
                    index = find(ys == min(ys));
                    Ymin(k) = min(ys);
                    RotModelY(k2-1,:) = RotModelY(k2-1,:) - Ymin(k);
                    RotModelCoMY(k2-1) = RotModelCoMY(k2-1) - Ymin(k);
                    
                    for i = 1:sides(k)
                        num = uint32(1000*round(RotModelX(k2-1,i),3));
                        
                        CheckPos = num;%str2double(sprintf('%d',1000*round(RotModelX(k2-1,i),3))); % Gets rid of sci. notation.
                        if CheckPos > length(Track) || CheckPos < 1
                            break;
                        end
                        ys(i) = round(RotModelY(k2-1,i),3) - (Track(CheckPos));
                        A(i) = abs(round(RotModelY(k2-1,i),3) - (Track(CheckPos)- Ymin(k))) < 0.1/poly.Scale;
                        
                    end
                else
                    L(2) = L(1);
                    locations = find(A==1);
                    if length(locations) > 0
                        touchMax = max(RotModelX(k2-1,locations));
                        L(1) = locations(find(touchMax == RotModelX(k2-1,locations)));
                    end
                end
                
                
                
                thetaoffset = -1; % -atand(abs((RotModelY(k2,L(1)))/(RotModelX(k2,L(1)))));
               
                YminNew(k) = RotModelY(k2-1,L(1));
                XminNew(k) = RotModelX(k2-1,L(1));
                
                
                 perimX = RotModelX(k2-1,1:max(find(RotModelX(k2-1,:) ~= 0)));
                 perimY = RotModelY(k2-1,1:max(find(RotModelY(k2-1,:) ~= 0)));
                 perime = perimeter(perimX,perimY);
                 perime = 1;
                 for sid = 2:sides(k)+1

                     xSpace = linspace(RotModelX(k2 - 1,sid-1), RotModelX(k2 - 1,sid), 150);%round(10 * lengthOfSide / perime,0));
                     xSpace = xSpace(1:max(find(xSpace > 0)));
                     
                     x1 = RotModelX(k2 - 1,sid);
                     x2 = RotModelX(k2 - 1,sid - 1);
                     y1 = RotModelY(k2 - 1,sid);
                     y2 = RotModelY(k2 - 1,sid - 1);
                     
                     m = (y1 - y2) / (x1 - x2);
                     b = RotModelY(k2 - 1,sid) - (m * RotModelX(k2 - 1,sid));
%                      for zz = 1:
%                     f(q) = (m * q) + b;
%                     
                     lengthOfSide = sqrt((RotModelY(k2 - 1,sid) - RotModelY(k2 - 1,sid-1))^2 +  (RotModelX(k2 - 1,sid) - RotModelX(k2 - 1,sid-1))^2);
%                     
                     under(1) = -999;
                     for pp = 1:length(xSpace)
                         num = uint32(1000*round(xSpace(pp),3));
                         
                        CheckPosit(pp) = num;%str2double(sprintf('%d',1000*round(xSpace(pp),3))); % Gets rid of sci. notation.
                        
                        if CheckPosit(pp) > length(Track) || CheckPosit(pp) < 1
                            break;
                        end
                        under(pp) = round(m * (xSpace(pp)) + b,3) - (Track(CheckPosit(pp)));
                        
                     end
                     pp = 1;
                    if length(under) ~= length(xSpace)
                        break;
                    end
                     if min(under) < 0
                         minX = xSpace(min(find(under < 0)));
                         
                         RotModelY(k2-1,:) = RotModelY(k2-1,:) - under(min(find(under < 0)));
                         RotModelCoMY(k2-1) = RotModelCoMY(k2-1) - under(min(find(under < 0)));
                         
                         YminNew(k) = m * (xSpace(pp)) + b - under(min(find(under < 0)));
                         XminNew(k) = minX;
%                          L(1) = locations(find(touchMin == xSpace(locations)))
                     end
                 end
                
                
                %L(1) = find(RotModelX(k2-1,())));
                
                % Compares values of each of the points (excluding
                % that of the center of rotation) with the y-values
                % of the track.
                
                
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
            trackXX = 0:1/1000:(x-1)/1000;
%             plot(RotModelX(k2,1:sides(k)+1),RotModelY(k2,1:sides(k)+1),'k',RotModelCoMX(k2),RotModelCoMY(k2),'rx',trackXX,Track,'g');
%             axis([-2.5*4 + max(RotModelX(k2,:)) 0.676*4 + max(RotModelX(k2,:)) -1*4 + min(RotModelY(k2,:)) 2.5*3 + max(RotModelY(k2,:))]);
           
%             title('Title'); xlabel('Distance traveled'); pause(0.001);
        end
        
        % Ideal path is the overall slope of the track plus the center of
        % mass of the species.
        IdealPath = zeros(1,x);
        if PathType == 2 || PathType == 5
            for t = 1:x
                IdealPath(t) = CoMY(k) + (t/1000)*((Track(end) - Track(1))/(x/1000));
            end
        else
            for t = 1:x
                IdealPath(t) = CoMY(k);
            end
        end
        
                t = CoMX(k):(max(RotModelCoMX)+CoMX(k))/1000:max(RotModelCoMX);
                IdealPath = zeros(1,length(t));
                for t2 = 1:length(t)
                    IdealPath(t2) = CoMY(k);
                end                
%                 plot(RotModelCoMX,RotModelCoMY,'r',t,IdealPath,'b--')
%                 plot(PosX(k,1:sides+1),PosY(k,1:sides+1),'k',RotModelCoMX,RotModelCoMY,'r',t,IdealPath,'b--');
        
%          plot(PosX(k,1:sides(k)+1),PosY(k,1:sides(k)+1),'k',RotModelCoMX,RotModelCoMY,'r',t,IdealPath,'b--',0:1/1000:(x-1)/1000,Track,'g');%,0:1/1000:(x-1)/1000,IdealPath
%          axis([-1.588*4 + min(PosX(k,:)) 1.588*4 + max(RotModelX(k2,:)) -0.75*4 + Ymin(k) 1.75*4 + Ymin(k)]);
%          title('Title'); legend('Shape','Movement of CoM','Ideal Path'); pause(0.01);
        % clc; clearvars RotModel* t t2 thetaoffset; clearvars -REGEXP New$; pause(0.01); disp('Viewing complete! View another species or type "0" if done.')
        Track = Track - VerticalOffset; % Used to set track back into place after graphing.
        
        diff(1) = 0;
        for i = 1:length(RotModelCoMY)
                        
                        ind = mod(i,length(IdealPath));
%                         fprintf('i: %i, Ind: %i', i ,ind);
                  if i > length(IdealPath)
                      break;
                  end
                        diff(i) = abs(RotModelCoMY(i) - IdealPath(i));
                  
        end
                
        residual = (sum(diff)/length(RotModelCoMY)) * (1 / poly.Scale);
    else
        error('Invalid species selected. Must be an integer from 1 to 100.');
    end


end