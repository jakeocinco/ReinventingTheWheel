% Polygon Track Function
% This code does the following:
%   - Allows the user to choose from a variety of different track types to
%   choose from.
%   - Allows the user to customize preset tracks by modifying specific
%   sections.

function [PathType,maxcount,Track,DerTrack] = PolygonTrackFunction
clc; clearvars; clear;
PathType = 3; %menu('Which type of path will be drawn?','Parabolic','Linear','Zigzag','Cycloidal','Freehand');
AcceptPreview = 2;
switch PathType
    case 1 % Parabolic pathtype
        while AcceptPreview == 2
            OriginalTrackPeriod = input('Specify a distance between minimum points (m): ');
            OriginalTrackAmplitude = input('Specify a maximum height between minimum points (m): ' );
            % Ax = B
            % y(x) = ax^2 + bx + c
            % y = Amplitude; x = Period
            A = [0 0 1; (OriginalTrackPeriod/2)^2 (OriginalTrackPeriod/2) 1; OriginalTrackPeriod^2 OriginalTrackPeriod 1];
            B = [0; OriginalTrackAmplitude; 0];
            OriginalTrackCoefficients = A\B;
            Track = zeros(1,5*OriginalTrackPeriod+1);
            offset = 0;
            count = 0;
            while count/1000 < 5*OriginalTrackPeriod+1
                count = count + 1;
                if rem(count/1000,OriginalTrackPeriod) == 0
                    offset = offset + OriginalTrackPeriod;
                end
                Track(count) = OriginalTrackCoefficients(1)*(count/1000 - offset)^2 + OriginalTrackCoefficients(2)*(count/1000 - offset) + OriginalTrackCoefficients(3);
            end
            maxcount = count;
            %plot(0:1/1000:(maxcount-1)/1000,Track,'k');
            axis([0 maxcount/1000 0 maxcount/5000]); title('Track Preview');
            while AcceptPreview == 2
                AcceptPreview = menu('Is this track acceptable?','Yes','No');
                if AcceptPreview == 1
                    break;
                end
                Action = menu('What do you want to do?','Set new specifications','Modify a specific section','Copy a specific section');
                switch Action
                    case 0
                        clc; clearvars -EXCEPT AcceptPreview; break;
                    case 1
                        clc; clearvars -EXCEPT AcceptPreview; break;
                    case 2
                        LeftEndpoint = input('Specify a left endpoint of track to modify (m): ');
%                         plot(0:1/1000:(maxcount-1)/1000,Track,'k',LeftEndpoint,Track(1000*LeftEndpoint+1),'rx');
                        RightEndpoint = input('Specify a right endpoint of track to modify (m): ');
                        if RightEndpoint <= LeftEndpoint
                            error('Invalid range selected')
                        end
                        HighlightedTrack = Track(LeftEndpoint*1000 + 1:1:RightEndpoint*1000 + 1);
%                         plot(0:1/1000:(maxcount-1)/1000,Track,'k',(LeftEndpoint:1/1000:RightEndpoint),HighlightedTrack,'r');
%                         axis([0 maxcount/1000 0 maxcount/5000]); title('Track Editor');
                        TrackPeriod = input('Specify a distance between these endpoints (m): ');
                        TrackAmplitude = input('Specify a maximum height between these endpoints (m): ' );
                        % Ax = B
                        % y(x) = ax^2 + bx + c
                        % y = Amplitude; x = Period
                        A = [0 0 1; (TrackPeriod/2)^2 (TrackPeriod/2) 1; TrackPeriod^2 TrackPeriod 1];
                        B = [Track(1000*LeftEndpoint+1); TrackAmplitude; Track(1000*RightEndpoint+1)];
                        TrackCoefficients = A\B;
                        
                        if TrackPeriod ~= OriginalTrackPeriod
                            maxcount = length(Track) + 1000*(TrackPeriod - OriginalTrackPeriod) + 1;
                        end
                        offset = 0;
                        trackoffset = LeftEndpoint + 1/1000;
                        for count = 1:maxcount % Offset modifications.
                            if rem(count/1000,OriginalTrackPeriod) == 0
                                offset = offset + OriginalTrackPeriod;
                            end
                            
                            if count < 1000*LeftEndpoint + 1 % Unchanged left-hand portion of track
                                Track(count) = Track(count);
                            elseif count >= 1000*LeftEndpoint + 1 && count <= 1000*(LeftEndpoint+TrackPeriod) + 1 % Altered portion of selected track
                                Track(count) = TrackCoefficients(1)*(count/1000 - trackoffset)^2 + TrackCoefficients(2)*(count/1000 - trackoffset) + TrackCoefficients(3);
                            else % Offset portion of right-hand track
                                Track(count) = OriginalTrackCoefficients(1)*(count/1000 - offset)^2 + OriginalTrackCoefficients(2)*(count/1000 - offset) + OriginalTrackCoefficients(3);
                            end
                        end
                        
                        offset = LeftEndpoint;
                        for count = 1000*LeftEndpoint+1:1:1000*(LeftEndpoint+TrackPeriod)+1 % Alteration of selected track.
                            Track(count) = TrackCoefficients(1)*(count/1000 - offset)^2 + TrackCoefficients(2)*(count/1000 - offset) + TrackCoefficients(3);
                        end
                        
%                         plot(0:1/1000:(maxcount-1)/1000,Track,'k');
%                         axis([0 maxcount/1000 0 maxcount/5000]); title('Track Preview');
                        
                    case 3
                        LeftEndpoint = input('Specify a left endpoint of track to modify (m): ');
                        plot(0:1/1000:(maxcount-1)/1000,Track,'k',LeftEndpoint,Track(1000*LeftEndpoint+1),'rx');
                        RightEndpoint = input('Specify a right endpoint of track to modify (m): ');
                        if RightEndpoint <= LeftEndpoint
                            error('Invalid range selected')
                        end
                        HighlightedTrack = Track(LeftEndpoint*1000 + 1:1:RightEndpoint*1000 + 1);
                        plot(0:1/1000:(maxcount-1)/1000,Track,'k',(LeftEndpoint:1/1000:RightEndpoint),HighlightedTrack,'r');
                        axis([0 maxcount/1000 0 maxcount/5000]); title('Track Editor');
                end
            end
        end
        
    case 2 % Linear pathtype
        while AcceptPreview == 2
            OriginalTrackSlope = input('Specify the slope of the track: ');
            OriginalTrackLength = input('Specify the length of the track (m): ');
            % y(x) = mx + b
            Track = zeros(1,1000*OriginalTrackLength+1);
            count = 0;
            for count = 1:1:1000*OriginalTrackLength+1
                Track(count) = OriginalTrackSlope*(count/1000);
            end
            maxcount = count;
            plot(0:1/1000:(maxcount-1)/1000,Track,'k'); title('Track Preview')
            
            while AcceptPreview == 2
                AcceptPreview = menu('Is this track acceptable?','Yes','No');
                if AcceptPreview == 1
                    break;
                end
                Action = menu('What do you want to do?','Set new specifications','Modify a specific section');
                switch Action
                    case 0
                        break;
                    case 1
                        break;
                    case 2
                        LeftEndpoint = input('Specify a left endpoint of track to modify (m): ');
                        plot(0:1/1000:(maxcount-1)/1000,Track,'k',LeftEndpoint,Track(1000*LeftEndpoint+1),'rx');
                        RightEndpoint = input('Specify a right endpoint of track to modify (m): ');
                        if RightEndpoint <= LeftEndpoint
                            error('Invalid range selected')
                        end
                        HighlightedTrack = Track(LeftEndpoint*1000 + 1:1:RightEndpoint*1000 + 1);
                        plot(0:1/1000:(maxcount-1)/1000,Track,'k',(LeftEndpoint:1/1000:RightEndpoint),HighlightedTrack,'r');
                        
                        TrackSlope = input('Specify the slope of this segment: ');
                        TrackLength = input('Specify the length of this segment: ');
                        
                        if TrackLength ~= RightEndpoint - LeftEndpoint
                            maxcount = length(Track) + 1000*(TrackLength - (RightEndpoint-LeftEndpoint)) + 1;
                        end
                        
                        rightoffset = (LeftEndpoint+TrackLength)+(1/1000);
                        leftoffset = LeftEndpoint+(1/1000);
                        for count = 1:maxcount
                            if count <= 1000*LeftEndpoint + 1 % Unmodified left-hand portion of track.
                                Track(count) = Track(count);
                            elseif count > 1000*LeftEndpoint + 1 && count <= Track(1000*(LeftEndpoint + TrackLength)+1) % Modified selected portion of track.
                                Track(count) = TrackSlope*(count/1000-leftoffset) + Track(1000*LeftEndpoint+1);
                            else % Offset right-hand portion of track.
                                Track(count) = OriginalTrackSlope*(count/1000-rightoffset) + Track(1000*(LeftEndpoint + TrackLength)+1);
                            end
                        end
                        plot(0:1/1000:(maxcount-1)/1000,Track,'k');
                end
            end
        end
    case 3 % Zig-Zag pathtype
        OriginalTrackPeriod = 2;%input('Specify a distance between minimum points (m): ');
        OriginalTrackAmplitude = .75;%input('Specify the maximum height between minimum points (m): ');
        % y(x) = -|(2b/a)(x-(a/2))| + b
        % Period = a
        % Amplitude = b
        Track = zeros(1,5*OriginalTrackPeriod+1);
        offset = OriginalTrackPeriod/2;
        count = 0;
        while count/1000 < 5*OriginalTrackPeriod+1
            count = count + 1;
            if rem(count/1000,OriginalTrackPeriod) == 0
                offset = offset + OriginalTrackPeriod;
            end
%             if rem(count,2) == 0
                Track(count) = -abs(((2*OriginalTrackAmplitude)/OriginalTrackPeriod)*((count/1000)-offset)) + OriginalTrackAmplitude;
%             end
%             Track2(count) = -abs(((2*OriginalTrackAmplitude)/OriginalTrackPeriod)*((count/1000)-offset)) + OriginalTrackAmplitude;
        end
        maxcount = count;
         plot(1:length(Track),Track,'k');
         axis([0 length(Track) 0 3]); title('Track Preview');
        
        while AcceptPreview == 2
            AcceptPreview = 1;%menu('Is this track acceptable?','Yes','No');
            if AcceptPreview == 1
                break;
            end
            Action = menu('What do you want to do?','Set new specifications','Modify a specific section');
            switch Action
                case 0
                    clc; clearvars -EXCEPT AcceptPreview; break;
                case 1
                    clc; clearvars -EXCEPT AcceptPreview; break;
                case 2
                    LeftEndpoint = input('Specify a left endpoint of track to modify (m): ');
%                     plot(0:1/1000:(maxcount-1)/1000,Track,'k',LeftEndpoint,Track(1000*LeftEndpoint+1),'rx');
                    axis([0 maxcount/1000 0 maxcount/5000]); title('Track Editor');
                    RightEndpoint = input('Specify a right endpoint of track to modify (m): ');
                    if RightEndpoint <= LeftEndpoint
                        error('Invalid range selected')
                    end
                    HighlightedTrack = Track(LeftEndpoint*1000 + 1:1:RightEndpoint*1000 + 1);
                    plot(0:1/1000:(maxcount-1)/1000,Track,'k',(LeftEndpoint:1/1000:RightEndpoint),HighlightedTrack,'r');
                    axis([0 maxcount/1000 0 maxcount/5000]); title('Track Editor');
                    TrackPeriod = input('Specify a distance between these endpoints (m): ');
                    TrackAmplitude = input('Specify a maximum height between these endpoints (m): ' );
                    
                    if RightEndpoint - LeftEndpoint ~= TrackPeriod
                        maxcount = length(Track) + 1000*(TrackPeriod - (RightEndpoint-LeftEndpoint)) + 1;
                    end
                    
                    offset = TrackPeriod/2;
                    trackoffset = OriginalTrackPeriod/2;
                    for count = 1:maxcount
                        if rem(count/1000,OriginalTrackPeriod) == 0
                            trackoffset = trackoffset + OriginalTrackPeriod;
                        end
                        if count < 1000*LeftEndpoint + 1 % Unchanged left-hand portion of track
                            Track(count) = Track(count);
                        elseif count >= 1000*LeftEndpoint + 1 && count <= 1000*(LeftEndpoint+TrackPeriod) + 1 % Altered portion of selected track
                            Track(count) = -abs(((2*TrackAmplitude)/TrackPeriod)*((count/1000)-offset-LeftEndpoint)) + TrackAmplitude + Track(1000*LeftEndpoint);
                        else % Offset portion of right-hand track
                            Track(count) = -abs(((2*OriginalTrackAmplitude)/OriginalTrackPeriod)*((count/1000)-trackoffset)) + OriginalTrackAmplitude;
                        end
                    end
                    plot(0:1/1000:(maxcount-1)/1000,Track,'k');
                    axis([0 maxcount/1000 0 maxcount/5000]); title('Track Editor');
            end
        end
    case 4 % Cycloidal pathtype
        sides = 4;%input('Specify the number of sides the polygon needs to ride ideally: ');
        if rem(sides,1) ~= 0 || sides < 3
            error('Number of sides must be an integer greater than or equal to 3');
        end
        sidelength = 2;%input('Enter one of the side lengths of the polygon (m): ');
        if sidelength <= 0
            error('Radius must be greater than zero');
        end
        % y(x) = d - a*cosh(x/a); d = a/(cosd(180/n))
        a = sidelength/(2*tand(180/sides)); % Apothem of the polygon
        d = a/(cosd(180/sides)); % Hypotenuse of the apothem-side triangle
        Track = zeros(1,sidelength*sides);
        x = zeros(1,sidelength*sides);
        x(1) = 0;
        y = -acosh(d);
        for count = 1:1000*sidelength*sides
            x(count+1) = x(count) + acosh(d)/500;
            y = y + acosh(d)/500;
            if rem(count+1,1000) == 0
                y = -acosh(d);
            end
            Track(count) = d - (a*cosh(y/a));
        end
        maxcount = count;
        %plot(x(1:1000*sidelength*sides),Track,'k');
        axis([0 maxcount/1000 0 0.78*maxcount/1000]); title('Track Preview');
        
    case 5 % User-drawn pathtype
        
end

DerTrack = zeros(1,maxcount);
% for count = 1:maxcount
%     if count == 1
%         DerTrack(count) = (Track(count+1) - Track(count)) / 0.001;
%     elseif count == maxcount
%         DerTrack(count) = (Track(count) - Track(count-1)) / 0.001;
%     else
%         DerTrack(count) = (Track(count+1) - Track(count-1)) / (2*0.001);
%     end
% end
end