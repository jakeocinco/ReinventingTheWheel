function p = perimeter(XPos, YPos)
    
    for i = 2:length(XPos)
   
            z(i) = sqrt((XPos(i) - XPos(i-1))^2 + (YPos(i) - YPos(i-1))^2);
  
        
    end
    p = sum(z);
end