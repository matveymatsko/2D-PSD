function [elevation] = elevationRemover(elevation, unwantedElevation)

x = size(elevation,1); 
y = size(elevation,2); 

for i = 1:x    
    for j = 1:y  
        if elevation(i,j) == unwantedElevation
            elevation(i,j) = 0; 
        end 
    end
end

end