function [ scaledData ] = ScaleMotionSensorData( data, numberOfBits, maxValue )
%SCALEMOTIONSENSORDATA scales data to be in range
   offset = 2 ^ (numberOfBits - 1);
   scaledData = data ./ offset  * maxValue; 
end



