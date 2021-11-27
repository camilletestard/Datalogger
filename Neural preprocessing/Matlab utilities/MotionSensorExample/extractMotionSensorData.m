function [extractedData wordsOfData] = extractMotionSensorData(startBlockIdx, motionSensorData, motionDataCell, MetaData)
% Extracts accelerometer, gyroscope, and magnetometer data and scales them into
% proper range and units 
    
  
    dataStartOffset = num2cell(motionSensorData(startBlockIdx ...
        + MetaData.offsetIdx) + 1); % where this type of data starts relative to block
    wordsOfData = num2cell(motionSensorData(startBlockIdx ...
        + MetaData.numberOfWordsIndex)); % length of data in 2-byte words

    % extract data in each direction (x, y and z)
    extractedData.x = cell2mat(cellfun(@(x, y, z) x(y:MetaData.numberOfAxes:y+z-1), ...
        motionDataCell, dataStartOffset, wordsOfData, 'un', 0)); 
    
    extractedData.y = cell2mat(cellfun(@(x, y, z) x(y + 1:MetaData.numberOfAxes:y+z-1), ...
        motionDataCell, dataStartOffset, wordsOfData, 'un', 0));
    
    extractedData.z = cell2mat(cellfun(@(x, y, z) x(y + 2:MetaData.numberOfAxes:y+z-1), ...
        motionDataCell, dataStartOffset, wordsOfData, 'un', 0));
    
    
    %  scale data based on range
    extractedData = structfun(@(x) x.*((MetaData.maxVal)/2^(MetaData.numberOfBits - 1)), extractedData, 'UniformOutput', 0);



end