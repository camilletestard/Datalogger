function motionSensorData = readMSData(fileName, consts, numChannels, motionSensorIndex)
% Reads in data as 16 bit words and extracts motion sensor column.

    if ~(exist(fileName, 'file') == 2) 
        error('The file does not exist');
    end
    fid = fopen(fileName);
    myData = fread(fid, (consts.FILE_SIZE_IN_BYTES/2), 'int16');
    fclose(fid);
    [motionSensorData] = extractMSColumn(myData, numChannels, motionSensorIndex);
    
end

function [motionSensorData] = extractMSColumn(myData, numChannels, motionSensorIndex)

    myDataMatrix = reshape(myData, numChannels, length(myData)/numChannels)'; 
    motionSensorData = myDataMatrix(:, motionSensorIndex + 1);
    
end
