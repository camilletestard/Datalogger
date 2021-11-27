function [startBlockIdx, motionDataCell] = divideMSDataIntoBlocks(consts, motionSensorData)
% Divides motion sensor data into the blocks in which they were written
% The beggining of each block is recognized by a constant identifier

    ConstIDArray = [consts.CONST_ID1; consts.CONST_ID2]; % block start constant identifiers
    startBlockIdx = findstr(motionSensorData', ConstIDArray'); % index of where each block starts
    if isempty(startBlockIdx)
        error('No motion sensor data found. Check that your motion sensor channel index has been entered correctly.');
    end
    numBlocks = length(startBlockIdx);
    
    % change form of data to cell array of blocks
    motionDataMat = reshape(motionSensorData, length(motionSensorData)/numBlocks, numBlocks);
    motionDataCell = num2cell(motionDataMat, 1)';
    
end
