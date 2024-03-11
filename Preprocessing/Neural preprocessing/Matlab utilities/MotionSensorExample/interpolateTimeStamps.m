function timeStampsInterp = interpolateTimeStamps(timeStamps, motionDataCell, consts, dataType)
% Interpolation of timestamps to assign a time stamp for each data point
% To minimize the effects of jitter, we interpolate between multiple blocks
% The timestamps in the last block are extrapolated based on the expected
% frequency of each signal (i.e. accelerometer, gyroscope, magnetometer)

    switch dataType
        case 'Accl'
            numberOfWordsIndex = consts.ACCL_NUMBER_OF_WORDS_INDEX;
            frequency = consts.FREQUENCY_ACCL;
        case 'Gyro'
            numberOfWordsIndex = consts.GYRO_NUMBER_OF_WORDS_INDEX;
            frequency = consts.FREQUENCY_GYRO;
        case 'Mag'
            numberOfWordsIndex = consts.MAG_NUMBER_OF_WORDS_INDEX;
            frequency = consts.FREQUENCY_MAG;
    end
    
    numBlocksForInterp = 10;% number of blocks to combine for interpolation, this averaging should help reduce effects of jitter

    numWordsPerAxis = cell2mat(cellfun(@(x) x(numberOfWordsIndex + 1)./consts.NUMBER_OF_AXES, ...
        motionDataCell, 'un', 0)); % get number of words in each block for each axis

    % check if this is first file
    isFirstFile = numWordsPerAxis(1) == 0; % first block of first file in a recording has no motion sensor data
    
     % Interpolate data between every numBlocksForInterp'th timestamp
    
     if isFirstFile 
        startBatchIndex = [2:numBlocksForInterp:length(numWordsPerAxis)];
    else
        startBatchIndex = [1:numBlocksForInterp:length(numWordsPerAxis)];
    end


    endBatchIndex = [startBatchIndex(1:end - 1)+numBlocksForInterp-1 length(numWordsPerAxis)-1];
    
    numWordsBatched = arrayfun(@(x, y) sum(numWordsPerAxis(x:y)), startBatchIndex,...
        endBatchIndex, 'un', 1); % total number of words in each timestamp range "batch"
    
     timeStampsInterp = arrayfun(@(x1, x2, n)  double(timeStamps(x1)) + (0:n-1).*(double(timeStamps(x2 + 1)) - double(timeStamps(x1)))/(n),...
         startBatchIndex, endBatchIndex, double(numWordsBatched), 'un', 0);  % actual interpolation
     timeStampsInterp = cell2mat(timeStampsInterp);
     
     % extrapolate for last block
      lastBlockTimeStamps = timeStampsInterp(end) + [1:numWordsPerAxis(end)]*frequency;
      
      timeStampsInterp = [timeStampsInterp lastBlockTimeStamps]; 
     
     
end