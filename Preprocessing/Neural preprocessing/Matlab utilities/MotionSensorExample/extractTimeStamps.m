function [timeStampsInterp] = extractTimeStamps(motionDataCell, consts)
% Time stamps are given for each block and must be interpolated to acquire time stamps for each data point within a block
% To minimize the effects of jitter, we interpolate timestamps across multiple blocks in interpolateTimeStamps function 


    timeStamps = cellfun(@(x) x(consts.TIMESTAMP_INDEX + 1), motionDataCell, 'un', 0); % extract time stamps

    % interpret as unsigned 32 bit integer
    timeStamps = cellfun(@(x) int16(x), timeStamps, 'un', 0); % change into 16 bit integers from double
    timeStamps = cell2mat(cellfun(@(x) (typecast(x, 'uint32')*consts.RESOLUTION), timeStamps, 'un', 0)); % read as 32 bit unsigned integer
    
    % interpolate data points
    timeStampsInterp.Accl = interpolateTimeStamps(timeStamps, motionDataCell, consts, 'Accl');
    timeStampsInterp.Gyro =  interpolateTimeStamps(timeStamps, motionDataCell, consts, 'Gyro');
    timeStampsInterp.Mag =  interpolateTimeStamps(timeStamps, motionDataCell, consts, 'Mag');



    

end






