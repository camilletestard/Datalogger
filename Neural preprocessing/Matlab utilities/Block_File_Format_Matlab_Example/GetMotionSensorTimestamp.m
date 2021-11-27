function [ timestamps ] = GetMotionSensorTimestamp( data, blockStartIndices, motionSensorDataType, frequency )
%GETMOTIONSENSORTIMESTAMP gets timestamp of a motion sensor point
    switch motionSensorDataType
        case MotionSensorEnum.Accelerometer
            sizeOffset = MotionSensorConstants.AccelerometerLengthPosition;
        case MotionSensorEnum.Gyroscope
            sizeOffset = MotionSensorConstants.GyroscopeLengthPosition;
        case MotionSensorEnum.Magnetometer
            sizeOffset = MotionSensorConstants.MagnetometerLengthPosition;

    end
    
    numPointsInBlocks = data(blockStartIndices + sizeOffset) ./ MotionSensorConstants.NumberOfAxes;
    
    % get timestamp starting each block
    
    timestamps = cell(length(numPointsInBlocks), 1);
    for i = 1:length(blockStartIndices)
        timestampWords = data(blockStartIndices(i) + MotionSensorConstants.TimestampIndex);
        firstTimestampInBlock = double(typecast(int16(timestampWords), 'uint32') * MotionSensorConstants.TimeResolution);
        lastTimestampInBlock = firstTimestampInBlock + double(frequency / 1000 * (numPointsInBlocks(i) - 1));
        timestamps{i} = firstTimestampInBlock : frequency / 1000:lastTimestampInBlock;
    end
    
    timestamps = cell2mat(timestamps');
    



end

