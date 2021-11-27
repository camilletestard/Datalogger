function [ extractedData ] = ExtractMotionSensorDataByType(motionSensorData, blockStartIndices, motionSensorDataType)

%EXTRACTMOTIONSENSORDATABYTYPE extract motion sensor data based on type
%(accelerometer, gyroscope, motion sensor)
    numberOfBlocks = length(blockStartIndices);
        switch motionSensorDataType
            case MotionSensorEnum.Accelerometer
                offset = blockStartIndices + motionSensorData(blockStartIndices + MotionSensorConstants.AccelerometerOffsetPosition);
                size = motionSensorData(blockStartIndices + MotionSensorConstants.AccelerometerLengthPosition);
            case MotionSensorEnum.Gyroscope
                offset = blockStartIndices + motionSensorData(blockStartIndices + MotionSensorConstants.GyroscopeOffsetPosition);
                size = motionSensorData(blockStartIndices + MotionSensorConstants.GyroscopeLengthPosition);
            case MotionSensorEnum.Magnetometer
                offset = blockStartIndices + motionSensorData(blockStartIndices + MotionSensorConstants.MagnetometerOffsetPosition);
                size = motionSensorData(blockStartIndices + MotionSensorConstants.MagnetometerLengthPosition);
        end
        
    extractedData = zeros(sum(size), 1);

    destPosition = 1;
    for blockIdx = 1:numberOfBlocks
        sourceStart = offset(blockIdx);
        extractedData(destPosition:destPosition + size(blockIdx) - 1) = motionSensorData(sourceStart:sourceStart + size(blockIdx) -1);
        destPosition = destPosition + size(blockIdx);
    end

    


end

