% This is an example demonstrating how to correctly extract data from block
% file format files.

clc; clear;

%% ===============Set up for file parse============= %%

folder = 'G:\';
file = 'NEUR0000.DF1'
fileName = fullfile(folder, file);
fid = fopen(fileName, 'r');
rawData = uint8(fread(fid, Inf, 'uint8'));
fclose(fid);


%% ================Extract metadata from block header =============== %%


constId = (hex2num(HeaderConstants.HexConstId));
constIdBytes = typecast(constId, 'uint8');
blockStartIndices = FindDataBlockStart(rawData, constIdBytes);
numberOfBlocks = length(blockStartIndices);
startOfFirstHeader = blockStartIndices(1);
endOfFirstHeader = startOfFirstHeader + HeaderConstants.HeaderTotalBytes;
firstHeader = rawData(startOfFirstHeader:endOfFirstHeader);
HeaderStruct = ExtractHeaderData(firstHeader);

%% ================Extract neural data from blocks =============== %%

%check where each type data is in partition info
neuralIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.NeuralData), HeaderStruct.PartitionInfo, 'un', 0)));
isNeuralPresent = ~isempty(neuralIndex);
if isNeuralPresent
    neuralDataAsBytes = ExtractDataByType(rawData, HeaderStruct, neuralIndex, blockStartIndices, numberOfBlocks);
end


%% ================Convert neural data to physical units (V) =============== %%
% get meta data from file start event using event file reader
if isNeuralPresent
    numberOfAdcBits = 16;
    offset = 2 ^ (numberOfAdcBits - 1);
    voltageResolution = 1.95e-7;
    numberOfChannels = 64;
    frequency = 32000;

    %cast bytes to unsigned 16 bit integers and store as float
    neuralData= single(typecast(neuralDataAsBytes, 'uint16'));

    %convert uint16 values to voltage values
    for dataPointIndex = 1:length(neuralData)
        neuralData(dataPointIndex) = voltageResolution * (neuralData(dataPointIndex) - offset); 
    end

    % sort by channels
    neuralDataMat = reshape(neuralData, numberOfChannels, []);

    % get timestamps of neural data
    timestampsNeural = GetTimestamps(HeaderStruct.Timestamp, frequency, size(neuralDataMat, 2));
end
%% ================Extract audio data from blocks =============== %%

%check where each type data is in partition info
audioIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.Audio), HeaderStruct.PartitionInfo, 'un', 0)));
isAudioPresent = ~isempty(audioIndex);
if isAudioPresent
    audioDataAsBytes = ExtractDataByType(rawData, HeaderStruct, audioIndex, blockStartIndices, numberOfBlocks);
end

%% ================Scale audio data and save as wav =============== %%
% get meta data from file start event using event file reader

if isAudioPresent
    numberOfAudioBits = 15;
    frequency = 200000;
    isAudioSigned = true;

    audioData = ScaleAudioData(isAudioSigned, audioDataAsBytes, numberOfAudioBits); 
    % save audio data as wav file
    % audiowrite('C:\example200kHz.wav', audioData, frequency, 'BitsPerSample', 16);
    % get timestamps of audio data
    timestampsAudio = GetTimestamps(HeaderStruct.Timestamp, frequency, length(audioData));
end

%% ================Extract motion sensor data from blocks =============== %%

%check where each type data is in partition info
motionSensorIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.MotionSensor), HeaderStruct.PartitionInfo, 'un', 0)));
isMotionSensorPresent = ~isempty(motionSensorIndex);
if isMotionSensorPresent
    motionSensorAsBytes = ExtractDataByType(rawData, HeaderStruct, motionSensorIndex, blockStartIndices, numberOfBlocks);

end

%% ================Sort motion sensor data by data type =============== %%
% The values for acclMax and gyroMax are chosen by the user. They can be found using the Event
% File Viewer in the file started event.
acclMax = 2*MotionSensorConstants.G; % m/s^2, maximum value of selected range
gyroMax = 250; %  degrees per second,  maximum value of selected range
magMax = MotionSensorConstants.Magnetometer9250Range; %  Teslas,  maximum value of selected range

% extract motion sensor data from inner block

if isMotionSensorPresent

    
    motionSensorData = single(typecast(motionSensorAsBytes, 'int16'));    
    blockStartIndices = FindMotionSensorBlockStart(motionSensorData, MotionSensorConstants.ConstId);    
    accelerometerDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndices, MotionSensorEnum.Accelerometer);
    gyroscopeDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndices, MotionSensorEnum.Gyroscope);
    magnetometerDataTemp = ExtractMotionSensorDataByType(motionSensorData, blockStartIndices, MotionSensorEnum.Magnetometer);
    
    scaledAccelerometer = ScaleMotionSensorData(accelerometerDataTemp, MotionSensorConstants.AccelerometerNumberOfBits, acclMax);
    scaledGyroscope = ScaleMotionSensorData(gyroscopeDataTemp, MotionSensorConstants.GyroscopeNumberOfBits, gyroMax);
    scaledMagnetometer = ScaleMotionSensorData(magnetometerDataTemp, MotionSensorConstants.Magnetometer9250NumberOfBits, magMax);

    AccelerometerData = SortDataByAxis(scaledAccelerometer);
    GyroscopeData = SortDataByAxis(scaledGyroscope);
    MagnetometerData = SortDataByAxis(scaledMagnetometer);
    
    % get timestamp of a particular point
    index = 1000; % want the timestamp of the 200000th point
    timestampsAccelerometer = GetMotionSensorTimestamp(motionSensorData, blockStartIndices, MotionSensorEnum.Accelerometer, MotionSensorConstants.AccelerometerFrequency);
    timestampsGyroscope = GetMotionSensorTimestamp(motionSensorData, blockStartIndices, MotionSensorEnum.Gyroscope, MotionSensorConstants.GyroscopeFrequency);
    timestampsMagnetometer = GetMotionSensorTimestamp(motionSensorData, blockStartIndices, MotionSensorEnum.Magnetometer, MotionSensorConstants.MagnetometerFrequency);

    % example: plot accelerometer
    plot(timestampsAccelerometer, AccelerometerData.X);
    hold on;
    plot(timestampsAccelerometer, AccelerometerData.Y);
    hold on;
    plot(timestampsAccelerometer, AccelerometerData.Z);
    ylim([-acclMax, acclMax]);
end

%% ================Extract RTK data if present =========== %%
RTKIndex = find(cell2mat(arrayfun(@(x) x.DataType == uint32(DataTypeEnum.RTK), HeaderStruct.PartitionInfo, 'un', 0)));
isRTKPresent = ~isempty(RTKIndex);
if isRTKPresent
    RTKDataPartitions = ExtractDataByType(rawData, HeaderStruct, RTKIndex, blockStartIndices, numberOfBlocks); % extracts allocated bytes for RTK, not all of these bytes will contain data
    [RTKData, RTKTimestamps] = ExtractRTKData(RTKDataPartitions, HeaderStruct, rawData, blockStartIndices);
end
