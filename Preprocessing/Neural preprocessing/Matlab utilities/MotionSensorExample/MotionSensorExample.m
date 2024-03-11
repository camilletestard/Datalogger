% This script reads in a single file, extracts the motion sensor data, and plots the accelerometer,
% gyroscope and magnetometer data over time. The timestamps are in ms since
% midnight

clc; clear;


%% =================================== Initialize ========================================= %%

consts = Constants();
fileName = 'C:\NEUR0000.DT2';%your file name including full path
numChannels = 32; % enter the number of neural channels in your data set
motionSensorIndex = 3; % zero indexed channel for motion sensor data
acclMax = 2*consts.G; % m/s^2, maximum value of selected range
gyroMax = 250; %  degrees per second,  maximum value of selected range
magMax = 1200e-6; %  Teslas,  maximum value of selected range


%% ==================== Prepare for data extraction ========================== %%

motionSensorData = readMSData(fileName, consts, numChannels, motionSensorIndex); 


[startBlockIdx, motionDataCell] = divideMSDataIntoBlocks(consts, motionSensorData);

%% ==================== Time stamps ======================= %%


[timeStamps] = extractTimeStamps(motionDataCell, consts);



%% ==================== Extract accelerometer data ========================== %%

  

AcclMetadata = struct('offsetIdx', consts.ACCL_OFFSET_INDEX,...
    'numberOfWordsIndex', consts.ACCL_NUMBER_OF_WORDS_INDEX, ...
    'numberOfBits', consts.NUMBER_OF_BITS_ACCL,...
    'numberOfAxes', consts.NUMBER_OF_AXES, 'maxVal', acclMax); 



motionSensorStruct.AccelerometerData = extractMotionSensorData(startBlockIdx,...
    motionSensorData, motionDataCell, AcclMetadata);


%% ==================== Extract gyroscope data ========================== %%


GyroMetadata = struct('offsetIdx', consts.GYRO_OFFSET_INDEX,...
    'numberOfWordsIndex', consts.GYRO_NUMBER_OF_WORDS_INDEX,...
    'numberOfBits', consts.NUMBER_OF_BITS_GYRO,...
    'numberOfAxes', consts.NUMBER_OF_AXES,  'maxVal', gyroMax);  

motionSensorStruct.GyroscopeData = extractMotionSensorData(startBlockIdx,...
    motionSensorData, motionDataCell, GyroMetadata);


%% ==================== Extract magnetometer data ========================== %%


MagMetadata = struct('offsetIdx', consts.MAG_OFFSET_INDEX,...
    'numberOfWordsIndex', consts.MAG_NUMBER_OF_WORDS_INDEX,...
    'numberOfBits', consts.NUMBER_OF_BITS_MAG,...
    'numberOfAxes', consts.NUMBER_OF_AXES, 'maxVal', magMax);  

motionSensorStruct.MagnetometerData = extractMotionSensorData(startBlockIdx,...
    motionSensorData, motionDataCell, MagMetadata);




%% ==================== Plot data ========================== %%

subplot(3, 1, 1);
x = timeStamps.Accl;
plot(x, motionSensorStruct.AccelerometerData.x, 'b', x, motionSensorStruct.AccelerometerData.y, 'r', x, motionSensorStruct.AccelerometerData.z, 'g');
xlabel('Time (ms since midnight)');
ylabel('Acceleration (m/s^2)')
legend('x', 'y', 'z');

subplot(3, 1, 2);
x = timeStamps.Gyro;
plot(x, motionSensorStruct.GyroscopeData.x, 'b', x, motionSensorStruct.GyroscopeData.y, 'r', x, motionSensorStruct.GyroscopeData.z, 'g');
xlabel('Time (ms since midnight)');
ylabel('Angular velocity (deg/s)')
legend('x', 'y', 'z');

subplot(3, 1, 3);
x = timeStamps.Mag;
plot(x, motionSensorStruct.MagnetometerData.x, 'b', x, motionSensorStruct.MagnetometerData.y, 'r', x, motionSensorStruct.MagnetometerData.z, 'g');
xlabel('Time (ms since midnight)');
ylabel('Magnetic flux density (Tesla)')
legend('x', 'y', 'z');






