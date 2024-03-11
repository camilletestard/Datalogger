classdef MotionSensorConstants
    properties (Constant = true)
        
        
        % =================== General info =================== %
%         FILE_SIZE_IN_BYTES = 2^24;
          NumberOfAxes = 3; %x, y, z



        % =================== Header indices and values =================== %
        
%         HexConstId =  
        AccelerometerOffsetPosition = 2;
        GyroscopeOffsetPosition = 3;
        MagnetometerOffsetPosition = 4;
        AccelerometerLengthPosition = 6;
        GyroscopeLengthPosition = 7;
        MagnetometerLengthPosition = 8;
        TimestampIndex = [10 11];
        
        ConstId = [13579 24680];

        
        
        
        
        % ================ Accelerometer info ================= %         
        G = 9.81;
        AccelerometerNumberOfBits = 16;
        AccelerometerFrequency = 1000; %kHz
        
        
        
        
        % =================== Gyroscope info ==================== % 
        GyroscopeNumberOfBits = 16;
        GyroscopeFrequency = 1000; % kHz
        

        
        % ================ Magnetometer info ================= %
        % use 9150 values for Spikelog16 and Ratlog64
        Magnetometer9150NumberOfBits = 13; 
        Magnetometer9150Range = 1200e-6;
        
        % for all other loggers use 9250 values
        Magnetometer9250NumberOfBits = 14; 
        Magnetometer9250Range = 4800e-6;
        
        MagnetometerFrequency = 1000; % kHz;
        MagnetometerNumberOfBits = 13;
        
        % ================ Timestamp info ================= % 
        TimeResolution = 62.5e-3; %ms
        
        
    end
end