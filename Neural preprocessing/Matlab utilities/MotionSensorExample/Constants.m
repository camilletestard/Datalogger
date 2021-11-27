classdef Constants
    properties (Constant = true)
        
        
        % =================== General info =================== %
        FILE_SIZE_IN_BYTES = 2^24;
        NUMBER_OF_AXES = 3; %x, y, z



        % =================== Header indices and values =================== %
        
        CONST_ID1_INDEX = 0;
        CONST_ID2_INDEX = 1;
        ACCL_OFFSET_INDEX = 2;
        GYRO_OFFSET_INDEX = 3;
        MAG_OFFSET_INDEX = 4;
        RESERVED1_INDEX = 5;
        ACCL_NUMBER_OF_WORDS_INDEX = 6;
        GYRO_NUMBER_OF_WORDS_INDEX = 7;
        MAG_NUMBER_OF_WORDS_INDEX = 8;
        RESERVED2_INDEX = 9;
        TIMESTAMP_INDEX = [10 11];
        
        CONST_ID1 = 13579;
        CONST_ID2 = 24680;
        RESERVED1_VALUE = 0;
        RESERVED2_VALUE = 0;
        
        
        
        
        
        % ================ Accelerometer info ================= %         
        G = 9.81;
        NUMBER_OF_BITS_ACCL = 16;
        FREQUENCY_ACCL = 1; %kHz
        
        
        
        
        % =================== Gyroscope info ==================== % 
        NUMBER_OF_BITS_GYRO = 16;
        FREQUENCY_GYRO = 1; % kHz
        

        
        % ================ Magnetometer info ================= % 
        NUMBER_OF_BITS_MAG = 13;
        FREQUENCY_MAG = 1/9; % kHz;
        
        
        
        % ================ Timestamp info ================= % 
        RESOLUTION = 62.5e-3; %ms
        
        
    end
end