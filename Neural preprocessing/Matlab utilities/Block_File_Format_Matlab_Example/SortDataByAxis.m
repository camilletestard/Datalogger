function [ AxisData ] = SortDataByAxis( data )
%SORTDATABYAXIS sorts motion sensor data by axis
    AxisData.X = data(1:MotionSensorConstants.NumberOfAxes:length(data));
    AxisData.Y = data(2:MotionSensorConstants.NumberOfAxes:length(data));    AxisData.X = data(1:MotionSensorConstants.NumberOfAxes:length(data));
    AxisData.Z = data(3:MotionSensorConstants.NumberOfAxes:length(data));
end

