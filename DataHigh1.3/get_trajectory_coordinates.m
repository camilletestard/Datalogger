%get trajectory coordinates
%Find figure handles from data high Projection3D (GUI) 
cd('C:\Users\Camille Testard\Desktop\Grants_Fellowships\R37_Grant_renewal\pilot_data\session2')
load('trajectory_data.mat'); close all

%get new epoch starts
D.epochStarts = handles.D.epochStarts; 
%Get trajectory coordinates
p = handles.orths * handles.D(1).data; close all
