# Script to pull together ME data

library(xlsx)

#Load data
file = file.choose() # choose the formatted behavior file
log <- xlsx::read.xlsx(file, 1)

LeftVideo_file = file.choose() # choose the left video ME file
ME_left_file = read.delim(LeftVideo_file, header = TRUE, sep = " ")
ME_left = ME_left_file[,1]

RightVideo_file = file.choose() # choose the right video ME file
ME_right_file = read.delim(RightVideo_file, header = TRUE, sep = " ")
ME_right = ME_right_file[,1]

#First important note: a few frames difference between the two videos.
diff_vids = nrow(ME_right_file) - nrow(ME_left_file)
print(paste('Difference in frames between videos is:', abs(diff_vids), 'frames.'))
if (diff_vids <0){
  ME_right[(length(ME_right)+1):(length(ME_right)+abs(diff_vids))]=0
} else {
  ME_left[(length(ME_left)+1):(length(ME_left)+abs(diff_vids))]=0
}

MotionEnergy = cbind(ME_left, ME_right)

mac =0
if (mac == 1) { home = '~'} else
{home = 'C:/Users/GENERAL'}

setwd(paste(home,'/Dropbox (Penn)/Datalogger/Deuteron_Data_Backup/Ready to analyze output/', sep="" ))
write.csv(MotionEnergy, 'Amos_2021-07-29/ME.csv',row.names=F)

# #Visualize ME
# plot(1:3000, ME_right[1:3000], type = "S", frame = FALSE, 
#      col = "red", xlab = "x", ylab = "y")
# lines(1:3000, ME_left[1:3000], col = "blue", type = "S", lty = 2)

#Get recording length in sec
length_recording = log[nrow(log),"Time"]/1000
#Important note: the length of the recording session should not be more than 2sec different between the
#deuteron log and "unit rasters". That is because in the script "log_read_sorted_NEX_file_CombinedChannels.m"
#Unit_rasters is initialized by the line: Unit_rasters(unit,:) = zeros(1,round(length_recording*temp_resolution)+1);
#round and +1 can lead to up to 2sec difference with the length of recording from the log.

#Get video length in sec
length_video = (log[nrow(log)-1,"Time"] - log[2,"Time"])/1000
frame_rate = 29.973
expected_num_frames = frame_rate*length_video
diff_frames = abs(expected_num_frames-length(ME_right))
