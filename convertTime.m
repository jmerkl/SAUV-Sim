function t = convertTime(time)
[hour, minutes] = strtok(time, ':');
minutes = str2num(minutes(2:end));
t = str2num(hour)*3600 + minutes*60;
end