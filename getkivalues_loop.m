projdir = '/scratch/kg98/kristina/Projects/MMH002/';
addpath([projdir,'code/'])
fileID = fopen([projdir,'MMH002_PET_subjects.txt']);
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
% compute numsubs
numSubs = length(ParticipantIDs);

for i = 41:41
    subject = ParticipantIDs{i}
    getkivalues_02(projdir, subject)
end