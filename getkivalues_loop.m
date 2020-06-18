projdir = '/scratch/kg98/kristina/Projects/MMH002/';
fileID = fopen([projdir,'MMH002_PET_subjects.txt']);
ParticipantIDs = textscan(fileID,'%s');
ParticipantIDs = ParticipantIDs{1};
% compute numsubs
numSubs = length(ParticipantIDs);

for i = 1:numSubs
    subject = ParticipantIDs{i}
    getkivalues_02(projdir, subject)
end