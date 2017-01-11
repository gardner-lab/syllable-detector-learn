load('song.mat');

files = dir('mat/*.mat');

song = NaN*zeros(size(song, 1), length(files));

w = waitbar(0, 'Loading...');

for i = 1:length(files)
    waitbar(i / length(files), w);
    load(strcat(files(i).folder, filesep, files(i).name));
    song(:, i) = data.audio;
end

[~, foo] = find(~isnan(song(1, :)));

song = song(:, foo);

close(w);

save('song.mat', 'song', 'nonsong', 'fs');
