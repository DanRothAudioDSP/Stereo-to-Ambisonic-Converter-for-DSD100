%%Stereo to Ambisonic
%parameters
n = 1; %ambisonic order
azi = [0, 90, 180, 270]; %azimuths
azi = azi .* pi/180;
elev = zeros(size(azi)); %elevations, 0 for no height
elev = elev .* pi/180;
Fs = 44100; %sample rate

%calculate Spherical Harmonics
clear front
clear right
clear back
clear left
front = zeros((n+1)^2);
left = zeros((n+1)^2);
back = zeros((n+1)^2);
right = zeros((n+1)^2);
i = 1;
for nn = 0:n
    for m = -nn:nn
            front(i) = SphericalHarmonics(azi(1), 0, nn, m);
            left(i) = SphericalHarmonics(azi(2), 0, nn, m);
            back(i) = SphericalHarmonics(azi(3), 0, nn, m);
            right(i) = SphericalHarmonics(azi(4), 0, nn, m);
            i = i + 1;
    end
end

%positions
vocalPos = left;
bassPos = back;
drumsPos = front;
otherPos = right;

%load audio
%Input and output paths
DSD100_path = ['DSD100  input path']
ambisonic_DSD100_path = ['output path']

dev_test = ["Dev","Test"]

%It takes the names of each song folder
%folder_list = {s.name};
%Make sure it takes only the correct folders (avoid '.', '..' and other
%hidden folders) by properly modifying the start number of the array
%folder_list = string(folder_list(3:end));

%For loop for all DSD100 dataset
%for idx = 1:numel(folder_list)
for g = 1:2

    dataset_path = [char(DSD100_path) filesep 'Sources' filesep char(dev_test(g))]
    
    s = dir(dataset_path);
    s = s([s.isdir]);
    s(strncmp({s.name},'.',1)) = []

    folder_list = string({s.name})

    file_ang = fopen(['ADSD100_crossangles_' char(dev_test(g)) '.txt'],'w');

    for idx = 1:50

        %Input and output path for each song
        inputpath = [char(dataset_path) filesep char(folder_list(idx)) filesep]
        outputpath_mix = [char(ambisonic_DSD100_path) filesep 'Mixtures' filesep char(dev_test(g)) filesep char(folder_list(idx)) filesep]
        outputpath_sources = [char(ambisonic_DSD100_path) filesep 'Sources' filesep char(dev_test(g)) filesep char(folder_list(idx)) filesep]
        
        clear shHrirs
        shHrirs = audioread('Google resonance spherical harmonic HRIRs');
        
        fprintf(file_ang,'%6.4f %6.4f %6.4f %6.4f\n',azi);
        
        clear vocal
        clear bass
        clear drums
        clear other
        vocal = audioread([inputpath, 'vocals.wav']);
        bass = audioread([inputpath, 'bass.wav']);
        drums = audioread([inputpath, 'drums.wav']);
        other = audioread([inputpath, 'other.wav']);
        vocal = sum(vocal,2); %make mono
        bass = sum(bass,2);
        drums = sum(drums,2);
        other = sum(other,2);

        %create B-format matrices for audio
        clear v
        clear b
        clear d
        clear o
        v = zeros(size(vocal,1),length(vocalPos));
        b = zeros(size(bass,1),length(bassPos));
        d = zeros(size(drums,1),length(drumsPos));
        o = zeros(size(other,1),length(otherPos));

            %encode to B-format
            for c = 1:length(vocalPos)
                v(:,c) = vocalPos(c) .* vocal;
                b(:,c) = bassPos(c) .* bass;
                d(:,c) = drumsPos(c) .* drums;
                o(:,c) = otherPos(c) .* other;
            end
            
            %create matrix to access sources
            sources = ["drums.wav", "vocal.wav", "bass.wav", "other.wav"];
            clear ambisources
            ambisources = [d, v, b, o];
            
            %Checking if the output directory exists. If not, creates it
            if ~exist(outputpath_sources, 'dir')
                mkdir(outputpath_sources)
                disp(outputpath_sources)
                disp('THE OUTPUT FOLDER DID NOT EXIST, CREATING IT...')
            end
            
            for k = 1:4
                %binaural decode for sources
                audioSOut = shbinauralrendersymmetric(ambisources(:,(4*k-3):(4*k)), shHrirs);

                %sum to stereo file with normalization
                norm = max([max(abs(audioSOut(:,1))),max(abs(audioSOut(:,2)))]);
                audioSOut(:,1) = audioSOut(:,1)/norm;
                audioSOut(:,2) = audioSOut(:,2)/norm;
                audiowrite([outputpath_sources char(sources(k))], audioSOut, Fs);
                clear audioSOut
            end

            %sum to ambisonic mix
            ambisoundField = zeros(size(vocal,1), length(vocalPos));
            for h = 1:length(front)
                ambisoundField(:,h) = v(:,h) + b(:,h) + d(:,h) + o(:,h);
            end

        %binaural decode

        % Check if the |audioIn| and |shHrirs| have the same number of channels.
        numAudioInChannels = size(ambisoundField, 2);
        numShHrirsChannels = size(shHrirs, 2);
        if numAudioInChannels ~= numShHrirsChannels
            error('Number of channels in input signal and HRIRs must be the same.');
        end

        audioInLength = size(ambisoundField, 1);
        shHrirLength = size(shHrirs, 1);
        audioOutLength = audioInLength + shHrirLength - 1;

        % Pre-allocate data matrices for speed.
        audioOut = single(zeros(audioOutLength, 2));

        % Zero-pad |audioIn| and |shHrirs| matrices for frequency domain convolution.
        audioIn = [ambisoundField; zeros(shHrirLength - 1, numAudioInChannels)];
        shHrirs = [shHrirs; zeros(audioInLength - 1, numShHrirsChannels)];

            for channel = 1:numShHrirsChannels
                harmonic = channel - 1; % because of 1-based indexing used in MATLAB.
                [~, m] = getnm(harmonic);
                filteredAudioIn = fftfilt(ambisoundField(:, channel), shHrirs(:, channel));
                if m < 0
                    % Asymetric spherical harmonic case.
                    audioOut(:, 1) = audioOut(:, 1) + filteredAudioIn;
                    audioOut(:, 2) = audioOut(:, 2) - filteredAudioIn;
                else
                    % Symetric spherical harmonic case.
                    audioOut(:, 1) = audioOut(:, 1) + filteredAudioIn;
                    audioOut(:, 2) = audioOut(:, 2) + filteredAudioIn;
                end
            end

        %sum to stereo file
        norm = max([max(abs(audioOut(:,1))),max(abs(audioOut(:,2)))]);
        audioOut(:,1) = audioOut(:,1)/norm;
        audioOut(:,2) = audioOut(:,2)/norm;
        
        %Checking if the output directory exists. If not, creates it
        if ~exist(outputpath_mix, 'dir')
            mkdir(outputpath_mix)
            disp(outputpath_mix)
            disp('THE OUTPUT FOLDER DID NOT EXIST, CREATING IT...')
        end
        
        audiowrite([outputpath_mix 'ambisonic.wav'], audioOut, Fs);
        clear audioOut
        clear audioInLength
        clear shHrirlength
        clear audioOutLength
        clear ambisoundField
        
    end
end
