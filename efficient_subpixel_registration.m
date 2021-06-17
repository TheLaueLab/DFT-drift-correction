%DFT-drift-correction - script for microscope drift correction by cross-correlating the brightfield images of cells.

%Authors: Aleks Ponjavic, University of Leeds; Aleksandra Jartseva, University of Cambridge
%Based on algorithm from Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
%Function for reading TIFF files >4gb from: Tristan Ursell (2020). imread_big -- read in TIFF stacks larger than 4GB (https://www.mathworks.com/matlabcentral/fileexchange/61376-imread_big-read-in-tiff-stacks-larger-than-4gb), MATLAB Central File Exchange. Retrieved October 9, 2020.
%% Disclaimer
% I have made every effort to evaluate the proper working of this code
% under many different conditions. However, it is the responsibility of
% the user to ensure that this registration code is adequate and working 
% correcntly for their application.
% Feel free to e-mail me with questions or comments at A.Ponjavic@leeds.ac.uk.

%%Inputs - files
% - Main folder with data.
% - Reference brightfield images. Since reading large files slows the scripts down a lot, there is an option to divide the single brightfield series into several files. The user needs to specify the last frame number in each file.
% - Localisations file for correction

%%Inputs - parameters
% - Pixel size in nm.
% - Upsampling factor: image cross-correlation will be with accuracy 1/upFactor pixels.
% - Are the images flipped?
% - Occasional correction with averaging: used to remove the need for continuous BF imaging and also to improve cross-correlation precision. Imaging might be done in two modes: taking a burst of images once per cycle (i.e. WL_image_WL_image_WL_image_WL), or at the beginning and end of each cycle (i.e. WL_image_WL; WL_image_WL; WL_image_WL). NB! Even if images are taken once per cycle, there needs to be a WL image both at the very beginning and end of the series! Two parameters that need to be specified in occasional mode are the number of super-res images taken in one cycle, and the number of BF images taken in a burst for averaging purposes.

%%Output
% - Output is written in the same folder as input, with an appended user-specified string.
% - If desired, whitelight images can be plotted while the script runs.
% - If input is a single bead, SD before and after correction will be calculated nd displayed for comparison.

clear all

%% Parameters from user %%

folder = ''; %main folder
fileNamesRef = {'ref1.tif', 'ref2.tif'}; % the reference brightfield image(s) of cells
lastFrameInFile = [10000 20000]; % the last frame number in each file (global)

fileNameLocs = ''; % the localisations file for correction
fileExtLocs = '.xls'; % the localisation file's extension
xcol = 10; % column with the x coordinate; 1 for ViSP, 10 for PeakFit format
ycol = 11; % column with the y coordinate; 2 for ViSP, 11 for PeakFit format
zcol = 12; % column with the z coordinate; 3 for ViSP, 12 for PeakFit format
fcol = 1; % column with the frame number; 5 for ViSP, 1 for PeakFit format
readOffset = 9; % offset loc file reading to the line where data starts; 0 for ViSP, 9 for PeakFit format

pixelSize = 271; %nm
upFactor = 100; %upsampling factor: image cross-correlation will be with accuracy 1/upFactor pixels

flip_x = false; %is the BF image flipped wrt super-res image? (x)
flip_y = true; %is the BF image flipped wrt super-res image? (y)

is.occasional = true; %are WL images taken occasionally, rather than continuously?
twice = true; %if occasionally, are they taken once per cycle or twice?
cycle = 500; %num super-res images per one cycle
burst = 20; %num WL images taken together (to increase accuracy)

isSingleBead = false; %is this a control imaging of a single bead? If so, SD will be calculated and displayed
plotWL = false; %plot every 100th WL frame?
outputAppend = '_corr'; %string to append to the output

%% Initial setup %%
fileRefs = strcat(folder, fileNamesRef);

% Read in the localisations file
dat = dlmread([folder fileNameLocs fileExtLocs], '	', readOffset, 0);

% Find indices of localisations in the first frame
locs = find(dat(:,fcol)==min(dat(:,fcol)));

% Initialise the oldDat and newDat objects, and write localisations from the
% first frame into them.
oldDat(locs,:) = dat(locs,:); %contains the old coordinates of localisations
newDat(locs,:) = dat(locs,:); %contains the corrected coordinates of localisations

%% Align every subsequent frame to the first frame
if ~is.occasional
    % Read in the first frame of the pseudoWL image
    f = im2double(imread_big(fileRefs{1}, [1 1]));
    for i = (min(dat(:,fcol))+1):max(dat(:,fcol))
        %% Read in the i'th frame of the pseudoWL image
        if(i/100==floor(i/100))
            i %display count for every 100th frame
        end
        locs = find(dat(:,fcol)==i); %indices of localisations in the i'th frame
        if(~isempty(locs)) %only perform FT if there is a localisation
            fileIdx = min(find(lastFrameInFile>=i));
            if (i <= lastFrameInFile(1))
                frame = i;
            else
                frame = i - lastFrameInFile(fileIdx-1);
            end
            g = im2double(imread_big(fileRefs{fileIdx}, [frame frame]));

            %% Sample image registration and shift correction
            % dftregistration.m from Guizar et al, Opt Lett 2008.
            % The function receives the FT of f and g and the upsampling factor.
            % Finds the peak of the cross-correlation function, using DFT and
            % cross-correlation theorem. Upsampling the cross-correlation function
            % allows to have finer than pixel precision. It is achieved by assuming
            % that all higher frequencies in the product of the DFTs are 0 before
            % doing the reverse FT.
            % The unique algorithm from Guizar et al in Aleks's test with
            % upFactor = 100 yielded the pixel shift errors of 0.0016 and 0.0043
            % pixels (no idea how Aleks measured this!) in x and y, respectively,
            % which is well within the expected accuracy of 0.01.
            % Using the conventional zero-padded FFT approach with the same
            % accuracy would have required computation of a 25,600x25,600 FFT,
            % which would have required >19GB RAM and a very comfortable chair.
            %
            % NB The code expects DC of the FTs at (1,1), so don't use fftshift.
            [output, Greg] = dftregistration(fft2(f),fft2(g),upFactor); %fft2 = 2D fast FT
            xx = output(4);
            yy = output(3);
    
            % Shift the localisations in the i'th frame, write the result into
            % newDat
            oldDat(locs,:) = dat(locs,:);
            newDat(locs,:) = dat(locs,:);
            %shift
            if flip_x
                newDat(locs,xcol)=dat(locs,xcol)-xx*pixelSize;
            else
                newDat(locs,xcol)=dat(locs,xcol)+xx*pixelSize;
            end
            if flip_y
                newDat(locs,ycol)=dat(locs,ycol)-yy*pixelSize;
            else
                newDat(locs,ycol)=dat(locs,ycol)+yy*pixelSize;
            end
            %% Plot the initial and the registered pseudoWL images in real time for every 100th frame
            if(i/100==floor(i/100) & plotWL)
                figure(1);
                subplot(1,2,1);
                imshow(abs(g));
                title(['Initial image, i = ', num2str(i)])
                subplot(1,2,2);
                imshow(abs(ifft2(Greg))); % registered image = initial -> FT -> shift -> FT
                title('Registered image')
                drawnow
            end
        end
    end
%% If drift correction is occasional, align the WL frames to the previous set of WL images and interpolate in between
elseif (is.occasional & ~twice)
    % Read in the first burst of the pseudoWL images, calculate error
    f1 = im2double(imread_big(fileRefs{1}, [1 1]));
    Fxx = zeros(burst, 1);
    Fyy = zeros(burst, 1);
    for i = 2:burst
        f = im2double(imread_big(fileRefs{1}, [i i]));
        [Foutput, Freg] = dftregistration(fft2(f1),fft2(f),upFactor);
        Fxx(i) = Foutput(4);
        Fyy(i) = Foutput(3);
    end
    AvFxx = mean(Fxx);
    AvFyy = mean(Fyy);
    disp(['Mean shift in first burst: ', num2str(AvFxx), ' in x, ', num2str(AvFyy), ' in y']);
    % Match WL images to first frame, average for each burst
    AvGxx = zeros(ceil(max(dat(:,fcol))/cycle)+1, 1);
    AvGxx(1) = AvFxx;
    AvGyy = zeros(ceil(max(dat(:,fcol))/cycle)+1, 1);
    AvGyy(1) = AvFyy;
    for c = (1:ceil(max(dat(:,fcol))/cycle))
        (c-1)*cycle+1 %display count for first frame of every cycle
        go = ~isempty(find(dat(:,fcol)>((c-1)*cycle) & dat(:,fcol)<=((c+1)*cycle)));
        if(go) %do not perform FT if there are no localisations to be corrected
            Gxx = zeros(burst, 1);
            Gyy = zeros(burst, 1);
            for i = 1:burst
                frameAbs = c*burst+i;
                fileIdx = min(find(lastFrameInFile>=frameAbs));
                if (frameAbs <= lastFrameInFile(1))
                    frame = frameAbs;
                else
                    frame = frameAbs - lastFrameInFile(fileIdx-1);
                end
                g = im2double(imread_big(fileRefs{fileIdx}, [frame frame]));
                [Goutput, Greg] = dftregistration(fft2(f1),fft2(g),upFactor);
                Gxx(i) = Goutput(4);
                Gyy(i) = Goutput(3);
            end
            AvGxx(c+1) = mean(Gxx);
            AvGyy(c+1) = mean(Gyy);
            %Interpolate
            xx = interp1([1 cycle], [AvGxx(c) AvGxx(c+1)], 1:cycle);
            yy = interp1([1 cycle], [AvGyy(c) AvGyy(c+1)], 1:cycle);

            % Shift the localisations in the corresponding frames, write the result into
            % newDat
            for i = 1:cycle
                frame = (c-1)*cycle+i;
                locs = find(dat(:,fcol)==frame); %indices of localisations in the frame
                oldDat(locs,:) = dat(locs,:);
                newDat(locs,:) = dat(locs,:);
                %shift
                if flip_x
                    newDat(locs,xcol)=dat(locs,xcol)-xx(i)*pixelSize;
                else
                    newDat(locs,xcol)=dat(locs,xcol)+xx(i)*pixelSize;
                end
                if flip_y
                    newDat(locs,ycol)=dat(locs,ycol)-yy(i)*pixelSize;
                else
                    newDat(locs,ycol)=dat(locs,ycol)+yy(i)*pixelSize;
                end
            end
            %% Plot the initial and the registered pseudoWL images in real time for last frame in each cycle
            if plotWL
                figure(1);
                subplot(1,2,1);
                imshow(abs(g));
                title(['Initial image, frame ', num2str(c*burst-1)])
                subplot(1,2,2);
                imshow(abs(ifft2(Greg))); % registered image = initial -> FT -> shift -> FT
                title('Registered image')
                drawnow
            end
        end
    end
elseif (is.occasional & twice)
    % Read in the first burst of the pseudoWL images, calculate error
    f1 = im2double(imread_big(fileRefs{1}, [1 1]));
    Fxx = zeros(burst, 1);
    Fyy = zeros(burst, 1);
    for i = 2:burst
        f = im2double(imread_big(fileRefs{1}, [i i]));
        [Foutput, Freg] = dftregistration(fft2(f1),fft2(f),upFactor);
        Fxx(i) = Foutput(4);
        Fyy(i) = Foutput(3);
    end
    AvFxx = mean(Fxx);
    AvFyy = mean(Fyy);
    disp(['Mean shift in first burst: ', num2str(AvFxx), ' in x, ', num2str(AvFyy), ' in y']);
    % Match WL images to first frame, average for each burst
    for c = (1:ceil(max(dat(:,fcol))/cycle))
        (c-1)*cycle+1 %display count for first frame of every cycle
        go = ~isempty(find(dat(:,fcol)>((c-1)*cycle) & dat(:,fcol)<=((c+1)*cycle)));
        if(go) %do not perform FT if there are no localisations to be corrected
            Gxx1 = zeros(burst, 1);
            Gyy1 = zeros(burst, 1);
            Gxx2 = zeros(burst, 1);
            Gyy2 = zeros(burst, 1);
            for i = 1:burst
                %start
                if (c ~= 1) %to avoid reading the first burst twice
                    frameAbs1 = 2*(c-1)*burst+i;
                    fileIdx1 = min(find(lastFrameInFile>=frameAbs1));
                    if (frameAbs1 <= lastFrameInFile(1))
                        frame1 = frameAbs1;
                    else
                        frame1 = frameAbs1 - lastFrameInFile(fileIdx1-1);
                    end
                    g1 = im2double(imread_big(fileRefs{fileIdx1}, [frame1 frame1]));
                    [Goutput1, Greg1] = dftregistration(fft2(f1),fft2(g1),upFactor);
                    Gxx1(i) = Goutput1(4);
                    Gyy1(i) = Goutput1(3);
                end
                %end
                frameAbs2 = (2*c-1)*burst+i;
                fileIdx2 = min(find(lastFrameInFile>=frameAbs2));
                if (frameAbs2 <= lastFrameInFile(1))
                    frame2 = frameAbs2;
                else
                    frame2 = frameAbs2 - lastFrameInFile(fileIdx2-1);
                end
                g2 = im2double(imread_big(fileRefs{fileIdx2}, [frame2 frame2]));
                [Goutput2, Greg2] = dftregistration(fft2(f1),fft2(g2),upFactor);
                Gxx2(i) = Goutput2(4);
                Gyy2(i) = Goutput2(3);
            end
            if (c == 1)
                AvGxx1 = AvFxx;
                AvGyy1 = AvFyy;
            else
                AvGxx1 = mean(Gxx1);
                AvGyy1 = mean(Gyy1);
            end
            AvGxx2 = mean(Gxx2);
            AvGyy2 = mean(Gyy2);
            %Interpolate
            xx = interp1([1 cycle], [AvGxx1 AvGxx2], 1:cycle);
            yy = interp1([1 cycle], [AvGyy1 AvGyy2], 1:cycle);
            %xx = xx - AvFxx; %matching was done to f1; AvFxx corrects for inaccuracy of f1 wrt rest of the first burst
            %yy = yy - AvFyy;

            % Shift the localisations in the corresponding frames, write the result into
            % newDat
            for i = 1:cycle
                frame = (c-1)*cycle+i;
                locs = find(dat(:,fcol)==frame); %indices of localisations in the frame
                oldDat(locs,:) = dat(locs,:);
                newDat(locs,:) = dat(locs,:);
                %shift
                if flip_x
                    newDat(locs,xcol)=dat(locs,xcol)-xx(i)*pixelSize;
                else
                    newDat(locs,xcol)=dat(locs,xcol)+xx(i)*pixelSize;
                end
                if flip_y
                    newDat(locs,ycol)=dat(locs,ycol)-yy(i)*pixelSize;
                else
                    newDat(locs,ycol)=dat(locs,ycol)+yy(i)*pixelSize;
                end
            end
            %% Plot the initial and the registered pseudoWL images in real time for last frame in each cycle
            if plotWL
                figure(1);
                subplot(1,2,1);
                imshow(abs(g));
                title(['Initial image, frame ', num2str(c*burst-1)])
                subplot(1,2,2);
                imshow(abs(ifft2(Greg))); % registered image = initial -> FT -> shift -> FT
                title('Registered image')
                drawnow
            end
        end
    end
end
%% Plot the result, show stats and write the corrected localisations file
figure(2);
plot(oldDat(:,xcol),oldDat(:,ycol),'.r')
hold on
plot(newDat(:,xcol)+200,newDat(:,ycol)+200,'.b') %displace the localisations in newDat by 200 nm for visibility
axis equal
title('Localisations')
legend('Old', 'New')
% Calculate the standard deviation of localisations - only applicable to
% a single bead.
if (isSingleBead)
    disp(['Standard deviation of corrected bead position is ', ...
        num2str(std(newDat(:,xcol))), ' in x, and ', ...
        num2str(std(newDat(:,ycol))), ' in y.'])
end
% Write the corrected localisations file
dlmwrite([folder fileNameLocs outputAppend fileExtLocs], newDat, 'delimiter', '\t','precision',7);

