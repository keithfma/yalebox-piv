function [] = yaleboxMovie(which,streak)
% function [] = yaleboxMovie(which,streak)
%
% Create movies of sandbox experiments by cropping, resizing, annotating
% and appending images using external ImageMagick commands, and compiling
% into a movie using an external ffmepg command. The movie parameters are
% coded within the function, while the type of movie to make is defined by
% the input arguments.
%
% which = String. Choose the type of movie to make: 'side', 'top', or 'all'
% streak = Logical. Choose whether to make a streak movie (1) or not (0).
%
% Keith Ma, Yale University, 11/12

nversion = 1.1;
fprintf('\nYALEBOXMOVIE version %.1f\n',nversion);
% vNext: Mark double s-points!
% v1.1: Additional cropping for rotated side images, to avoid excessive
% black space
% v1.0: Consolidated software into a single function, other minor improvements.
% v<1.0: Single programs for each type of movie, slightly worse method for making a scalebar

warning off all

%% Define parameters

%... Side parameters
sideFileStub = 'pRetroshear03SideF_*.png';
sideCoords = 'pRetroshear03SideF_Coords.mat';
sideLevel = 0.5;
sideEfoldNum = 3;

sideAngle = -10; % clockwise angle from horizontal to table top (used to remove blackspace)

sideScale.fontsize = 80; % in points
sideScale.xBL = 0.2; % % in m (woco)
sideScale.yBL = -0.05; % in m (woco)

sideCounter.fontsize = 80; % in points
sideCounter.xw = 0.42; % in m (woco)
sideCounter.yw = -0.11; % in m (woco)

sideQuality = 10; % quantizer scale (quality and file size) of final video, ranges from 1 (best, largest) to 31 (worst, smallest)

%... Top parameters
topFileStub = 'pRetroshear03Top_*.png';
topCoords = 'pRetroshear03Top_Coords.mat';
topLevel = 0.5;
topEfoldNum = 3;

topScale.fontsize = 60; % in points
topScale.xBL = -0.53; % % in m (woco)
topScale.yBL = 0.22; % in m (woco)

topCounter.fontsize = 70; % in points
topCounter.xw = -0.5; % in m (woco)
topCounter.yw = 0.03; % in m (woco)

topQuality = 30; % quantizer scale (quality and file size) of final video, ranges from 1 (best, largest) to 31 (worst, smallest)

%... Video constants
tri.b = 1; % in cm
tri.h = sqrt(3)/2*tri.b;

fps = 5;
ffmpeg_max_width = 1920; % ffmpeg maximum image width (1080p HD has dimensions of 1920*1080)
ffmpeg_max_height = 1080; % ffmpeg maximum image height

%% basic initialization steps

% select which move to make from input argument 'which'
makeSide = strcmp(which,'side');
makeTop = strcmp(which,'top');
makeAll = strcmp(which,'all');
if sum([makeTop makeSide makeAll])~=1; error('Bad value for input argument ''which''\n'); end

% unit convertions
tri.b = tri.b/100; % cm->m
tri.h = tri.h/100; % cm->m


%% initialize side, if needed
if makeSide || makeAll
    
    % get image file names
    sideFiles = dir(sideFileStub);
    if isempty(sideFiles); error('No side images found using the side file stub %s\n',sideFileStub); end
    nImage = numel(sideFiles);
    
    % get original image size
    iminfo = imfinfo(sideFiles(1).name);
    side.width = iminfo.Width;
    side.height = iminfo.Height;
    
    % get coordinate system
    if ~exist(sideCoords,'file'); error('Side coordinate file ''%s'' not found\n',sideCoords); end
    load(sideCoords,'Jw2p','Jp2w','xw0','yw0');
    side.xw0 = xw0;
    side.yw0 = yw0;
    side.Jw2p = Jw2p;
    side.Jp2w = Jp2w;
    side.world2pixel = @(XW,YW) [side.Jw2p(1,1).*(XW-side.xw0)+side.Jw2p(1,2).*(YW-side.yw0), side.Jw2p(2,1).*(XW-side.xw0)+side.Jw2p(2,2).*(YW-side.yw0)];
    side.pixel2world = @(XP,YP) [side.Jp2w(1,1).*XP+side.Jp2w(1,2).*YP+side.xw0, side.Jp2w(2,1).*XP+side.Jp2w(2,2).*YP+side.yw0];
    
    % get left side world coordinate limit for the side image (tabletop may be tilted relative to the pixel axes)
    temp = side.pixel2world(ones(side.height,1),(1:side.height)'); % get world coorindates for all pixels at the left image edge
    [~,ind] = min(abs(temp(:,2))); % find the pixel coordinate nearest to the table top (where y_world = 0)
    side.xL = temp(ind,1); % extract the x-ccordinate at the edge of the image, at the table top
    
    % get right side world coordinate limit for the side image (tabletop may be tilted relative to the pixel axes)
    temp = side.pixel2world(side.width*ones(side.height,1),(1:side.height)'); % get world coorindates for all pixels at the right image edge
    [~,ind] = min(abs(temp(:,2))); % find the pixel coordinate nearest to the table top (where y_world = 0)
    side.xR = temp(ind,1); % extract the x-ccordinate at the edge of the image, at the table top
    
    % get counter position
    temp = side.world2pixel(sideCounter.xw,sideCounter.yw); % convert counter position to pixel coords
    sideCounter.xp = temp(1);
    sideCounter.yp = temp(2);
    
    % build s-point triangle vertices in pixel coords
    temp = side.world2pixel(-tri.b/2,0); % bottom left vertex
    sideTri.xBL = temp(1);
    sideTri.yBL = temp(2);
    temp = side.world2pixel(tri.b/2,0); % bottom right vertex
    sideTri.xBR = temp(1);
    sideTri.yBR = temp(2);
    temp = side.world2pixel(0,tri.h); % peak vertex
    sideTri.xPK = temp(1);
    sideTri.yPK = temp(2);
    
    % build scale variables
    temp = side.world2pixel(0.1,0)-side.world2pixel(0,0); % 0.1m scale bar length in pixel units
    sideScale.width = round(temp(1));
    temp = side.world2pixel(sideScale.xBL,sideScale.yBL); % convert scale position to pixel coords
    sideScale.xBL = temp(1);
    sideScale.yBL = temp(2);
    
    % build counter variables
    temp = side.world2pixel(sideCounter.xw,sideCounter.yw); % convert counter position to pixel coords
    sideCounter.xp = temp(1);
    sideCounter.yp = temp(2);
    
    % setup cropping/resize for a side-only movie: no cropping, resize so that dimensions are <= ffmpegs limitations
    if makeSide
        side.cropWidth = side.width;
        side.cropHeight = side.height;
        side.cropOffsetX = 0;
        
        if side.width>ffmpeg_max_width % side image is too wide for ffmpeg, resize to fit
            side.resizeWidth = ffmpeg_max_width;
        else
            side.resizeWidth = side.width;
        end
        
        side.resizeHeight = ffmpeg_max_height;
    end
    
    % test/set threshold value if needed
    if streak
        I = imread(sideFiles(end).name);
        while 1
            f = figure;
            imshow(im2bw(I,sideLevel));
            button = questdlg('Is this threshold level acceptable?');
            if strcmp(button,'Yes')
                break
            else
                sideLevel = input('Threshold value [0 1]: ');
            end
            close(f)
        end
        
        % prepare decay background
        sideDecayFact = exp(-1/sideEfoldNum);
        SIDELAST = zeros([size(I,1), size(I,2)]); % black background for first image
        
    end
    
end

%% initialize top, if needed
if makeTop || makeAll
    
    % get image file names
    topFiles = dir(topFileStub);
    if isempty(topFiles); error('No top images found using the top file stub %s\n',topFileStub); end
    nImage = numel(topFiles); % overwrites nImage when making an All movie, but this is fine since there is a test later to make sure the sereis are of the same length
    
    % get coordinate system
    if ~exist(topCoords,'file'); error('Top coordinate file ''%s'' not found\n',topCoords); end
    load(topCoords,'Jw2p','Jp2w','xw0','yw0');
    top.xw0 = xw0;
    top.yw0 = yw0;
    top.Jw2p = Jw2p;
    top.Jp2w = Jp2w;
    top.world2pixel = @(XW,YW) [top.Jw2p(1,1).*(XW-top.xw0)+top.Jw2p(1,2).*(YW-top.yw0), top.Jw2p(2,1).*(XW-top.xw0)+top.Jw2p(2,2).*(YW-top.yw0)];
    top.pixel2world = @(XP,YP) [top.Jp2w(1,1).*XP+top.Jp2w(1,2).*YP+top.xw0, top.Jp2w(2,1).*XP+top.Jp2w(2,2).*YP+top.yw0];
    
    % get original image size
    iminfo = imfinfo(topFiles(1).name);
    top.width = iminfo.Width;
    top.height = iminfo.Height;
    
    % get left side world coordinate limit for the top image (never tilted)
    temp = top.pixel2world(1,1); % get world coorindates for the upper left corner
    top.xL = temp(1); % extract the x-ccordinate at the upper left corner
    
    % get right side world coordinate limit for the top image (never tilted)
    temp = top.pixel2world(top.width,1); % get world coorindates for the upper left corner
    top.xR = temp(1); % extract the x-ccordinate at the upper left corner
    
    % build spoint triangles in pixel coords
    temp = top.world2pixel(-tri.b/2,0); % bottom left vertex
    topTri.xBL = temp(1);
    topTri.yBL = temp(2);
    temp = top.world2pixel(tri.b/2,0); % bottom right vertex
    topTri.xBR = temp(1);
    topTri.yBR = temp(2);
    temp = top.world2pixel(0,tri.h); % bottom peak vertex
    topTri.xBP = temp(1);
    topTri.yBP = temp(2);
    % adjust position of the bottom triangle to the image edge
    offset = topTri.yBR;
    topTri.yBR = topTri.yBR-offset;
    topTri.yBL = topTri.yBL-offset;
    topTri.yBP = topTri.yBP-offset;
    
    temp = top.world2pixel(-tri.b/2,0.25); % top left vertex (ASSUMES STANDARD BOX WIDTH)
    topTri.xTL = temp(1);
    topTri.yTL = temp(2);
    temp = top.world2pixel(tri.b/2,0.25); % top right vertex
    topTri.xTR = temp(1);
    topTri.yTR = temp(2);
    temp = top.world2pixel(0,0.25-tri.h); % top peak vertex
    topTri.xTP = temp(1);
    topTri.yTP = temp(2);
    % adjust position of the top triangle to the image edge
    offset = top.height-topTri.yTR;
    topTri.yTR = topTri.yTR+offset;
    topTri.yTL = topTri.yTL+offset;
    topTri.yTP = topTri.yTP+offset;
    
    % build scale variables
    temp = top.world2pixel(0.1,0)-top.world2pixel(0,0); % 0.1m scale bar length in pixel units
    topScale.width = round(temp(1));
    temp = top.world2pixel(topScale.xBL,topScale.yBL); % convert scale position to pixel coords
    topScale.xBL = temp(1);
    topScale.yBL = temp(2);
    
    % get counter position
    temp = top.world2pixel(topCounter.xw,topCounter.yw); % convert counter position to pixel coords
    topCounter.xp = temp(1);
    topCounter.yp = temp(2);
    
    % setup cropping/resize for a top-only movie: no cropping, resize so that dimensions are <= ffmpegs limitations
    if makeTop
        top.cropWidth = top.width;
        top.cropHeight = top.height;
        top.cropOffsetX = 0;
        
        if top.width>ffmpeg_max_width % top image is too wide for ffmpeg, resize to fit
            top.resizeWidth = ffmpeg_max_width;
        else
            top.resizeWidth = top.width;
        end
        
        top.resizeHeight = ffmpeg_max_height;
    end
    
    % test/set threshold value if needed
    if streak
        
        I = imread(topFiles(end).name);
        while 1
            f = figure;
            imshow(im2bw(I,topLevel));
            button = questdlg('Is this threshold level acceptable?');
            if strcmp(button,'Yes')
                break
            else
                topLevel = input('Threshold value [0 1]: ');
            end
            close(f)
        end
        
        % prepare decay background
        topDecayFact = exp(-1/topEfoldNum);
        TOPLAST = zeros([size(I,1), size(I,2)]); % black background for first image
        
    end
    
end

%% initialize all, if needed
if makeAll % crop the smaller of the two images on both sides, resize to the larger image
    
    % check that image series make sense
    if numel(sideFiles)~=numel(topFiles); error('Different number of side and top images, resolve this issue and retry.\n'); end
    
    % crop to the smaller of the two images on the left side
    if abs(side.xL)>abs(top.xL) % side extends farther than top, crop side
        temp = side.world2pixel(top.xL,0); % get pixel position in the side image of the farthest left position in the top image
        side.cropOffsetX = ceil(temp(1,1));        
        side.cropHeight =   ceil( side.height-max(tand(sideAngle),0)*side.cropOffsetX ); % remove unneeded blackspace after cropping width    
        top.cropOffsetX = 0;
        
    elseif abs(side.xL)<abs(top.xL) % top extends farther than side, crop top
        side.cropOffsetX = 0;         
        side.cropHeight = side.height;
        temp = top.world2pixel(side.xL,0); % get pixel position in the top image of the farthest left position in the side image
        top.cropOffsetX = ceil(temp(1,1));
        
    else % miraculously the same, crop neither
        side.cropOffsetX = 0;
        side.cropHeight = side.height;
        top.cropOffsetX = 0;
        
    end
    
    % crop to the smaller of the two images on the right side
    if abs(side.xR)>abs(top.xR) % side extends farther than top, crop side
        temp = side.world2pixel(top.xR,0); % get pixel position in the side image of the farthest right position in the top image
        side.cropWidth = ceil(temp(1,1))-side.cropOffsetX;
        side.cropOffsetY = ceil( (side.width-ceil(temp(1,1)))*max(tand(-sideAngle),0) );   % remove unneeded blackspace after cropping width    
        top.cropWidth = top.width-top.cropOffsetX;
        
    elseif abs(side.xR)<abs(top.xR) % top extends farther than side, crop top
        side.cropWidth = side.width-side.cropOffsetX;
        side.cropOffsetY = 0;
        temp = top.world2pixel(side.xR,0); % get pixel position in the top image of the farthest left position in the side image
        top.cropWidth = ceil(temp(1,1))-top.cropOffsetX;
        
    else % miraculously the same, crop neither
        side.cropWidth = side.width-side.cropOffsetX;
        side.cropOffsetY = 0;
        top.cropWidth = top.width-top.cropOffsetX;
        
    end
    
    % don't crop (top) in the vertical (arbitrarily large number)
%     side.cropHeight = 1e5;
    top.cropHeight = 1e5;
    
    % resize smaller image to the width of the larger image
    if side.cropWidth>top.cropWidth        
        side.resizeWidth = side.cropWidth;
        top.resizeWidth = side.cropWidth;        
    else          
        side.resizeWidth = top.cropWidth;
        top.resizeWidth = top.cropWidth;        
    end
    
    
    
end


%% Create images


% initialize extra trimming variables
trim.width = ffmpeg_max_width; % arbitrarily high
trim.height = ffmpeg_max_height; % arbitrarily high
finalTrim = sprintf(' -resize %ix%i -crop %ix%i+0+0 ',ffmpeg_max_width,ffmpeg_max_height,trim.width,trim.height);

for i = [1 1:nImage] % repeat the first image to deal with trimming
    
    if makeSide || makeAll % prepare side image
        
        if streak % create side streak image
            IM = imread(sideFiles(i).name); % load current image
            IM = double(im2bw(IM,sideLevel)); % threshold current image
            IM = IM+sideDecayFact.*SIDELAST; % add streaks
            IM(IM>1) = 1; % limit maximum values
            SIDELAST = IM; % update streaks
            imwrite(IM,'streakSide.png');
            thisSideImage = 'streakSide.png';
            
        else % don't
            thisSideImage = sideFiles(i).name;
        end
        
        % prepare commands
        side.addScalebar = sprintf(' -verbose -gravity northwest -fill yellow -draw "polygon %.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f" -pointsize %i -draw "text %.1f,%.1f '' 10 cm''" ',...
            sideScale.xBL,sideScale.yBL,sideScale.xBL+sideScale.width,sideScale.yBL,sideScale.xBL+sideScale.width,sideScale.yBL+sideScale.fontsize,sideScale.xBL,sideScale.yBL+sideScale.fontsize,sideScale.fontsize,sideScale.xBL+sideScale.width,sideScale.yBL);
        
        side.addSptTri = sprintf(' -verbose -gravity southwest -fill yellow -draw "polygon %.1f,%.1f %.1f,%.1f %.1f,%.1f" ',...
            sideTri.xBL,sideTri.yBL,sideTri.xBR,sideTri.yBR,sideTri.xPK,sideTri.yPK);
        
        side.cropAndResize = sprintf(' -crop %ix%i+%i+%i +repage -resize %i ',...
            side.cropWidth,side.cropHeight,side.cropOffsetX,side.cropOffsetY,side.resizeWidth);
        
        side.addCounter = sprintf(' -fill yellow -pointsize %i -draw "text %.1f,%.1f ''%i''" ',...
            sideCounter.fontsize,sideCounter.xp,sideCounter.yp,i);
        
        quality = sideQuality;
        
    end
    
    if makeTop || makeAll % prepare top image
        
        if streak % create top streak image
            IM = imread(topFiles(i).name); % load current image
            IM = double(im2bw(IM,topLevel)); % threshold current image
            IM = IM+topDecayFact.*TOPLAST; % add streaks
            IM(IM>1) = 1; % limit maximum values
            TOPLAST = IM; % update streaks
            imwrite(IM,'streakTop.png');
            thisTopImage = 'streakTop.png';
            
        else % don't
            thisTopImage = topFiles(i).name;
        end
        
        top.addScalebar = sprintf(' -verbose -gravity northwest -fill yellow -draw "polygon %.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f" -pointsize %i -draw "text %.1f,%.1f '' 10 cm''" ',...
            topScale.xBL,topScale.yBL,topScale.xBL+topScale.width,topScale.yBL,topScale.xBL+topScale.width,topScale.yBL+topScale.fontsize,topScale.xBL,topScale.yBL+topScale.fontsize,topScale.fontsize,topScale.xBL+topScale.width,topScale.yBL);
        
        top.addSptTri = sprintf(' -verbose -gravity southwest -fill yellow -draw "polygon %.1f,%.1f %.1f,%.1f %.1f,%.1f" -draw "polygon %.1f,%.1f %.1f,%.1f %.1f,%.1f" ',...
            topTri.xBL,topTri.yBL,topTri.xBR,topTri.yBR,topTri.xBP,topTri.yBP,topTri.xTL,topTri.yTL,topTri.xTR,topTri.yTR,topTri.xTP,topTri.yTP);
        
        top.cropAndResize = sprintf(' -crop %ix%i+%i+0  +repage -resize %i ',...
            top.cropWidth,top.cropHeight,top.cropOffsetX,top.resizeWidth);
        
        top.addCounter = sprintf(' -fill yellow -pointsize %i -draw "text %.1f,%.1f ''%i''" ',...
            topCounter.fontsize,topCounter.xp,topCounter.yp,i);
        
        quality = topQuality;
        
    end
    
    if makeSide % make side image
        
        runCommand = ['convert ' thisSideImage side.addScalebar '- | convert -' side.addSptTri '- | convert -' side.addCounter,...
            '- | convert -' side.cropAndResize '- | convert -' finalTrim sprintf('m%03i.png',i)];
        
    end
    
    if makeTop % make top image
        
        runCommand = ['convert ' thisTopImage ' -flip ' top.addScalebar '- | convert -' top.addSptTri '- | convert -' top.addCounter,...
            '- | convert -' top.cropAndResize '- | convert -' finalTrim sprintf('m%03i.png',i)];
        
    end
    
    
    if makeAll
        
        makeTempImage = ['convert ' thisSideImage side.addScalebar '- | convert -' side.addSptTri '- | convert -' side.addCounter,...
            '- | convert -' side.cropAndResize '- | convert -' finalTrim 'temp.png'];
        
        fprintf('\nCreate temporary side image: %s\n',makeTempImage);
        system(makeTempImage,'-echo');
        
        runCommand = ['convert ' thisTopImage ' -flip ' top.addScalebar '- | convert -' top.addSptTri '- | convert -' top.addCounter,...
            '- | convert -' top.cropAndResize '- | convert -' finalTrim sprintf('- | convert - temp.png -append m%03i.png',i)];
        
    end
    
    % make image!
    fprintf('\nCreate image: %s\n',runCommand);
    system(runCommand,'-echo');
    
    % set trim on first step
    if i==1
        
        % check image
        I = imread('m001.png');
        f = figure;
        imshow(I);
        button = questdlg('Is this image acceptable?');
        close(f); 
        if ~strcmp('Yes',button); keyboard; error('Image not acceptable, program terminated'); end
        
        iminfo = imfinfo('m001.png');
        imwidth = iminfo.Width;
        imheight = iminfo.Height;
        if round(imwidth/2)~=imwidth/2
            trim.width = imwidth-1;
        else
            trim.width = imwidth;
        end
        
        if round(imheight/2)~=imheight/2
            trim.height = imheight-1;
        else
            trim.height = imheight;
        end
        
        finalTrim = sprintf(' -resize %ix%i -crop %ix%i+0+0 ',ffmpeg_max_width,ffmpeg_max_height,trim.width,trim.height);
    end
    
end

%% Create video
ffmpegCommand = sprintf('ffmpeg -r %i -i m%%03d.png -q %i -y movie.avi',fps,quality); 
fprintf(1,'\nCreate videos: %s\n',ffmpegCommand);
system(ffmpegCommand,'-echo');


