function [Jp2w,Jw2p,xw0,yw0] = getwoco(imfile,xspt)
% function [Jp2w,Jw2p,xw0,yw0] = getwoco(imfile,xspt)
%
% Compute transformation from pixel to world coordinates.  Done by
% interactively inputing reference points and computing a least-squares
% solution for the transformation parameters. 
%
% Transformation from pixel to world coordinates is computed as:
%
%   [xW,yW] = Jp2w*[xP;yP]+[xw0;yw0]; (for coordinates)
%   [uW,vW] = Jp2w*[uP;vP]; (for displacements)
%
% Inverse transformation from world to pixel coordinates is computed as:
%
%   [xP,yP] = Jw2p*([xP;yP]-[xw0;yw0]); (for coordinates)
%   [uP,vP] = Jw2p*[uP;vP]; (for displacements)
%
% Arguments:
% 
% imfile = String. Filename of woco image.
%
% xspt = Scalar.  Known x-position of the true s-point in world
% x-coordinate.  Useful for cases where the woco grid is centered rather
% than lined up with the front of the centerpeice as it should be.  Default
% is 0.

% Returns:
%
% Jp2w, Jw2P = Jacobian matrix for the transformation from pixel to world
% coordinates and vice versa.
%
% xw0, yw0 = Offsets in the x- and y-directions
%
% Saves:
% 
% All returned variables in the generic MAT file "getwoco_output.mat"
%
% Keith Ma, Yale Univerisity 2012

try
    
    if nargin==1
        xspt = 0;
    end
    
    wocoim = imread(imfile);
    
    % click in all points
    f = figure;
    imshow(wocoim)
    title('Navigate with w-a-s-d-i-o, New segment with n, ENTER to finish or retry')
    [xp,yp] = myginput;
    
    % adjust any points that lie outside the image
    [nr,nc] = size(wocoim);
    outside = xp<1 | xp>nc | yp<1 | yp>nr;
    m = msgbox(sprintf('Adjusted %i point(s) that were outside the image\n',sum(outside)),'!','warn');
    uiwait(m);
    xp(xp<1)=1; yp(yp<1)=1;
    xp(xp>nc)=nc; yp(yp>nr)=nr;
    
    % prepare
    npts = numel(xp);
    xw = nan(size(xp));
    yw = nan(size(xp));
    
    % enter woco for all points
    close(f);
    figure
    imshow(wocoim)
    title('Enter world coordinates for each point (m)')
    for i = 1:npts
        hold on
        p = plot(xp(i),yp(i),'*r','MarkerSize',15);
        inputwoco = inputdlg({'X (m):','Y (m):'},'Input woco for the red point');
        xw(i) = str2double(inputwoco{1});
        yw(i) = str2double(inputwoco{2});
        set(p,'Color','k');
    end
    
    % adjust for known offset
    xw = xw-xspt;
    
    % compute pixel-to-world transformations
    
    Jp2w = nan(2,2); % empty jacobian matrix, filled below
    
    % x-dir
    A = [ones(npts,1), xp(:), yp(:)];
    param = A\xw;
    xw0 = param(1);
    Jp2w(1,1) = param(2);
    Jp2w(1,2) = param(3);
    
    % y-dir
    param = A\yw;
    yw0 = param(1);
    Jp2w(2,1) = param(2);
    Jp2w(2,2) = param(3);
    
    
    % compute world-to-pixel transformations (inverse of the above, inverse of a diagnoal matrix is 1/elements)
    Jw2p = inv(Jp2w);
    
    % test transform functions by displaying a downsampled woco image in true coordinates
    
    % create pixel coord grids
    im = double(wocoim(:,:,1));
    [xpi ypi] = meshgrid(1:size(im,2),1:size(im,1));
    
    % downsample so that display is possible!
    d = 5;
    im = downsample(im,d); im = downsample(im',d)';
    xpi = downsample(xpi,d); xpi = downsample(xpi',d)';
    ypi = downsample(ypi,d); ypi = downsample(ypi',d)';
    
    % create woco grids
    n = numel(xpi);
    [nr nc] = size(xpi);
    xvec = reshape(xpi,n,1);
    yvec = reshape(ypi,n,1);
    xW = reshape(Jp2w(1,1).*xvec+Jp2w(1,2).*yvec+xw0, nr, nc);
    yW = reshape(Jp2w(2,1).*xvec+Jp2w(2,2).*yvec+yw0, nr, nc);
    
    % recover pixel grids
    n = numel(xW);
    [nr nc] = size(xW);
    xvec = reshape(xW,n,1);
    yvec = reshape(yW,n,1);
    xpf = reshape(Jw2p(1,1).*(xvec-xw0)+Jw2p(1,2).*(yvec-yw0), nr, nc);
    ypf = reshape(Jw2p(2,1).*(xvec-xw0)+Jw2p(2,2).*(yvec-yw0), nr, nc);
    
    % display image with all 3 grids
    close all
    figure
    p = pcolor(xpi,ypi,im); set(p,'EdgeColor','none');
    colormap('gray')
    title('Initial pixel coords');
    axis equal
    axis tight
    
    figure
    p = pcolor(xW,yW,im); set(p,'EdgeColor','none');
    colormap('gray')
    title('World coords derived from pixel coords');
    axis equal
    axis tight
    
    figure
    p = pcolor(xpf,ypf,im); set(p,'EdgeColor','none');
    colormap('gray')
    title('Pixel coords recovered from world coords');
    axis equal
    axis tight
    
    % Save results
    save('getwoco_output.mat','Jp2w','Jw2p','xw0','yw0');
    
catch err
    fprintf(2,'ERROR, entering debug mode\n');
    keyboard
end





% SUPPORT CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out1,out2,out3] = myginput(arg1)
%GINPUT Graphical input from mouse.
%   [X,Y] = GINPUT(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse.  Data points are entered by pressing
%   a mouse button or any key on the keyboard except carriage return,
%   which terminates the input before N points are entered.
%
%   [X,Y] = GINPUT gathers an unlimited number of points until the
%   return key is pressed.
%
%   [X,Y,BUTTON] = GINPUT(N) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   Examples:
%       [x,y] = ginput;
%
%       [x,y] = ginput(5);
%
%       [x, y, button] = ginput(1);
%
%   See also GTEXT, WAITFORBUTTONPRESS.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 5.32.4.15 $  $Date: 2010/08/16 21:04:38 $

out1 = []; out2 = []; out3 = []; y = [];
c = computer;
if ~strcmp(c(1:2),'PC')
    tp = get(0,'TerminalProtocol');
else
    tp = 'micro';
end

if ~strcmp(tp,'none') && ~strcmp(tp,'x') && ~strcmp(tp,'micro'),
    if nargout == 1,
        if nargin == 1,
            out1 = trmginput(arg1);
        else
            out1 = trmginput;
        end
    elseif nargout == 2 || nargout == 0,
        if nargin == 1,
            [out1,out2] = trmginput(arg1);
        else
            [out1,out2] = trmginput;
        end
        if  nargout == 0
            out1 = [ out1 out2 ];
        end
    elseif nargout == 3,
        if nargin == 1,
            [out1,out2,out3] = trmginput(arg1);
        else
            [out1,out2,out3] = trmginput;
        end
    end
else
    
    fig = gcf;
    figure(gcf);
    
    if nargin == 0
        how_many = -1;
        b = [];
    else
        how_many = arg1;
        b = [];
        if  ischar(how_many) ...
                || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
                || ~(fix(how_many) == how_many) ...
                || how_many < 0
            error('MATLAB:ginput:NeedPositiveInt', 'Requires a positive integer.')
        end
        if how_many == 0
            % If input argument is equal to zero points,
            % give a warning and return empty for the outputs.
            
            warning ('MATLAB:ginput:InputArgumentZero',...
                ['GINPUT(N) requires a positive, integer, scalar input. ',...
                'Returning without collecting any user input.\n',...
                'This will become an error in a future release']);
        end
    end
    
    % Suspend figure functions
    state = uisuspend(fig);
    
    toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
    if ~isempty(toolbar)
        ptButtons = [uigettool(toolbar,'Plottools.PlottoolsOff'), ...
            uigettool(toolbar,'Plottools.PlottoolsOn')];
        ptState = get (ptButtons,'Enable');
        set (ptButtons,'Enable','off');
    end
    
    oldwarnstate = warning('off', 'MATLAB:hg:Figure:Pointer');
    set(fig,'Pointer','crosshair');
    warning(oldwarnstate);
    
    fig_units = get(fig,'Units');
    char = 0;
    
    % We need to pump the event queue on unix
    % before calling WAITFORBUTTONPRESS
    drawnow
    
    while how_many ~= 0
        % Use no-side effect WAITFORBUTTONPRESS
        waserr = 0;
        try
            keydown = wfbp;
        catch %#ok<CTCH>
            waserr = 1;
        end
        if(waserr == 1)
            if(ishghandle(fig))
                set(fig,'Units',fig_units);
                uirestore(state);
                error('MATLAB:ginput:Interrupted', 'Interrupted');
            else
                error('MATLAB:ginput:FigureDeletionPause', 'Interrupted by figure deletion');
            end
        end
        % g467403 - ginput failed to discern clicks/keypresses on the figure it was
        % registered to operate on and any other open figures whose handle
        % visibility were set to off
        figchildren = allchild(0);
        if ~isempty(figchildren)
            ptr_fig = figchildren(1);
        else
            error('MATLAB:ginput:FigureUnavailable','No figure available to process a mouse/key event');
        end
        %         old code -> ptr_fig = get(0,'CurrentFigure'); Fails when the
        %         clicked figure has handlevisibility set to callback
        if(ptr_fig == fig)
            if keydown % OTHER KEY PRESS
                char = get(fig, 'CurrentCharacter');
                button = abs(get(fig, 'CurrentCharacter'));
                scrn_pt = get(0, 'PointerLocation');
                set(fig,'Units','pixels');
                loc = get(fig, 'Position');
                % We need to compensate for an off-by-one error:
                pt = [scrn_pt(1) - loc(1) + 1, scrn_pt(2) - loc(2) + 1];
                set(fig,'CurrentPoint',pt);
                
                % SEGMENT CODE, called if 'n' is pressed
                if strcmp('n',get(fig, 'CurrentCharacter'))
                    out1 = [out1;NaN]; %#ok<AGROW>
                    y = [y;NaN]; %#ok<AGROW>
                    disp('NaN inserted')
                    continue
                end
                % SEGMENT CODE

            else % MOUSE CLICK
                
                button = get(fig, 'SelectionType');
                if strcmp(button,'open')
                    button = 1;
                elseif strcmp(button,'normal')
                    button = 1;
                elseif strcmp(button,'extend')
                    button = 2;
                elseif strcmp(button,'alt')
                    button = 3;
                else
                    error('MATLAB:ginput:InvalidSelection', 'Invalid mouse selection.')
                end
                

                    
                
            end
            axes_handle = gca;
            drawnow;
            pt = get(axes_handle, 'CurrentPoint');
            
            how_many = how_many - 1;
            
            % PLOT POINT ON MOUSE CLICK
            hold on
            plot(pt(1,1),pt(1,2),'*r'); 
            
            if(char == 13) % & how_many ~= 0)
                % if the return key was pressed, char will == 13,
                % and that's our signal to break out of here whether
                % or not we have collected all the requested data
                % points.
                % If this was an early breakout, don't include
                % the <Return> key info in the return arrays.
                % We will no longer count it if it's the last input.
                break;
            end
            
            % INSET NAN FOR N BUTTON PRESS
%             if strcmp(char,'n') 
%                 out1 = [out1;NaN]; %#ok<AGROW>
%                 y = [y;NaN]; %#ok<AGROW>
%                 disp('NaN inserted')
                % END UPDATED CODE
%             else
                out1 = [out1;pt(1,1)]; %#ok<AGROW>
                y = [y;pt(1,2)]; %#ok<AGROW>
                b = [b;button]; %#ok<AGROW>
%             end
        end
    end
    
    uirestore(state);
    if ~isempty(toolbar) && ~isempty(ptButtons)
        set (ptButtons(1),'Enable',ptState{1});
        set (ptButtons(2),'Enable',ptState{2});
    end
    set(fig,'Units',fig_units);
    
%     if ~strcmp(char,'n')
        if nargout > 1
            out2 = y;
            if nargout > 2
                out3 = b;
            end
        else
            out1 = [out1 y];
        end
%     end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = []; %#ok!

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
  h=findall(fig,'type','uimenu','accel','C');   % Disabling ^C for edit menu so the only ^C is for
  set(h,'accel','');                            % interrupting the function.
  
  %%%%%%Keep waiting for button press until user hits 'Enter' or clicks the
  %%%%%%mouse or Cntrl+C. Continue the while loop and zoom in/out if user
  %%%%%%presses r/v
  while(1)
      keydown = waitforbuttonpress;
      current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
      if~isempty(current_char) && (keydown == 1)           % If the character was generated by the
          if(current_char == 3)                       % current keypress AND is ^C, set 'waserr'to 1
              waserr = 1;                             % so that it errors out.
              break;
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %Your keypressfcn functionality here % UPDATE IN v0.2
          panfrac = 0.5;
          if(current_char == 'a') % pan left
              xlim = get(gca,'XLim');
              xrng = abs(diff(xlim));
              shift = panfrac*xrng;
              set(gca,'XLim',[xlim(1)-shift,xlim(2)-shift])
          elseif(current_char == 'd') % pan right 25%
              xlim = get(gca,'XLim');
              xrng = abs(diff(xlim));
              shift = panfrac*xrng;
              set(gca,'XLim',[xlim(1)+shift,xlim(2)+shift])
          elseif(current_char == 'w') % pan up 10%
              ylim = get(gca,'YLim');
              yrng = abs(diff(ylim));
              shift = panfrac*yrng;
              set(gca,'YLim',[ylim(1)+shift,ylim(2)+shift])
          elseif(current_char == 's') % pan down 25%
              ylim = get(gca,'YLim');
              yrng = abs(diff(ylim));
              shift = panfrac*yrng;
              set(gca,'YLim',[ylim(1)-shift,ylim(2)-shift])
          elseif(current_char == 'i') % zoom in x2
              zoom(2);
          elseif(current_char == 'o') % zoom out x2
              zoom(0.5);
          elseif(current_char == 13 || current_char == 'n') % added 'n' to escape from this loop for segment creation
              break;
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
      else
          break;
      end
      
  end
  %%%%%%%%%End of While loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
set(h,'accel','C');                                 % Set back the accelerator for edit menu.
catch %#ok!
  waserr = 1;
end
drawnow;
if(waserr == 1)
   set(h,'accel','C');                                % Set back the accelerator if it errored out.
   error('MATLAB:ginput:Interrupted', 'Interrupted');
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%