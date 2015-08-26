function [displacement,spin,Dv,Dd,D2x,D2y,WkStar,AkStar] = yalebox_decompose_step(x,y,uX,uY,iUXY)
%
% Compute derived quantities (deformation parameters) from PIV data.
% Required that the S-point is at x,y = 0, that base of wedge is located in
% the last row of the matrices, and that units are in meters.
%
% x,y = Double 2D matrices. World coordinates of displacement vectors (piv window centers), in m from s-point (x) or table (y).
% uX,uY = Double 2D matrices. Displacement vector components in m.
% iUXY = Logical 2D matrix.  Logical mask identifying where data was computed using piv (1) and where this is impossible (0).
% 
% displacement = 
% spin = 
% Dv = 
% Dd = 
% D2x = 
% D2y = 
% WkStar = 
% AkStar = 
%
% Mark Brandon, Keith Ma, Sept 2012

% v1.1: Minor modifications, now used in routine pic post-processing.
%
% v1.0:  Separated from Mark's decompose.m program.  Replaced the padding
% routine with a new one based on extrapolation with inpaint_nans.

%% Initialize...

% choose methods
padSize = 4; % set pad size for the number of rows added to the top and bottom of the grid. The padding is currently set to four rows, to account for the size of the derivative filters and the Gaussian smoothing filter.
useMatlabGradient = 0; % Choose which algorithm to use for the gradient calculation. 1 = MATLABs gradient.m, 0 = 3rd-party option, currently derivitivesByFilters.m
useGaussianSmoothing = 0; % Turn gaussian smoothing filter on (1) or off (0).  Acts only on derived quantities (NOT uX and uY).

% define constants
rad = pi/180;

% define a few vars
[m,n] = size(uX); % record starting size of grids
dX = x(1,2)-x(1,1); % calculate grid spacing. Note that negative values are allowed to account for the transition from matrix coordinates to Cartesian coordinates.
dY = y(2,1)-y(1,1);

% CHANGE: keep units in m, change only for plotting
% CHANGE: don't normalized velocities, do this for plotting only

%% Treat data boundaries...

% set regions with no data to 0
uX(~iUXY) = 0;
uY(~iUXY) = 0;

% zero-pad data out to the required distance
uX = padarray(uX, [padSize, padSize]);
uY = padarray(uY, [padSize, padSize]);
%uX = zeroPadArray(uX,padSize,padSize,padSize,padSize);
%uY = zeroPadArray(uY,padSize,padSize,padSize,padSize);

% replace 0's with NaNs
uX(uX==0) = NaN;
uY(uY==0) = NaN;

% use snazzy inpaint_nans.m software to extrapolate from data.  This method preserves gradients at the margins!
uX = inpaint_nans(uX,1);
uY = inpaint_nans(uY,1);

[mPad,nPad]=size(uX); % padded matrix dimensions

% WHY CHANGE PADDING ROUTINES:  3 problems with Mark's padding, 1) repeating
% displacement values from the upper row does NOT preserve gradients, 2)
% there can be gradients in uY at the boundaries, so uY must be padded
% properly as well, and 3) the sides must be padded too, since these data
% are used to compute errors as well as to produce plots.


%% Calculate deformation parameters

%... Calculate deformation-gradient tensor F = E + 1,
% where E is the displacement-gradient tensor.  Two options for the
% gradient calcuation.

if useMatlabGradient
    [F11,F12] = gradient(uX,dX,dY);
    [F21,F22] = gradient(uY,dX,dY); 
else
    [F11,F12] = derivativesByFilters(uX,'x','y',dX,dY,'seven');
    [F21,F22] = derivativesByFilters(uY,'x','y',dX,dY,'seven');
end

F11 = F11 + 1;
F22 = F22 + 1;

%... Calculate the spin tensor W, and the stretch tensor D in its principal form
% D1 and D2 are the principal extension and shortening rates, respectively
% theta contains the principal shortening rate directions in radians
D1 = nan(mPad,nPad);
D2 = nan(mPad,nPad);
D2x = nan(mPad,nPad);
D2y = nan(mPad,nPad);
spin = nan(mPad,nPad);
for i = 1:mPad
    for j = 1:nPad
        %... Skip if F contains one or more nans
        if ~isnan(F11(i,j)*F12(i,j)*F21(i,j)*F22(i,j))
            F = [F11(i,j), F12(i,j); F21(i,j),  F22(i,j)];
            %... Polar decomposition, F = VR, B = V^2 = FF', and R = (V^-1)F,
            % where B is the Left Cauchy-Green tensor.
            B = F*F';
            %... Calculate eigen solution and sort with eigenvalues in descending order
            [T,lambda] = eig(B);
            [lambda order] = sort(diag(lambda),'descend');
            S = sqrt(lambda);
            D1(i,j) = log(S(1));
            D2(i,j) = log(S(2));
            T = T(:,order);
            %... Calculate unit vector for maximum shortening rate direction
            theta = acos(T(1,2));
            D2x(i,j) = cos(theta);
            D2y(i,j) = sin(theta);
            %... Calculate spin
            invV = T*diag(lambda.^(-1/2))*T';
            R = invV*F;
            spin(i,j) = atan2(R(2,1),R(1,1));
        end
    end
end
clear F11 F12 F21 F22

%... Calculate invariants, assuming plane strain (D3=0)
% Dt = sqrt(D1.^2 + D2.^2);
Dd = sqrt(((D1-D2).^2 + D1.^2 + D2.^2)/3);
Dv = sqrt(1/3)*(D1 + D2);

%... Smooth displacement field using Gaussian kernel with sigma = 1.2
% pixels, as recommended on p. 512 in Adrian and Westerweel, 2011, Particle
% Image Velocimetry. The filter in imgaussian defaults to a size equal to
% +/-3sigma. Must be enabled in the "choose methods" section.

if useGaussianSmoothing
    spin = imgaussian(spin,1.2,3);
    Dv = imgaussian(Dv,1.2,3);
    Dd = imgaussian(Dd,1.2,3);
    % Dt = imgaussian(Dt,1.2,3);
    D2x = imgaussian(D2x,1.2,3);
    D2y = imgaussian(D2y,1.2,3);
end

%... Remove padding
uX = removePad(uX,padSize); % function is simple and included below, separated for visual simplicity.
uY = removePad(uY,padSize);
spin = removePad(spin,padSize);
Dv = removePad(Dv,padSize);
Dd = removePad(Dd,padSize);
% Dt = removePad(Dt,padSize);
D2x = removePad(D2x,padSize);
D2y = removePad(D2y,padSize);

%.. restore mask
uX(~iUXY) = NaN;
uY(~iUXY) = NaN;
spin(~iUXY) = NaN;
Dv(~iUXY) = NaN;
Dd(~iUXY) = NaN;
% Dt(~iUXY) = NaN;
D2x(~iUXY) = NaN;
D2y(~iUXY) = NaN;

%... Calculate displacement magnitude
displacement = sqrt(uX.^2 + uY.^2);

%... Calculate kinematic numbers
% Wk = 2*spin./(sqrt(2)*Dt);
% Ak = Dv./(sqrt(2)*Dt);
WkStar = 2*spin./(sqrt(2)*Dd);
AkStar = Dv./(sqrt(2)*Dd);



function trimmed = removePad(A,pad)
% reused code to trim off padding

trimmed = A(pad+1:end-pad,pad+1:end-pad);


