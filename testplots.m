function [] = testplots(x,y,u,v,cval,mask)
% make a few standard plots to test the sandboxpiv code in progress




%% Compute strain

% NaN mask
mask = double(mask);
mask(mask==0) = NaN;
u = u.*mask;
v = v.*mask;

% % add a few rows at the base (TEMP)
% u = [min(min(u))*ones(5,size(u,2)); u];
% v = [zeros(5,size(u,2)); v];
% 
% % smooth
% u = nansmooth2(u,1);
% v = nansmooth2(v,1);

% compute gradients
[L11,L12]=nangrad2(u);
[L21,L22]=nangrad2(v);

% Calculate rate of deformation tensor (D) and decompose it into its principal form
%       D1 and D2 are the principal extension and shortening rates, respectively
%       D2angle contains the principal shortening rate directions in radians
[m n] = size(u);
D1=NaN(m,n);
D2=NaN(m,n);
D2angle = NaN(m,n);
for i=1:m
    for j=1:n
        %... Skip if L contains one or more NaNs
        if ~isnan(L11(i,j)*L12(i,j)*L21(i,j)*L22(i,j))
            D=[L11(i,j) (L12(i,j)+L21(i,j))/2; (L12(i,j)+L21(i,j))/2  L22(i,j)];
            [eigvect,eigval] = eig(D);
            %... Check to see that both principal strain rates are signficant (>eps)
            if abs(diag(eigval)) > eps(diag(eigval))
                if eigval(1,1) > eigval(2,2)
                    D1(i,j)= eigval(1,1);
                    D2(i,j)= eigval(2,2);
                    D2angle(i,j) = acos(eigvect(1,2));
                else
                    D1(i,j)= eigval(2,2);
                    D2(i,j)= eigval(1,1);
                    D2angle(i,j) = acos(eigvect(1,1));
                end
            else
                %... insignificant principal rates, set eigen solution to NaNs
                D1(i,j) = NaN;
                D2(i,j)= NaN;
                D2angle(i,j) = NaN;
            end
        end
    end
end

% Calculate invariants and kinematic numbers
omega = (L12-L21)/2; % vorticity
Dt = sqrt(D1.^2 + D2.^2); % 2nd Invariant = measure of strain rate
Dv = D1 + D2; % Trace = dilatancy


%% Plots

figure
subplot(3,1,1)
M = sqrt(u.^2+v.^2);
p = pcolor(x,y,M); set(p,'EdgeColor','None');
psigma = nanstd(M(:));
pmean = nanmean(M(:));
nstd = 3;
caxis([pmean-nstd*psigma pmean+nstd*psigma])
colorbar
title('Velocity magnitude')
axis equal
axis tight
hold on

dfact = 5;
du = downsample(u,dfact); du = downsample(du',dfact)';
dv = downsample(v,dfact); dv = downsample(dv',dfact)';
dx = downsample(x,dfact); dx = downsample(dx',dfact)';
dy = downsample(y,dfact); dy = downsample(dy',dfact)';
quiver(dx,dy,du,dv,1,'Color','k')

subplot(3,1,2)
p = pcolor(x,y,Dt); set(p,'EdgeColor','None');  
psigma = nanstd(Dt(:));
pmean = nanmean(Dt(:));
nstd = 2;
caxis([pmean-nstd*psigma pmean+nstd*psigma])
colorbar
title('Strain (2nd invariant)')
axis equal
axis tight

subplot(3,1,3)
p = pcolor(x,y,omega);  set(p,'EdgeColor','None');
psigma = nanstd(omega(:));
pmean = nanmean(omega(:));
nstd = 2;
caxis([pmean-nstd*psigma pmean+nstd*psigma])
colorbar
title('Vorticity')
axis equal
axis tight

