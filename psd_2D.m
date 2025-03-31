% Calculates radially averaged 2D power spectrum. 
% It is better to remove mean and tilt in your topography before applying
% this function
% Here, same pixelwidth in x and y direction is assumed which is typical
% of measurement instruments.
% ===
% inputs:
% z : height topography, a matrix of n*m size (SI units, i.e. meters)
% PixelWidth: size of each Pixel in topography/image (SI units, i.e.
% meters). If you don't know your pixelwidth just devide topography length by
% number of pixels in length.
% ===
% outputs:
% q: wavevectors, which is 2pi/lambda. lambda is wavelength of your
% roughness components.
% C: Radially averaged power spectrum (2D PSD)
% PSD: this structure contain any other interesting information you may
% need.
% in order to plot final result just use:
% loglog(q,C)
% 2D FFT of surface topography could be seen by:
% imagesc(1+log10(abs(PSD.Hm)))

% The original version can be found in the following link: 
% https://www.mathworks.com/matlabcentral/fileexchange/54297-radially-averaged-surface-roughness-topography-power-spectrum-psd
% =========================================================================
% make the size an even number (For PSD)


PixelWidth = 0.5;


%% read the tif file for the Wonderland of Rocks and its joint sets 
% z = imread('small_portion_PSD.tif'); % first figure 
% z = imread('mountains_part_WoR.tif');
% z = imread('random_part_WoR.tif'); % third figure
% z = imread('another2_random_part_WoR.tif'); 
% z = imread('another_small_portion.tif'); 
% z = imread('another_random_part_WoR.tif'); 
% z = imread('mountains2_part_WoR.tif');
z = imread('joint_sets2_WoR.tif'); 
% z = imread('WoR_Joint_Set.tif'); 
% z = imread('distinctRotatedCross.tif'); 
% z = imread('emptyLand.tif'); 

%% read the tif file for the Canyonlands and its grabens 

% z = imread('CanyonLandDEM1.tif'); 
% z = imread('CroppedCanyonLandDEM1.tif'); % second figure
% z = imread('CanyonLandDEM2.tif'); 
% z = imread('CroppedCanyonLandDEM2.tif');
% z = imread('CanyonLandDEM3.tif'); 
% z = imread('CroppedCanyonLandDEM3.tif'); 

% this is for generating two diffraction grates in each corner
% z = diffraction_2D_Ver2; 

% z = tifDEMConverter("viz.USGS10m_hillshade.tif");

% z = tifDEMConverter("grabensCanyonLand1.tif");
% z = tifDEMConverter("grabensCanyonLand2.tif");
% z = single(z); 


z = elevationRemover(z, -9999); % ONLY do this when you want
                                % to remove certain elevation

% z = subMatrix; %REMEMBER TO RUN THE rotatingDEM SCRIPT BEFOREHAND! 

% z = circularWindow(z); % applying the circular padding
% z = rectangularWindow(z); % applying the rectangular padding 


% x = size(z,1); 
% y = size(z,2); 
% 
% for i = 1:x    
%     for j = 1:y  
%         if z(i,j) > 1355
%             z(i,j) = 1355; 
%         elseif z(i,j) < 1350
%             z(i,j) = 1350; 
%         end 
%     end
% end

figure(1) 
imagesc(z)
colormap gray
clim([round(min(z(z > 0))) round(max(z, [], "all"))])
daspect([1 1 1]);

[n,m] = size(z);
if mod(n,2)
    z = z(2:end,:);
    n = n -1;
end
if mod(m,2)
    z = z(:,2:end);
    m = m -1;
end
Ln = n * PixelWidth; % width of image
Lm = m * PixelWidth; % length of image
a = PixelWidth; % lattice spacing in meter



% =========================================================================
% Window function (up to user)
% win = tukeywin(n,0.25)*tukeywin(m,0.25)'; % this made a stronger signal
% win = hann(n)*hann(m)';
% win = kaiser(n,10)*kaiser(m,10)';
% win = bartlett(n) * bartlett(m)';
win = (1 - (((0:n-1)-((n-1)/2))/((n+1)/2)).^2)'*(1 - (((0:m-1)-((m-1)/2))/((m+1)/2)).^2); % Welch window
recInt = trapz(0:n-1,trapz(0:m-1,((rectwin(n)* rectwin(m)').^2),2),1); % integral of squared rectangular window
winInt = trapz(0:n-1,trapz(0:m-1,(win.^2),2),1); % integral of square of selected window function

% U = (winInt/recInt); % Normalization constant
% z_win = z.*win; % Height applied with the window(s)

% U = 1/(2*pi); 

U = 1/sqrt(2*pi); % BASIC Normalization constant
z_win = z; 


% =========================================================================
% Calculate 2D PSD
Hm = fftshift(fft2(z_win,n,m));
Cq = (1/U)*(a^2/((n*m)*((2*pi)^2)).*((abs((Hm))).^2));
Cq(n/2+1,m/2+1) = 0; % remove the mean 


% =========================================================================
% corresponding wavevectors to Cq values after fftshift has been applied

qx_1=zeros(m,1);
for k=0:m-1
    qx_1(k+1)=(2*pi/m)*(k);
end
qx_2 = fftshift(qx_1);
qx_3 = unwrap(qx_2-2*pi);
qx = qx_3/a;

qy_1=zeros(n,1);
for k=0:n-1
    qy_1(k+1)=(2*pi/n)*(k);
end
qy_2 = fftshift(qy_1);
qy_3 = unwrap(qy_2-2*pi);

qy = qy_3/a;



% =========================================================================
% Radial Averaging

[qxx,qyy] = meshgrid(qx,qy);
[~,rho] = cart2pol(qxx,qyy);
% rho = floor(rho);
J = 100; % resolution in q space (increase if you want)
qrmin = log10(sqrt((((2*pi)/Lm)^2+((2*pi)/Ln)^2)));
qrmax = log10(sqrt(qx(end).^2 + qy(end).^2)); % Nyquist
% lambda = (2*pi)./(10.^linspace(qrmin,qrmax,J)); %wavelengths 
% q = floor(10.^linspace(qrmin,qrmax,J));
q = 10.^linspace(qrmin,qrmax,J);



% =========================================================================
% Averaging Cq values
C_AVE = zeros(1,length(q));
ind = cell(length(q)-1 , 1);
for j = 1:length(q)-1
    ind{j} = find( rho > q(j) & rho <=(q(j+1)));
    C_AVE(j) = mean(Cq(ind{j}));
end
ind = ~isnan(C_AVE);
C = C_AVE(ind);
q = q(ind);
lambda = (2*pi)./q; 



% =========================================================================
PSD.Hm = Hm; % 2D FFT
PSD.C = C;
PSD.q = q;
PSD.Cq = Cq;
PSD.qx = qx;
PSD.qy = qy;
PSD.z_win = z_win; % z profile after window function applied



% =========================================================================
% This displays the 2D PSD computed directly from FFT 
figure(2)
loglog(q,C);
xlabel('Wavevector: log q (m^{-1})')
ylabel('Power Spectrum: log C (m^4)')
title('Radially Averaged 2D PSD')



% =========================================================================
% This displays the window that represents artificial symmetry of the
% wavevectors in x and y direction. Type of window can be adjusted in the
% code
TwoDimPSD = 1+log10(abs(PSD.Hm));
% TwoDIMPSD = smoothdata2(TwoDimPSD, "gaussian", 20);

% this controls how much noise you want the 2D PSD to have the higher the 
% threshold, less pixel intensities will show up.

% threshold = 4.9; % use for joint sets 
% threshold = 5.7; % use for grabens 
% threshold = 5.1; % use for grabens (shows more wavevectors)
% threshold = 5.4; 
% threshold = 0;
% threshold = 5.5;
threshold = 6; % use for padding
% threshold = 6.3;

TwoDimPSD(TwoDimPSD < threshold) = 0; % removes the pixel intensity ratio

figure(3)
imagesc(TwoDimPSD)
xlabel('Wavevector, q_{x} (m^{-1})')
ylabel('Wavevector, q_{y} (m^{-1})')
colormap turbo
colorbar
% clim([0 9]);
daspect([1 1 1])

figure(4) 
surf(TwoDimPSD)
shading interp



% =========================================================================
% This displays the local maxima of the PSD, which are the most frequent 
% wavelengths in the PSD 
figure(5)
[pks,locs] = findpeaks(log(C));
newQ = log(q);
plot(log(q),log(C),newQ(locs),pks,"o")



% =========================================================================
% This outputs the frequent wavelengths in a separate table 

frequentWavelength = zeros(1, size(locs, 2));
for i = 1:size(locs, 2)
    frequentWavelength(i) = (2 * pi) ./ q(locs(i));
end

waveVectorandLength = zeros(2, size(q, 2));
waveVectorandLength(1, :) = q;
waveVectorandLength(2, :) = lambda;

% Create a table with vertical data
results = table(q(:), lambda(:), 'VariableNames', {'Wavevector', 'Wavelength'});



% =========================================================================
% Table for the local maxima (peaks) with corresponding wavelengths

peakWavelengths = (2 * pi) ./ q(locs);
maximaTable = table(q(locs)', pks', peakWavelengths', 'VariableNames', {'Wavevector at Maximum', 'Max PSD Value', 'Wavelength at Maximum'});



% % =========================================================================
% % This section outputs a cross section of the 2D PSD 
% smoothHorizontalCrossSection = TwoDimPSD(541,:); 
% smoothVerticalCrossSection = TwoDimPSD(:,541); 
% figure(5)
% subplot(2,2,1)
% plot(smoothHorizontalCrossSection);
% title('Horizontal Cross Section')
% subplot(2,2,2)
% plot(smoothVerticalCrossSection); 
% title('Vertical Cross Section') 
% subplot(2,2,3)
% inverse1 = ifft(smoothHorizontalCrossSection);
% plot(inverse1);
% subplot(2,2,4) 
% inverse2 = ifft(smoothVerticalCrossSection);
% plot(inverse2);
% 
% % figure(6)
% % subplot(1,2,1)
% % plot(fft(inverse1));
% % subplot(1,2,2)
% % plot(fft(inverse2));
% 
% 
% 
% % =========================================================================
% % Find the PSD of a grid that contains one pixel 
% onePxGrid = zeros(size(z), 'single');
% onePxGrid((size(z,1)/2),(size(z,2)/2)) = 1; 
% 
% halfSquareLength = 199;
% onePxGrid((size(z,1)/2)-halfSquareLength:(size(z,1)/2)+halfSquareLength, (size(z,2)/2)-halfSquareLength:(size(z,2)/2)+halfSquareLength) = 2;
% 
% figure(7)
% imagesc(onePxGrid)
% colormap turbo
% colorbar
% % clim([0 3]);
% daspect([1 1 1])
% 
% % onePxPSD = abs(fftshift(fft2(onePxGrid)));
% onePxPSD = psd_2D_function(onePxGrid); 
% 
% figure(8)
% imagesc(onePxPSD)
% xlabel('Wavevector, q_{x} (m^{-1})')
% ylabel('Wavevector, q_{y} (m^{-1})')
% colormap turbo
% colorbar
% % clim([0 3]);
% daspect([1 1 1])
% 
% onePxCross = zeros(size(onePxPSD)); 
% onePxCross(size(onePxCross)/2,:) = onePxPSD(size(onePxPSD)/2,:); 
% onePxCross(:,size(onePxCross)/2) = onePxPSD(:,size(onePxPSD)/2); 
% magnitude = (max(TwoDimPSD,[],"all")/max(onePxPSD,[],"all"));
% onePxCross = onePxCross .* magnitude; 
% 
% 
% 
% 
% 
% % =========================================================================
% % Smooth the TwoDimPSD with some smooth function
% % smoothTwoDimPSD = smoothdata2(TwoDimPSD); 
% % smoothTwoDimPSD = smoothTwoDimPSD .* (9.844/6.95); 
% smoothTwoDimPSD = smoothdata2(TwoDimPSD,"movmean", 9);
% % try reducing the window to 9, so smoothdata2(TwoDimPSD,movmedian,9)
% 
% % "movmean" — Average over each 2-D window of A. This method is useful for 
% %  reducing periodic trends in data. THIS IS DEFAULT FOR THE FUNCTION 
% % smoothdata2
% 
% % "movmedian" — Median over each 2-D window of A. This method is useful for
% % reducing periodic trends in data when outliers are present.
% 
% % "gaussian" — Gaussian-weighted average over each 2-D window of A.
% 
% % "lowess" — Linear regression over each 2-D window of A. This method can 
% % be computationally expensive but results in fewer discontinuities.
% 
% % "loess" — Quadratic regression over each 2-D window of A. This method is 
% % slightly more computationally expensive than "lowess".
% 
% % "sgolay" — Savitzky-Golay filter, which smooths according to a quadratic 
% % polynomial that is fitted over each 2-D window of A. This method can be
% % more effective than other methods when the data varies rapidly.
% 
% figure(9)
% imagesc(smoothTwoDimPSD)
% colorbar
% clim([0 9]);
% daspect([1 1 1])
% 
% % Take the two cross sections
% % smoothHorizontalCrossSection = smoothTwoDimPSD(541,:); 
% % smoothVerticalCrossSection = smoothTwoDimPSD(:,541); 
% % figure(10)
% % subplot(1,2,1)
% % plot(smoothHorizontalCrossSection);
% % title('Horizontal Cross Section')
% % subplot(1,2,2)
% % plot(smoothVerticalCrossSection); 
% % title('Vertical Cross Section') 
% 
% 
% 
% % =========================================================================
% % Apply basic arithmetic principles 
% % smoothFinalResult = TwoDimPSD - smoothTwoDimPSD; 
% smoothFinalResult = TwoDimPSD - onePxCross; 
% 
% threshold2 = 0;
% smoothFinalResult(smoothFinalResult <= threshold2) = NaN;
% 
% figure(11)
% imagesc(smoothFinalResult)
% colorbar
% % clim([0 9]);
% daspect([1 1 1])
% 
% 
% smoothHorizontalCrossSection = smoothFinalResult(541,:); 
% smoothVerticalCrossSection = smoothFinalResult(:,541); 
% figure(12)
% subplot(1,2,1)
% plot(smoothHorizontalCrossSection);
% title('Horizontal Cross Section')
% subplot(1,2,2)
% plot(smoothVerticalCrossSection); 
% title('Vertical Cross Section') 
% 
% 
% 
% % =========================================================================
% % Graph the PSD in a 3D manner
% figure(13) 
% surf(smoothFinalResult) 
% shading interp % used for removing the gridlines on the graph
% 
% % =========================================================================
% % Using lowpass function to filter the frequency 
% % 
% % lowPassPSD = lowpass(TwoDimPSD,0.01);
% % 
% % figure(14) 
% % imagesc(lowPassPSD)
% % xlabel('Wavevector, q_{x} (m^{-1})')
% % ylabel('Wavevector, q_{y} (m^{-1})')
% % colormap turbo
% % colorbar
% % clim([0 9]);
% % daspect([1 1 1])