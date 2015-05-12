function [rasterX rasterY] = createRaster(S, barheight, fS, colr, varargin)
% createRaster(S, barheight, binSize, offset) creates a 
% raster plot using the spike sequences from S 
% S         is expected to be a cell array of sequences of 1s and 0s.
% barheight dictates how long the bars will be
% fS        Sampling frequency (ephys). 
% offset    sets the limits on the minimum x-axis value and is
%           the time before stimulus presentation that the spike
%           train recording begins from
%
% Sample of how this function might be called:
% >> createRaster(timestamps, 0.4, 3, 200)
%

% Vivek Jayaraman, 4th January 2002
%

% Validate input args
error(nargchk(3,6,nargin));

if (nargin < 6)
    offset = 0;
    xmin = 0;
else
    offset = varargin{2};
    xmin = -offset;
end

binSize = 1/fS;

% First, figure out how much padding's needed for various spike
% trains. Take the longest one, pad the others with zeros.
N = size(S,2);

rasterX = [];
rasterY = [];

y0 = N; 

scaledOffset = round(offset/binSize);
% Simulation loop
%for (j=N:-1:1)        
for (j = 1:N)
    % Set up raster plot
    tndx = find(S(:,j));
    if ~isempty(tndx),
        % Make every third guy NaN so MATLAB won't plot it :-)
        x = zeros(length(tndx)*3,1)+NaN; 
        y = x;
        n = length(x);
        x(1:3:n) = tndx-scaledOffset;
        x(2:3:n) = tndx-scaledOffset;
        y(1:3:n) = y0;
        y(2:3:n) = y0-barheight;
        
        rasterX = [rasterX;x];
        rasterY = [rasterY;y];  
    end;
    y0 = y0 - 1;
end

% figure;
%h = patch([1000 1000 2000 2000], [0 110 110 0], [0.9 0.9 0.9]);
%set(h, 'edgecolor', 'none');
%hold on;
% Use cool trick above now to get raster plot! :-)
plot(rasterX*binSize,rasterY, colr);
%lh = line([0 0], [0 N]);
%set(lh, 'color', 'g');
% xlabel('time(msec)');
% ylabel('spike trains for different cells');
%axis([xmin max(rasterX*binSize) 0 max(rasterY)]);
% axis([0 1000 0 110]);
