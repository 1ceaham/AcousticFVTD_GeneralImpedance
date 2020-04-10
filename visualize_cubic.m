function visualize_cubic(vP, nElDims, usedVec, space, mode)
% visualize_cubic plots anti-aliased time domain pressure data from a
% simulation made with a cubic grid of sample points.
%
%   INPUTS
%   vP: Velocity potential.
%   nElDims: The number of elements along each dimension of the rectangular
%       solid that bounds the problem domain.
%   usedVec: A vector describing which of the elements in the grid
%       corresponding to nElDims actually has data associated with it, as
%       stored (in order) in vP.
%   space: An identifier used to determine default slice locations.
%       Options: 'box', 'theater', 'hallway'.
%   mode: An identifier describing how to transform the data to make the
%       visualization more clear. Options: 'db', 'compand'.
%       

%% Error checking

implementedSpaces = {'box','theater','hallway'};
implementedModes = {'db','compand'};

if ~any(strcmpi(space, implementedSpaces))
    error(['Specified space must be ''' ...
        strjoin(implementedSpaces(1:end-1), ''', ''') ''', or ''' ...
        implementedSpaces{end} '''.'])
end

if ~any(strcmpi(mode, implementedModes))
    error(['Specified mode must be ''' ...
        strjoin(implementedModes(1:end-1), ''', ''') ''', or ''' ...
        implementedModes{end} '''.'])
end

%% Reorganize linear vP into cubic grid

nVisSteps = size(vP,2)-1;
vPrectangular = nan(prod(nElDims), size(vP,2)); % Start with NaNs so that unused elements are invisible
vPrectangular(usedVec,:) = vP;
pMesh = reshape(diff(vPrectangular,[],2),[nElDims nVisSteps]); % Reshape into original grid

%% Anti-aliasing
% Assuming 7x oversampling, we can low-pass to recover the part where the
% dispersion error is not too high

[b,a] = butter(9,2/7);
pMesh = filter(b,a,pMesh,[],4);

%% Determine how to visualize data

switch mode
    case 'db'
        % dB scale, absolute pressure
        pMesh = db(pMesh);
        pMesh(pMesh<-90) = -90;
    case 'compand'
        % Companded version, pos/neg pressure
        % https://www.allaboutcircuits.com/technical-articles/companding-logarithmic-laws-implementation-and-consequences/
        mu = 10000; % Companding factor, picked arbitrarily
        pMesh = sign(pMesh).*log(1+mu*abs(pMesh))./log(1+mu);
end

pMax = max(pMesh(:));
pMin = min(pMesh(:));

%% Determine slices

% Things here are a little weird since it seems like the slice function
% switches X and Y axes?

switch space
    case 'box'
        xSlice = round(nElDims(2)/2);
        ySlice = round(nElDims(1)/2);
        zSlice = round(nElDims(3)/2);
    case 'theater'
        xSlice = round([0.25 0.75]*nElDims(2));
        ySlice = round(nElDims(1)/2);
        zSlice = round(2/3*nElDims(3));
    case 'hallway'
        xSlice = round(nElDims(2)/2);
        ySlice = round(nElDims(1)/16):round(nElDims(1)/8):nElDims(1);
        zSlice = round(nElDims(3)/2);
end

%% Create figure and scrollbar

figure
hmvP = slice(pMesh(:,:,:,1),xSlice,ySlice,zSlice);

% Draggable scrollbar to visualize the time series
hscrollbar1 = uicontrol('style','slider','units','normalized', ...
    'position',[0 0 1 .05],'min',1,'max',nVisSteps,'Value',1);
addlistener(hscrollbar1,'ContinuousValueChange', ...
    @(src, evt) updateSurfaces(src, evt, hmvP, xSlice, ySlice, zSlice, pMesh));

% Other options for limits, default viewing angle
% hl1 = addlistener(hscrollbar1,'ContinuousValueChange',@(src,evt) set(hmvP, 'ZData', vPmesh(:,:,slicingPlane,round(evt.Source.Value))));
% hl2 = addlistener(hscrollbar1,'ContinuousValueChange',@(src,evt) set(hmvP.Parent, 'CLim', max(abs(get(hmvP.Parent, 'CLim'))) * [-1 1]));
% hl2 = addlistener(hscrollbar1,'ContinuousValueChange',@(src,evt) set(hmvP.Parent, 'CLim', [min(vPmesh(:,:,slicingPlane,round(evt.Source.Value)),[],'all') vPmax]));
% hl2 = addlistener(hscrollbar1,'ContinuousValueChange',@(src,evt) set(hmvP.Parent, 'CLim', [vPmin vPmax]));
% view(90,90)

% Make it look good
axis image
set(gca,'PlotBoxAspectRatioMode','auto');
set(hmvP, 'EdgeAlpha', 0.3)
% xlabel('X'); ylabel('Y'); zlabel('Z');

%% Colorbar options depending on mode

colorbar
switch mode
    case 'db'
        set(gca, 'CLim', [pMin pMax])
    case 'compand'
        symMax = max(abs([pMin pMax]));
        set(gca, 'CLim', [-symMax symMax])
        intX = [1 1.5 2.2 3.7 4.5 5];
        colormap(flipud([interp1(intX, [0.5 1 0.9 0.9 0 0], linspace(1,5,64)).' ...
            interp1(intX, [0 0 0.9 0.9 1 0.5], linspace(1,5,64)).' ...
            interp1(intX, [0.5 1 0.9 0.9 0 0], linspace(1,5,64)).']))
end

end

function updateSurfaces(~, evt, hmvP, xSlice, ySlice, zSlice, vPmesh)
    % Update the color data in place so we don't have to replot the patches
    tIdx = round(evt.Source.Value);
    curSurf = 1;
    for xSidx = 1:length(xSlice)
        set(hmvP(curSurf), 'CData', squeeze(vPmesh(:,xSlice(xSidx),:,tIdx)))
        curSurf = curSurf + 1;
    end
    for ySidx = 1:length(ySlice)
        set(hmvP(curSurf), 'CData', squeeze(vPmesh(ySlice(ySidx),:,:,tIdx)))
        curSurf = curSurf + 1;
    end
    for zSidx = 1:length(zSlice)
        set(hmvP(curSurf), 'CData', squeeze(vPmesh(:,:,zSlice(zSidx),tIdx)))
        curSurf = curSurf + 1;
    end
end