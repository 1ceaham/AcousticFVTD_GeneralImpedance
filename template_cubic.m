% function vP = template_cubic(space, boundary, varargin)
% template_cubic sets up and runs an acoustic FVTD simulation using a cubic
% mesh with staircase boundaries.
%
% By default, this file is set up as a script as a demonstration of how to
% prepare data for the finite volume time domain implementation in fvtd.m.
% To use it as a function that returns the velocity potential, simply
% uncomment the first line of the file. Furthermore, while the differing
% approaches for each space occurred in the course of development, they
% have purposefully been left non-homogeneous to demonstrate possible ways
% to parameterize and simulate spaces.
%
%   INPUTS
%   space: One of the default spaces configured in this file. Options:
%       'box', 'theater', and 'hallway'. The first two are based on the
%       examples given in Bilbao et al.
%   boundary: One of the default boundary conditions configured in this
%       file. Options: 'theater', 'resistiveBox', 'uniform', and
%       'curtainCeiling', where the first two correspond to boundary
%       conditions given in Bilbao et al.
%
%   OPTIONAL INPUTS
%   'nofig': Suppresses the visualization of the resulting pressure data.
%

%% Parse inputs

if exist('varargin','var')
    % Running as a function
    if any(strcmpi(varargin, 'nofig')), visualize = 0; else, visualize = 1; end
else
    clear
    close all
    
    visualize = 1;
    space = 'theater';
    boundary = 'theater';
end

%% Error checking

implementedSpaces = {'box','theater','hallway'};
implementedBCs = {'everywhere','curtainCeiling','resistiveBox','theater', ...
    'reflective'};

if ~any(strcmpi(space, implementedSpaces))
    error(['Specified space must be ''' ...
        strjoin(implementedSpaces(1:end-1), ''', ''') ''', or ''' ...
        implementedSpaces{end} '''.'])
end

if ~any(strcmpi(boundary, implementedBCs))
    error(['Specified boundary condition must be ''' ...
        strjoin(implementedBCs(1:end-1), ''', ''') ''', or ''' ...
        implementedBCs{end} '''.'])
end

%% Geometry

% Specify a bounding box for the problem domain to be meshed with cubic
% elements. For simplicity, assume the space is in the first octant.

switch space
    case 'box'
        Lx = 4*sqrt(5);
        Ly = 4*sqrt(3);
        Lz = 4;

        boundingBox = [Lx Ly Lz];
    case 'theater'
        backWidth = 7.6;
        backHeight = 2.7;
        frontWidth = 10;
        frontHeight = 4.5;
        backFront = 10;

        % Call the screen a cylinder centered at the back wall tangent at the front
        % two corners - not exact, but this isn't CAD / meshing software
        screenAngleMax = atan(frontWidth/2 / backFront);
        radius = sqrt((frontWidth/2)^2 + backFront^2);

        % Linearize into 10 elements
        screenAngles = linspace(-screenAngleMax, screenAngleMax, 11).';
        screenWidths = sin(screenAngles)*radius;
        screenDepths = cos(screenAngles)*radius;

        % Now generate the convex hull
        points = [-backWidth/2 0 frontHeight-0.4;
            backWidth/2 0 frontHeight-0.4;
            -backWidth/2 0 frontHeight-backHeight;
            backWidth/2 0 frontHeight-backHeight];
        points = [points; [screenWidths screenDepths repmat(frontHeight,length(screenWidths),1)]];
        points = [points; [screenWidths screenDepths zeros(length(screenWidths),1)]];
        points(:,1) = points(:,1) + frontWidth/2;
        
        % Get the Delaunay triangulation for fast geometry checks
        DT = delaunayTriangulation(points);
        % tetramesh(DT)
        
        % Can check the volume for agreement with paper
        [~,theaterVolume] = convhull(points);
        
        boundingBox = max(points);
    case 'hallway'
        hallWidth = 1.59;
        hallHeight = 2.375;

        alcoveWidth = 0.8;
        alcoveHeight = 2.2;

        isAlcove = [0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
        nAlcoves = sum(isAlcove);

        sectionLengths = fliplr([3 1.265 3.37 2.629 3.37 2.629 3.37 2.629 3.37 2.629 3.37 2.629 3.37 2.629 4.635]);
        orig = [0 cumsum(sectionLengths(1:end-1))];

        boundingBox = [sum(sectionLengths) hallWidth+alcoveWidth hallHeight];
end

%% Set up mesh in bounding box

fs = 2000; % Should be relatively quick; dispersion error ok up to 250 Hz
% 7x oversampling corresponds to cutoff given in Southern et al. 2011
% fs = 700*7; % Ok through 500 Hz octave band, (0.12 - 0.15 * fs)
c = 340;
T = 1/fs;
d = T*c*sqrt(3); % Call it d here to differentiate from h, the matrix version

nElDims = round(boundingBox./d); % How many elements along each primary axis
nElBB = prod(nElDims); % Total number of bounding box elements

xVec = d/2:d:d*nElDims(1);
yVec = d/2:d:d*nElDims(2);
zVec = d/2:d:d*nElDims(3);
[x,y,z] = ndgrid(xVec,yVec,zVec);

%% Connectivity

conn = connectivity_cubic(nElDims);

%% Cull unused cells

switch space
    case 'box'
        % Assume they're all valid
        isInVolumeVec = ones(nElDims);
    case 'theater'
        % Use the Delaunay triangulation
        isInVolumeVec = ~isnan(pointLocation(DT,[x(:) y(:) z(:)]));
    case 'hallway'
        hallVolLengths = [sectionLengths sectionLengths(isAlcove>0)];
        hallVolOrigins = [orig orig(isAlcove>0)];
        volumes = [hallVolOrigins.' ...
            [zeros(length(sectionLengths),1);
            repmat(hallWidth,nAlcoves,1)] ...
            zeros(length(hallVolLengths),1) ...
            (hallVolLengths+hallVolOrigins).' ...
            [repmat(hallWidth,length(sectionLengths),1); repmat(hallWidth+alcoveWidth,nAlcoves,1)] ...
            [repmat(hallHeight,length(sectionLengths),1); repmat(alcoveHeight,nAlcoves,1)]];

        % % Original method
        % isInVolume = zeros(1,nEl);
        % for i = 1:nEl
        %     isInVolume(i) = any(all([x(i)>volumes(:,1) y(i)>volumes(:,2) z(i)>volumes(:,3) ...
        %         x(i)<volumes(:,4) y(i)<volumes(:,5) z(i)<volumes(:,6)],2));
        % end
        
        % Vectorized method
        xIn = xVec > volumes(:,1) & xVec < volumes(:,4);
        yIn = yVec > volumes(:,2) & yVec < volumes(:,5);
        zIn = zVec > volumes(:,3) & zVec < volumes(:,6);

        isInVolumeVec = zeros(nElDims);
        for volIdx = 1:size(volumes,1)
            [mx, my, mz] = ndgrid(xIn(volIdx,:), yIn(volIdx,:), zIn(volIdx,:));
            inVolume = mx & my & mz;
            isInVolumeVec = isInVolumeVec | inVolume;
        end
        
        % all(isInVolumeVec(1:nEl) == isInVolume)
        
        % % Delaunay method, slower but easier to get which volume
        % % Generate points from volumes and find triangulations
        % ptIndices = [1 2 3; 1 2 6; 1 5 3; 1 5 6; 4 2 3; 4 2 6; 4 5 3; 4 5 6];
        % for volIdx = 1:size(volumes,1)
        %     for ptIdx = 1:size(ptIndices,1)
        %         points((volIdx-1)*8 + ptIdx, :) = volumes(volIdx, ptIndices(ptIdx,:));
        %     end
        % 
        %     DT{volIdx} = delaunayTriangulation(points(1+(volIdx-1)*8 : volIdx*8, :));
        % end
        % 
        % figure
        % hold on
        % cellfun(@tetramesh, DT)
        % axis equal
        % 
        % inVolumeDT = cell2mat(cellfun(@(t) pointLocation(t, [x(:) y(:) z(:)]), DT, 'uni', 0));
        % isInVolumeVecDT = any(~isnan(inVolumeDT),2);
end

% TODO: Calculate volume / surface area errors
% dimError = d*(nElDims)-boundingBox; % pct error (d*(nElDims)-boundingBox)./boundingBox

% scatter3(x(isInVolumeVec),y(isInVolumeVec),z(isInVolumeVec),'b.')
% axis equal
% axis vis3d

%% Generate reduced versions

usedVec = find(isInVolumeVec);
unusedVec = find(~isInVolumeVec);

conn(ismember(conn, unusedVec)) = 0; % Set references to eliminated cells to zero

conn = conn(usedVec,:); % Eliminate useless cells

[~, conn] = ismember(conn, usedVec); % Change indices for remaining cells to their order in the reduced set

% Cubic cells without fitting mean h, S, and V are the same everywhere
nEl = length(usedVec);
h = repmat(d, nEl, 6);
S = repmat(d.^2, nEl, 6);
V = repmat(d.^3, nEl, 1);

%% Boundary conditions

% Boundary constants
% Material #1 (Double hung curtain in front of wall)
L1 = [24.4272 23.2179 9.6630e-2 9.4285e-2 9.4588e-2 9.8014e-2 1.0520e-1 1.0114e-1 8.8486e-2];
R1 = [1672.6264 1210.5657 39.3922 35.8868 55.0224 76.1128 54.2331 58.0354 99.5043];
C1 = [8.4126e-7 1.8279e-7 8.4863e-7 9.3606e-8 3.3616e-8 1.6608e-8 9.2427e-9 6.4229e-9 5.3150e-9];

% Material #2 (Porous material on wall; acoustic treatment)
L2 = [26.0297 26.6955 4.1593e-2 3.2549e-2 4.2137e-2];
R2 = [1749.91 1583.26 330.129 457.115 629.982];
C2 = [7.49563e-7 1.57614e-7 3.83857e-7 4.47654e-8 1.15506e-8];

% Material #3 (Porous material on thick panel; seat cushions)
L3 = [22.7431 8.1499e-2 6.768e-2 6.4599e-2 1.6493e-2 7.1131e-3];
R3 = [750.5833 336.6442 491.0421 676.5088 493.5131 60.0617];
C3 = [2.1842e-7 7.3503e-7 7.8999e-8 2.6469e-8 3.4710e-8 2.3215e-8];

switch boundary
    case 'everywhere'
        % Repeat for every side of every node for simplicity
        Ljlm = repmat(permute(L2,[3 1 2]),nEl,6);
        Rjlm = repmat(permute(R2,[3 1 2]),nEl,6);
        Cjlm = repmat(permute(C2,[3 1 2]),nEl,6);
    case 'curtainCeiling'
        % Repeat for every side of every node for simplicity
        % We also have to pad Material #2 so it's the same size as #1
        L2(9) = 0; R2(9) = 0; C2(9) = 0;
        Ljlm = repmat(permute(L2,[3 1 2]),nEl,6);
        Rjlm = repmat(permute(R2,[3 1 2]),nEl,6);
        Cjlm = repmat(permute(C2,[3 1 2]),nEl,6);

        % Replace the "ceiling" side with curtains
        Ljlm(:,6,:) = repmat(permute(L1,[3 1 2]),nEl,1);
        Rjlm(:,6,:) = repmat(permute(R1,[3 1 2]),nEl,1);
        Cjlm(:,6,:) = repmat(permute(C1,[3 1 2]),nEl,1);
    case 'theater'
        % Repeat #2 for every side of every node for simplicity
        % (The paper specifies #2 on the "floor," but that may be a typo
        % since that should be the seat cushions, so we'll switch them.)
        % We also have to pad Material #2 and #3 to match the size of #1
        L2(9) = 0; R2(9) = 0; C2(9) = 0;
        L3(9) = 0; R3(9) = 0; C3(9) = 0;
        
        Ljlm = repmat(permute(L2,[3 1 2]),nEl,6);
        Rjlm = repmat(permute(R2,[3 1 2]),nEl,6);
        Cjlm = repmat(permute(C2,[3 1 2]),nEl,6);

        % Replace the Y+ side with curtains
        Ljlm(:,5,:) = repmat(permute(L1,[3 1 2]),nEl,1);
        Rjlm(:,5,:) = repmat(permute(R1,[3 1 2]),nEl,1);
        Cjlm(:,5,:) = repmat(permute(C1,[3 1 2]),nEl,1);
        
        % Replace the Z- side with cushions
        Ljlm(:,5,:) = repmat(permute(L3,[3 1 2]),nEl,1);
        Rjlm(:,5,:) = repmat(permute(R3,[3 1 2]),nEl,1);
        Cjlm(:,5,:) = repmat(permute(C3,[3 1 2]),nEl,1);
    case 'resistiveBox'
        if ~strcmpi(space, 'box'), error('This BC needs to be rewritten to work with other spaces.'), end
        
        alpha = 0.05;
        rho = 1.21;

        admittance = alpha / (rho*c);

        closeWall = find(x == min(x));
        farWall = find(x == max(x));

        % Initialize
        Ljlm = zeros(nEl,6);
        Rjlm = zeros(nEl,6);
        Cjlm = zeros(nEl,6);

        Rjlm(closeWall,1) = 1/admittance;
        Rjlm(farWall,4) = 1/admittance;

        Cjlm(closeWall,1) = inf;
        Cjlm(farWall,4) = inf;
    case 'reflective'
        Ljlm = zeros(nEl,6);
        Rjlm = zeros(nEl,6);
        Cjlm = zeros(nEl,6);
end

%% Initial conditions

vP = zeros(nEl,fs/2); % Velocity potential, sampled at each interior node

switch space
    case 'box'
        rImp = boundingBox/2;
    case 'theater'
        rImp = boundingBox/2;
    case 'hallway'
        rImp = [1 hallWidth/2 1.5];
end

variance2 = 0.4;
dists = sum(([x(usedVec) y(usedVec) z(usedVec)] - repmat(rImp,nEl,1)).^2,2);
vP(:,1:7) = repmat(1/(variance2*sqrt(2*pi))*exp(-1*dists.^2/(2*variance2)),1,7) ...
    .* repmat(gausswin(7).',length(dists),1);
vP(vP<1e-20) = 0; % Clamp gaussian

% % Other options
% vP(1,1:7) = gausswin(7).';
% vP(:,2) = 1/(variance2*sqrt(2*pi))*exp(-1*dists.^2/(2*variance2));
% vP(:,2) = rand(size(vP(:,1)))*2-1;
% vP(:,1) = vP(:,2);

%% Run simulation

vP = fvtd(T,h,S,V,conn,Ljlm,Rjlm,Cjlm,vP);

%% Visualize results

if visualize
    visualize_cubic(vP, nElDims, usedVec, space, 'compand')
end
