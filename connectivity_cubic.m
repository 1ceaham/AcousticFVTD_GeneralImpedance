function conn = connectivity_cubic(nElDims)
% connectivity_cubic assembles a connectivity matrix for a cubic grid.
%
%   INPUTS
%   nElDims: The number of cubic elements along each axis of the grid.
%
%   OUTPUTS
%   conn: The connectivity matrix for each element. The columns of the
%       matrix correspond to neighbors in the following directions, in
%       order: -X, -Y, -Z, +X, +Y, +Z.
%

nElBB = prod(nElDims); % Total number of elements in the bounding box

%% Original method (slower but easy to understand)

% elIdx = reshape(1:nElBB,nElDims);
% connLoop = zeros(nElBB,6);
% for i = 1:nElBB
%     [ix,iy,iz] = ind2sub(size(elIdx),i);
%     
%     if (ix-1 > 0), connLoop(i,1) = elIdx(ix-1,iy,iz); end
%     if (iy-1 > 0), connLoop(i,2) = elIdx(ix,iy-1,iz); end
%     if (iz-1 > 0), connLoop(i,3) = elIdx(ix,iy,iz-1); end
%     if (ix+1 <= nElDims(1)), connLoop(i,4) = elIdx(ix+1,iy,iz); end
%     if (iy+1 <= nElDims(2)), connLoop(i,5) = elIdx(ix,iy+1,iz); end
%     if (iz+1 <= nElDims(3)), connLoop(i,6) = elIdx(ix,iy,iz+1); end
% end

%% Analytic method (faster but more complicated)

conn(:,1) = 0:nElBB-1;
conn(nElDims(1)+1:nElDims(1):end,1) = 0;

conn(:,2) = (1:nElBB) - nElDims(1);
conn(mod(0:nElBB-1, nElDims(1)*nElDims(2)) < nElDims(1),2) = 0;

conn(:,3) = (1:nElBB) - nElDims(1)*nElDims(2);
conn(1:nElDims(1)*nElDims(2),3) = 0;

conn(:,4) = 2:nElBB+1;
conn(nElDims(1):nElDims(1):end,4) = 0;

conn(:,5) = (1:nElBB) + nElDims(1);
conn(mod((0:nElBB-1) + nElDims(1), nElDims(1)*nElDims(2)) < nElDims(1),5) = 0;

conn(:,6) = (1:nElBB) + nElDims(1)*nElDims(2);
conn(end-nElDims(1)*nElDims(2)+1:end,6) = 0;

%% Check to make sure the two are equivalent

% all(conn == connLoop, 'all')
