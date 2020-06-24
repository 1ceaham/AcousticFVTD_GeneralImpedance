function vP = fvtd(T,h,S,V,conn,Ljlm,Rjlm,Cjlm,vP,varargin)
% FVTD Compute acoustic finite volume time domain simulation.
%
%   FVTD is the computational core of running an acoustic finite volume
%   time domain (FVTD) simulation, an implementation of the method given in
%   Bilbao et al. 2016 [1]. Given a mesh and corresponding materials
%   properties, this function computes the acoustic velocity potential
%   resulting from an initial perturbation. If requested, it can also
%   verify the numerical conservation of energy in terms of the interior,
%   boundaries, and dissipation.
%
%   INPUTS
%   T: Sampling period, the inverse of the sample rate.
%   h: Inter-cell distances. Given a collection of N cells where the
%       greatest number of faces on a cell is maxFaces, this is a matrix of
%       size N-by-maxFaces. For the cubic case, this would be N-by-6.
%   S: Inter-cell surface areas. Same size as h.
%   V: Cell volumes. Matrix of size N-by-1.
%   conn: For each cell, the indices of the cells it shares a face with.
%       Same size as h, where 0 implies a connection to a boundary. If a
%       row is all zeros, the cell is not connected to anything, but will
%       be computed nonetheless.
%   Ljlm: Inductive branch for each face of each cell. Given a maximum
%       branch length M, this is a matrix of size N-by-maxFaces-by-M.
%       Impedances are only defined at boundaries, but for convenience,
%       these are extracted in the function (according to conn) from these
%       larger matrices to make them easier to specify.
%   Rjlm: Resistive branch for each face of each cell. See Ljlm.
%   Cjlm: Capacitative branch for each face of each cell. See Ljlm.
%   vP: The initial conditions for the velocity potential (Psi). Also 
%       determines the number of temporal steps to simulate by the number
%       of samples preallocated. Given a desired length of nSteps, this is
%       a matrix of size N-by-nSteps. T*nSteps gives the simulation length
%       in seconds. Computation begins with the 3rd sample, vP(:,3).
%
%   OUTPUTS
%   vP: Computed response to the initial velocity potential (Psi).
%
%   OPTIONAL INPUTS
%   'verifyEnergy': Compute and plot the stored and dissipated energy to
%       ensure that the energy balance is preserved to machine precision.
%   'noProgress': Disable the display of progress and the total time taken.
%
%   REFERENCES
%   [1] S. Bilbao, B. Hamilton, J. Botts, and L. Savioja, “Finite Volume
%       Time Domain Room Acoustics Simulation under General Impedance
%       Boundary Conditions,” IEEE/ACM Transactions on Audio, Speech, and
%       Language Processing, vol. 24, no. 1, pp. 161–173, Jan. 2016, doi:
%       10.1109/TASLP.2015.2500018.
%

%% Parse varargin

if any(strcmpi(varargin, 'verifyEnergy')), verifyEnergy = 1; else, verifyEnergy = 0; end
if any(strcmpi(varargin, 'noProgress')), displayProgress = 0; else, displayProgress = 1; end

%% Constants

c = 340;
rho = 1.21;

nEl = size(conn,1);
nFaces = size(conn,2);
M = size(Ljlm, 3);

%% Stability

quotient = h; % If there are any zero entries, we need to ignore them
quotient(logical(h)) = c^2*T^2*S(logical(h))./(2*h(logical(h)));
criteria = sum(quotient, 2)./V;

if any(criteria > 1), error('Stability criteria not satisfied.'); end

%% Set up simulation

% The indices j, l, and m refer to interior nodes, faces of each node, and
% RLC branches respectively.

etajlm = 1./(Ljlm/T + Rjlm/2 + T./(4*Cjlm)); % (eq. 40)

% Boundary connection masks
gamma = conn == 0; 
gMask = repmat(gamma,1,1,M)>0;

% Now pull out the linear indices for each boundary
Llm = Ljlm(gMask); 
Rlm = Rjlm(gMask);
Clm = Cjlm(gMask);
etalm = etajlm(gMask);

% Set zero capacitance entries to inf
Clm(Clm == 0) = Inf;

clearvars('Ljlm','Rjlm','Cjlm') % Only necessary to get actual boundaries

%% Generate boundary-indexed S

Sjlm = repmat(S, 1, 1, M);
Slm = Sjlm(gMask);

%% Storage

nSteps = size(vP,2);
% TODO: Introduce spatially downsampled vP for saving, and only have vPcur and vPnext for processing
vlb = zeros(size(Llm,1),3); % v_l bar, average velocity at each boundary for each RLC branch; 3 time samples
gl = vlb; % g_l, the internal variable at each boundary for each RLC branch; 3 time samples
% TODO: We might only need two samples for gl?

if verifyEnergy
    hi = zeros(1,nSteps-1); % Energy conservation variables
    hb = hi; q = hi;
    vjk = zeros([size(conn) nSteps]); % Inter-node velocities
    p = zeros(size(vP)); % Pressures
end

bValue = zeros(size(gMask)); % Intermediate computation for vP update

%% Summation indexing generation

% In order to implement the FVTD update scheme in a vectorized fashion, we
% need two operations: subtracting each node's velocity potential from its
% neighbors, and finding the velocity potential corresponding to each
% boundary. Thus, we use indexing to form the desired computation with
% elementwise matrix operations and then sum across the proper dimensions
% in order to recover the desired value.

interRowMask = repmat((1:nEl).',1,nFaces); % Each node referring to itself
interSubtractionIdx = conn + interRowMask.*gamma; % Replace empty sides with reference to self

% Thus, the one row of the sum of a node's neighbors' velocity potential
% minus its own looks something like:
%   sum(vP([1 3 5 1 7 9]) - vP([1 1 1 1 1 1])) 
%   = (vP(3)-vP(1)) + (vP(5)-vP(1)) + (vP(7)-vP(1)) + (vP(9)-vP(1))
% (In this case, the specified node has 2 sides that are boundaries)

interBoundConn = repmat(interRowMask.*gamma,1,1,M); % Sides corresponding to boundaries, repeated for each RLC branch
interBoundIdx = interBoundConn(interBoundConn>0); % List of nodes corresponding to each boundary

% In general, boundary-related variables (like vlb or or any of the *lm
% constants) are in order by the "side" of the cube, meaning all boundaries
% on side 1 first, then side 2, etc., and then repeated for the number of
% branches in the generalized impedance, but should all be indexed in the
% same way such that when it is re-indexed with gMask it should go back
% exactly into the spots to be manipulated in a similar manner as interior
% node related variables.

%% Setup for progress

if displayProgress
    curPctIdx = 1;
    pctPts = (1:99)./100;
    fprintf('Computation progress: 00%%\n')
    tic
end

%% Vectorized implementation

% TODO: Finish getting rid of jlm variables by indexing the result of
% Slm.*etalm with gMask
denominator = (1+(rho*c^2*T./(4*V)).*sum(sum( Sjlm.*etajlm.*gMask ,3),2));

for step = 2:nSteps-1
    modIdx = mod([step-1 step step+1] - 1, 3) + 1;
    % The slowest part of the computation when there are relatively many
    % boundaries is indexing into vP for the boundary summations:
    % vP(interBoundIdx,step[+/-]1). Trivially, the step-1 term is used
    % twice, so can be pulled out into its own variable, but this should be
    % looked at in detail to further accelerate the computation.
    vPinterBoundMinus = vP(interBoundIdx,step-1);
    
    % Intermediate steps for clarity (and to avoid more complicated indexing)
    bValue(gMask) = Slm.*etalm.*( ...
        2*Llm.*vlb(:,modIdx(2))./T ...
        - gl(:,modIdx(2))./Clm ...
        - rho*vPinterBoundMinus./(2*T) ); % Do in separate step so that it can be indexed into a form that can be summed properly below
    iValue = S./h.*reshape(vP(interSubtractionIdx,step) - vP(interRowMask,step), size(conn)); % Just for convenience

    % Compute velocity potential (eq. 39)
    vP(:,step+1) = ( 2*vP(:,step) - vP(:,step-1) ...
        + (c^2*T^2./V).*sum( iValue ,2) ...
        - (c^2*T^2./(2*V)).*sum(sum( bValue ,3),2) ) ...
        ./ denominator;
    
    % Compute boundary quantities (eq. 41a-b)
    vlb(:,modIdx(3)) = etalm.*((Llm/T - Rlm/2 - T./(4*Clm)).*vlb(:,modIdx(2)) ...
        - gl(:,modIdx(2))./Clm ...
        + rho/(2*T)*(vP(interBoundIdx,step+1)-vPinterBoundMinus));
    gl(:,modIdx(3)) = (T/2)*(vlb(:,modIdx(3))+vlb(:,modIdx(2))) + gl(:,modIdx(2));

    %% Energy balance validation
    
    % As suggested in the text, conservation of energy can be used to debug
    % and verify that the method is performing as expected (see Fig 13).
    % Since the boundary and dissipated energies refer to the boundary
    % quantities that are not stored, we have to compute them in the loop.
    
    if verifyEnergy
        vjk(:,:,step) = reshape(vP(interRowMask,step)-vP(interSubtractionIdx,step),size(conn))./h;
        p(:,step) = rho*(vP(:,step)-vP(:,step-1))/T;

        hi(step) = sum(V.*(p(:,step).^2)/(2*rho*c^2) + sum(rho.*S.*h.*vjk(:,:,step).*vjk(:,:,step-1)/4,2)); % (eq. 27a)
        hb(step) = sum(Slm.*(Llm.*vlb(:,modIdx(2)).^2 + 1./Clm.*gl(:,modIdx(2)).^2))/2; % (eq. 33)
        q(step) = sum(Slm.*Rlm.*((vlb(:,modIdx(2))+vlb(:,modIdx(1)))/2).^2); % (eq. 32)
    end

    %% Display some progress
    if displayProgress
        if (curPctIdx <= length(pctPts)) && (step/nSteps > pctPts(curPctIdx))
            fprintf('\b\b\b\b%02d%%\n',round(pctPts(curPctIdx)*100))
            curPctIdx = curPctIdx + 1;
        end
    end
end

if displayProgress, toc; end
% whos % Display sizes of all variables to see what takes up the most

%% Post-hoc computation

% p2 = rho*diff(vP,1,2)/T; % Ok for all except first and last value: all(p(:,2:end-1) == p2(:,1:end-1), 'all')
% vjk2 = reshape(vP(interRowMask,:)-vP(interSubtractionIdx,:),[size(conn) nSteps])./h; % Same as above: all(vjk(:,:,2:end-1) == vjk2(:,:,2:end-1), 'all')

%% Visualize results

if verifyEnergy
    % Energy balance
    figure
    plot(hi(2:end))
    hold on
    plot(hb(2:end))
    plot(cumsum(q(2:end)*T))
    htot = hi(2:end)+hb(2:end)+cumsum(q(2:end)*T);
    plot(htot)
    legend('Internal energy','Boundary energy','Cumulative dissipated energy','Sum')
end


