%% double check the orientation of xrms
% Select sequence file
function [xrms,x0] = check_orient_xrms(baseFolder, matrix_size)
    if nargin<2||isempty(matrix_size)
        matrix_size=48;
    end
    
disp('Select main meas sequence file:');
[measSeqName, seqFolder] = uigetfile('*.seq', 'Select main meas sequence file', baseFolder);
if measSeqName == 0
    error('Sequence file selection was cancelled');
end

% Select raw data file
disp('Select main raw data file:');
[measFilename, rawDir] = uigetfile('*.dat', 'Select main raw data file', baseFolder);

if measFilename == 0
    error('measurement data file selection was cancelled');
end


measureFile = fullfile(rawDir, measFilename);

seqFile = fullfile(seqFolder, measSeqName);
measSeqParams = extract_seq_params(seqFile);

reader = createRawDataReader(measureFile, 1);
if isfield(measSeqParams, 'nshot')
    reader.acquisitionParams.nShot = measSeqParams.nshot;
end
if isfield(measSeqParams, 'nseg')
    reader.acquisitionParams.nSeg = measSeqParams.nseg;
end

reader.acquisitionParams.nShot_off = 20;
reader.acquisitionParams.traj_type = 'pulseq';
%
reader.acquisitionParams.pulseqTrajFile_name = fullfile(seqFolder,measSeqName);
% check if the hash from pulseq sequence and from twix match each other
isMatch = check_hash(measureFile,reader.acquisitionParams.pulseqTrajFile_name);

% Ensure consistency in number of shot-off points
nShotOff = reader.acquisitionParams.nShot_off;
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI
t_tot = bmTraj(reader.acquisitionParams);                       % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements
% ==============================================
% Warning: due to the memory limit, make sure the matrix size <=240

N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./240;

% ------
nCh = size(y_tot, 1);
nCh
nFr = 1;
x0 = cell(nCh, 1);
for i = 1:nFr
    for iCh = 1:nCh
    x0{iCh} = bmMathilda(y_tot(iCh,:), t_tot, ve_tot, [], N_u, N_u, dK_u, [], [], [], []);
    disp(['Processing channel: ', num2str(iCh),'/', num2str(nCh)])
   
    end
end

%
bmImage(x0);


% Root mean square across the channels
% Initialize an array to store sum of squared images
[nx, ny, nz] = size(x0{1});  % Get the dimensions (240,240,240)
numCoils = numel(x0);  % Number of coils (20)

sum_of_squares = zeros(nx, ny, nz, 'single');  % Preallocate in single precision

% Compute sum of squared images

for coil = 1:numCoils
    % straightforward
    % sum_of_squares = sum_of_squares + abs(x0{coil}).^2;
    % eliminate extra square-root step
    sum_of_squares = sum_of_squares + real(x0{coil}.*conj(x0{coil}));
end

% Compute the root mean square (RMS)
xrms = sqrt(sum_of_squares / numCoils);  % Normalize by the number of coils
bmImage(xrms)

end