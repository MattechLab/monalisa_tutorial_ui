function info = read_twix_system_info(rawFile)
%READ_TWIX_SYSTEM_INFO Read scanner/version fields without loading raw samples.
% Uses MONALISA's header parser. For RAID containers, inspect the last scan,
% matching the tutorial's image-data reader. Return selected metadata only.
assert(isfile(rawFile),'TWIX file not found: %s',rawFile);
stat = dir(rawFile);
fid = fopen(rawFile,'rb','ieee-le');
assert(fid>=0,'Cannot open TWIX file: %s',rawFile);
fileCleanup = onCleanup(@() fclose(fid));
first = fread(fid,1,'uint32'); second = fread(fid,1,'uint32');
assert(~isempty(first) && ~isempty(second),'Truncated TWIX file.');
offset = 0; measurementLength = stat.bytes;
containerFormat = 'VB';
if first<10000 && second>=1 && second<=64
    containerFormat = 'RAID';
    % Directory entries are 152 bytes: measID, fileID, offset, length, names.
    assert(fseek(fid,16+152*(second-1),'bof')==0,'Cannot seek to RAID entry.');
    offset = fread(fid,1,'uint64'); measurementLength = fread(fid,1,'uint64');
end
assert(isscalar(offset) && isscalar(measurementLength) && offset>=0 && ...
    offset+measurementLength<=stat.bytes,'Invalid TWIX measurement bounds.');
assert(fseek(fid,offset,'bof')==0,'Cannot seek to measurement header.');
headerBytes = fread(fid,1,'uint32');
assert(isscalar(headerBytes) && headerBytes>=8 && headerBytes<=measurementLength && ...
    headerBytes<=32*1024*1024,'Invalid or oversized TWIX header.');
hdr = read_twix_hdr_JH_for_monalisa(fid);
assert(ftell(fid)<=offset+headerBytes,'Header parser exceeded the declared header size.');
info = struct('rawFile',rawFile,'scannerModel',hdr.Dicom.ManufacturersModelName, ...
    'softwareVersion',hdr.Dicom.SoftwareVersions, ...
    'measuredBaseline',hdr.Meas.tMeasuredBaselineString, ...
    'fieldStrengthT',hdr.Dicom.flMagneticFieldStrength, ...
    'headerSignature',hdr.Dicom.tSequenceVariant, ...
    'containerFormat',containerFormat,'headerBytes',headerBytes);
end
