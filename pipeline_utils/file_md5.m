function digest = file_md5(path)
%FILE_MD5 Hash file bytes without requiring a MATLAB JVM.
assert(isfile(path),'File not found: %s',path);
if usejava('jvm')
    fid = fopen(path,'rb'); assert(fid>=0,'Cannot read %s',path);
    cleanup = onCleanup(@() fclose(fid));
    md = java.security.MessageDigest.getInstance('MD5');
    while ~feof(fid)
        bytes = fread(fid,1024*1024,'*uint8');
        md.update(bytes);
    end
    digest = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
else
    if ismac
        command = ['/sbin/md5 -q ' shell_quote(path)];
    elseif isunix
        command = ['md5sum ' shell_quote(path)];
    else
        command = ['certutil -hashfile ' shell_quote(path) ' MD5'];
    end
    [status,result] = system(command);
    assert(status==0,'MD5 command failed: %s',result);
    found = regexp(result,'(?mi)^\s*([0-9a-f]{32})(?:\s|$)','tokens','once');
    assert(~isempty(found),'Could not read an MD5 digest from: %s',result);
    digest = lower(found{1});
end
end
