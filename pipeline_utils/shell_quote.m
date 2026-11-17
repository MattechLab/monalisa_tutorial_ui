function quoted = shell_quote(path)
%SHELL_QUOTE Quote one path argument for the platform's system shell.
path = char(path);
if ispc
    assert(~contains(path,'"') && ~contains(path,'%'), 'Unsupported character in shell path.');
    quoted = ['"' path '"'];
else
    quoted = ['''' strrep(path,'''','''"''"''') ''''];
end
end
