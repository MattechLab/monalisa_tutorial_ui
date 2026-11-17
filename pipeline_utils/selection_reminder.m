function selection_reminder(titleText,message)
%SELECTION_REMINDER Show instructions before native file/folder selection.
% macOS pickers can hide their title; this modal reminder remains visible
% until OK is clicked, then the caller opens its picker.
fprintf('\n%s\n%s\n',titleText,message);
if ismac
    uiwait(msgbox(message,titleText,'help','modal'));
end
end
