function writefigs(H,pdfname)
% Write figures to pdf
 
 try
% Append to pdf
if isstruct(H)  % H is either a struct with figure handles, or a single fig. handle
    fls = fieldnames(H);
    for fs = 1:length(fls)
        h = H.(fls{fs});
        if ishandle(h)
            export_fig(h,'-dpdf',pdfname,'-painters','-append');  % write to pdf
        end
    end
else
    export_fig(H,'-dpdf',pdfname,'-painters','-append');  % write to pdf
end

catch
fprintf('writefigs failed. export_fig package needed.\n');
end