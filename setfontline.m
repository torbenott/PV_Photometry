function setfontline(fs, ls, fn)
    
% formatfig(fs,ls)
%   fs: font size
%   ls: line size
%   fn: font name

if nargin == 0
    fs = 14; % paper 14.  summary 7
    ls = 2;
    fn = 'Arial';
elseif nargin == 1
    ls = 1;
    fn = 'Arial';
elseif nargin == 2
    fn = 'Arial';
elseif nargin > 3
%     error('Numbers of input argument is incorrect');
end

h = get(gcf,'children');
n = size(h);
for i=1:n
    t = get(h(i),'type');
    if strcmp(t,'axes')
        set(h(i),'fontsize',fs,'fontname',fn);
        set(get(h(i),'Xlabel'),'fontsize',fs,'fontname',fn);
        set(get(h(i),'Ylabel'),'fontsize',fs,'fontname',fn);
        set(get(h(i),'Title'),'fontsize',fs,'fontname',fn);
        ch = get(h(i),'children');
        nch = length(ch);
        for j=1:nch
           try set(ch(j),'linewidth',ls);
           catch error
           end
        end
    end
end

