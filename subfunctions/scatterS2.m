function h = scatterS2(varargin) 
% scatterS2 plot scatters of data on unit sphere
% scatterS2(A), A: 3 x n matrix, consisting of n column direction vectors.
% scatterS2(A,paramstruct), see below for paramstruct.
% 
% Input: paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%     paramstruct = struct('field1',values1,...
%                          'field2',values2,...
%                          'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%    fields           values
%
%    titlestring      string for title
%
%    isize            numeric (default 10)
%
%    icolor           0  (default) fully black and white version (everywhere)
%                     string (any of 'r', 'g', 'b', etc.) that single color
%                     1  (default)  color version (Matlab 7 color default)
%                     2  time series version (ordered spectrum of colors)
%                     nx3 color matrix:  a color label for each data point
%                             (to be used everywhere, except SiZer & QQ
%                              useful for comparing classes)
%
%
% Last modified by Sungkyu Jung, Feb 2013

if nargin==0, help plotshapes, return, end     %   Display help
data = varargin{1};
[d n]=size(data);
if d ~= 3 ; disp('Warning from scatterS2.m: input data set is not on S^2'); return,end;

isize = 10;
icolor = repmat([0 0 0],n,1);
titlestring = 'Scatter on unit 2-sphere';

if nargin > 1
    paramstruct = varargin{2};
    if isfield(paramstruct,'isize') ;    %  then change to input value
        isize = paramstruct.isize;
    end ;
    if isfield(paramstruct,'titlestring') ;    %  then change to input value
        titlestring = paramstruct.titlestring;
    end ;
    if isfield(paramstruct,'icolor') ;    %  then change to input value
        icolor = paramstruct.icolor;
        if length(icolor)==1 && icolor <  4
            icolor = repmat('b',n,1);
        end
        if ischar(icolor) || size(icolor,1) == 1
            icolor = repmat(icolor,n,1);
        end
    end ;
end

scatter3(data(1,:),data(2,:),data(3,:),isize,icolor,'filled');
hold on;
view(geodmeanSk(data,0.01))
% draw sphere
[q w e] = sphere(100);
surface(q,w,e,'Edgecolor','none','Facecolor',[0.9 0.9 0.9],'facealpha',0.2);
light('position',[1 0 1]);
xlabel('x');ylabel('y');zlabel('z');
title(titlestring)
axis equal;
hold off;

