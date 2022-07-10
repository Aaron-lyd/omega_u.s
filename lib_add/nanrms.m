function y = nanrms(varargin)

narginchk(1,2);
y = rms(varargin{:},'omitnan');
