function [MaxDD, MaxDDIndex] = maxdrawdown(Data, Format)
%MAXDRAWDOWN Calculate maximum drawdown for one or more price series.
%	Given a T x N matrix of Data with T observations of N total return price series (also known as
%	total equity), compute maximum drawdown for each series in an N vector MaxDD and identify start
%	and end indexes of maximum drawdown period for each series in a 2 x N matrix MaxDDIndex.
%
%		MaxDD = maxdrawdown(Data);
%		MaxDD = maxdrawdown(Data, Format);
%		[MaxDD, MaxDDIndex] = maxdrawdown(Data, Format);
%
%	Inputs:
%		Data - T x N matrix with T samples of N total equity time series with earliest data in row
%			T(1,:) and most recent data in row T(end,:).
%
%	Optional Inputs:
%		Format - String to indicate format of data. Options are:
%			'return' (default) - Compute maximum drawdown as a maximum percentage drop from a peak.
%			'arithmetic' - Compute maximum drawdown of a Brownian motion with drift (differences of
%				data from peak to trough).
%			'geometric' - Compute maximum drawdown of a geometric Brownian motion with drift
%				(differences of log of data from peak to trough).
%
%	Outputs:
%		MaxDD - 1 x N vector with maximum drawdown for each of N time series.
%		MaxDDIndex - 2 x N vector of latest start and earliest end indexes for each maximum drawdown
%			period for each total equity time series, where the first row contains the start indexes
%			and the second row contains the end indexes of each maximum drawdown period.
%
%	Notes:
%		Drawdown is the drop in total returns from the start to the end of a period. If the total
%		equity time series is increasing over an entire period, drawdown is zero. Otherwise, it is a
%		positive number. Maximum drawdown is an ex-ante proxy for "downside risk" that computes the
%		largest drawdown over all intervals of time that can be formed within a specified interval
%		of time.
%
%		Maximum drawdown is sensitive to quantization error, i.e., daily and monthly data over
%		identical time periods will usually have different values for maximum drawdown.
%
%		If raw data are in return form, convert to a price series in the function call with:
%			MaxDD = maxdrawdown(ret2tick(Data), Format);
%
%		If a series never has a drawdown, the value MaxDD is 0 and MaxDDIndex has NaNs for index
%		values.
%
%		Maximum drawdown requires positive input data unless the format is 'arithmetic'.
%
%		If two or more periods exist with the same maximum drawdown, the indexes
%		for the earliest period are returned.
%
%	See also emaxdrawdown

%	Copyright 1995-2018 The MathWorks, Inc.

% Step 1 - check arguments

if nargin < 1 || isempty(Data)
	disp('Error 1');
end

if ~isscalar(Data) && isvector(Data) && isa(Data,'double')
	Data = Data(:);
	[~, N] = size(Data);
elseif ismatrix(Data) && size(Data, 1) > 1 && isa(Data,'double')
	[~, N] = size(Data);
else
	disp('Error 2');
end

if nargin < 2 || isempty(Format)
	choice = 1;
else
    Format = convertStringsToChars(Format);
    if ~ischar(Format) || size(Format,1) ~= 1
       disp('Error 3');
    else
        choice = find(strncmpi(Format,{'return','arithmetic','geometric'},length(Format)));
        if isempty(choice)
         disp('Error 4');
        end
    end
end

% Step 2 - compute maximum drawdown

%if choice == 1 || choice == 3
%	if any(any(Data <= 0))
%		disp('Error 5');
%	end
%end

if choice == 3
    Data = log(Data);
end

MaxDDIndex = ones(2,N);
% each peak can potentially be the start of a max drawdown
peaks = cummax(Data);       
if choice == 1      % 'return' format
    drdn = (peaks - Data)./peaks;
else                % 'arithmetic' or 'geometric' formats
    drdn = (peaks - Data);
end
[MaxDD, MaxDDIndex(2, :)] = max(drdn);
MaxDD(isnan(MaxDD)) = 0;

% find the peak value that result in MaxDD
linearIdx = sub2ind(size(Data), MaxDDIndex(2, :), 1:N);  
peakLeadsToMaxDD = peaks(linearIdx);   

% find the first idx with the peak value
mask = peakLeadsToMaxDD == Data;
[~, MaxDDIndex(1, :)] = max(mask);

k = MaxDDIndex(1,:) == MaxDDIndex(2,:);
MaxDDIndex(:,k) = NaN;




