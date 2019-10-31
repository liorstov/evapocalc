function [scale, shape] = wblfit_LR_data_portion(data, data_portion)
% estimates parameters of the Weibull distribution underlying the sample in
% 'data' using the method of the linear regression in the transformed
% coordinates: log(log(1/(1-ecdf))); log(data); (where ecdf is the
% empirical cumulative distribution function). Only the quantitative
% information in the data included in 'data_portion' is used; the weight in
% probability of the rest of the data is cosidered - see Marra & al (2019)  
%   
%   Marra & al 2019: Marra F, D Zoccatelli, M Armon, E Morin, 2019. A
%       simplified MEV formulation to model extremes emerging from multiple
%       nonstationary underlying processes. Adv. Water Resour., 127,
%       280-290, doi:10.1016/j.advwatres.2019.04.002 
% 
% Input
% - data: array of data to be fitted
% - data_portion: 2-elements array with the limits in probability of the
%       data to be used for the parameters estimation
%       e.g. data_portion = [0.75, 1] uses the largest 25% of the data
% 

% computes empirical cdf
X = sort(data,'ascend');
ECDF = (1:numel(X))'/(numel(X)+1); % ecdf computed with Weibull plotting positions

% uses data in data portion
fidx = max(1,fix(numel(X)*data_portion(1)));
tidx = ceil(numel(X)*data_portion(2));

% independent variable is the frequency, dependent variable the amount
Results = weightedfit([log(log(1./(1-ECDF(fidx:tidx)))), log(X(fidx:tidx)), ones(size(X(fidx:tidx)))]);
shape = 1/Results.slope;
scale = exp(-Results.Intercept/Results.slope).^(-1/shape);

end

function Results=weightedfit(data)
% 
%Coded by Ebo Ewusi-Annan (University of Florida, 2011)
% modified by F. Marra, 2018
% 
% This code fits makes a linear fit to a data set (using y =bx+a) where
% each data point has a different or constant standard deviation. Your data
% should have three or two columns. The first column should be the
% independent variable(x) and the second column should be the dependent
% variable(y). Column three should contain your standard deviations for
% each datapoint. In the situations where you do not specify a column
% three, the code assigns a weight of one to all data points and this
% corresponds to the regular linear fits. 
% INPUTS: 
%data = 3 columns; column 1 = x, column2 = y and column 3 = standard dev.
% 
% OUTPUTS
%Result.slope= b; Fitted slope
%Result.Intercept = a; Fitted intercept
% 
%REFERENCES
%1. Willam H. Press, Saul A. Teukolsky and Willan T. Vetterling (1997).
%Numerical Recipes in Fortran.
%2. Philip R. Bevington and D. Keith Robinson (2003). Data Reduction and
%Error Analysis for the Physical Sciences.
   
x= data(:,1);
y=data(:,2);
[s t]= size(data);
stdv=ones(s,1);
if t==3; stdv=data(:,3); end
w = 1./stdv.^2;
S = sum(w);
Sx = sum(w.*x);
Sy = sum(w.*y);
Sxx= sum(w.*x.^2);
Sxy= sum(w.*x.*y);
Delta = S*Sxx - (Sx)^2;
a = (Sxx*Sy - Sx*Sxy)./Delta;
b = (S*Sxy - Sx*Sy)./Delta;
Results.slope=b;

end