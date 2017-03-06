% Load raw data (.csv file)
% Location of raw input file
source_raw = '\\ad.uws.edu.au\dfshare\HomesHWK$\90929058\My Documents\GitHub\Alexis et al 2017 input\CUP_EddyPro_QC_170207.csv';
% Create datastore to access collection of data
ds_raw = datastore(source_raw);
% Select variable of interest
ds_raw.SelectedVariableNames = {'DateTime','NEE_c'}; % for example
% Read selected variables, save it in the workspace as a table
Raw_Data = readall(ds_raw);
% Get data from the table, change format if necessary
NEE = Raw_Data.NEE_c;
DateTime_CUP_cell = Raw_Data.DateTime;
DateTime_CUP = datetime(DateTime_CUP_cell,'InputFormat','dd/MM/yyyy HH:mm');
% clear unused variables
clearvars DateTime_CUP_cell Raw_Data ds_raw source_raw; 

% Scripts will use steps presented above to load the necessary raw data input to
% produce output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load processed data (netCDF file)
% Location of processed input file 
source_processed = '\\ad.uws.edu.au\dfshare\HomesHWK$\90929058\My Documents\GitHub\Alexis et al 2017 input\CumberlandPlain_2014_L6_EP_moderate.nc';
% Information on file, including variable name
finfo = ncinfo(source_processed);
% Read variable from .cd file
vardata = ncread(source_processed,'ER'); % for example
ER = reshape(vardata,[],1); % to get a 1 column vector instead of 3 dim matrix
% clear unused variables
clearvars source_processed finfo vardata;




