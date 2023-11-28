%Tiffany Anderson

%This script calculates rates and intercepts using ST
% 
% Modified from NKauai_ST_corr.m
% 10 Feb 2015
% Run on Oahu data - Sept 8, 2015
% 
% Modified: 4/3/2017 Modified TAnderson to include boundaries in plots.  
%           Set up to run on updated Kauai shorelines. 
% Modified: Sept 7, 2018 TAnderson Run on Oahu updated dataset (including 2014? imagery)
% Modified: Sept 11, 2018 TAnderson. Changed name from
%           'AA_ST_smooth_MAIN_Kauai2017.m' to 'AA_ST_smooth_MAIN_Oahu2018.m'
% Modified: Sept 17, 2019 TAnderson.  Run on Maui updated dataset (including 
%           2016 imagery).  Also, changed name from
%           'AA_ST_smooth_MAIN_Oahu2018.m' to
%           'AA_ST_smooth_MAIN_Maui2019.m'
%           Also added flag for adding an additional degree of freedom for
%           the rare case of only 3 shorelines. (see note below for
%           justiication) Will output to screen the transect #s where an 
%           additional degree of freedom was added. (So far, only Maui)
% Modified: Dec 14, 2020, TAnderson.  Changed name to
%           'AA_ST_smooth_MAIN_Oahu2020.m' to run on updated Oahu data 
%           collected with drone. 
% Modified: Feb 24, 2023, TAnderson.  Changed name to AA_ST_smooth_MAIN.m
%           to reduce redundancy of filenames.
%           ALSO, included reading in vegetation data, and saving that data
%           in the matlab structures and workspaces.  

clear all

% ----------------------------------------------------------------------
% BEGINNING OF USER INPUT AREA

% flag that, if true, adds one degree of freedom to the t-distribution. 
% After truncation, there may be only 3 shorelines; this results in one 
% degree of freedom, which would create an unusably large uncertainty. 
% Thus, we assume that the limited data, and results, match those of the 
% surrounding area. It is sensible that the shoreline there, would behave 
% similarly to the surrounding area. Looking at the rates, this appears 
% true. Also, the uncertainty is larger than the surrounding areas and  
% does not predict rates that are less well-behaved than the surrounding 
% areas. 
flag_df_one_add_one = true;
% flag_df_one_add_one = false; % use this for Kauai, where the only transects 
                             % with less than 4 shorelines is Kapaa,
                             % transects 89-98.  This is near the Kapaa
                             % Pool, where rock revetment was built
                             % pre-1975, to stop erosion in response to 
                             % nearshore dredging.  Shoreline here deviates
                             % from bearby, so hard to justify adding an
                             % additional degree of freedom.  

% Input region name (NKauai, EKauai, SKauai, WKauai)
% region = 'NKauai';
% region = 'EKauai';
% region = 'SKauai';
% region = 'WKauai';
 
% Input region name (NOahu, EOahu, SOahu, WOahu)
%  region = 'NOahu';
region = 'EOahu';
% region = 'SOahu';
% region = 'WOahu';

% Input region name (Kihei, WMaui, NMaui)
%  region = 'Kihei';
% region = 'WMaui';
% region = 'NMaui';

savefig = true; % figures are saved if this is set to true. 
% savefig = false;

% END OF USER INPUT AREA
% ----------------------------------------------------------------------


m2ft = 3.28084; % conversion factor for meters to feet

% make directory called '<region>_results' if it does not exist
if ~isdir([region '_results'])
    mkdir([region '_results'])
end
% make directory called 'Figs' in existing dir '<region>_results' 
% if it does not exist
if ~isdir([region '_results/Figs'])
    mkdir([region '_results/Figs'])
end
    
% Define directory name and full filename where section file is located
regiondirname = [pwd '/' region '_dat/'];
sectfilename = [region '_sections.txt'];
sectfile = [regiondirname sectfilename];
sectid = fopen(sectfile);
sects = textscan(sectid,'%s');
fclose(sectid); % close file
sects = sects{1}; % convert to cell of strings

results_all = cell(length(sects),1); 

for sectnum = 1:length(sects)
    
    sectname = sects{sectnum};
    fprintf('Processing %s...\n',sectname)


    % Define directory name and full filename where data file is located
    dirname = [regiondirname sectname '_dat/'];
    modelname = [sectname 'toedist.txt'];
    datfile = [dirname modelname];

    % Define full filename of truncation file (identifies hardened shorelines)
    truncname = [sectname 'truncation.txt'];
    truncfile = [dirname truncname];
    
    % Define full filename of boundary file (identifies alongshore breaks in continuity)
    boundname = [sectname 'boundary.txt'];
    boundfile = [dirname boundname];
    
    % Define directory name and full filename where veg line data file is located
    vegname = [sectname 'vegdist.txt'];
    vegfile = [dirname vegname];


    
    % ---------------------------  
    % load data, veg, boundary, and truncation files
    data_raw = load(datfile); % loads data matrix (data before any trucation)
    trunc = load(truncfile); % loads truncation info into matrix
    bounds = load(boundfile); % load boundary file
    veg_raw = load(vegfile); % loads veg matrix
    
    if isempty(data_raw)
        error('Data file is empty.  Please check file format. Text MS-DOS works.')
    end
    if isempty(bounds)
        error('Boundary file is empty.  Please do not use "NaN" - define entire boundary.')
    end
    if isnan(bounds)
        error('Please do not use "NaN" - define one boundary instead.')
    end
    if isempty(trunc)
        error('Truncation file is empty. If no trunction, please use "NaN"')
    end
    if isempty(veg_raw)
        error('Veg dist file is empty.  Please check file format. Text MS-DOS works.')
    end
    if ~(length(veg_raw(1,2:end)) == length(data_raw(1,2:end)))
        error('Different number of survey dates in vegdist and toedist file.  Please check data.')
    end
    if ~(length(veg_raw(3:end,1)) == length(data_raw(3:end,1)))
        error('Different number of transects in vegdist and toedist file.  Please check data.')
    end
    if ~(sum( veg_raw(1,2:end) - data_raw(1,2:end) ) < 0.0001)
        error('Dates in veg dist file do not match those in toedist.  Please check data.')
    end
    if ~(sum( veg_raw(3:end,1) - data_raw(3:end,1) ) <= 0)
        error('Transect numbers in veg dist file do not match those in toedist.  Please check data.')
    end


    
    
    
    % -----------------------------------------------
    % define time series segments 
    t_data = data_raw(1,2:end);
    
    % Truncate full data series so that duplicate seawall positions are set to NaN
    data_trunc = data_raw; % will contain data after truncation.
    
    % put truncated values in matrix for plotting later
    truncated_data = nan(size(data_raw));
    % put hard shoreline values in matrix for plotting later
    hardshore_data = nan(size(data_raw));  % hard shoreline only if it is NOT the last data point in the time series. 
    hardshore_data_flat = nan(size(data_raw)); % all hard shorelines after the first one, have same position as first instance of hard shoreline
    hardshore_data_flat_truncated = nan(size(data_raw)); % flat hard shoreline positions that were truncated only (not including the first instance of the hard shoreline)

    % Based on trunctaion file, remove duplicate data due to seawall 
    % construction or other hard structure by setting those data to NaN 
    % NOTE: This does not change the size of the matrix (doesn't remove 
    %       any rows or columns), it only replaces the duplicate data with NaNs
    if isnan(trunc(1,1))
        % do nothing - no transect data to truncate 
    else
        for ind = 1:size(trunc,1)
            %finds the START row index (start transect)
            [row_start, ~] = find(data_trunc(:,1)==trunc(ind,1),1,'last'); % used "last" becuase sometimes the first two rows have the number "0", which could be erroneously taken as the transect number "0". 
            %finds the END row index (end transect)
            [row_end, ~] = find(data_trunc(:,1)==trunc(ind,2),1,'last');
            %keep columns based on columns 3 and 4 in trunc file
            cols_keep = trunc(ind,3):trunc(ind,4);
            cols_tmp = 1:size(data_trunc,2)-1;
            cols_tmp(cols_keep) = NaN;
            %find all columns that are NOT nanned. 
            cols_remove = cols_tmp(~isnan(cols_tmp)) + 1; % add 1 because indexing into matrix data_orig (date indices are offset by one due to transect numbers being in the first column)
            
            % only consider data as hard shoreline if it is NOT the last data in the time series 
            %   (truncating early data can lead to last data point to be end of truncated series)
            if any(sum(~isnan(data_raw(row_start:row_end,(trunc(ind,4)+1):end)),2) > 1)
                hard_dat = data_raw(row_start:row_end,(trunc(ind,4)+1):end); %the data that will be removed...?
                hard_dat_nan = isnan(hard_dat); % boolean matrix 
                hardshore_data(row_start:row_end,(trunc(ind,4)+1):end) = hard_dat; %everything is nan except for hard_data...?
                hard_dat_flat = repmat(hard_dat(:,1),1,size(hard_dat,2));  % helps when plotting so hard shorelines look stable (instead of moving a little)
                hard_dat_flat(hard_dat_nan) = NaN; % put nans in places where original data was nan
                hardshore_data_flat(row_start:row_end,(trunc(ind,4)+1):end) = hard_dat_flat;
                hardshore_data_flat_truncated(row_start:row_end,(trunc(ind,4)+2):end) = hard_dat_flat(:,2:end); % only inlcudes data that were removed from analysis
 
            end
            clear hard_dat hard_dat_flat
            
            truncated_data(row_start:row_end,cols_remove) = data_raw(row_start:row_end,cols_remove); % all data points that were removed in the truncation process
            data_trunc(row_start:row_end,cols_remove) = NaN;
        end
    end
    clear i row_start row_end cols_keep cols_tmp cols_remove
    truncated_data = truncated_data(3:end,2:end); % only keep the data portion (no tr num or dates or errors)
    hardshore_data = hardshore_data(3:end,2:end); % only keep the data portion (no tr num or dates or errors)
    hardshore_data_flat = hardshore_data_flat(3:end,2:end); % only keep the data portion (no tr num or dates or errors)
    hardshore_data_flat_truncated = hardshore_data_flat_truncated(3:end,2:end); % only keep the data portion (no tr num or dates or errors)
    
    % perform ST regression
    % NOTE: Rates in ST structure are negative (retreat) and positive (advance)
    [ST_sect, ST_sect_figs] = AAA_ST_regress(data_trunc, bounds, flag_df_one_add_one);
    ST_sect.truncated_data_m = truncated_data; % add truncated data matrix to structure (meters)
    ST_sect.truncated_data_ft = truncated_data*m2ft; % add truncated data matrix to structure (feet)
    ST_sect.hardshore_data_m = hardshore_data; % add hard shoreline data matrix to structure (meters)
    ST_sect.hardshore_data_ft = hardshore_data*m2ft; % add hard shoreline data matrix to structure (feet)
    ST_sect.hardshore_data_flat_m = hardshore_data_flat;
    ST_sect.hardshore_data_flat_ft = hardshore_data_flat*m2ft;
    ST_sect.hardshore_data_flat_truncated_m = hardshore_data_flat_truncated; % only the hard shore data that was removed in analysis (m)
    ST_sect.hardshore_data_flat_truncated_ft = hardshore_data_flat_truncated*m2ft; % only the hard shore data that was removed in analysis (ft)
    
%     % create matrix of errors only at hard shorelines, and add to ST_sect struct
%     hardshore_err = ST_sect.m_err;
%     hardshore_err(isnan(ST_sect.hardshore_data_m)) = NaN; % error corresponding to hard shorelines
%     hardshore_err_truncated = hardshore_err;
%     hardshore_err_truncated(isnan(ST_sect.hardshore_data_flat_truncated_m)) = NaN; % remove first instance of hard shoreline
%     ST_sect.hardshore_err_m = hardshore_err; % add hard shoreline err matrix to structure (meters)
%     ST_sect.hardshore_err_ft = hardshore_err*m2ft; % add hard shoreline err matrix to structure (feet)
%     ST_sect.hardshore_err_truncated_m = hardshore_err_truncated; % add truncated hard shoreline err matrix to structure (meters)
%     ST_sect.hardshore_err_truncated_ft = hardshore_err_truncated*m2ft; % add truncated hard shoreline err matrix to structure (feet)
    
    
    if (savefig)
        % save figures of data and parameters 
        saveas(ST_sect_figs.data_METERS,...    % save data plot (m)
            [region '_results/Figs/' sectname '_' ...
            int2str(ST_sect.t_data(1)) '_' ...
            int2str(ST_sect.t_data(end)) '_data_METERS.fig'])
        saveas(ST_sect_figs.params_METERS,...  % save parameter plot (m)
            [region '_results/Figs/' sectname '_' ...
            int2str(ST_sect.t_data(1)) '_' ...
            int2str(ST_sect.t_data(end)) '_params_METERS.fig'])
        saveas(ST_sect_figs.data_FEET,...    % save data plot (ft)
            [region '_results/Figs/' sectname '_' ...
            int2str(ST_sect.t_data(1)) '_' ...
            int2str(ST_sect.t_data(end)) '_data_FEET.fig'])
        saveas(ST_sect_figs.params_FEET,...  % save parameter plot (ft)
            [region '_results/Figs/' sectname '_' ...
            int2str(ST_sect.t_data(1)) '_' ...
            int2str(ST_sect.t_data(end)) '_params_FEET.fig'])
    end
    close(ST_sect_figs.data_METERS)   % close data figure (m)
    close(ST_sect_figs.params_METERS) % close parameter figure (m)
    close(ST_sect_figs.data_FEET)   % close data figure (ft)
    close(ST_sect_figs.params_FEET) % close parameter figure (ft)

    results_area.name = sectname;
    results_area.ST = ST_sect;
    results_area.dates = ST_sect.t_data;
    results_area.data_raw = data_raw;
    results_area.veg_raw = veg_raw;
    clear res_tmp ST_sect ST_sect_figs
    
    % put area results into cell array with all other areas
    results_all{sectnum} = results_area;
    
    % save workspace
    save([region '_results/workspace_' sectname])

end

save([region '_results/workspace_ALL']) % save workspace
save([region '_results/AA_results_all'],'results_all') % save "results_all" struct

fprintf(1,'Pau: %s rates.\n', region);


