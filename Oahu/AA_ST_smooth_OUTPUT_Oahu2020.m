% Tiffany Anderson
% 
% This script formats variables for output (tables)
% Also, makes tiled plots of data and trend at each transect
%
% 4 April 2017 Created as "AA_ST_smooth_OUTPUT_Kauai2017.m"
% Modified: April 5, 2018  Added single plots for ArcOnline
% Modified: Sept 7, 2018  Used this script to run Mokuleia, Oahu
% Modified: Sept 10, 2018 Changed name to "AA_ST_smooth_OUTPUT_Oahu2018.m"
%           Edited so that table output has header names as requested by
%           Jade for easy integration into ArcONLINE system.  Also, tables
%           include all areas for the island (not separate files per area).  
% 
%           Here are the headings for the m and ft historical shoreline rates:
%           Transect SRate_ft Suncert_ft Sstd_ft Num DegF AreaName Island
%			Transect SRate_m  Suncert_m  Sstd_m  Num DegF AreaName Island
% Modified: Jan 3, 2019. T. Anderson.  Removed extra zero in front of
%           single digit transects in the filename of single plots.  So, 
%           instead of transect # 00 in the filename, it is 0.  
% Modified: Dec 5, 2020 T. Anderson.  Edited table column headers to match 
%           the attribute names in the ArcONLINE shapefiles.  These were 
%           given by Kammie T. in late 2020. 
%           New header names are:
%           Transect sRate_m sUncert_m sStd_m Num DegF AreaName Island GISID           GISID
%           Transect sRate_ft sUncert_ft sStd_ft Num DegF AreaName Island GISID
%           Also, added the island name to the csv file of rates 
%             (ex. Oahu_rates_meters.csv)



clear all

% ----------------------------------------------------------------------
% BEGINNING OF USER INPUT AREA

% % Input region name (NKauai, EKauai, SKauai, WKauai)
% % regions = {'NKauai','EKauai','SKauai','WKauai'};
% regions = {'WKauai'};
% island = 'Kauai';

% % Input region name (NOahu, EOahu, SOahu, WOahu)
% regions = {'NOahu','EOahu','SOahu','WOahu'};
% regions = {'EOahu','SOahu'};
regions = {'EOahu'};
island = 'Oahu';

% Input region name (Kihei, WMaui, NMaui)
% regions = {'Kihei','WMaui','NMaui'};
% regions = {'Kihei'};
% regions = {'WMaui'};
% regions = {'NMaui'};
% island = 'Maui';

% write_tables = false;
write_tables = true;

transect_plots_multiple_meters = false;   % (meters) multiple transects per page, for printing
transect_plots_multiple_feet = false;     % (feet) multiple transects per page, for printing 
transect_plots_single_meters = false;     % (meters) singel transect per page, for ArcOnline
transect_plots_single_feet = false;       % (feet) singel transect per page, for ArcOnline

transect_plots_multiple_meters = true;   % (meters) multiple transects per page, for printing
transect_plots_multiple_feet = true;     % (feet) multiple transects per page, for printing 
transect_plots_single_meters = true;     % (meters) singel transect per page, for ArcOnline
transect_plots_single_feet = true;       % (feet) singel transect per page, for ArcOnline


% END OF USER INPUT AREA
% ----------------------------------------------------------------------

nmap_all = nan(length(regions),1);
map_hawn_names = [];
map_names = [];
map_sectnames = [];
map_trs = [];
map_regions = [];
map_island = [];
map_GISID = [];

for rind = 1:length(regions)
    
    region = regions{rind};

    % load file containing map areas and corresponding transect ranges
    areafile = [pwd '/' island '_map_areas.xlsx'];
    [NUM,TXT,RAW] = xlsread(areafile,region);
    map_hawn_names = [map_hawn_names; RAW(2:end,3)]; % Hawaiian map area names (used for plotting)
    map_names = [map_names; RAW(2:end,4)];   % map area names
    map_sectnames = [map_sectnames; RAW(2:end,5)];       % section names 
    map_trs = [map_trs; cell2mat(RAW(2:end,6:7))]; % transect numbers: [start end]
    nmap = length(RAW(2:end,1)); % number of map areas in this region
    nmap_all(rind) = nmap;    % put number of map areas for this region in a vector    
    map_regions = [map_regions; RAW(2:end,2)];
    map_island = [map_island; RAW(2:end,1)];
    clear areafile NUM TXT RAW nmap region
    
    
%     % load file containing map areas and corresponding transect ranges
%     regiondirname = [pwd '/' region '_dat/'];
%     areafilename = [region '_map_areas.txt'];
%     areafile = [regiondirname areafilename];
%     areaid = fopen(areafile);
%     datscan = textscan(areaid,'%s%s%d%d'); % read contents (4 columns) of area id file
%     fclose(areaid); % close file
%     map_names = [map_names; datscan{1}];   % map area names
%     map_sectnames = [map_sectnames; datscan{2}];       % section names 
%     map_trs = [map_trs; [datscan{3} datscan{4}]]; % transect numbers: [start end]
%     nmap = length(datscan{1}); % number of map areas in this region
%     nmap_all(rind) = nmap;    % put number of map areas for this region in a vector    
%     tmp_regions(1:nmap,1) = deal({region});  % region for each map area (same for all)
%     map_regions = [map_regions; tmp_regions];
%     tmp_island(1:nmap,1) = deal({island});   % island for each map area (same for all)
%     map_island = [map_island; tmp_island];
%     clear datscan areaid tmp_regions tmp_island areafile ...
%         areafilename regiondirname nmap region
      
end % read in all map areas, section names, and transects for each map area

% initialize cell array that will hold table and figure data; and a helper vector
out_results_all = cell(sum(nmap_all),1); % output matrices stored in cell array
region_ind_base = [0; cumsum(nmap_all(1:end-1))];
    
for rind = 1:length(regions)

    region = regions{rind};

    % load results from ST regression
    load([region '_results/AA_results_all.mat'])
    resmat = cell2mat(results_all); % convert cell array to structure array
    clear results_all

    nmap = nmap_all(rind); % number of map areas in this region

    for ind_area = 1:nmap % make output matrix for each map_area 

        % find index in resmat, where section for this map area is located.
        sect_ind = find(strcmp({resmat.name}, map_sectnames(region_ind_base(rind)+ind_area)));

    %     fprintf('map name: %s; sect ind: %d; sect name: %s\n', ...
    %         map_names{ind_area}, sect_ind, resmat(sect_ind).name)

        if isempty(sect_ind)
            error('Check section and map area names. Cannont find section.\n')
        end


        ST = resmat(sect_ind).ST;
        x_data = ST.x_data;
        tr_start_ind = find(x_data==map_trs(region_ind_base(rind)+ind_area,1),1);
        tr_end_ind = find(x_data==map_trs(region_ind_base(rind)+ind_area,2),1);

        % if can't find start or end index in the section, then quit
        if isempty(tr_start_ind) || isempty(tr_end_ind)   
            error('Cannot find map start/end transect index in section.\n')
        end


        Ncols = 6; % columns are: Tr num, sm rate, sm rate 95% unc, sm rate std, num shorelines, df
        Ntr = length(x_data(tr_start_ind:tr_end_ind));
        Ntime = length(ST.t_data);

        % output matrix (meters)
        out_mat_m = nan(Ntr, Ncols);               % initialize output matrix
        out_mat_m(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        out_mat_m(:,2) = ST.rate_sm(tr_start_ind:tr_end_ind);               % smoothed rate (meters)
        out_mat_m(:,3) = ST.rate_var_sm(tr_start_ind:tr_end_ind).^(1/2)...
            .*tinv(0.975,ST.df(tr_start_ind:tr_end_ind));                   % uncertainty (smoothed 95% conf interval (+/-), m)
        out_mat_m(:,4) = ST.rate_var_sm(tr_start_ind:tr_end_ind).^(1/2);    % smoothed standard deviation of rate (m)
    %     out_mat_m(:,5) = sum(~isnan(ST.y_data(tr_start_ind:tr_end_ind,:)),2); % number of shorelines
        out_mat_m(:,5) = ST.ndat(tr_start_ind:tr_end_ind);                                           % number of shorelines
        out_mat_m(:,6) = ST.df(tr_start_ind:tr_end_ind);                    % degrees of freedom (num shorelines - 2) 
        nrow = size(out_mat_m,1);                                           % number of rows (transects)
        area_name = cell2mat(map_names(region_ind_base(rind)+ind_area));              % area name for this area in for loop
        mat_col_area_name(1:nrow,1) = deal({area_name});                      % make cell array containing area name (size = #transects x 1)
        mat_col_island(1:nrow,1) = deal({island});                          % make cell array containing island name (size = #transects x 1)
        tr_num_str_trimmed = cellfun(@strtrim, cellstr(num2str(x_data(tr_start_ind:tr_end_ind))),'UniformOutput',false); % transect numbers without whitespace
        mat_col_GISID(1:nrow,1) = cellfun(@(x) [island '_' area_name '_' x],tr_num_str_trimmed,'UniformOutput',false);  % make cell array containing GISID (<IslandName>_<AreaName>_<TransectID> (size = #transects x 1)
        out_cell_m = cat(2,num2cell(out_mat_m),mat_col_area_name,mat_col_island,mat_col_GISID); % cell array of output for writing to table
        
        % output matrix (feet)
        out_mat_ft = nan(Ntr, Ncols);               % initialize output matrix
        out_mat_ft(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        out_mat_ft(:,2) = ST.rate_sm_ft(tr_start_ind:tr_end_ind);            % smoothed rate (meters)
        out_mat_ft(:,3) = ST.rate_var_sm_ft(tr_start_ind:tr_end_ind).^(1/2)...
            .*tinv(0.975,ST.df(tr_start_ind:tr_end_ind));                    % uncertainty (smoothed 95% conf interval (+/-), m)
        out_mat_ft(:,4) = ST.rate_var_sm_ft(tr_start_ind:tr_end_ind).^(1/2); % smoothed standard deviation of rate (m)
        out_mat_ft(:,5) = sum(~isnan(ST.y_data(tr_start_ind:tr_end_ind,:)),2); % number of shorelines
        out_mat_ft(:,6) = ST.df(tr_start_ind:tr_end_ind);                    % degrees of freedom (num shorelines - 2) 
        out_cell_ft = cat(2,num2cell(out_mat_ft),mat_col_area_name,mat_col_island,mat_col_GISID); % cell array of output for writing to table
        
        % transect plot matrix. Data values, (meters)
        % NOTE: positions are flipped, so they match trends (negative==erosion)
        plt_dat_m = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_dat_m(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_dat_m(:,2:Ntime+1) = -ST.y_data(tr_start_ind:tr_end_ind,:);      % truncated data, and trs with <3 shorelines removed (meters)

        % transect plot matrix. Data values, (feet)
        plt_dat_ft = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_dat_ft(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_dat_ft(:,2:Ntime+1) = -ST.y_data_ft(tr_start_ind:tr_end_ind,:);   % truncated data, and trs with <3 shorelines removed (feet)

        % transect plot matrix. Modeled values (trend), (meters)
        plt_mod_m = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_mod_m(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_mod_m(:,2:Ntime+1) = -ST.pos_sm(tr_start_ind:tr_end_ind,:);      % smoothed predicted positions (meters)

        % transect plot matrix. Modeled values (trend), (feet)
        plt_mod_ft = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_mod_ft(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_mod_ft(:,2:Ntime+1) = -ST.pos_sm_ft(tr_start_ind:tr_end_ind,:);   % smoothed predicted positions (feet)

        % Data errors matrix, (meters)
        plt_err_m = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_err_m(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_err_m(:,2:Ntime+1) = ST.m_err(tr_start_ind:tr_end_ind,:);      % data errors (meters)

        % Data errors matrix, (feet)
        plt_err_ft = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_err_ft(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_err_ft(:,2:Ntime+1) = ST.m_err_ft(tr_start_ind:tr_end_ind,:);   % data errors (meters)

        % Truncated data matrix, (meters)
        plt_trunc_m = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_trunc_m(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_trunc_m(:,2:Ntime+1) = -ST.truncated_data_m(tr_start_ind:tr_end_ind,:);      % data that were truncated (meters)

        % Truncated data matrix, (feet)
        plt_trunc_ft = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_trunc_ft(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_trunc_ft(:,2:Ntime+1) = -ST.truncated_data_ft(tr_start_ind:tr_end_ind,:);   % data that were truncated (feet)

        % Hard shoreline (flat) data matrix, (meters) 
        plt_hard_m = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_hard_m(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_hard_m(:,2:Ntime+1) = -ST.hardshore_data_flat_m(tr_start_ind:tr_end_ind,:);      % hardshore data points (meters)

        % Hard shoreline (flat) data matrix, (feet)
        plt_hard_ft = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_hard_ft(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_hard_ft(:,2:Ntime+1) = -ST.hardshore_data_flat_ft(tr_start_ind:tr_end_ind,:);   % hardshore data points (feet)

        % Hard shoreline truncated (flat) data matrix, (meters) 
        plt_hard_trunc_m = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_hard_trunc_m(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_hard_trunc_m(:,2:Ntime+1) = -ST.hardshore_data_flat_truncated_m(tr_start_ind:tr_end_ind,:);      % hardshore data points (meters)

        % Hard shoreline truncated (flat) data matrix, (feet)
        plt_hard_trunc_ft = nan(Ntr, Ntime+1);               % initialize output matrix
        plt_hard_trunc_ft(:,1) = x_data(tr_start_ind:tr_end_ind);                   % transect numbers
        plt_hard_trunc_ft(:,2:Ntime+1) = -ST.hardshore_data_flat_truncated_ft(tr_start_ind:tr_end_ind,:);   % hardshore data points (feet)


        % put matrices in structure, and store the structure in a cell array
        out_results_area.island = island;
        out_results_area.region = region;
        out_results_area.area_name = map_names{region_ind_base(rind)+ind_area};
        out_results_area.hawaiian_name = map_hawn_names{region_ind_base(rind)+ind_area};
        out_results_area.mat_meters = out_cell_m;
        out_results_area.mat_feet = out_cell_ft;
        out_results_area.plt_dat_m = plt_dat_m;
        out_results_area.plt_dat_ft = plt_dat_ft;
        out_results_area.plt_mod_m = plt_mod_m;
        out_results_area.plt_mod_ft = plt_mod_ft;
        out_results_area.plt_err_m = plt_err_m;
        out_results_area.plt_err_ft = plt_err_ft;
        out_results_area.plt_trunc_m = plt_trunc_m;
        out_results_area.plt_trunc_ft = plt_trunc_ft;
        out_results_area.plt_hard_m = plt_hard_m;
        out_results_area.plt_hard_ft = plt_hard_ft;
        out_results_area.plt_hard_trunc_m = plt_hard_trunc_m;
        out_results_area.plt_hard_trunc_ft = plt_hard_trunc_ft;
        out_results_area.t_data = ST.t_data;
        out_results_all{region_ind_base(rind)+ind_area} = out_results_area;

        clear tr_start_ind tr_end_ind sect_ind ST x_data Ntr ...
            out_mat_m out_mat_ft plt_dat_m plt_dat_ft...
            plt_mod_m plt_mod_ft plt_err_m plt_err_ft...
            plt_trunc_m plt_trunc_ft plt_hard_m plt_hard_ft...
            plt_hard_trunc_m plt_hard_trunc_ft out_results_area...
            nrow mat_col_area_name mat_col_island mat_col_GISID

    end
    clear resmat

end % will make output matrices and plot data for all map areas and regions





% ---------------------------------------------------------------------------
% Write results to csv files. 
% Headers have the following format (as requested by Jade for integration with ArcONLINE
%   Transect SRate_ft Suncert_ft Sstd_ft Num DegF AreaName Island
%	Transect SRate_m  Suncert_m  Sstd_m  Num DegF AreaName Island
if write_tables

    out_path = ['AA_' island '_results/tables/'];
    if ~exist(out_path,'dir')
        mkdir(out_path)
    end
    
    fname_ft = [out_path island '_rates_feet.csv'];
    fname_m = [out_path island '_rates_meters.csv'];

    
    hdr_m = {'Transect','sRate_m','sUncert_m','sStd_m', ...
        'Num','DegF','AreaName','Island','GISID'};
    hdr_ft = {'Transect','sRate_ft','sUncert_ft','sStd_ft', ...
        'Num','DegF','AreaName','Island','GISID'};
    fmt_hdr = [repmat('%s,', 1, length(hdr_m)) '\n']; % file formatting for header
    fmt_dat = '%u,%4.6f,%4.6f,%4.6f,%u,%u,%s,%s,%s\n'; % file formatting for data results

    
    % write files for all map areas on entire island
    tmp = cellfun(@(x) x.('mat_meters'), out_results_all,'UniformOutput',false);
    out_mat_meters = vertcat(tmp{:})'; clear tmp
    tmp = cellfun(@(x) x.('mat_feet'), out_results_all,'UniformOutput',false);
    out_mat_feet = vertcat(tmp{:})'; clear tmp

    fid_m = fopen(fname_m, 'w');
    fprintf(fid_m, fmt_hdr, hdr_m{:});
    fprintf(fid_m, fmt_dat, out_mat_meters{:});
    fclose(fid_m);

    fid_ft = fopen(fname_ft, 'w');
    fprintf(fid_ft, fmt_hdr, hdr_ft{:});
    fprintf(fid_ft, fmt_dat, out_mat_feet{:}); 
    fclose(fid_ft);

end




if transect_plots_multiple_meters

    % -----------------------------------
    % METERS
    % plot data and trend at each transect 
    % multiple plots per page

    plotdims = [5,3]; % dimensions of tiled plots [row,col]
    

    for ind = 1:sum(nmap_all)
   
        region = out_results_all{ind}.region;
        area_name = out_results_all{ind}.area_name;
        t_data = out_results_all{ind}.t_data;
        x_data = out_results_all{ind}.plt_dat_m(:,1);
        Ntr = length(x_data);
        Ntime = length(t_data);
        y_data_m = out_results_all{ind}.plt_dat_m(:,2:end);
        y_mod_m = out_results_all{ind}.plt_mod_m(:,2:end);
        y_err_m = out_results_all{ind}.plt_err_m(:,2:end);
        y_hard_m = out_results_all{ind}.plt_hard_m(:,2:end);
        y_hard_trunc_m = out_results_all{ind}.plt_hard_trunc_m(:,2:end);

        % calculate weighted mean y-positions (weighted by squared error)
        nan_inds = isnan(y_data_m);     % logical matrix = 1 where nan (includes truncation)
        nan_mat = ones(size(nan_inds));
        nan_mat(nan_inds) = NaN;        % matrix of ones, except nan where truncated or missing data.
        y_mean_m = nansum(y_data_m./y_err_m.^2,2)./nansum(y_err_m.^(-2).*nan_mat,2);

        % shift y-data so that it is relative to the weighted mean y-positions
        y_data_ctr_m = y_data_m - repmat(y_mean_m,1,Ntime);
        y_mod_ctr_m = y_mod_m - repmat(y_mean_m,1,Ntime);
        y_hard_ctr_m = y_hard_m - repmat(y_mean_m,1,Ntime);
        y_hard_trunc_ctr_m = y_hard_trunc_m - repmat(y_mean_m,1,Ntime);

        % set modeled points to zero where truncated. 
        trunc_keep_mat_m = nan(size(y_data_m));
        for ii = 1:length(x_data)
            istart = find(~isnan(y_data_m(ii,:)),1,'first');
            iend = find(~isnan(y_data_m(ii,:)),1,'last');
            trunc_keep_mat_m(ii,istart:iend) = 1;
        end
        y_mod_ctr_m_trunc = y_mod_ctr_m.*trunc_keep_mat_m;


        % set time min/max and position(including errorbar) min/max for this study area (axis limits)
        tmin = min(t_data);
        tmax = max(t_data);

        % all data values +/- error in this section
        tmp_dat1_m = y_data_ctr_m + y_err_m;
        tmp_dat2_m = y_data_ctr_m - y_err_m;

        ymin = min([nanmin(nanmin(tmp_dat1_m)) nanmin(nanmin(tmp_dat2_m))]);
        ymax = max([nanmax(nanmax(tmp_dat1_m)) nanmax(nanmax(tmp_dat2_m))]);
        clear tmp_dat1_m tmp_dat2_m 

        
        % create folders to store plots in (if not already made)
        if ~exist(sprintf('%s_results/tr_plots_METERS_%s/png/',region,area_name),'dir')
            mkdir(sprintf('%s_results/tr_plots_METERS_%s/png/',region,area_name))
        end
        if ~exist(sprintf('%s_results/tr_plots_METERS_%s/eps/',region,area_name),'dir')
            mkdir(sprintf('%s_results/tr_plots_METERS_%s/eps/',region,area_name))
        end
        if ~exist(sprintf('%s_results/tr_plots_METERS_%s/pdf/',region,area_name),'dir')
            mkdir(sprintf('%s_results/tr_plots_METERS_%s/pdf/',region,area_name))
        end


        num_pages = 1;
        fignum = 1;


        % loop plots data and trend at each transect for this map area
        for tr_ind = 1:Ntr

            % output to screen
            fprintf(1,'Plotting (multiple) in meters: %s (section %d/%d), Transect# %d/%d.\n',...
                area_name, ind, size(out_results_all,1), x_data(tr_ind,1), x_data(Ntr,1));


            if fignum ==1 % create new page
                figure('position',[.25*72 0*72 8.5*72 11*72],'Visible','Off')
                set(gcf,'Paperpositionmode','manual')
                set(gcf,'PaperUnits','inches')
                set(gcf,'Paperposition',[.25 0 8.5 11])

                % display map area name on page
                % normalized units [left,bottom,width,height]
                axname = axes('Position',[0.1,0.98,0.7,0.1], 'Visible','off'); 
                text(0,0,[area_name ' - smoothed shoreline change rates'])
            end



            subplot(plotdims(1),plotdims(2),fignum)
            hold on


            % plot data relative to weighted mean y-positions
            % METERS
            plot(t_data, y_data_ctr_m(tr_ind,:), ...  % plot data
                'r+','MarkerSize',6,'linewidth',1)
            plot(repmat(t_data',2,1), ...             % plot error bars
                [y_data_ctr_m(tr_ind,:) + y_err_m(tr_ind,:);...
                 y_data_ctr_m(tr_ind,:) - y_err_m(tr_ind,:)],...
                'c','linewidth',1)
            plot(t_data, y_hard_trunc_ctr_m(tr_ind,:), ... % plot hard shorelines (only those not already plotted; so, truncated)
                'r+','MarkerSize',6,'linewidth',1)   
            plot(repmat(t_data',2,1), ...             % errorbars for hard shorelines (only those not already plotted; so, truncated)
                [y_hard_trunc_ctr_m(tr_ind,:) + y_err_m(tr_ind,:);...
                 y_hard_trunc_ctr_m(tr_ind,:) - y_err_m(tr_ind,:)],...
                 'c','linewidth',1)
            plot(t_data, y_hard_ctr_m(tr_ind,:), ...  % plot squares around hardened shoreline data (all hard shorelines - including first)
                'sb','MarkerSize',7,'linewidth',0.75)
            plot(t_data, y_mod_ctr_m_trunc(tr_ind,:), ...   % plot model (trend)
                'b','linewidth',1)



            % set axis limits 
            ymult = 5; % multiple of desired axis (example, if ymult = 5, and ymin = -37, then set y-axis min to closet (lower) multiple of 5 (-40))
            tmult = 10; % time (x-axis) will be closest multiple of 10. 
            axis([floor(tmin/tmult)*tmult,ceil(tmax/tmult)*tmult,...
                floor(ymin/ymult)*ymult,ceil(ymax/ymult)*ymult])
            v=axis;
            xlabel('yr','fontsize',8)
            ylabel('m','fontsize',8)
            pm = char(177); % plus/minus symbol
            text(0.03, 0.92,...
                ['\color{blue}' sprintf('Smoothed rate = %4.2f %c %4.2f m/yr',...
                cell2mat(out_results_all{ind}.mat_meters(tr_ind,2)),pm,...
                cell2mat(out_results_all{ind}.mat_meters(tr_ind,3)))],...
                'fontsize',6,'Units','normalized')

            tr_number = x_data(tr_ind,1);
            title_str = ['Transect #',num2str(tr_number)];
            if sum(~isnan(y_hard_ctr_m(tr_ind,:)))>0
                title([title_str '*'],'fontsize',10);
            else
                title(title_str,'fontsize',10);
            end



            if fignum==prod(plotdims) || tr_ind==Ntr
                eval(sprintf('print -dpng %s_results/tr_plots_METERS_%s/png/%s_METERS_%02.0f',region,area_name,area_name,num_pages));
                eval(sprintf('print -depsc %s_results/tr_plots_METERS_%s/eps/%s_METERS_%02.0f',region,area_name,area_name,num_pages));
                eval(sprintf('print -dpdf %s_results/tr_plots_METERS_%s/pdf/%s_METERS_%02.0f',region,area_name,area_name,num_pages));
                close
                % save and close figure
                if tr_ind==Ntr
                    break
                else
                    fignum = 1;
                    num_pages = num_pages + 1;
                end
            else
                fignum = fignum + 1;
            end

        end



    end

    fprintf(1,'Done plotting (multiple) meters.\n');
end





if transect_plots_multiple_feet

    % -----------------------------------
    % FEET
    % plot data and trend at each transect 
    % multiple plots per page

    plotdims = [5,3]; % dimensions of tiled plots [row,col]


    for ind = 1:sum(nmap_all)
        
        region = out_results_all{ind}.region;
        area_name = out_results_all{ind}.area_name;
        t_data = out_results_all{ind}.t_data;
        x_data = out_results_all{ind}.plt_dat_ft(:,1);
        Ntr = length(x_data);
        Ntime = length(t_data);
        y_data_ft = out_results_all{ind}.plt_dat_ft(:,2:end);
        y_mod_ft = out_results_all{ind}.plt_mod_ft(:,2:end);
        y_err_ft = out_results_all{ind}.plt_err_ft(:,2:end);
        y_hard_ft = out_results_all{ind}.plt_hard_ft(:,2:end);
        y_hard_trunc_ft = out_results_all{ind}.plt_hard_trunc_ft(:,2:end);


        % calculate weighted mean y-positions (weighted by squared error)
        nan_inds = isnan(y_data_ft);    % logical matrix = 1 where nan (includes truncation)
        nan_mat = ones(size(nan_inds)); 
        nan_mat(nan_inds) = NaN;        % matrix of ones, except nan where truncated or missing data.
        y_mean_ft = nansum(y_data_ft./y_err_ft.^2,2)./nansum(y_err_ft.^(-2).*nan_mat,2);

        % shift y-data so that it is relative to the weighted mean y-positions
        y_data_ctr_ft = y_data_ft - repmat(y_mean_ft,1,Ntime);
        y_mod_ctr_ft = y_mod_ft - repmat(y_mean_ft,1,Ntime);
        y_hard_ctr_ft = y_hard_ft - repmat(y_mean_ft,1,Ntime);
        y_hard_trunc_ctr_ft = y_hard_trunc_ft - repmat(y_mean_ft,1,Ntime);
        
        % set modeled points to zero where truncated. 
        trunc_keep_mat_ft = nan(size(y_data_ft));
        for ii = 1:length(x_data)
            istart = find(~isnan(y_data_ft(ii,:)),1,'first');
            iend = find(~isnan(y_data_ft(ii,:)),1,'last');
            trunc_keep_mat_ft(ii,istart:iend) = 1;
        end
        y_mod_ctr_ft_trunc = y_mod_ctr_ft.*trunc_keep_mat_ft;


        % set time min/max and position(including errorbar) min/max for this study area (axis limits)
        tmin = min(t_data);
        tmax = max(t_data);

        % all data values +/- error in this section
        tmp_dat1_ft = y_data_ctr_ft + y_err_ft;
        tmp_dat2_ft = y_data_ctr_ft - y_err_ft;

        ymin = min([nanmin(nanmin(tmp_dat1_ft)) nanmin(nanmin(tmp_dat2_ft))]);
        ymax = max([nanmax(nanmax(tmp_dat1_ft)) nanmax(nanmax(tmp_dat2_ft))]);
        clear tmp_dat1_ft tmp_dat2_ft

        
        % create folders to store plots in (if not already made)
        if ~exist(sprintf('%s_results/tr_plots_FEET_%s/png/',region,area_name),'dir')
            mkdir(sprintf('%s_results/tr_plots_FEET_%s/png/',region,area_name))
        end
        if ~exist(sprintf('%s_results/tr_plots_FEET_%s/eps/',region,area_name),'dir')
            mkdir(sprintf('%s_results/tr_plots_FEET_%s/eps/',region,area_name))
        end
        if ~exist(sprintf('%s_results/tr_plots_FEET_%s/pdf/',region,area_name),'dir')
            mkdir(sprintf('%s_results/tr_plots_FEET_%s/pdf/',region,area_name))
        end



        num_pages = 1;
        fignum = 1;


        % loop plots data and trend at each transect for this map area
        for tr_ind = 1:Ntr

    %         tr_ind = 117; % moloaa transect #119 (hard shoreline)

            % output to screen
            fprintf(1,'Plotting (multiple) in feet: %s (section %d/%d), Transect# %d/%d.\n',...
                area_name, ind, size(out_results_all,1), x_data(tr_ind,1), x_data(Ntr,1));


            if fignum ==1 % create new page
                figure('position',[.25*72 0*72 8.5*72 11*72],'Visible','Off')
                set(gcf,'Paperpositionmode','manual')
                set(gcf,'PaperUnits','inches')
                set(gcf,'Paperposition',[.25 0 8.5 11])

                % display map area name on page
                % normalized units [left,bottom,width,height]
                axname = axes('Position',[0.1,0.98,0.7,0.1], 'Visible','off'); 
                text(0,0,[area_name ' - smoothed shoreline change rates'])
            end


            subplot(plotdims(1),plotdims(2),fignum)
            hold on


            % plot data relative to weighted mean y-positions
            % FEET
            plot(t_data, y_data_ctr_ft(tr_ind,:), ...  % plot data
                'r+','MarkerSize',6,'linewidth',1)
            plot(repmat(t_data',2,1), ...              % plot error bars
                [y_data_ctr_ft(tr_ind,:) + y_err_ft(tr_ind,:);...
                 y_data_ctr_ft(tr_ind,:) - y_err_ft(tr_ind,:)],...
                'c','linewidth',1)
            plot(t_data, y_hard_trunc_ctr_ft(tr_ind,:), ... % plot hard shorelines (only those not already plotted; so, truncated)
                'r+','MarkerSize',6,'linewidth',1)   
            plot(repmat(t_data',2,1), ...              % errorbars for hard shorelines (only those not already plotted; so, truncated)
                [y_hard_trunc_ctr_ft(tr_ind,:) + y_err_ft(tr_ind,:);...
                 y_hard_trunc_ctr_ft(tr_ind,:) - y_err_ft(tr_ind,:)],...
                 'c','linewidth',1)
            plot(t_data, y_hard_ctr_ft(tr_ind,:), ...  % plot squares around hardened shoreline data (all hard shorelines - including first)
                'sb','MarkerSize',7,'linewidth',0.75)
            plot(t_data, y_mod_ctr_ft_trunc(tr_ind,:), ...   % plot model (trend)
                'b','linewidth',1)



            % set axis limits 
            ymult = 5; % multiple of desired axis (example, if ymult = 5, and ymin = -37, then set y-axis min to closet (lower) multiple of 5 (-40))
            tmult = 10; % time (x-axis) will be closest multiple of 10. 
            axis([floor(tmin/tmult)*tmult,ceil(tmax/tmult)*tmult,...
                floor(ymin/ymult)*ymult,ceil(ymax/ymult)*ymult])
            v=axis;
            xlabel('yr','fontsize',8)
            ylabel('ft','fontsize',8)
            pm = char(177); % plus/minus symbol
            text(0.03, 0.92,...
                ['\color{blue}' sprintf('Smoothed rate = %4.2f %c %4.2f ft/yr',...
                cell2mat(out_results_all{ind}.mat_feet(tr_ind,2)),pm,...
                cell2mat(out_results_all{ind}.mat_feet(tr_ind,3)))],...
                'fontsize',6,'Units','normalized')

            tr_number = x_data(tr_ind,1);
            title_str = ['Transect #',num2str(tr_number)];
            if sum(~isnan(y_hard_ctr_ft(tr_ind,:)))>0
                title([title_str '*'],'fontsize',10);
            else
                title(title_str,'fontsize',10);
            end



            if fignum==prod(plotdims) || tr_ind==Ntr
                eval(sprintf('print -dpng %s_results/tr_plots_FEET_%s/png/%s_FEET_%02.0f',region,area_name,area_name,num_pages));
                eval(sprintf('print -depsc %s_results/tr_plots_FEET_%s/eps/%s_FEET_%02.0f',region,area_name,area_name,num_pages));
                eval(sprintf('print -dpdf %s_results/tr_plots_FEET_%s/pdf/%s_FEET_%02.0f',region,area_name,area_name,num_pages));
                close
                % save and close figure
                if tr_ind==Ntr
                    break
                else
                    fignum = 1;
                    num_pages = num_pages + 1;
                end
            else
                fignum = fignum + 1;
            end

        end



    end

    fprintf(1,'Done plotting (multiple) meters.\n');
end


% create single plots for each transect; all plots for entire island are 
% put in one folder
if transect_plots_single_meters

    % -----------------------------------
    % METERS
    % plot data and trend at each transect 
    % ONE plot per page

    out_dir = ['AA_' island '_results/plots_single/METERS/'];
    if ~exist(out_dir,'dir')
        mkdir(out_dir)
    end


    for ind = 1:sum(nmap_all)    
  
        region = out_results_all{ind}.region;
        area_name = out_results_all{ind}.area_name;
        hawaiian_name = out_results_all{ind}.hawaiian_name;
        t_data = out_results_all{ind}.t_data;
        x_data = out_results_all{ind}.plt_dat_m(:,1);
        Ntr = length(x_data);
        Ntime = length(t_data);
        y_data_m = out_results_all{ind}.plt_dat_m(:,2:end);
        y_mod_m = out_results_all{ind}.plt_mod_m(:,2:end);
        y_err_m = out_results_all{ind}.plt_err_m(:,2:end);
        y_hard_m = out_results_all{ind}.plt_hard_m(:,2:end);
        y_hard_trunc_m = out_results_all{ind}.plt_hard_trunc_m(:,2:end);

        % calculate weighted mean y-positions (weighted by squared error)
        nan_inds = isnan(y_data_m);     % logical matrix = 1 where nan (includes truncation)
        nan_mat = ones(size(nan_inds));
        nan_mat(nan_inds) = NaN;        % matrix of ones, except nan where truncated or missing data.
        y_mean_m = nansum(y_data_m./y_err_m.^2,2)./nansum(y_err_m.^(-2).*nan_mat,2);

        % shift y-data so that it is relative to the weighted mean y-positions
        y_data_ctr_m = y_data_m - repmat(y_mean_m,1,Ntime);
        y_mod_ctr_m = y_mod_m - repmat(y_mean_m,1,Ntime);
        y_hard_ctr_m = y_hard_m - repmat(y_mean_m,1,Ntime);
        y_hard_trunc_ctr_m = y_hard_trunc_m - repmat(y_mean_m,1,Ntime);

        % set modeled points to zero where truncated. 
        trunc_keep_mat_m = nan(size(y_data_m));
        for ii = 1:length(x_data)
            istart = find(~isnan(y_data_m(ii,:)),1,'first');
            iend = find(~isnan(y_data_m(ii,:)),1,'last');
            trunc_keep_mat_m(ii,istart:iend) = 1;
        end
        y_mod_ctr_m_trunc = y_mod_ctr_m.*trunc_keep_mat_m;


        % set time min/max and position(including errorbar) min/max for this study area (axis limits)
        tmin = min(t_data);
        tmax = max(t_data);

        % all data values +/- error in this section
        tmp_dat1_m = y_data_ctr_m + y_err_m;
        tmp_dat2_m = y_data_ctr_m - y_err_m;

        ymin = min([nanmin(nanmin(tmp_dat1_m)) nanmin(nanmin(tmp_dat2_m))]);
        ymax = max([nanmax(nanmax(tmp_dat1_m)) nanmax(nanmax(tmp_dat2_m))]);
        clear tmp_dat1_m tmp_dat2_m 


%             % create folders to store plots in (if not already made)
%             if ~exist(sprintf('%s_results/tr_plots_single_METERS_%s/png/',region,area_name),'dir')
%                 mkdir(sprintf('%s_results/tr_plots_single_METERS_%s/png/',region,area_name))
%             end



        num_pages = 1;
        fignum = 1;


        % loop plots data and trend at each transect for this map area
        for tr_ind = 1:Ntr

            % output to screen
            fprintf(1,'Plotting (single) in meters: %s (section %d/%d), Transect# %d/%d.\n',...
                area_name, ind, size(out_results_all,1), x_data(tr_ind,1), x_data(Ntr,1));


            % On Windows systems, a pixel is 1/96th of an inch
            fig = figure('OuterPosition',[0,0,600,450],'Visible','Off'); % size in pixels on screen
%             fig = figure('OuterPosition',[0,0,600,450]); % size in pixels on screen
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 600/96 450/96]; % 96 pixels per inch on windows machine
            fig.PaperPositionMode = 'manual';


            hold on

            % plot data relative to weighted mean y-positions
            % METERS
            hline_dat = plot(t_data, y_data_ctr_m(tr_ind,:), ...  % plot data
                'r+','MarkerSize',10,'linewidth',2);
            hline_unc = plot(repmat(t_data',2,1), ...             % plot error bars
                [y_data_ctr_m(tr_ind,:) + y_err_m(tr_ind,:);...
                 y_data_ctr_m(tr_ind,:) - y_err_m(tr_ind,:)],...
                'c','linewidth',2);
            plot(t_data, y_hard_trunc_ctr_m(tr_ind,:), ... % plot hard shorelines (only those not already plotted; so, truncated)
                'r+','MarkerSize',10,'linewidth',2);   
            plot(repmat(t_data',2,1), ...             % errorbars for hard shorelines (only those not already plotted; so, truncated)
                [y_hard_trunc_ctr_m(tr_ind,:) + y_err_m(tr_ind,:);...
                 y_hard_trunc_ctr_m(tr_ind,:) - y_err_m(tr_ind,:)],...
                 'c','linewidth',2)
            hline_box = plot(t_data, y_hard_ctr_m(tr_ind,:), ...  % plot squares around hardened shoreline data (all hard shorelines - including first)
                'sb','MarkerSize',12,'linewidth',1);
            hline_rate = plot(t_data, y_mod_ctr_m_trunc(tr_ind,:), ...   % plot model (trend)
                'b','linewidth',2);



            % set axis limits 
            ymult = 5; % multiple of desired axis (example, if ymult = 5, and ymin = -37, then set y-axis min to closet (lower) multiple of 5 (-40))
            tmult = 10; % time (x-axis) will be closest multiple of 10. 
            axis([floor(tmin/tmult)*tmult,ceil(tmax/tmult)*tmult+tmult,... % add tmult to y-axis to allow room for legend
                floor(ymin/ymult)*ymult,ceil(ymax/ymult)*ymult])
            v=axis;
            xlabel('Date (years)','fontsize',12,'FontWeight','bold')
            ylabel({'  \bf\fontsize{12}Position (m)  ' '\fontsize{12}\bf\leftarrow \fontsize{11}landward    seaward \fontsize{12}\rightarrow'})
            pm = char(177); % plus/minus symbol
%             text(0.08, 0.93,...
%                 ['\color{blue}' sprintf('Smoothed rate = %4.2f %c %4.2f m/yr (%c 95%% confidence)',...
%                 out_results_all{ind}.mat_meters(tr_ind,2),pm,...
%                 out_results_all{ind}.mat_meters(tr_ind,3),pm)],...
%                 'fontsize',12,'Units','normalized','FontWeight','bold')
            xcoord = 0.194; ycoord = 0.88; xlen = 0.06; 
            annotation('line',[xcoord,xcoord+xlen],[ycoord,ycoord],'LineWidth',2,'Color','b') % rate line
            text(0.18, 0.95,...                                                               % rate line text
                sprintf('Smoothed rate = %4.2f %c %4.2f m/yr (95%% conf.)',...
                cell2mat(out_results_all{ind}.mat_meters(tr_ind,2)),pm,...
                cell2mat(out_results_all{ind}.mat_meters(tr_ind,3))),...
                'fontsize',12,'Units','normalized','FontWeight','bold')
            yoff = 0.06;
            annotation('line',[xcoord,xcoord+xlen],[ycoord-yoff,ycoord-yoff],... % data line
                'LineWidth',2,'Color','r')
            text(0.18, 0.95-0.072,'LWMs',...                                      % data line text
                'fontsize',12,'Units','normalized','FontWeight','bold')
            xoff1 = 0.185;
            annotation('line',[xcoord+xoff1,xcoord+xoff1+xlen],[ycoord-yoff,ycoord-yoff],... % data error line
                'LineWidth',2,'Color','c')
            text(0.41, 0.95-0.072,'LWM conf. (Etp)',...                                      % data error line text
                'fontsize',12,'Units','normalized','FontWeight','bold')


            tr_number = x_data(tr_ind,1);
            title_str = [area_name ', Transect #',num2str(tr_number)];
            title_str_hawaiian = [hawaiian_name ', Transect #',num2str(tr_number)];
%             title(title_str,'fontsize',16);
            title(title_str_hawaiian,'fontsize',16);

            if sum(~isnan(y_hard_ctr_m(tr_ind,:)))>0
                text(0.72, 0.872 ,'   ','FontSize',2,...                     % hard shoreline box
                    'LineWidth',1,'EdgeColor','b','Units','normalized')
                text(0.745, 0.95-0.07,'Hard shoreline',...                   % hard shoreline box text
                    'fontsize',12,'Units','normalized','FontWeight','bold')
            end

%             legend([hline_dat,hline_unc(1),hline_box],...
%                 {'data','data errors','hard shoreline'},...
%                 'Orientation','horizontal',...
%                 'FontSize',12,'FontWeight','bold',...
%                 'Units','normalized','Position',[0.34,0.395,0.3,0.8])
%             legend('boxoff')


            text(0.78, -0.113,['Plot created: ' date],...
                'Color', [.5 .5 .5], 'Units','normalized')


%                 eval(sprintf('print -dpng -r0 %s_results/tr_plots_single_METERS_%s/png/%s_METERS_%02.0f',region,area_name,area_name,tr_number));
            eval(sprintf('print -dpng -r0 %s%s_METERS_%d',out_dir,area_name,tr_number));
            close

%             % png, add metadata for area_name, transect_num
% %                 plt_filename = sprintf('%s_results/tr_plots_single_METERS_%s/png/%s_METERS_%02.0f.png',region,area_name,area_name,tr_number);
%             plt_filename = sprintf('%s%s_METERS_%02.0f.png',out_dir,area_name,tr_number);
%             plt_info = imfinfo(plt_filename);
%             plt = imread(plt_filename);
%             imwrite(plt,plt_filename,...
%                 'AreaName',area_name,...
%                 'TransectNumber',num2str(tr_number),...
%                 'XResolution',plt_info.XResolution,...
%                 'YResolution',plt_info.YResolution,...
%                 'ResolutionUnit',plt_info.ResolutionUnit);
%             clear plt_filename plt_info plt


        end



    end    

    clear out_dir
    fprintf(1,'Done plotting (single) meters.\n');
end




if transect_plots_single_feet

    % -----------------------------------
    % FEET
    % plot data and trend at each transect 
    % ONE plots per page
    
    out_dir = ['AA_' island '_results/plots_single/FEET/'];
    if ~exist(out_dir,'dir')
        mkdir(out_dir)
    end


    for ind = 1:sum(nmap_all)   
    
        region = out_results_all{ind}.region;
        area_name = out_results_all{ind}.area_name;
        hawaiian_name = out_results_all{ind}.hawaiian_name;
        t_data = out_results_all{ind}.t_data;
        x_data = out_results_all{ind}.plt_dat_ft(:,1);
        Ntr = length(x_data);
        Ntime = length(t_data);
        y_data_ft = out_results_all{ind}.plt_dat_ft(:,2:end);
        y_mod_ft = out_results_all{ind}.plt_mod_ft(:,2:end);
        y_err_ft = out_results_all{ind}.plt_err_ft(:,2:end);
        y_hard_ft = out_results_all{ind}.plt_hard_ft(:,2:end);
        y_hard_trunc_ft = out_results_all{ind}.plt_hard_trunc_ft(:,2:end);


        % calculate weighted mean y-positions (weighted by squared error)
        nan_inds = isnan(y_data_ft);    % logical matrix = 1 where nan (includes truncation)
        nan_mat = ones(size(nan_inds)); 
        nan_mat(nan_inds) = NaN;        % matrix of ones, except nan where truncated or missing data.
        y_mean_ft = nansum(y_data_ft./y_err_ft.^2,2)./nansum(y_err_ft.^(-2).*nan_mat,2);

        % shift y-data so that it is relative to the weighted mean y-positions
        y_data_ctr_ft = y_data_ft - repmat(y_mean_ft,1,Ntime);
        y_mod_ctr_ft = y_mod_ft - repmat(y_mean_ft,1,Ntime);
        y_hard_ctr_ft = y_hard_ft - repmat(y_mean_ft,1,Ntime);
        y_hard_trunc_ctr_ft = y_hard_trunc_ft - repmat(y_mean_ft,1,Ntime);
        

        % set modeled points to zero where truncated. 
        trunc_keep_mat_ft = nan(size(y_data_ft));
        for ii = 1:length(x_data)
            istart = find(~isnan(y_data_ft(ii,:)),1,'first');
            iend = find(~isnan(y_data_ft(ii,:)),1,'last');
            trunc_keep_mat_ft(ii,istart:iend) = 1;
        end
        y_mod_ctr_ft_trunc = y_mod_ctr_ft.*trunc_keep_mat_ft;


        % set time min/max and position(including errorbar) min/max for this study area (axis limits)
        tmin = min(t_data);
        tmax = max(t_data);

        % all data values +/- error in this section
        tmp_dat1_ft = y_data_ctr_ft + y_err_ft;
        tmp_dat2_ft = y_data_ctr_ft - y_err_ft;

        ymin = min([nanmin(nanmin(tmp_dat1_ft)) nanmin(nanmin(tmp_dat2_ft))]);
        ymax = max([nanmax(nanmax(tmp_dat1_ft)) nanmax(nanmax(tmp_dat2_ft))]);
        clear tmp_dat1_ft tmp_dat2_ft


%             % create folders to store plots in (if not already made)
%             if ~exist(sprintf('%s_results/tr_plots_single_FEET_%s/png/',region,area_name),'dir')
%                 mkdir(sprintf('%s_results/tr_plots_single_FEET_%s/png/',region,area_name))
%             end



        num_pages = 1;
        fignum = 1;



        % loop plots data and trend at each transect for this map area
        for tr_ind = 1:Ntr

            % output to screen
            fprintf(1,'Plotting (single) in feet: %s (section %d/%d), Transect# %d/%d.\n',...
                area_name, ind, size(out_results_all,1), x_data(tr_ind,1), x_data(Ntr,1));


            % On Windows systems, a pixel is 1/96th of an inch
            fig = figure('OuterPosition',[0,0,600,450],'Visible','Off'); % size in pixels on screen
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 600/96 450/96]; % 96 pixels per inch on windows machine
            fig.PaperPositionMode = 'manual';

            hold on


            % plot data relative to weighted mean y-positions
            % FEET
            hline_dat = plot(t_data, y_data_ctr_ft(tr_ind,:), ...  % plot data (note: line handles were for legend, that I didn't use because I couldn't format it well enough)
                'r+','MarkerSize',10,'linewidth',2);
            hline_unc = plot(repmat(t_data',2,1), ...              % plot error bars
                [y_data_ctr_ft(tr_ind,:) + y_err_ft(tr_ind,:);...
                 y_data_ctr_ft(tr_ind,:) - y_err_ft(tr_ind,:)],...
                'c','linewidth',2);
            plot(t_data, y_hard_trunc_ctr_ft(tr_ind,:), ... % plot hard shorelines (only those not already plotted; so, truncated)
                'r+','MarkerSize',10,'linewidth',2);   
            plot(repmat(t_data',2,1), ...             % errorbars for hard shorelines (only those not already plotted; so, truncated)
                [y_hard_trunc_ctr_ft(tr_ind,:) + y_err_ft(tr_ind,:);...
                 y_hard_trunc_ctr_ft(tr_ind,:) - y_err_ft(tr_ind,:)],...
                 'c','linewidth',2)
            hline_box = plot(t_data, y_hard_ctr_ft(tr_ind,:), ...  % plot squares around hardened shoreline data (all hard shorelines - including first)
                'sb','MarkerSize',12,'linewidth',1);
            hline_rate = plot(t_data, y_mod_ctr_ft_trunc(tr_ind,:), ...   % plot model (trend)
                'b','linewidth',2);


            % set axis limits 
            ymult = 5; % multiple of desired axis (example, if ymult = 5, and ymin = -37, then set y-axis min to closet (lower) multiple of 5 (-40))
            tmult = 10; % time (x-axis) will be closest multiple of 10. 
            axis([floor(tmin/tmult)*tmult,ceil(tmax/tmult)*tmult+tmult,... % add tmult to y-axis to allow room for legend
                floor(ymin/ymult)*ymult,ceil(ymax/ymult)*ymult])
            v=axis;
            xlabel('Date (years)','fontsize',12,'FontWeight','bold')
            ylabel({'  \bf\fontsize{12}Position (ft)  ' '\fontsize{12}\bf\leftarrow \fontsize{11}landward    seaward \fontsize{12}\rightarrow'})
            pm = char(177); % plus/minus symbol
%             text(0.08, 0.93,...
%                 ['\color{blue}' sprintf('Smoothed rate = %4.1f %c %4.1f ft/yr (%c 95%% confidence)',...
%                 out_results_all{ind}.mat_feet(tr_ind,2),pm,...
%                 out_results_all{ind}.mat_feet(tr_ind,3),pm)],...
%                 'fontsize',12,'Units','normalized','FontWeight','bold')
            xcoord = 0.194; ycoord = 0.88; xlen = 0.06; 
            annotation('line',[xcoord,xcoord+xlen],[ycoord,ycoord],'LineWidth',2,'Color','b') % rate line
            text(0.18, 0.95,...                                                               % rate line text
                sprintf('Smoothed rate = %4.1f %c %4.1f ft/yr (95%% conf.)',...
                cell2mat(out_results_all{ind}.mat_feet(tr_ind,2)),pm,...
                cell2mat(out_results_all{ind}.mat_feet(tr_ind,3))),...
                'fontsize',12,'Units','normalized','FontWeight','bold')
            yoff = 0.06;
            annotation('line',[xcoord,xcoord+xlen],[ycoord-yoff,ycoord-yoff],... % data line
                'LineWidth',2,'Color','r')
            text(0.18, 0.95-0.072,'LWMs',...                                      % data line text
                'fontsize',12,'Units','normalized','FontWeight','bold')
            xoff1 = 0.185;
            annotation('line',[xcoord+xoff1,xcoord+xoff1+xlen],[ycoord-yoff,ycoord-yoff],... % data error line
                'LineWidth',2,'Color','c')
            text(0.41, 0.95-0.072,'LWM conf. (Etp)',...                                      % data error line text
                'fontsize',12,'Units','normalized','FontWeight','bold')


            tr_number = x_data(tr_ind,1);
            title_str = [area_name ', Transect #',num2str(tr_number)];
            title_str_hawaiian = [hawaiian_name ', Transect #',num2str(tr_number)];
%             title(title_str,'fontsize',16);
            title(title_str_hawaiian,'fontsize',16);
            
            if sum(~isnan(y_hard_ctr_ft(tr_ind,:)))>0
                text(0.72, 0.872 ,'   ','FontSize',2,...                     % hard shoreline box
                    'LineWidth',1,'EdgeColor','b','Units','normalized')
                text(0.745, 0.95-0.07,'Hard shoreline',...                   % hard shoreline box text
                    'fontsize',12,'Units','normalized','FontWeight','bold')
            end

%             legend([hline_dat,hline_unc(1),hline_box],...
%                 {'data','data errors','hard shoreline'},...
%                 'Orientation','horizontal',...
%                 'FontSize',12,'FontWeight','bold',...
%                 'Units','normalized','Position',[0.34,0.395,0.3,0.8])
%             legend('boxoff')


            text(0.78, -0.113,['Plot created: ' date],...
                'Color', [.5 .5 .5], 'Units','normalized')


%                 eval(sprintf('print -dpng -r0 %s_results/tr_plots_single_FEET_%s/png/%s_FEET_%02.0f',region,area_name,area_name,tr_number));
            eval(sprintf('print -dpng -r0 %s%s_FEET_%d',out_dir,area_name,tr_number));
            close

%             % png, add metadata for area_name, transect_num
%             plt_filename = sprintf('%s%s_FEET_%02.0f.png',out_dir,area_name,tr_number);
%             plt_info = imfinfo(plt_filename);
%             plt = imread(plt_filename);
%             imwrite(plt,plt_filename,...
%                 'AreaName',area_name,...
%                 'TransectNumber',num2str(tr_number),...
%                 'XResolution',plt_info.XResolution,...
%                 'YResolution',plt_info.YResolution,...
%                 'ResolutionUnit',plt_info.ResolutionUnit);
%             clear plt_filename plt_info plt


        end



    end

    clear out_dir
    fprintf(1,'Done plotting (single) feet.\n');

end



     



