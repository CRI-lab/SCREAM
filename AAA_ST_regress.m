function [ST, figs_out] = AAA_ST_regress(data_orig, bounds, flag_df_one_add_one)

% 4/3/2017 T.Anderson. Modified so that rates are POSITIVE if advancing,
% and NEGATIVE if retreating (inland).  This is the opposite of what was 
% done in the past. 
% Also added results in feet (in addition to meters) and plots in feet. 
%
% 3/9/2018 T.Anderson.  Modified so that the weighted time mean is
% subtracted from each time, for each transect.  This makes the intercept
% and rate variances independent (no cross-correlation).  It doesn't really
% change the variance of the rate, but it does change the variance of the
% intercept.  
%
% 9/17/2018 T.Anderson.  Modified to include flag 'flag_df_one_add_one',
% which, if true, adds one degree of freedom to the t-distribution 
% when there are only three shorelines (which would give one degree of 
% freedom).  This is done to keep the uncertainty reasonable and usable 
% (not so huge) by assuming that the uncertainty at these transects is 
% similar to that of neighboring transects. It is sensible that the 
% shoreline there would behave similarly to the surrounding area. Looking 
% at the rates, this appears true. Also, notice that sparse data (only 
% three) is rare, does not span extended alongshore distances, and data are 
% highly correlated in the alongshore direction (so uncertainty between 
% neighboring transects is likely to be silimar). 


    m2ft = 3.28084; % conversion factor for meters to feet

    % Determine transects with < three data points (indices are for y_data,
    % not the entire data_orig matrix). 
    crap_tr_inds = find(sum(~isnan(data_orig(3:end,2:end)),2)<3);
    good_tr_inds = (1:size(data_orig,1)-2)';
    good_tr_inds(crap_tr_inds) = [];
    
    if isempty(crap_tr_inds)
        data_orig_cond = data_orig;
    else
        data_orig_cond = data_orig;
        data_orig_cond(crap_tr_inds+2,:) = []; % remove data at this transect (removes entire tranect)
    end
    
  
    % Define distances between transects
    dx_tr = 20; % meters

    
    % Define some variables from data input
    x_data = data_orig_cond(3:end,1); %transect numbers
    x_dist = (x_data-x_data(1)+1)*dx_tr;  %convert transect numbers to distances;first transect is 20m, 2nd is 40m, etc
    y_data = data_orig_cond(3:end,2:end); %shoreline positions (m) (distance from baseline)
    t_data = data_orig_cond(1,2:end); %survey times (year)
    I = length(x_data); %number of transects (spatial)
    J = length(t_data); %number of surveys (time)
    m_err = repmat(data_orig_cond(2,2:end),I,1); % put measurement errors in matrix m_err

    % create initial data covariance matrix
    % (diagonal sparse matrix of squared measurement errors)
    Cd_orig = spdiags(reshape(m_err.^2,J*I,1),0,J*I,J*I);
    Cd_inv_orig = inv(Cd_orig);
    Cd_orig_half = Cd_orig.^(1/2);

    
      

    % **********************************************
    % Calculate rates and parameters using ST method
    % **********************************************

    % determine weighted mean of times for each transect
    m_err2_mat_zero = m_err.^-2; % calculate each err^(-2)
    m_err2_mat_zero(isnan(y_data)) = 0; % set to zero where there is no data
    t_weighted_means  = m_err2_mat_zero*t_data'./sum(m_err2_mat_zero,2);
    
    % define time and data arrays (with and without NaN data values)
    t_data_mat = repmat(t_data,I,1);
    t_data_mat_demean = t_data_mat - repmat(t_weighted_means,1,J);
%     t_data_demean = t_data - mean(t_data);
    d_orig = reshape(y_data,J*I,1);
    nan_inds = find(isnan(d_orig));
    d_mod = d_orig; % set d_mod (d model) by removing data gaps (zeros)
    d_mod(nan_inds) = [];
    N = length(d_mod); % number of total data points 


    % define system matrix (includes NaN data)
    i = (1:length(d_orig))';
    j = repmat((1:I)',J,1);
%     sp = reshape(repmat(t_data_demean,I,1),J*I,1);
    sp = reshape(t_data_mat_demean,J*I,1);
    G_r = sparse(i,j,sp); % system matrix (portion) for rate only
    G_i = sparse(i,j,ones(size(sp))); % system matrix (portion) for intercept only
    G_orig = [G_r G_i]; % complete system matrix for rate and intercept
    clear i j sp G_r G_i 


    % remove rows in system matrix and covariance matrix, where data is nan
    G_mod = G_orig;  % set G_mod (G model) by removing rows/cols for missing data
    G_mod(nan_inds,:) = [];
    Cd0_mod = Cd_orig;
    Cd0_mod(nan_inds,:) = [];
    Cd0_mod(:,nan_inds) = [];
    Cd0_mod_inv = inv(Cd0_mod);


    % invert for ST rates and intercepts
    m_ST = (G_mod'*Cd0_mod_inv*G_mod)\G_mod'*Cd0_mod_inv*d_mod; %#ok<MINV> %top half are rates, bottom half are intercepts
    ST_r = m_ST(1:I);
    ST_int = m_ST(I+1:end);    
    ST_y = reshape(G_orig*m_ST,I,J); % make predictions for dates where no data exist (NaNs)
    ST_res = y_data - ST_y; % when subtracting NaNs, will be NaN.


    % compute covariance scaling factor (alpha)
    STi = reshape(repmat(1:I,J,1),I*J,1);
    STj = (1:I*J)';
    STs = reshape(ST_res',I*J,1);
    ST_tr = sparse(STi,STj,STs);
    clear STi STj STs
    CSTtr = spdiags(reshape(m_err.^2',I*J,1),0,I*J,I*J);
    remove_nan = isnan(reshape(ST_res',I*J,1));
    ST_tr(:,remove_nan) = [];
    CSTtr(remove_nan,:) = [];
    CSTtr(:,remove_nan) = [];
    alpha_df = (sum(~isnan(ST_res),2)-2); % degrees of freedom used to calculate alpha
    
    if flag_df_one_add_one % if alpha == 1 (three shorelines), add one degree of freedom to assume that error is close to neighboring transects with slightly more data
        alpha_one = (alpha_df == 1);
        alpha_df(alpha_one) = 2;
        if sum(alpha_one)>0
            sprintf('Adding degree of freedom to Transects: %d \n',x_data(alpha_one))
        end
    end
    
    ST_alpha = full(diag(ST_tr*inv(CSTtr)*ST_tr')./alpha_df);
    alpha_zero = find(alpha_df<1); % find df == 0 (only 2 shorelines) This shouldn't happen - checked all data files
%     ST_alpha(alpha_zero) = ST_alpha(alpha_zero-1); % Set this alpha to the one next to it. 
%     ST_alpha(alpha_zero) = nan; % set the alpha with zero df to nan because it will be infinity (can't divide by zero)
%     ST_alpha(alpha_zero) = nanmean(ST_alpha); % use the mean of alphas in this section as the alpha value for statistical calculations    
    clear CSTtr remove_nan ST_tr


    % define new scaled data covariance matrix
    Cd_hat = spdiags(repmat(ST_alpha,J,1),0,J*I,J*I)*Cd_orig;
    Cd_hat(nan_inds,:) = [];
    Cd_hat(:,nan_inds) = [];
    Cd_hat_inv = inv(Cd_hat);

    % ST model covariance matrix
    ST_Cm = inv(G_mod'*Cd_hat_inv*G_mod); %#ok<MINV>

    % ST rate and intercept variances
    ST_r_var = diag(ST_Cm(1:I,1:I));
    ST_int_var = diag(ST_Cm(I+1:end,I+1:end));

    % ST position prediction variance
    ST_var = diag(G_orig*ST_Cm*G_orig'); %#ok<MINV>
    ST_y_var = reshape(ST_var,I,J);
    
    % Number of data (shoreline positions) per transect, there were used in
    % regression
    ST_ndat = sum(~isnan(y_data),2);

    % Student's t-distribution
%     ST_df = sum(~isnan(ST_res),2)-2;
    ST_df = alpha_df;  % use the same degrees of freedome for the T-distribution later, and for that used in determining the uncertainty
    ST_df_less2 = x_data(ST_df<2);
%     ST_df(ST_df<2) = 2; % set minimum df to 2 because tinv can blow up for 
%                         % transects where there are only a few shorelines
%                         % due to seawall truncation. 
%                         % This should never happen becuase data was checked to ensure only transects with four or more data are used. 
                        
%     ST_tconf95two = tinv(percent95two,ST_df);
%     ST_tconf90two95one = tinv(percent90two95one,ST_df);
%     ST_tconf80two90one = tinv(percent80two90one,ST_df);
%     ST_tconf80one = tinv(percent80one,ST_df);

    
    % Set the rates, etc. to NaN for the transects with <3 data points, and
    % also set the rates, etc. to NaN for "smoothed" data results.
    % Note: No smoothing when transects are removed. 
    if ~isempty(crap_tr_inds)
    
        sprintf(['Rates will not be smoothed. %d out of %d transects '...
            'have less than 3 data points.\n'], ...
            length(crap_tr_inds), I+length(crap_tr_inds))
        
        % set non-smooth rates, etc. to NaN for transects with sparse data,
        % and insert those values into the rate, etc. arrays.
        x_data = data_orig(3:end,1); %transect numbers
        y_data = data_orig(3:end,2:end); %shoreline positions (m) (distance from baseline)
        t_data = data_orig(1,2:end); %survey times (year)
        m_err = repmat(data_orig(2,2:end),I,1); % put measurement errors in matrix m_err
        [I,J] = size(y_data);
        nbounds = size(bounds,1)/2;
        
        ST_r_alt = nan(I,1);
        ST_r_alt(good_tr_inds) = ST_r;
        ST_r_var_alt = nan(I,1);
        ST_r_var_alt(good_tr_inds) = ST_r_var;
        ST_int_alt = nan(I,1);
        ST_int_alt(good_tr_inds) = ST_int;
        ST_int_var_alt = nan(I,1);
        ST_int_var_alt(good_tr_inds) = ST_int_var;
        
        ST_y_alt = nan(I,J);
        ST_y_var_alt = nan(I,J);
        for ind = 1:length(good_tr_inds)
            ST_y_alt(good_tr_inds(ind),:) = ST_y(ind,:);
            ST_y_var_alt(good_tr_inds(ind),:) = ST_y_var(ind,:);
        end
        
        ST_alpha_alt = nan(I,1);
        ST_alpha_alt(good_tr_inds) = ST_alpha;
        ST_res_alt = y_data - ST_y_alt; % when subtracting NaNs, will be NaN.
        alpha_df_alt = sum(~isnan(ST_res_alt),2)-2;
        alpha_df_alt(crap_tr_inds) = NaN;
        alpha_zero_alt = crap_tr_inds; % find df <= 0 (< 3 shorelines) 
        ST_df_alt = alpha_df_alt;
        ST_df_less2_alt = sort([x_data(ST_df_alt<2); x_data(crap_tr_inds)]);
        
        % rename variables for insertion into ST structure
        ST_r = ST_r_alt;
        ST_r_var = ST_r_var_alt;
        ST_int = ST_int_alt;
        ST_int_var = ST_int_var_alt;
        ST_y = ST_y_alt;
        ST_y_var = ST_y_var_alt;
        ST_alpha = ST_alpha_alt;
        ST_res = ST_res_alt;
        alpha_df = alpha_df_alt;
        alpha_zero = alpha_zero_alt; 
        ST_df = ST_df_alt;
        ST_df_less2 = ST_df_less2_alt;
        clear ST_r_alt ST_r_var_alt ST_int_alt ST_int_var_alt ST_y_alt...
            ST_y_var_alt ST_alpha_alt ST_res_alt alpha_df_alt...
            alpha_zero_alt ST_df_alt ST_df_less2_alt
    
        % set "smooth" variables to NaN
        ST_r_sm = NaN;
        ST_r_var_sm_corr  = NaN;
        ST_y_sm = NaN;
        ST_y_var_sm = NaN;
        ST_Cm_sm = NaN;
        S = NaN;
        rvar_ac_damp = NaN;
        
    else  % only smooth if there were no transects removed
        
        % Create smoothing matrix according to boundaries defined in file
        nbounds = size(bounds,1)/2;
        Snb = cell(nbounds,1);
        mat_str = [];
        for nb = 1:nbounds
            trb = bounds(nb*2-1:nb*2,1);
            numtr = find(x_data==trb(2))-find(x_data==trb(1))+1;
            Snb{nb} = make_smooth_mat(numtr);
            if isempty(mat_str)
                msg = sprintf('Snb{%d}',nb);
            else
                msg = sprintf(',Snb{%d}',nb);
            end
            mat_str = [mat_str msg];
            clear msg
        end
        S = eval(['blkdiag(' mat_str ')']);
        clear Snb mat_str nb trb numtr


        % ---------------------------------------------------------------------
        % Calculate autocorrelation of rate erros for use in calculating variance of
        % smoothed rates (to include alongshore correlation of rates so as not
        % to unfairly reduce variance in a smoothed rate). A smoothed rate is a
        % weighted average. 

        % define number of transects in the largest boundary area
        bounds_inds = nan(size(bounds,1),1);
        for k = 1:size(bounds,1)
            bounds_inds(k) = find(x_data==bounds(k));
        end
        bmax_trn = max(bounds_inds(2:2:size(bounds,1))-bounds_inds(1:2:size(bounds,1))+1);
        clear k

        % separate rate variances in to boundary areas. Put in matrix where each col is 
        % a boundary area.
        rvar_bounds = nan(bmax_trn,nbounds); % initialize structure that will hold rate variances, seperated by boundary area. (One column for each boundary area). 
        for nb = 1:nbounds
            trb = bounds(nb*2-1:nb*2,1);    % [start end] transect numbers for boundary
            bind1 = find(x_data==trb(1));   % start index of boundary
            bind2 = find(x_data==trb(2));   % end index of boundary
            numtr = bind2-bind1+1;          % number of transects in this boundary
            rvar_bounds(1:numtr,nb) = ST_r_var(bind1:bind2);
        end
        clear nb trb bind1 bind2 numtr 

        % calculate mean rate variance of each boundary area
        rates_bounds_mean = nanmean(rvar_bounds);

        % remove mean of each boundary section
        rates_bounds_nm = rvar_bounds-repmat(rates_bounds_mean,bmax_trn,1);

        % replace nans with zeros
        rates_bounds_nm(isnan(rvar_bounds)) = 0;

        % initialize vector to hold autocorrelation to zeros so the large lags with no data will be zero
        rvar_ac = zeros(I,1); 

        % calculate autocovariance of all shoreline residuals
        for k = 1:bmax_trn 

            resshift = [rates_bounds_nm(k:end,:); zeros(k-1,nbounds)];
            rvar_ac(k) = trace(rates_bounds_nm'*resshift)/sum(sum(resshift~=0))...
                *(bmax_trn-(k-1))/bmax_trn;

        end
        clear resshift k bmax_trn

        % re-weight res_ac so it's zero lag is 1
        rvar_ac = rvar_ac/rvar_ac(1);

        % damp autocovariance with cosine function so that autocorr function 
        % goes to zero at about 3/4 the total number of lags (governed by l)
        l = 6;
        damp = cos(pi*(0:I-1)/(2*(I-1))).^l;
        rvar_ac_damp = rvar_ac.*damp';
        clear damp
        % ---------------------------------------------------------------------


        % construct correlated rate var matrix. 
        r_corr_mat = toeplitz(rvar_ac_damp);
        r_var_mat = sqrt(ST_r_var)*sqrt(ST_r_var)';
        r_corr_var_mat = r_corr_mat.*r_var_mat;
        clear r_corr_mat r_var_mat 

        % create a model covariance matrix that includes the smoothed 
        % correalted rate errors.
        ST_Cm_sm = ST_Cm;
        ST_Cm_sm(1:I,1:I) = S*r_corr_var_mat*S'; % replace uncorrelated rate 
                                                 % vars with correlated vars


        % Calculate smoothed ST rates, rate variancess, ints, int var, and positions 
        % using boundaries in the boundary file.
        ST_r_sm = S*ST_r;
        ST_r_var_sm_corr = diag(S*r_corr_var_mat*S');
        ST_y_sm = reshape(G_orig*[S*ST_r;ST_int],I,J); % only rates smoothed, not intercepts
        ST_var_sm = diag(G_orig*ST_Cm_sm*G_orig'); % only rate variance correlated and smoothed 
        ST_y_var_sm = reshape(ST_var_sm,I,J);
        clear ST_var_sm
    
    end % end if ~isempty(crap_tr_inds) % smooth only if all transects have >=3 data points
   
    
    
    % %%%%%%%% LOOK HERE.  NOTICE THAT RATES ARE FLIPPED (NEGATIVE
    % %%%%%%%%             INDICATES RETREAT)
    % put results into a structure called ST
    ST = struct('data_orig',data_orig,'bounds',bounds,...
        'x_data',x_data,'t_data',t_data',...
        't_weighted_means',t_weighted_means,...
        't_data_mat_demean',t_data_mat_demean,...
        'y_data',y_data,'m_err',m_err,...
        'm_err2_mat_zero',m_err2_mat_zero,... % m_err.^(-2), and zero where there's no data
        'ndat',ST_ndat,...   % number of data used to calculate each rate
        'rate',-ST_r,'rate_var',ST_r_var,...
        'rate_sm',-ST_r_sm,'rate_var_sm',ST_r_var_sm_corr,...
        'int',ST_int,'int_var',ST_int_var,...
        'pos',ST_y,'pos_var',ST_y_var,...
        'pos_sm',ST_y_sm,'pos_var_sm',ST_y_var_sm,...
        'alpha',ST_alpha,'alpha_zero',alpha_zero,...
        'alpha_df',alpha_df,...
        'Cd',Cd_hat,'Cm',ST_Cm,'Cm_sm',ST_Cm_sm,...
        'res',ST_res,'df',ST_df,'df_less2',ST_df_less2,'S',S,...
        'rvar_ac_damp',rvar_ac_damp,...
        'y_data_ft',y_data*m2ft,'m_err_ft',m_err*m2ft, ...
        'rate_ft',-ST_r*m2ft,'rate_var_ft',ST_r_var*(m2ft^2),...
        'rate_sm_ft',-ST_r_sm*m2ft,'rate_var_sm_ft',ST_r_var_sm_corr*(m2ft^2),...
        'int_ft',ST_int*m2ft,'int_var_ft',ST_int_var*(m2ft^2),...
        'pos_ft',ST_y*m2ft,'pos_var_ft',ST_y_var*(m2ft^2),...
        'pos_sm_ft',ST_y_sm*m2ft,'pos_var_sm_ft',ST_y_var_sm*(m2ft^2));
    
    
    % -----------------------------------------------------------
    % plot data 
    
    % -----------------
    % (METERS)
    % ------------------
    fig1 = figure(1);
    clf; grid on
    plot(ST.x_data,ST.y_data)
    legend(num2str(ST.t_data))
    xlabel('Transect number')
    ylabel('Distance from offshore baseline (m)')
    
    % -----------------
    % (FEET)
    % ------------------
    fig101 = figure(101);
    clf; grid on
    plot(ST.x_data,ST.y_data_ft)
    legend(num2str(ST.t_data))
    xlabel('Transect number')
    ylabel('Distance from offshore baseline (ft)')

    
    
    % ----------------------------------------------------------------------
    % plot smoothed ST rates, non-smoothed ST rates, and unsmoothed intercepts
    
    % --------------
    % METERS
    % --------------
    
    % determine y-limits (rates and intercepts) for plotting boundaries 
    max_rate_conf = max([ST.rate + tinv(0.975,ST.df).*full(sqrt(ST.rate_var)); ...
                         ST.rate - tinv(0.975,ST.df).*full(sqrt(ST.rate_var))]);
    min_rate_conf = min([ST.rate + tinv(0.975,ST.df).*full(sqrt(ST.rate_var)); ...
                         ST.rate - tinv(0.975,ST.df).*full(sqrt(ST.rate_var))]);
    max_int_conf = max([ST.int + tinv(0.975,ST.df).*full(sqrt(ST.int_var)); ...
                         ST.int - tinv(0.975,ST.df).*full(sqrt(ST.int_var))]);
    min_int_conf = min([ST.int + tinv(0.975,ST.df).*full(sqrt(ST.int_var)); ...
                         ST.int - tinv(0.975,ST.df).*full(sqrt(ST.int_var))]);
    
    % each cell holds an array of indices corresponding to the transects of
    % that boundary segment
    binds = cell(nbounds,1); 
    for nb = 1:nbounds
        trb = bounds(nb*2-1:nb*2,1);    % [start end] transect numbers for boundary
        bind1 = find(x_data==trb(1));   % start index of boundary
        bind2 = find(x_data==trb(2));   % end index of boundary
        binds{nb} = bind1:bind2; 
    end
    clear nb trb bind1 bind2
                     
    fig10 = figure(10);
    clf
    subplot(211)
    hold on
    grid on
    for nb = 1:nbounds
        plot(x_data(binds{nb}),ST.rate(binds{nb}),'k')  % non-smoothed rates
        plot(x_data(binds{nb}),ST.rate(binds{nb}) ...
            + tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}),ST.rate(binds{nb}) ...
            - tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}(1))*ones(2,1),[max_rate_conf; min_rate_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
        plot(x_data(binds{nb}(end))*ones(2,1),[max_rate_conf; min_rate_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
    end
    if isempty(crap_tr_inds) % plot smoothed rates, if calculated.
        for nb = 1:nbounds
            plot(x_data(binds{nb}),ST.rate_sm(binds{nb}),'b') % smoothed rates
            plot(x_data(binds{nb}),ST.rate_sm(binds{nb}) ...
                + tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var_sm(binds{nb})),'b--') % 95% two tail
            plot(x_data(binds{nb}),ST.rate_sm(binds{nb}) ...
                - tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var_sm(binds{nb})),'b--') % 95% two tail
        end
    end
    xlabel('Transect')
    ylabel('Rate (m/yr)')
    legend('non-smoothed','95% conf','95% conf','boundary',...
        'smoothed (BLUE ONLY)')
    legend off
    title([int2str(t_data(1)) ' to ' int2str(t_data(end)) ...
        ', min ' int2str(nanmin(ST_df+2)) ...
        ' to max ' int2str(nanmax(ST_df+2)) ' surveys'])
    subplot(212)
    hold on
    grid on
    for nb = 1:nbounds
        plot(x_data(binds{nb}),ST.int(binds{nb}),'k')  % non-smoothed intercepts
        plot(x_data(binds{nb}),ST.int(binds{nb}) ...
            + tinv(0.975,ST.df(binds{nb})).*sqrt(ST.int_var(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}),ST.int(binds{nb}) ...
            - tinv(0.975,ST.df(binds{nb})).*sqrt(ST.int_var(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}(1))*ones(2,1),[max_int_conf; min_int_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
        plot(x_data(binds{nb}(end))*ones(2,1),[max_int_conf; min_int_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
    end
    xlabel('Transect')
    ylabel('Intercept (m)')
    legend('non-smoothed','95% conf','95% conf','boundary')
    legend off
    
    
    % --------------
    % FEET             
    % --------------
    
    % determine y-limits (rates and intercepts) for plotting boundaries 
    max_rate_conf = max([ST.rate_ft + tinv(0.975,ST.df).*full(sqrt(ST.rate_var_ft)); ...
                         ST.rate_ft - tinv(0.975,ST.df).*full(sqrt(ST.rate_var_ft))]);
    min_rate_conf = min([ST.rate_ft + tinv(0.975,ST.df).*full(sqrt(ST.rate_var_ft)); ...
                         ST.rate_ft - tinv(0.975,ST.df).*full(sqrt(ST.rate_var_ft))]);
    max_int_conf = max([ST.int_ft + tinv(0.975,ST.df).*full(sqrt(ST.int_var_ft)); ...
                         ST.int_ft - tinv(0.975,ST.df).*full(sqrt(ST.int_var_ft))]);
    min_int_conf = min([ST.int_ft + tinv(0.975,ST.df).*full(sqrt(ST.int_var_ft)); ...
                         ST.int_ft - tinv(0.975,ST.df).*full(sqrt(ST.int_var_ft))]);
    
    fig1001 = figure(1001);
    clf
    subplot(211)
    hold on
    grid on
    for nb = 1:nbounds
        plot(x_data(binds{nb}),ST.rate_ft(binds{nb}),'k')  % non-smoothed rates
        plot(x_data(binds{nb}),ST.rate_ft(binds{nb}) ...
            + tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var_ft(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}),ST.rate_ft(binds{nb}) ...
            - tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var_ft(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}(1))*ones(2,1),[max_rate_conf; min_rate_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
        plot(x_data(binds{nb}(end))*ones(2,1),[max_rate_conf; min_rate_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
    end
    if isempty(crap_tr_inds) % plot smoothed rates, if calculated.
        for nb = 1:nbounds
            plot(x_data(binds{nb}),ST.rate_sm_ft(binds{nb}),'b') % smoothed rates
            plot(x_data(binds{nb}),ST.rate_sm_ft(binds{nb}) ...
                + tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var_sm_ft(binds{nb})),'b--') % 95% two tail
            plot(x_data(binds{nb}),ST.rate_sm_ft(binds{nb}) ...
                - tinv(0.975,ST.df(binds{nb})).*sqrt(ST.rate_var_sm_ft(binds{nb})),'b--') % 95% two tail
        end
    end
    xlabel('Transect')
    ylabel('Rate (ft/yr)')
    legend('non-smoothed','95% conf','95% conf','boundary',...
        'smoothed (BLUE ONLY)')
    legend off
    title([int2str(t_data(1)) ' to ' int2str(t_data(end)) ...
        ', min ' int2str(nanmin(ST_df+2)) ...
        ' to max ' int2str(nanmax(ST_df+2)) ' surveys'])
    subplot(212)
    hold on
    grid on
    for nb = 1:nbounds
        plot(x_data(binds{nb}),ST.int_ft(binds{nb}),'k')  % non-smoothed intercepts
        plot(x_data(binds{nb}),ST.int_ft(binds{nb}) ...
            + tinv(0.975,ST.df(binds{nb})).*sqrt(ST.int_var_ft(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}),ST.int_ft(binds{nb}) ...
            - tinv(0.975,ST.df(binds{nb})).*sqrt(ST.int_var_ft(binds{nb})),'k--') % 95% two tail
        plot(x_data(binds{nb}(1))*ones(2,1),[max_int_conf; min_int_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
        plot(x_data(binds{nb}(end))*ones(2,1),[max_int_conf; min_int_conf],'color',[0.5 0.5 0.5]) % beginning of boundary
    end
    xlabel('Transect')
    ylabel('Intercept (ft)')
    legend('non-smoothed','95% conf','95% conf','boundary')
    legend off
    
    
    figs_out.data_METERS = fig1;
    figs_out.params_METERS = fig10;
    figs_out.data_FEET = fig101;
    figs_out.params_FEET = fig1001;
    
end

   
    
    




