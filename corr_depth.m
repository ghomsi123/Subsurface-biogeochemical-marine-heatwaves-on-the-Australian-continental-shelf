% Julia Araujo - FEB2026
% Correlation between MHW depth extent, DCM depth, maximum stratification
% depth and thermocline depth

clear

%% Parameters
file='C:\Users\Julia\Documents\projects_personal\mhw_australia\data\gliders_clean_vf.mat';

% Regions
regions={'TAS','NSW','QLD','SW_WA'};

markers.TAS='v';
markers.NSW='hexagram';
markers.QLD='^';
markers.SW_WA='o';

% Seasons
seasons={'winter','spring','summer','autumn'};
seas_months=[6 7 8;
             9 10 11;
             12 1 2;
             3 4 5];

% colors.winter=[39 139 185];
% colors.spring=[52 173 108];
% colors.summer=[235 72 51];
% colors.autumn=[246 163 19];

colors.winter=[39 139 185];
colors.spring=[246 163 19];
colors.summer=[235 72 51];
colors.autumn=[200 200 200];

% Percentage of vertical cumulative heat
alfa=.90;

%% Load data
load(file)

%% Define time
time=datetime(1950,1,1)+days(data.TIME);
time_months=month(time);

%% Calculate baseline for 'seasonal mean' profile
baseline=data.SEVERITY<=1;

%% Calculate stratification (N²)
[N2,N2_p]=gsw_Nsquared(data.TEMP,data.PSAL,data.PRES);
N2_depth=-gsw_z_from_p(N2_p,data.LATITUDE);

% Interpolate
for i=1:size(N2_depth,2)
    is_nan=isnan(N2_depth(:,i));
    if sum(~is_nan)>1
        data.N2(:,i)=interp1(N2_depth(~is_nan,i),N2(~is_nan,i),data.DEPTH,'linear','extrap');
    else
        data.N2(:,i)=NaN;
    end
    data.N2(is_nan,i)=NaN;
end

data.N2=-data.N2*10e3;

%% Mask profiles according to region and respective depth range
data.TEMP(data.DEPTH>90,ismember(data.REGION,{'TAS','NSW'}))=NaN;
data.TEMP(data.DEPTH>40,ismember(data.REGION,{'QLD'}))=NaN;
data.TEMP(data.DEPTH>30,ismember(data.REGION,{'SW_WA'}))=NaN;

data.N2(data.DEPTH>85,ismember(data.REGION,{'TAS','NSW'}))=NaN;
data.N2(data.DEPTH>40,ismember(data.REGION,{'QLD'}))=NaN;
data.N2(data.DEPTH>30,ismember(data.REGION,{'SW_WA'}))=NaN;

%% Loop for regions and seasons
for i=1:length(regions)
    r=regions{i};
    i_region=ismember(data.REGION,r);

    % Stratified profiles
    N2_max=max(data.N2(:,i_region & data.SEVERITY>1),[],1,'omitnan');
    N2_p=prctile(N2_max,75);
    i_strat=(max(data.N2,[],1,'omitnan')>N2_p)';

    for ii=1:length(seasons)
        s=seasons{ii};
        i_season=ismember(time_months,seas_months(ii,:));

        idx=i_region & i_season & i_strat & data.SEVERITY>1;
        if sum(idx)<=1
            continue
        end

        %% Calculate baseline profiles
        prof_b=mean(data.TEMP(:,baseline & i_region & i_season),2,'omitnan');

        %% Calculate MHW depth extent
        prof_tanom=data.TEMP(:,idx)-prof_b;
        
        % Remove profiles with no positive anomalies
        idx_depth=NaN(1,size(prof_tanom,2));
        for c=1:size(prof_tanom,2)
            k=find(prof_tanom(:,c)>0,1,'last');
            if ~isempty(k)
                idx_depth(c)=k;
            end
        end
        i_neg=isnan(idx_depth);
        prof_tanom(:,i_neg)=[];
        disp(num2str(sum(i_neg)))

        % Set zero to negative and NaNs for cumulative heat
        prof_tanom(prof_tanom<0)=0;
        prof_tanom(isnan(prof_tanom))=0;

        % Cumulative heat
        cum_heat=cumsum(prof_tanom,1)./sum(prof_tanom,1,'omitnan');
            
        [~,idx_depth]=min(abs(cum_heat-alfa),[],1);
        depth_mhw.(r).(s)=data.DEPTH(idx_depth);

        % Filter MHWs shallower than 5m depth
        i_shallow=depth_mhw.(r).(s)<5;
        depth_mhw.(r).(s)(i_shallow)=[];

        %% Calculate DCM depth
        chl=data.CPHL(:,idx); chl(1:6,:)=NaN;
        [~,idx_depth]=max(chl,[],1,'omitnan');
        idx_depth(i_neg)=[];
        
        depth_dcm.(r).(s)=data.DEPTH(idx_depth); 
        depth_dcm.(r).(s)(depth_dcm.(r).(s)==0)=NaN;
        depth_dcm.(r).(s)(i_shallow)=[];

        %% Calculate maximum stratification depth
        n2=data.N2(:,idx); n2(1:6,:)=NaN;
        [~,idx_depth]=max(n2,[],1,'omitnan');
        idx_depth(i_neg)=[];
        
        depth_strat.(r).(s)=data.DEPTH(idx_depth); 
        depth_strat.(r).(s)(depth_strat.(r).(s)==0)=NaN;
        depth_strat.(r).(s)(i_shallow)=[];

        %% Calculate thermocline depth
        [~,prof_dTdz]=gradient(data.TEMP(:,idx));
        prof_dTdz(1:5,:)=NaN;
        [~,idx_depth]=min(prof_dTdz,[],1,'omitnan');
        idx_depth(i_neg)=[];
        
        depth_tcline.(r).(s)=data.DEPTH(idx_depth); 
        depth_tcline.(r).(s)(depth_tcline.(r).(s)==0)=NaN;
        depth_tcline.(r).(s)(i_shallow)=[];

    end
end

%% PLOT
figure('Units','centimeters','Position',[5 2 18 16])
tiledlayout(2,2,'TileSpacing','compact')

nexttile;
hold on
for i=1:length(regions)
    r=regions{i};

    for ii=[3 4 2 1]
        s=seasons{ii};

        if ~isfield(depth_strat.(r),s)
            continue
        end

        scatter(depth_strat.(r).(s),depth_mhw.(r).(s),15,colors.(s)/255,'filled',markers.(r), ...
                'MarkerFaceAlpha',.3,'MarkerEdgeColor','k','MarkerEdgeAlpha',.3)
    end
end

axis equal
xticks(0:20:100)
yticks(0:20:100)
xlim([0 90])
ylim([0 90])

xlabel('Maximum stratification depth (m)','FontWeight','bold')
ylabel('MHW depth extent (m)','FontWeight','bold')

box on
grid on

fontsize(gcf,8,'points')
fontname(gcf,'Helvetica')

nexttile;
hold on
for i=1:length(regions)
    r=regions{i};

    for ii=[3 4 2 1]
        s=seasons{ii};

        if ~isfield(depth_strat.(r),s)
            continue
        end

        scatter(depth_dcm.(r).(s),depth_strat.(r).(s),15,colors.(s)/255,'filled',markers.(r), ...
                'MarkerFaceAlpha',.3,'MarkerEdgeColor','k','MarkerEdgeAlpha',.3)
    end
end

axis equal
xticks(0:20:100)
yticks(0:20:100)
xlim([0 90])
ylim([0 90])

xlabel('DCM depth (m)','FontWeight','bold')
ylabel('Maximum stratification depth (m)','FontWeight','bold')

box on
grid on

fontsize(gcf,8,'points')
fontname(gcf,'Helvetica')

nexttile;
hold on
for i=1:length(regions)
    r=regions{i};

    for ii=[3 4 2 1]
        s=seasons{ii};

        if ~isfield(depth_dcm.(r),s)
            continue
        end

        scatter(depth_dcm.(r).(s),depth_mhw.(r).(s),15,colors.(s)/255,'filled',markers.(r), ...
                'MarkerFaceAlpha',.3,'MarkerEdgeColor','k','MarkerEdgeAlpha',.3)
    end
end

axis equal
xticks(0:20:100)
yticks(0:20:100)
xlim([0 90])
ylim([0 90])

xlabel('DCM depth (m)','FontWeight','bold')
ylabel('MHW depth extent (m)','FontWeight','bold')

box on
grid on

fontsize(gcf,8,'points')
fontname(gcf,'Helvetica')

nexttile;
hold on
for i=1:length(regions)
    r=regions{i};

    for ii=[3 4 2 1]
        s=seasons{ii};

        if ~isfield(depth_dcm.(r),s)
            continue
        end

        scatter(depth_dcm.(r).(s),depth_tcline.(r).(s),15,colors.(s)/255,'filled',markers.(r), ...
                'MarkerFaceAlpha',.3,'MarkerEdgeColor','k','MarkerEdgeAlpha',.3)
    end
end

axis equal
xticks(0:20:100)
yticks(0:20:100)
xlim([0 90])
ylim([0 90])

xlabel('DCM depth (m)','FontWeight','bold')
ylabel('Thermocline depth (m)','FontWeight','bold')

box on
grid on

fontsize(gcf,8,'points')
fontname(gcf,'Helvetica')

% export_fig('C:\Users\Julia\Documents\projects_personal\mhw_australia\figs\depth_corr.png','-png','-transparent','-r360');





