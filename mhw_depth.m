% Julia Araujo - FEB2026
% Compare definitions of MHW depth extent

clear

%% Parameters
file='C:\Users\Julia\Documents\projects_personal\mhw_australia\data\gliders_clean_vf.mat';

% Regions
regions={'TAS','NSW','QLD','SW_WA'};

% Seasons
seasons={'winter','spring','summer','autumn'};
seas_months=[6 7 8;
             9 10 11;
             12 1 2;
             3 4 5];

% Percentage of vertical cumulative heat
alfa.p95=.95;
alfa.p90=.90;
alfa.p85=.85;

% Percentage of temperature difference
beta.p20=.2;
beta.p40=.4;

%% Load data
load(file)

%% Define time
time=datetime(1950,1,1)+days(data.TIME);
time_months=month(time);

%% Calculate baseline for 'seasonal mean' profile
baseline.all=true(size(time));
baseline.non=data.SEVERITY<=1;

%% Mask profiles according to region and respective depth range
data.TEMP(data.DEPTH>85,ismember(data.REGION,{'TAS','NSW'}))=NaN;
data.TEMP(data.DEPTH>40,ismember(data.REGION,{'QLD'}))=NaN;
data.TEMP(data.DEPTH>30,ismember(data.REGION,{'SW_WA'}))=NaN;

%% Loop for regions and seasons
for i=1:length(regions)
    r=regions{i};
    i_region=ismember(data.REGION,r);

    for ii=1:length(seasons)
        s=seasons{ii};
        i_season=ismember(time_months,seas_months(ii,:));

        idx=i_region & i_season & data.SEVERITY>1;

        %% Calculate baseline profiles
        prof_b.all=mean(data.TEMP(:,baseline.all & i_region & i_season),2,'omitnan');
        prof_b.non=mean(data.TEMP(:,baseline.non & i_region & i_season),2,'omitnan');

        %% Calculate MHW depth extent
        % METHOD A - Positive anomalies relative to seasonal mean (all profiles)
        prof_anom=data.TEMP(:,idx)-prof_b.all;
        
        idx_depth=NaN(1,size(prof_anom,2));
        for c=1:size(prof_anom,2)
            k=find(prof_anom(:,c)>0,1,'last');
            if ~isempty(k)
                idx_depth(c)=k;
            end
        end
        idx_depth(isnan(idx_depth))=[];

        mhw_depths.(r).(s).A=data.DEPTH(idx_depth);
        mhw_depth_mean.(r).(s).A=mean(mhw_depths.(r).(s).A);

        % METHOD B - Positive anomalies relative to seasonal mean (non-MHW profiles)
        prof_anom=data.TEMP(:,idx)-prof_b.non;
        
        idx_depth=NaN(1,size(prof_anom,2));
        for c=1:size(prof_anom,2)
            k=find(prof_anom(:,c)>0,1,'last');
            if ~isempty(k)
                idx_depth(c)=k;
            end
        end
        neg=isnan(idx_depth);
        idx_depth(neg)=[];

        mhw_depths.(r).(s).B=data.DEPTH(idx_depth);
        mhw_depth_mean.(r).(s).B=mean(mhw_depths.(r).(s).B);

        % METHODS C/D/E - % of vertical cumulative heat relative to seasonal mean (non-MHW profiles, alfa=0.85-0.95)
        prof_anom(:,neg)=[];
        prof_anom(prof_anom<0)=0;
        prof_anom(isnan(prof_anom))=0;
        cum_heat=cumsum(prof_anom,1)./sum(prof_anom,1,'omitnan');

        alfas=fieldnames(alfa);
        methods={'C','D','E'};
        for iii=1:length(alfas)
            a=alfas{iii};
            
            [~,idx_depth]=min(abs(cum_heat-alfa.(a)),[],1);
            mhw_depths.(r).(s).(methods{iii})=data.DEPTH(idx_depth);
            mhw_depth_mean.(r).(s).(methods{iii})=mean(mhw_depths.(r).(s).(methods{iii}));
        end
        
        % METHOD F - Shape-based metric of anomalies relative to seasonal mean (non-MHW profiles)
        prof_anom=data.TEMP(:,idx)-prof_b.non;

        prof_anom(prof_anom<0)=0;
        keep=any(prof_anom>0,1);
        prof_anom(:,keep);

        betas=fieldnames(beta);
        methods={'F','G'};
        for iii=1:length(betas)
            b=betas{iii};
            idx_depth=NaN(1,size(prof_anom,2));

            for c=1:size(prof_anom,2)
                tmax=max(prof_anom(:,c));
    
                if tmax==0 || isnan(tmax), continue; end
    
                thr=beta.(b)*tmax;        
                ok=prof_anom(:,c)>=thr;
            
                k=find(ok,1,'last');
                if ~isempty(k)
                    idx_depth(c)=k;
                end
            end
            idx_depth(isnan(idx_depth))=[];
    
            mhw_depths.(r).(s).(methods{iii})=data.DEPTH(idx_depth);
            mhw_depth_mean.(r).(s).(methods{iii})=mean(mhw_depths.(r).(s).(methods{iii}));
        end

    end
end

%% Transform structure into matrix
methods=fieldnames(mhw_depth_mean.(regions{1}).(seasons{1}));
for r=1:length(regions)
    for s=1:length(seasons)
        for m=1:length(methods)
            Z(s,m,r)=mhw_depth_mean.(regions{r}).(seasons{s}).(methods{m});
        end
    end
end

%% PLOT
%% Heatmap
figure('Units','centimeters','Position',[5 5 17 3])
tiledlayout(1,4,'TileSpacing','compact')

for r=1:length(regions)

    nexttile;
    pcolor([Z(:,:,r)-Z(:,4,r) zeros(size(Z,1),1); zeros(1,size(Z,2)+1)])
    colormap(brewermap(15,'-RdBu'))
    axis ij

    xticks(1.5:length(methods)+.5)
    xticklabels(methods)
    yticks(1:length(seasons))
    yticklabels('')

    if r==length(regions)
        colorbar
    end
    clim([-15 15])

end

fontsize(gcf,7,'points')
fontname(gcf,'Helvetica')

export_fig('C:\Users\Julia\Documents\projects_personal\mhw_australia\figs\mhw-depth_sens.png','-png','-transparent','-r360');





