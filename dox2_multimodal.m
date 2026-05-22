% Julia Araujo - FEV2026
% Clarify DOX2 multi-modal distribution

clear

%% Parameters 
file='C:\Users\Julia\Documents\projects_personal\mhw_australia\data\gliders_clean_vf.mat';

zmax.TAS=90;
zmax.NSW=90;
zmax.QLD=45;
zmax.SW_WA=35;

%% Load data
load(file)

%% Define some variables
time=datetime(1950,1,1)+days(data.TIME);
time_months=month(time);
depth=data.DEPTH;
eta_common=linspace(0,1,length(depth));

regions={'TAS','NSW','QLD','SW_WA'};

season_names={'winter','spring','summer','autumn'};
seasons.winter=[6 7 8];
seasons.spring=[9 10 11];
seasons.summer=[12 1 2];
seasons.autumn=[3 4 5];

mhw_type={'moderate','severe'};
mhw.moderate=data.SEVERITY>1 & data.SEVERITY<=2;
mhw.severe=data.SEVERITY>2 & data.SEVERITY<=3;

%% Mask profiles according to region and respective depth range
data.TEMP(data.DEPTH>90,ismember(data.REGION,{'TAS','NSW'}))=NaN;
data.DOX2(data.DEPTH>90,ismember(data.REGION,{'TAS','NSW'}))=NaN;
data.TEMP(data.DEPTH>45,ismember(data.REGION,{'QLD'}))=NaN;
data.DOX2(data.DEPTH>45,ismember(data.REGION,{'QLD'}))=NaN;
data.TEMP(data.DEPTH>35,ismember(data.REGION,{'SW_WA'}))=NaN;
data.DOX2(data.DEPTH>35,ismember(data.REGION,{'SW_WA'}))=NaN;

%% Take average profiles for regions and seasons
for i=1:4
    r=regions{i};
    s=season_names{i};

    for ii=1:2
        m=mhw_type{ii};

        %% Regions
        idx_regions=ismember(data.REGION,r) & mhw.(m);
        dummy_dox2=data.DOX2(:,idx_regions);
        dummy_temp=data.TEMP(:,idx_regions);

        % Calculate MLD and rescale depth
        dox2_i=NaN(size(dummy_dox2));

        for p=1:size(dummy_temp,2)
            dummy_prof=dummy_temp(:,p);
            if sum(~isnan(dummy_prof))>5
                
                dummy_mld=mld(depth(~isnan(dummy_prof)),dummy_prof(~isnan(dummy_prof)),'metric','threshold','refpres',0);
                
                idx_above=depth<=dummy_mld;
                idx_below=depth>dummy_mld;

                eta=NaN(size(depth));
                eta(idx_above)=0.5*depth(idx_above)/dummy_mld;
                eta(idx_below)=0.5+0.5*(depth(idx_below)-dummy_mld)/(zmax.(r)+1-dummy_mld);
                
                dox2_i(:,p)=interp1(eta,dummy_dox2(:,p),eta_common,'linear','extrap');
            else
                dox2_i(:,p)=NaN;
            end
        end

        profs.(r)(:,ii)=mean(dox2_i,2,'omitnan');

        clear dummy_temp dummy_dox2 dummy_prof dummy_mld idx_above idx_below p eta

        %% Seasons
        idx_seasons=ismember(time_months,seasons.(s)) & mhw.(m);
        dummy_dox2=data.DOX2(:,idx_seasons);
        dummy_temp=data.TEMP(:,idx_seasons);

        % Calculate MLD and rescale depth
        dox2_i=NaN(size(dummy_dox2));

        for p=1:size(dummy_temp,2)
            dummy_prof=dummy_temp(:,p);
            if sum(~isnan(dummy_prof))>5
                
                dummy_mld=mld(depth(~isnan(dummy_prof)),dummy_prof(~isnan(dummy_prof)),'metric','threshold','refpres',0);
                
                idx_above=depth<=dummy_mld;
                idx_below=depth>dummy_mld;

                eta=NaN(size(depth));
                eta(idx_above)=0.5*depth(idx_above)/dummy_mld;
                eta(idx_below)=0.5+0.5*(depth(idx_below)-dummy_mld)/(90+1-dummy_mld);
                
                dox2_i(:,p)=interp1(eta,dummy_dox2(:,p),eta_common,'linear','extrap');
            else
                dox2_i(:,p)=NaN;
            end
        end

        profs.(s)(:,ii)=mean(dox2_i,2,'omitnan');

        clear dummy_temp dummy_dox2 dummy_prof dummy_mld idx_above idx_below p eta

    end

end


%% PLOT
f=fieldnames(profs);

figure('Units','centimeters','Position',[5 5 16 3])
tiledlayout(1,4,'TileSpacing','compact')
for i=1:length(regions)
    r=regions{i};    
    nexttile

    plot(profs.(r),eta_common,'LineWidth',.8)
    vfill(175,180,[.9 .9 .9],'EdgeColor','none','bottom')
    vfill(215,220,[.9 .9 .9],'EdgeColor','none','bottom')
    colororder([hex2rgb('#e2a816');hex2rgb('#bf6610')])

    set(gca,'YDir','reverse')
    xlim([130 260])
    yticklabels('')
    hline(.5,'Color',[.7 .7 .7])
    
    grid on
    box on

    fontsize(gcf,7,'points')
    fontname(gcf,'Helvetica')

end
export_fig('C:\Users\Julia\Documents\projects_personal\mhw_australia\figs\dox2_regions_rescaled.png','-png','-transparent','-r360');

figure('Units','centimeters','Position',[5 5 16 3])
tiledlayout(1,4,'TileSpacing','compact')
for i=1:length(season_names)
    s=season_names{i};    
    nexttile

    plot(profs.(s),eta_common,'LineWidth',.8)
    vfill(175,180,[.9 .9 .9],'EdgeColor','none','bottom')
    vfill(215,220,[.9 .9 .9],'EdgeColor','none','bottom')
    colororder([hex2rgb('#e2a816');hex2rgb('#bf6610')])

    set(gca,'YDir','reverse')
    xlim([130 260])
    yticklabels('')
    hline(.5,'Color',[.7 .7 .7])

    grid on
    box on

    fontsize(gcf,7,'points')
    fontname(gcf,'Helvetica')

end
export_fig('C:\Users\Julia\Documents\projects_personal\mhw_australia\figs\dox2_seasons_rescaled.png','-png','-transparent','-r360');


