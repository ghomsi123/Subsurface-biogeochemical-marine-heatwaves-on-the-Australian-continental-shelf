% Julia Araujo - AUG2025
% Plot temperature, chlorophyll and dissolved oxygen profiles for non-MHW
% and MHW events of glider profiles off the Australian coast

clear

%% Parameters
file_data='C:\Users\Julia\Documents\projects_personal\mhw_australia\data\gliders_clean_vf.mat';
folder_fig='C:\Users\Julia\Documents\projects_personal\mhw_australia\figs\w-o_neg\';

vars={'TEMP','PSAL','CPHL','DOX2','N2','OSAT','AOU'};
regions={'QLD','SW_WA','NSW','TAS'};
seasons={'winter','spring','summer','autumn'};
mhw_flag={'no_mhw','mhw','all'};

x_lim.TEMP.QLD=[23 31];
x_lim.TEMP.SW_WA=[18 26];
x_lim.TEMP.NSW=[14 26];
x_lim.TEMP.TAS=[11 19];
x_lim.PSAL.QLD=[34 36];
x_lim.PSAL.SW_WA=[34 36];
x_lim.PSAL.NSW=[34 36];
x_lim.PSAL.TAS=[34 36];
x_lim.CPHL.QLD=[0 2.5];
x_lim.CPHL.SW_WA=[0 2.5];
x_lim.CPHL.NSW=[0 2.5];
x_lim.CPHL.TAS=[0 2.5];
x_lim.DOX2.QLD=[130 270];
x_lim.DOX2.SW_WA=[130 270];
x_lim.DOX2.NSW=[130 270];
x_lim.DOX2.TAS=[130 270];
x_lim.N2.QLD=[-3 13];
x_lim.N2.SW_WA=[-3 13];
x_lim.N2.NSW=[-3 13];
x_lim.N2.TAS=[-3 13];
x_lim.OSAT.QLD=[60 110];
x_lim.OSAT.SW_WA=[60 110];
x_lim.OSAT.NSW=[60 110];
x_lim.OSAT.TAS=[60 110];
x_lim.AOU.QLD=[-30 100];
x_lim.AOU.SW_WA=[-30 100];
x_lim.AOU.NSW=[-30 100];
x_lim.AOU.TAS=[-30 100];

dv.TEMP=2;
dv.PSAL=.6;
dv.CPHL=.5;
dv.DOX2=40;
dv.N2=3;
dv.OSAT=20;
dv.AOU=30;

%% Load data
load(file_data)
depth=data.DEPTH;

%% Define season months
time_months=month(time);
idx.winter=[6 7 8];
idx.spring=[9 10 11];
idx.summer=[12 1 2];
idx.autumn=[3 4 5];

%% Calculate stratification (N²)
[N2,N2_p]=gsw_Nsquared(data.TEMP,data.PSAL,data.PRES);
N2_depth=-gsw_z_from_p(N2_p,data.LATITUDE);

% Interpolate
for i=1:size(N2_depth,2)
    is_nan=isnan(N2_depth(:,i));
    if sum(~is_nan)>1
        data.N2(:,i)=interp1(N2_depth(~is_nan,i),N2(~is_nan,i),depth,'linear','extrap');
    else
        data.N2(:,i)=NaN;
    end
    data.N2(is_nan,i)=NaN;
end

data.N2=-data.N2*10e3;

%% Calculate oxygen variables
O2_sat=gsw_O2sol(data.PSAL,data.TEMP,data.PRES,data.LONGITUDE,data.LATITUDE);
data.OSAT=(data.DOX2./O2_sat)*100;
data.AOU=O2_sat-data.DOX2;
R=data.OSAT-100;

%% Calculate mean profiles and standard deviations
for j=1:length(vars)
    v=vars{j};
    disp(['Variable: ' v])

    for i=1:length(regions)
        r=regions{i};
    
        for ii=1:length(seasons)
            s=seasons{ii};
    
            for iii=1:length(mhw_flag)
                m=mhw_flag{iii};
    
                idx_region=ismember(data.REGION,r);
                idx_season=ismember(time_months,idx.(s));
    
                switch m
                    case 'no_mhw'
                        idx_mhw=data.SEVERITY<=1;
                    case 'mhw'
                        
                        % -----
                        % Calculate temperature anomalies (positive anomalies at surface)
                        % t_anom=data.TEMP-mean(data.TEMP(:,idx_region & idx_season & data.SEVERITY<1),2,'omitnan');
                        % t_anom=sum(t_anom(1:5,:)>0,1)>=1;
                        % -----
                        
                        idx_mhw=data.SEVERITY>1; % & t_anom';

                        % -----
                        % Calculate saturation anomaly
                        R_delta.(v).(r).(s)=mean(R(:,idx_region & idx_season & idx_mhw),2,'omitnan')-mean(R(:,idx_region & idx_season & ~idx_mhw),2,'omitnan');
                        % -----

                        not_nan=sum(~isnan(data.(v)),1,'omitnan')~=0; not_nan=not_nan';
                        disp([r ' in ' s ': ' num2str(sum(idx_region & idx_season & idx_mhw & not_nan))])
                    case 'all'
                        idx_mhw=true(size(data.SEVERITY));
                end
        
                profiles_mean.(v).(r).(s).(m)=mean(data.(v)(:,idx_region & idx_season & idx_mhw),2,'omitnan');
                profiles_stdv.(v).(r).(s).(m)=std(data.(v)(:,idx_region & idx_season & idx_mhw),[],2,'omitnan');
            
                % Shadded area for standard deviations
                val_max=profiles_mean.(v).(r).(s).(m)+profiles_stdv.(v).(r).(s).(m);
                val_min=profiles_mean.(v).(r).(s).(m)-profiles_stdv.(v).(r).(s).(m);
                std_x=[val_min;flipud(val_max)];
                std_y=[depth;flipud(depth)];
                profiles_stda.(v).(r).(s).(m)(:,2)=std_y(~isnan(std_x));
                profiles_stda.(v).(r).(s).(m)(:,1)=std_x(~isnan(std_x));

                % Calculate MLD
                if ismember(v,{'TEMP'})
                    dummy_temp=data.(v)(:,idx_region & idx_season & idx_mhw);
                    for p=1:size(dummy_temp,2)
                        dummy_prof=dummy_temp(:,p);
                        if sum(~isnan(dummy_prof))>5
                            dummy_mld(p)=mld(depth(~isnan(dummy_prof)),dummy_prof(~isnan(dummy_prof)),'metric','threshold','refpres',0);
                        else
                            dummy_mld(p)=NaN;
                        end
                    end
                    profiles_mld.(r).(s).(m)=mean(dummy_mld,'omitnan');
                    clear dummy_temp dummy_prof dummy_mld
                end
    
            end
        end
    end
end

clear std_x std_y val_min val_max ii i r s idx_region idx_season idx_mhw not_nan

%% PLOT
for j=1:length(vars)
    v=vars{j};

    figure('Units','centimeters','Position',[5 5 4 3])
    
    for i=1:length(regions)
        r=regions{i};
    
        for ii=1:length(seasons)
            s=seasons{ii};
    
            clf
            hold on
            
            % plot(profiles_mean.(v).(r).(s).all,depth,'Color',[.2 .2 .2],'LineWidth',1.2)
            plot(profiles_mean.(v).(r).(s).no_mhw,depth,'Color',[.4 .6 .8],'LineWidth',1.5)
            plot(profiles_mean.(v).(r).(s).mhw,depth,'Color',[.8 .1 .1],'LineWidth',1.5)
            
            % fill(profiles_stda.(v).(r).(s).all(:,1),profiles_stda.(v).(r).(s).all(:,2), ...
            %      [.2 .2 .2],'FaceAlpha',.1,'EdgeColor','none')
            fill(profiles_stda.(v).(r).(s).no_mhw(:,1),profiles_stda.(v).(r).(s).no_mhw(:,2), ...
                 [.4 .6 .8],'FaceAlpha',.2,'EdgeColor','none')
            fill(profiles_stda.(v).(r).(s).mhw(:,1),profiles_stda.(v).(r).(s).mhw(:,2), ...
                 [.8 .1 .1],'FaceAlpha',.2,'EdgeColor','none')

            if i==1 && ii==1
                ax=gca;
                pos=ax.Position;
            end
    
            xlim(x_lim.(v).(r))
            xticks(min(x_lim.(v).(r)):dv.(v):max(x_lim.(v).(r)))
            

            if ismember(r,{'TAS','NSW'})
                ylim([0 90])
                yticks(0:20:100)
            elseif ismember(r,{'QLD'})
                ylim([0 45])
                yticks(0:10:100)
            elseif ismember(r,{'SW_WA'})
                ylim([0 35])
                yticks(0:10:100)
            end
    
            if ~ismember(s,{'winter'})
                yticklabels(' ')
            end
    
            set(gca,'YDir','reverse','Position',pos)

            % hline(profiles_mld.(r).(s).all,'Color',[.2 .2 .2],'LineStyle','--')
            hline(profiles_mld.(r).(s).no_mhw,'Color',[.4 .6 .8],'LineStyle','--')
            hline(profiles_mld.(r).(s).mhw,'Color',[.8 .1 .1],'LineStyle','--')

            if ismember(v,{'N2'})
                vline(0,'Color',[.5 .5 .5])
            elseif ismember(v,{'OSAT'})
                vline(100,'Color',[.5 .5 .5])
            end
    
            box on
            grid on
    
            fontsize(gcf,7,'points')
            fontname(gcf,'Helvetica')
    
            export_fig([folder_fig 'profs_' v '_' r '_' s],'-png','-transparent','-r360');
        end
    end
end


%% R
figure('Units','centimeters','Position',[5 5 4 3])
for i=1:length(regions)
    r=regions{i};

    for ii=1:length(seasons)
        s=seasons{ii};

        clf
        hold on
        
        plot(R_delta.(v).(r).(s),depth,'Color',[0 0 0],'LineWidth',1.2)
        
        if i==1 && ii==1
            ax=gca;
            pos=ax.Position;
        end

        xlim([-12 18])
        xticks(-10:5:18)

        if ismember(r,{'TAS','NSW'})
            ylim([0 90])
            yticks(0:20:100)
        elseif ismember(r,{'QLD'})
            ylim([0 45])
            yticks(0:10:100)
        elseif ismember(r,{'SW_WA'})
            ylim([0 35])
            yticks(0:10:100)
        end

        if ~ismember(s,{'winter'})
            yticklabels(' ')
        end

        set(gca,'YDir','reverse','Position',pos)

        % hline(profiles_mld.(r).(s).all,'Color',[.2 .2 .2],'LineStyle','--')
        hline(profiles_mld.(r).(s).no_mhw,'Color',[.4 .6 .8],'LineStyle','--')
        hline(profiles_mld.(r).(s).mhw,'Color',[.8 .1 .1],'LineStyle','--')

        vline(0,'Color',[.5 .5 .5])
        vfill(-2,2,[.2 .2 .2],'FaceAlpha',.1,'EdgeColor','none')

        box on
        grid on

        fontsize(gcf,7,'points')
        fontname(gcf,'Helvetica')

        export_fig([folder_fig 'profs_O2A_' r '_' s],'-png','-transparent','-r360');
    end
end



