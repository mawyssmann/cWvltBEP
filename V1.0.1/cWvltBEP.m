classdef cWvltBEP
    %% *******************************************************************
    % FUNCTION CLASS DESCRIPTION
    % This function class performs wavelet-based analyses of individual bed
    % elevation profiles (BEPs) to determine dune length and height.

    % SOFTWARE DEPENDENCIES: 
    % 
    %  REQUIRED SOFTWARE
    %   1) "Cross wavelet and wavelet coherence toolbox for MATLAB" by
    %   Aslak Grinsted, John Moore, and Svetlana Jevrejeva. 
    %    Link for access: https://github.com/grinsted/wavelet-coherence
    %    Publication: Grinsted, A., J. C. Moore, S. Jevrejeva (2004),
    %    Application of the cross wavelet transform and wavelet coherence
    %    to geophysical time series, Nonlin. Process. Geophys., 11, 561566

    % SOFTWARE AUTHOR: Micah A. Wyssmann



    %% *******************************************************************
    % VERSION NOTES
    % Version 1.0.1 - 11/7/25
    %   First fully developed version made openly accessible

    %% *******************************************************************
    % PROPERTIES AND INITIALIZATION

    properties  % data container definition to define the entire class
        x           % x vector for BEP (MUST be in METERS)
        z           % z vector for BEP (MUST be in METERS)
        varz        % signal variance (m^2)
        dx          % signal delta_x (m)
        N           % number of samples in x/z

        % wavelet parameters independent of energy-normalized (L2 norm) vs. amplitude-normalized (L1 norm)
        scale       % scale for wavelet
        Nv          % # of voices per octave
        dj          % inverse of # voices per octave in wavelet scales (dj = 1/Nv)
        coi         % cone of influence
        period      % period for wavelet

        % energy-normalized wavelet data
        wave_L2      % wavelet coefficients (complex values) from L2 norm (TC98)
        P_L2         % wavelet power
        P95         % wavelet power associated with 95% significance level
        PR95        % power ratio relative to 95% significance level
        varw        % reconstructed variance from wavelet coefficients
        varRat      % ratio of wavelet variance to real variance
        
        % L1 norm wavelet data
        wave_L1      % wavelet coefficients (complex values) from L1 norm (cwt & J.M.Lilly)
        A_L1         % wavelet amplitude from L1 norm
        
        % Dune scale estimates and parameters
        lamwx       % length scale directly from Fourier Period
        AwxL1       % L1 amplitude (to use for estimating Hwx)
        lamwxs      % smoothed length scale directly from Fourier Period
        AwxL1s      % smoothed L1 amplitude (to use for estimating Hwx)
        CLm         % coefficient for Lw
        CHm         % coefficient for Hw
        CL_CIs      % configetLwxHwxdence interval coefficients for Lw
        CH_CIs      % confidence interval coefficients for Hw
        Lwx         % length scale
        Hwx         % height scale
        Lwx_CIs     % lower and upper confidence intervals for length scale
        Hwx_CIs     % lower and upper confidence intervals for length scale
    end
    
    methods     % class initialization methods
        function obj=cWvltBEP(x_m,z_m)
            obj.x = x_m(:);
            obj.z = detrend(z_m(:)); % adds a basic linear detrending
            obj.varz = var(obj.z);
            obj.dx = abs(obj.x(2) - obj.x(1));
            obj.N = length(x_m); 
        end
    end
    
    %% *******************************************************************
    % CLASS FUNCTIONS
    
    %% BASIC WAVELET CALCULATIONS

    methods 
        function obj = getWt(obj,Nv)
            if nargin == 1; Nv = 12; end % # voices per octave = 12 if not input
            
            % basics and L2 norm call
            obj.Nv = Nv; 
            obj.dj = 1/Nv;
            [obj.wave_L2,obj.period,obj.scale,obj.coi,sig95] = ...
                wt([obj.x(:),obj.z(:)],'MOTHER','MORLET','Dj',obj.dj);
            obj.P_L2 = (abs(obj.wave_L2)).^2; % wavelet power
            obj.P95 = mean((((sig95).^(-1)).*obj.P_L2)/obj.varz,2)'; % power at 95% significance
            obj.PR95 = sig95; % power ratio relative to 95% significance level 

            % L1 norm call
            kmin = min(obj.period.^(-1));
            kmax = max(obj.period.^(-1));
            klims = [kmin kmax];
            [obj.wave_L1,kL1,~] = cwt(obj.z(:),"amor",1/obj.dx,VoicesPerOctave=Nv,FrequencyLimits=klims);
            period_L1 = kL1(:).^(-1);
            obj.A_L1 = abs(obj.wave_L1);

            % compare and check L2 vs. L1 norm results in terms of period
            % compatability
            checkL2vsL1(obj,period_L1)

            % get variance associated with L2 norm wavelet
            obj = getWtL2Var(obj);
        end

        function obj = getWtL2Var(obj)
            prefac = repmat(obj.scale(:).^-1,1,length(obj.x));
            Pnorm = prefac.*obj.P_L2;
            C_delta = 0.776; % constant for a given wavelet, 0.776 for Morlet (see Torrence & Compo 1998)
            C = (obj.dj*obj.dx)/(C_delta);
            obj.varw = sum(C*Pnorm(:))/obj.N;
            obj.varRat = obj.varw/obj.varz;
        end
        
        function checkL2vsL1(obj,period_L1)
            period_PctDiff = (obj.period(:) - period_L1(:))./obj.period(:);
            
            if ~(length(obj.period) == length(period_L1))
                error('L1 and L2 norm approaches have different # of periods')
            end
            
            if max(period_PctDiff) > 0.001
                error('L1 and L2 norm period values differ by more than 0.1%')
            end
        end
    end
    
    %% BASIC WAVELET & GENERAL ANALYSIS PLOTTING
    
    methods 
        function plotWt(obj,plotver,log2C,plot95)
            if plotver == 1 % basic power plot defined by TC98
                C = obj.P_L2/obj.varz; % Power/sigma^2
                HCB = plotWt_basic(obj,C,log2C,plot95); 
                if log2C == 1
                    HCB.Label.String = '$log_2(|W_{x}(\lambda)|^2_{L2}/\sigma^2)$';
                else 
                    HCB.Label.String = '$|W_{x}(\lambda)|^2_{L2}/\sigma^2$';
                end
                HCB.Label.Interpreter = 'latex';
            elseif plotver == 2 % plot 95% significance
                C = obj.PR95;
                HCB = plotWt_basic(obj,C,log2C,plot95); 
                if log2C == 1
                    HCB.Label.String = '$log_2(PR_{95,L2})$';
                else
                    HCB.Label.String = '$PR_{95e}$';
                end
                HCB.Label.Interpreter = 'latex';
            elseif plotver == 3 % PLot amplitude from L1 norm
                C = obj.A_L1; 
                HCB = plotWt_basic(obj,C,log2C,plot95); 
                if log2C == 1
                    HCB.Label.String = '$log_2(|W_{a}(x,\lambda)|)$';
                else
                    HCB.Label.String = '$|W_{a}(x,\lambda)|$ (m)';
                end
                HCB.Label.Interpreter = 'latex';
            end
            set(gca,'FontName','Arial','FontSize',16)
            if exist('cbrewer') > 0
                cmp = cbrewer('seq','Blues',64,'linear');
            else
                cmp = parula(64);
            end
            colormap(cmp);
        end

        function plotWtandTS(obj,plotver,log2C,plot95,xlims,ylims,clims)
            % ----- HARDWIRED INPUTS (for visual figure presentation) -----
            dm = 0.7;       % max percentage of total window that figure takes up
            hfrac = 0.7;    % fraction of total height used by bottom pane
            % -------------------------------------------------------------
            
            figure('Position',[10 10 600 700]);

            shf = 0.5*(1-dm);
            x1 = shf; y1 = shf; w1 = dm; h1 = hfrac*dm;
            pos1 = [x1 y1 w1 h1];
            subplot('Position',pos1);
            plotWt(obj,plotver,log2C,plot95)
            moveColorbar(obj,'southoutside')
            % colorbar off
            cropWt(obj,xlims,ylims,clims);
            % set(gca,'XDir','reverse');
            
            x2 = shf; y2 = h1 + shf*2; w2 = dm; h2 = 1 - h1 - shf*2.5;
            pos2 = [x2 y2 w2 h2];
            subplot('Position',pos2);
            % x_RM = wtInfo.x/1609.34;
            plot(obj.x,obj.z,'k');
            xlabel('$x$ (m)','Interpreter','latex'); 
            ylabel('${\Delta}z$ (m)','Interpreter','latex');
            % set(gca,'XDir','reverse');
            set(gca,'FontName','Arial','FontSize',16);
            if isempty(xlims)
                xlim([min(obj.x) max(obj.x)]);
            else
                xlim(xlims);
            end
        end

        function HCB = plotWt_basic(obj,C,log2C,plot95)
            
            Yticks = 2.^(fix(log2(min(obj.period))):fix(log2(max(obj.period))));
            
            if log2C == 1
                C = log2(C);
            end
            
            H = imagesc(obj.x,log2(obj.period),C);
            
            if log2C == 1
                % Set colorbar limits and labels
                clim=get(gca,'clim'); %center color limits around log2(1)=0
                clim=[-1 1]*max(clim(2),3);
                set(gca,'clim',clim)
                
                cmax = floor(clim(2));
                if cmax > 8
                    cmax = 2*floor(cmax/2);
                    cvals = -cmax:2:cmax;
                else
                    cvals = -cmax:1:cmax;
                end
                
                HCB=colorbar;
                set(HCB,'ytick',cvals); % deleted -7:7
                barylbls=rats(2.^(get(HCB,'ytick')'),24);
                % barylbls([1 end],:)=' ';
                barylbls(:,all(barylbls==' ',1))=[];
                set(HCB,'yticklabel',barylbls);
            else
                HCB = colorbar;
            end
            
            % set y-axis labels and limits
            set(gca,'YLim',log2([min(obj.period),max(obj.period)]), ...
                'YDir','reverse', ...
                'YTick',log2(Yticks(:)), ...
                'YTickLabel',num2str(Yticks'), ...
                'layer','top')
            xlabel('$x$ (m)','Interpreter','latex')
            ylabel('$\lambda$ (m)','Interpreter','latex')
            hold on
            
            % Plot contours of 95% confidence
            if plot95 == 1
            [c,h] = contour(obj.x,log2(obj.period),obj.PR95,[1 1],'k'); %#ok
            set(h,'linewidth',2)
            end
            
            % Plot confidence interval
            % x2 = obj.x(:);
            tt=[obj.x([1 1])-obj.dx*.5; obj.x; obj.x([end end])+obj.dx*.5];
            hcoi=fill(tt,log2([obj.period([end 1]) obj.coi obj.period([1 end])]),'w');
            set(hcoi,'alphadatamapping','direct','facealpha',.5)
            
            hold off
            set(gca,'box','on','layer','top');
            % set(gca,'XDir','reverse');
        end

        function cropWt(obj,xlims,ylims,clims)
            if ~isempty(xlims); xlim(xlims); end
            if ~isempty(ylims); ylim(log2(ylims)); end
            if ~isempty(clims); clim(clims); end
        end

        function moveColorbar(obj,newLoc)
            ax = gca;
            HCB = ax.Colorbar;
            HCB.Location = newLoc;
        end
    end
    
    %% WAVELET-BASED CALCULATIONS FOR DUNE SCALE ESTIMATION

    methods % calculation methods
        function ys = getGausSmth(obj,x,y,ls,Plt)
            % Inputs: 1) x - irregular x locs; 2) y response variable at x (L/H);
            % 3) lengthscale for smoothing (2x sigma (std dev) of Gaus)

            if nargin < 4; Plt = 0; end

            x = x(:)'; % make sure horizontal vectors
            y = y(:)';
            [n,~] = size(x);
            
            Xi = repmat(x,n,1);
            X  = repmat(x(:),1,n);
            Y  = repmat(y(:),1,n);
            
            alpha = exp(-((X-Xi).^2)/(2*(ls/2)^2));
            
            s_alpha = sum(alpha,1);
            s_alphY = sum(alpha.*Y,1);
            
            ys = s_alphY./s_alpha; 

            if Plt == 1
                figure; plot(x,y,'.k',x,ys,'or');
            end
        end
        
        % function obj = getLwxHwx_OLD(obj,CLm,CL_CIs,CHm,CH_CIs,PeriodRng,PR95_th,mpow,ls)
        %     if nargin < 8; mpow = 1; end
        %     % if nargin < 9; uselg2 = 0; end
        % 
        %     Prat95 = obj.PR95; % pull PR95 for use
        % 
        %     if ~isempty(PR95_th) % zero anything < PR95_th
        %         Prat95(Prat95 < PR95_th) = 0;
        %     end
        % 
        %     periodMat = repmat(obj.period(:),1,obj.N); % matrix of periods
        %     if ~isempty(PeriodRng) % get rid of values outside of range
        %         Prat95(periodMat < PeriodRng(1)) = 0;
        %         Prat95(periodMat > PeriodRng(2)) = 0;
        %     end
        % 
        %     % DO CALCS WITH MAX PR95 METHOD - for amplitude (height)
        %     [M,I] = max(Prat95,[],1,'linear'); % get max
        % 
        %     Awxa = obj.A_L1(I);         % pull amplitude values associated with max PR95
        %     Awxa(M < PR95_th) = NaN;    % drop values with PR95 < PR95_th
        % 
        %     % DO CALCS WITH WEIGHTED AVERAGE PR95 METHOD - for length
        % 
        %     % Calculate weights and weighted average
        %     wts = Prat95.^mpow;  % initialize weights as sig95^m
        %     numer = sum(wts.*periodMat,1);
        %     denom = sum(wts,1);
        %     lamwxa = numer./denom;
        % 
        %     % replace points outside cone of influence (coi) with NaN
        %     fac = 0.5; % offset factor from coi where results are affected
        %     outsidecoi = (lamwxa > (fac*obj.coi));
        %     lamwxa(outsidecoi) = NaN; 
        %     Awxa(outsidecoi) = NaN;
        % 
        %     % smooth wavelet signal scale estimate data
        %     nw = round(ls/obj.dx);
        %     lamwxas  = smoothdata(lamwxa(:),'movmean',nw);
        %     Awxas  = smoothdata(Awxa(:),'movmean',nw);
        %     lamwxas(isnan(lamwxa))    = NaN;  % mask data after smoothing
        %     Awxas(isnan(Awxa))    = NaN; 
        % 
        %     % assign and save to structure
        %     obj.lamwx   = lamwxa(:);    % Raw lambda value
        %     obj.AwxL1   = Awxa(:);      % Raw amplitude value
        %     obj.lamwxs  = lamwxas(:);   % smoothed lambda value
        %     obj.AwxL1s  = Awxas(:);     % smoothed amplitude value
        %     obj.CLm     = CLm;          % median CL coefficient
        %     obj.CHm     = CHm;          % median CH coefficient
        %     obj.CL_CIs  = CL_CIs;       % 95% CI values for CL coefficient
        %     obj.CH_CIs  = CH_CIs;       % 95% CI values for CH coefficient
        %     obj.Lwx     = CLm*lamwxas(:); % Lwx value
        %     obj.Hwx     = 2*CHm*Awxas(:); % Hwx value
        %     obj.Lwx_CIs = [CL_CIs(1)*lamwxas(:), CL_CIs(2)*lamwxas(:)]; % Lwx 95% CIs
        %     obj.Hwx_CIs = [2*CH_CIs(1)*Awxas(:), 2*CH_CIs(2)*Awxas(:)]; % Hwx 95% CIs
        % end

        function obj = getLwxHwx(obj,PeriodRng,PR95_th,mpow,ls,CLm,CL_CIs,CHm,CH_CIs)            
            obj = getLamdaAmp(obj,PeriodRng,PR95_th,mpow); % Get weighted average lambda and max amplitude (saved in obj)
            obj = getSmoothedLamdaAmp(obj,ls); % Get smoothed version of lambda and max amplitude (saved in obj)
            obj = calcLwxHwx(obj,CLm,CL_CIs,CHm,CH_CIs); % Calculate Lwx and Hwx using smoothed lambda/amplitude and CL/CH coefficients
        end

        function obj = getLamdaAmp(obj,PeriodRng,PR95_th,mpow)
            Prat95 = obj.PR95; % pull PR95 for use
            
            if ~isempty(PR95_th) % zero anything < PR95_th
                Prat95(Prat95 < PR95_th) = 0;
            end

            periodMat = repmat(obj.period(:),1,obj.N); % matrix of periods
            if ~isempty(PeriodRng) % get rid of values outside of range
                Prat95(periodMat < PeriodRng(1)) = 0;
                Prat95(periodMat > PeriodRng(2)) = 0;
            end
            
            % DO CALCS WITH MAX PR95 METHOD - for amplitude (height)
            [M,I] = max(Prat95,[],1,'linear'); % get max

            Awxa = obj.A_L1(I);         % pull amplitude values associated with max PR95
            Awxa(M < PR95_th) = NaN;    % drop values with PR95 < PR95_th
            
            % DO CALCS WITH WEIGHTED AVERAGE PR95 METHOD - for length

            % Calculate weights and weighted average
            wts = Prat95.^mpow;  % initialize weights as sig95^m
            numer = sum(wts.*periodMat,1);
            denom = sum(wts,1);
            lamwxa = numer./denom;

            % replace points outside cone of influence (coi) with NaN
            fac = 0.5; % offset factor from coi where results are affected
            outsidecoi = (lamwxa > (fac*obj.coi));
            lamwxa(outsidecoi) = NaN; 
            Awxa(outsidecoi) = NaN;

            obj.lamwx   = lamwxa(:);    % Raw lambda value
            obj.AwxL1   = Awxa(:);      % Raw amplitude value
        end

        function obj = getSmoothedLamdaAmp(obj,ls)
            % smooth wavelet signal scale estimate data
            nw = round(ls/obj.dx);
            lamwxas  = smoothdata(obj.lamwx(:),'movmean',nw);
            Awxas  = smoothdata(obj.AwxL1(:),'movmean',nw);
            lamwxas(isnan(obj.lamwx))    = NaN;  % mask data after smoothing
            Awxas(isnan(obj.AwxL1))    = NaN; 

            % Assign and save to structure
            obj.lamwxs  = lamwxas(:);   % smoothed lambda value
            obj.AwxL1s  = Awxas(:);     % smoothed amplitude value
        end

        function obj = calcLwxHwx(obj,CLm,CL_CIs,CHm,CH_CIs) % This function updates calculations if coefficients are adjusted
            % Calculations and storate into structure are done simultaneously
            obj.CLm     = CLm;          % median CL coefficient
            obj.CHm     = CHm;          % median CH coefficient
            obj.CL_CIs  = CL_CIs;       % 95% CI values for CL coefficient
            obj.CH_CIs  = CH_CIs;       % 95% CI values for CH coefficient
            obj.Lwx     = CLm*obj.lamwxs(:); % Lwx value
            obj.Hwx     = 2*CHm*obj.AwxL1s(:); % Hwx value
            obj.Lwx_CIs = [CL_CIs(1)*obj.lamwxs(:), CL_CIs(2)*obj.lamwxs(:)]; % Lwx 95% CIs
            obj.Hwx_CIs = [2*CH_CIs(1)*obj.AwxL1s(:), 2*CH_CIs(2)*obj.AwxL1s(:)]; % Hwx 95% CIs
        end
        
        function [RMSE_L,RMSE_H] = getRMSE_LH(obj,xt,Lt,Ht)
            Lw_xt = interp1(obj.x,obj.Lwx,xt,'linear');
            RMSE_L = rmse(Lw_xt,Lt,'omitnan');

            Hw_xt = interp1(obj.x,obj.Hwx,xt,'linear');
            RMSE_H = rmse(Hw_xt,Ht,'omitnan');
        end

        function [RMSE_L,RMSE_H,MAE_L,MAE_H] = getRMSE_LH_Gsmth(obj,xt,Lt,Ht,ls,PLT)
            nw = round(ls/obj.dx); % number to smooth for wavelet data
            
            Lwxs = smoothdata(obj.Lwx,'movmean',nw);
            Lwxs_xt = interp1(obj.x,Lwxs,xt,'linear');
            Lts = getGausSmth(obj,xt,Lt,ls,0);
            RMSE_L = rmse(Lwxs_xt(:),Lts(:),'omitnan');
            MAE_L = mean(abs(Lwxs_xt(:)-Lts(:)),'omitnan');
            % MnPctErr_L = mean(abs((Lwxs_xt(:)-Lts(:))./Lts(:)),'omitnan');
            
            Hwxs = smoothdata(obj.Hwx,'movmean',nw);
            Hwxs_xt = interp1(obj.x,Hwxs,xt,'linear');
            Hts = getGausSmth(obj,xt,Ht,ls,0);
            RMSE_H = rmse(Hwxs_xt(:),Hts(:),'omitnan');
            MAE_H = mean(abs(Hwxs_xt(:)-Hts(:)),'omitnan');
            % MnPctErr_H = mean(abs((Hwxs_xt(:)-Hts(:))./Hts(:)),'omitnan');
            
            if PLT == 1
                figure('Position',[40 60 400 600])
                subplot(2,1,1)
                loglog(Lts,Lwxs_xt,'.k');
                xlabel('L (m)'); ylabel('Lw (m)')
                xl = xlim; hold on;
                loglog(xl,xl,'-k',xl,1.3*xl,'--k',xl,0.7*xl,'--k');

                subplot(2,1,2)
                loglog(Hts,Hwxs_xt,'.k');
                xlabel('H (m)'); ylabel('Hw (m)')
                xl = xlim; hold on;
                loglog(xl,xl,'-k',xl,1.3*xl,'--k',xl,0.7*xl,'--k');
            end
        end

        function [C_Lw_all,C_Hw_all] = getallCHLw(obj,xt,Lt,Ht)
            % get values interpolated onto grid for true values
            lamw_xt = interp1(obj.x,obj.lamwx,xt,'linear');
            Aw_xt = interp1(obj.x,obj.AwxL1,xt,'linear');
            
            C_Lw_all = Lt./lamw_xt;
            C_Hw_all = Ht./(2*Aw_xt);
        end

        function [C_Lw_all,C_Hw_all] = getallCHLw_smth(obj,xt,Lt,Ht,ls)
            % get values interpolated onto grid for true values
            nw = round(ls/obj.dx);
            
            lamws = smoothdata(obj.lamwx(:),'movmean',nw);
            Aws = smoothdata(obj.AwxL1(:),'movmean',nw);
            
            lamws_xt = interp1(obj.x,lamws,xt,'linear');
            Aws_xt = interp1(obj.x,Aws,xt,'linear');

            Lts = getGausSmth(obj,xt,Lt,ls,0);
            Hts = getGausSmth(obj,xt,Ht,ls,0);
            
            C_Lw_all = Lts(:)./lamws_xt;
            C_Hw_all = Hts(:)./(2*Aws_xt);
        end

        function [Lts,Hts,Lwm,Hwm] = getall_LwmHwm_LtHt_smth(obj,xt,Lt,Ht,ls,CLm,CHm)
            % get values interpolated onto grid for true values
            nw = round(ls/obj.dx);
            
            lamws = smoothdata(obj.lamwx(:),'movmean',nw);
            Aws = smoothdata(obj.AwxL1(:),'movmean',nw);
            
            lamws_xt = interp1(obj.x,lamws,xt,'linear');
            Aws_xt = interp1(obj.x,Aws,xt,'linear');

            Lts = getGausSmth(obj,xt,Lt,ls,0); Lts = Lts(:);
            Hts = getGausSmth(obj,xt,Ht,ls,0); Hts = Hts(:);
            
            Lwm = lamws_xt*CLm; Lwm = Lwm(:);
            Hwm = 2*Aws_xt*CHm; Hwm = Hwm(:);
        end

        function DataTable = getSummaryDataTable(obj)
            x_m = obj.x_m;

        end
    end

    %% POST-PROCESSING PLOTTING & VISUALIZATION

    methods
        function plot_EvalLwxHwx(obj,ver,Inputs)
        % versions
            % 1) raw (no normalization or comparison)
            % 2) raw plus dunes over top for compare
            % 3) normalized for comparison with simple signals
            %       inputs = [Lreal Hreal]

            if ver == 2  % comparison with known dune scales
                figure('Position',[40 60 950 500])
                p1 = subplot(3,2,[1 3 5]); % wavelet plot
                plotWt(obj,3,0,1)
                hold on; plot(obj.x,log2(obj.Lwx),'-r')
                c = get( ancestor(p1, 'axes'), 'Colorbar');
                c.Location = 'southoutside';
                plot(Inputs(:,1),log2(Inputs(:,2)),'or');
                
                subplot(3,2,2); % BEP series
                plot(obj.x,obj.z,'-k')
                xlabel('$x$ (m)','Interpreter','latex');
                ylabel('$z$ (m)','Interpreter','latex');
                set(gca,'FontName','Arial','FontSize',16);
                xlim([0 obj.x(end)])
                
                subplot(3,2,4) % 
                plot(obj.x,obj.Lwx,'-r',Inputs(:,1),Inputs(:,2),'or'); 
                xlabel('$x$ (m)','Interpreter','latex'); 
                ylabel('$L$ (m)','Interpreter','latex'); 
                xlim([0 obj.x(end)])
                set(gca,'FontName','Arial','FontSize',16);
                
                subplot(3,2,6)
                plot(obj.x,obj.Hwx,'-b',Inputs(:,1),Inputs(:,3),'ob'); 
                xlabel('$x$ (m)','Interpreter','latex'); 
                ylabel('$H$ (m)','Interpreter','latex'); 
                set(gca,'FontName','Arial','FontSize',16);
                xlim([0 obj.x(end)])
            elseif ver == 3 % normalization with a simple known scale
                Ltrue = Inputs(1);
                Htrue = Inputs(2); 
                
                figure('Position',[40 60 950 500])
                p1 = subplot(3,2,[1 3 5]); % wavelet plot
                plotWt(obj,3,0,1)
                hold on; plot(obj.x,log2(obj.Lwx),'.r')
                c = get( ancestor(p1, 'axes'), 'Colorbar');
                c.Location = 'southoutside';
                
                subplot(3,2,2); % BEP series
                plot(obj.x,obj.z,'-k')
                xlabel('$x$ (m)','Interpreter','latex');
                ylabel('$z$ (m)','Interpreter','latex');
                set(gca,'FontName','Arial','FontSize',16);
                xlim([0 obj.x(end)])
                
                subplot(3,2,4) % 
                plot(obj.x,Ltrue./(obj.Lwx),'.r',[0 obj.x(end)],[1 1],'--k'); 
                xlabel('$x$ (m)','Interpreter','latex'); 
                ylabel('$C_L$','Interpreter','latex'); 
                % ylabel('$L_{true}/L_{wx}$','Interpreter','latex'); 
                xlim([0 obj.x(end)])
                set(gca,'FontName','Arial','FontSize',16);
                
                subplot(3,2,6)
                plot(obj.x,Htrue./(obj.Hwx),'.b',[0 1],[1 1],'--k'); 
                xlabel('$x$ (m)','Interpreter','latex'); 
                ylabel('$C_H$','Interpreter','latex'); 
                % ylabel('$H_{true}/2H_{wx}$','Interpreter','latex'); 
                set(gca,'FontName','Arial','FontSize',16);
                xlim([0 obj.x(end)])
            end
            if ver == 4
                figure('Position',[50 60 600 500]);

                subplot(2,1,1) % 
                plot(obj.x,obj.Lwx,'-k',Inputs(:,1),Inputs(:,2),'ok'); 
                % xlabel('$x$ (m)','Interpreter','latex'); 
                ylabel('$L$ (m)','Interpreter','latex'); 
                xlim([0 obj.x(end)])
                set(gca,'FontName','Arial','FontSize',16);
                legend('Wavelet-based estimate','Ground truth')
                xlim([min(obj.x) max(obj.x)])
                
                subplot(2,1,2)
                plot(obj.x,obj.Hwx,'-k',Inputs(:,1),Inputs(:,3),'ok'); 
                xlabel('$x$ (m)','Interpreter','latex'); 
                ylabel('$H$ (m)','Interpreter','latex'); 
                set(gca,'FontName','Arial','FontSize',16);
                xlim([min(obj.x) max(obj.x)])
            end
        end

        function plot_EvalLwxHwx_CIs(obj,xref,Lref,Href,ls,GT) % CHm,CH_CIs
            % Inputs: 1-4) median and confidence intervals for CL CH; 5) ls
            % - smoothing window (in meters); 
            % INPUT NOTES: Confidence interval variables can be left empty
            % and it will plot without confidence intervals
            
            if nargin < 6; GT = 1; end

            % Harwired plotting variables
            color_value = [0 0.4470 0.7410];
            alpha_value = 0.15;
            line_width  = 0.75;
            
            figure('Position',[50 60 600 500]);

            % -- Lwx plot --
            subplot(2,1,1) % Lwx plot
            hold on;
            Bot = obj.Lwx_CIs(:,1);
            Top = obj.Lwx_CIs(:,2);
            X   = [obj.x(:); flipud(obj.x(:))];
            Y   = [Bot; flipud(Top)];
            X   = X(~isnan(Y));
            Y   = Y(~isnan(Y));
            
            fill(X,Y,color_value,'FaceAlpha',alpha_value,'EdgeColor','none');

            plot(obj.x,obj.Lwx,'LineWidth',line_width,'Color',color_value);
            if ~isempty(Lref)
                plot(xref,getGausSmth(obj,xref,Lref,ls,0),'.k','MarkerSize',8)
            end

            if ~isempty(Lref) % there is ground truth
                entries = cell(3,1);
                entries{1} = 'Wavelet estimate - 95% CI';
                entries{2} = 'Wavelet estimate - median';
                if GT == 1; entries{3} = 'Ground truth'; 
                else; entries{3} = 'Reference data'; end
            else % no ground truth/reference
                entries = cell(3,1);
                entries{1} = 'Wavelet estimate - 95% CI';
                entries{2} = 'Wavelet estimate - median';
            end

            legend(entries)
            ylabel('$L$ (m)','Interpreter','latex'); 
            hold off   
            box on
            set(gca,'FontName','Arial','FontSize',16);
            

            % -- Hwx plot --

            subplot(2,1,2) % Hwx plot
            hold on
            Bot = obj.Hwx_CIs(:,1);
            Top = obj.Hwx_CIs(:,2);
            X   = [obj.x(:); flipud(obj.x(:))];
            Y   = [Bot; flipud(Top)];
            X   = X(~isnan(Y));
            Y   = Y(~isnan(Y));

            fill(X,Y,color_value,'FaceAlpha',alpha_value,'EdgeColor','none');

            plot(obj.x,obj.Hwx,'LineWidth',line_width,'Color',color_value);
            if ~isempty(Href)
                plot(xref,getGausSmth(obj,xref,Href,ls,0),'.k','MarkerSize',8)
            end

            xlabel('$x$ (m)','Interpreter','latex'); 
            ylabel('$H$ (m)','Interpreter','latex'); 
            hold off   
            box on
            set(gca,'FontName','Arial','FontSize',16);
        end

    end 

end