%% INITIAL MODEL SETUP

% initialise grid coordinates
h = D./N;     % grid spacing
x = linspace(-h/2,D+h/2,N+2);
z = linspace(-h/2,D+h/2,N+2);
[X,Z] = meshgrid(x,z);

% initialise smooth random noise
rng(5);
rn = rand(N+2,N+2) - 0.5;
for i=1:smth
   rn(2:end-1,2:end-1) = rn(2:end-1,2:end-1) ...
                       + diff(rn(:,2:end-1),2,1)./8 ...
                       + diff(rn(2:end-1,:),2,2)./8;
    rn([1 end],:) = rn([2 end-1],:);
    rn(:,[1 end]) = rn(:,[2 end-1]);
end
rn = rn./max(abs(rn(:)));

% set initial condition
switch finit  % initial porosity
    case 'linear'
        f = f0 + (f1-f0) .* Z/D + df.*rn;
    case 'layer'
        f = f0 + (f1-f0) .* (1+erf(25*(Z/D-zlay)))/2 +   df.*rn;
end
switch Tinit  % initial temperature
    case 'linear'
        T = T0 + (T1-T0) .* Z/D + dT.*rn;
    case 'layer'
        T = T0 + (T1-T0) .* (1+erf(25*(Z/D-zlay)))/2 + dT.*rn;
end
switch Cinit  % initial concentration
    case 'linear'
        C = C0 + (C1-C0) .* Z/D + dC.*rn;
    case 'layer'
        C = C0 + (C1-C0) .* (1+erf(25*(Z/D-zlay)))/2 + dC.*rn;
end

% Set Layer and Fault porosity
f(abs(Z-LayerDepth) <= LayerWidth/2)= f_Layer;   % Create layer in porosity
ind = Z >= FaultDepth & abs((X-FaultPos)+ tand(FaultAngle)*(D-Z)...
    -0.5*(D-FaultDepth)) <= (0.5*FaultWidth/cosd(FaultAngle)); 
f(ind) = f_Fault;

for i=1:10
    f(2:end-1,2:end-1) = f(2:end-1,2:end-1) ...
                       + diff(f(:,2:end-1),2,1)./8 ...
                       + diff(f(2:end-1,:),2,2)./8;
    f([1 end],:) = f([2 end-1],:);
    f(:,[1 end]) = f(:,[2 end-1]);
end

DTDt = 0.*T(2:end-1,2:end-1);
DCDt = 0.*C(2:end-1,2:end-1);
dt   = 1e-6;


% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',14};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};

axh = 6.00; axw = 7.50; %   Height and width of axis
ahs = 1.50; avs = 1.00; %   Horzontal and vertial distance between axis
axb = 1.75; axt = 0.90; %   Bottom and top;Size of page relative to axis
axl = 1.75; axr = 0.90; %   Right and left; spacing of axis to page


% prepare and plot figure for mechanical solution fields
fh1 = figure(1); 
% clf; colormap(ocean);
fh = axb + 1*axh + 0*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh1,UN{:},'Position',[3 3 fw fh]);
set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh1,'Color','w','InvertHardcopy','off');
set(fh1,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

set(fh1, 'CurrentAxes', ax(1))
imagesc(x,z,f); axis equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Initial Porosity [vol]',TX{:},FS{:})
set(fh1, 'CurrentAxes', ax(2))
imagesc(x,z,T); axis equal tight;  box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});title('Initial Temperature [C]',TX{:},FS{:})
%text(0,0.9,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
set(fh1, 'CurrentAxes', ax(3))
imagesc(x,z,C); axis equal tight;  box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:});title('Initial Concentration',TX{:},FS{:})
% text(0,0.9,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
drawnow

% prepare solution & residual arrays for VP solver
w = zeros(N+1,N+2);  % vertical Darcy speed
u = zeros(N+2,N+1);  % horizontal Darcy speed
p = zeros(N+2,N+2);  % pore fluid pressure
F = zeros(N+2,N+2);  % residual for pressure equation


%% MAIN TIME STEPPING LOOP
m    = 0;
time = 0;
while time <= tend

    fprintf(1,'\n\n*****  step = %d,  time = %1.3e [s],  step = %1.3e \n\n',m,time,dt);

    % store previous rates of change
    To    = T;
    Co    = C;
    DTDto = DTDt;
    DCDto = DCDt;
    
    % calculate permeability [m2]
    k = a^2/b * f.^n;  % Kozeny-Carman relationship
    
    % calculate Darcy coefficient [m2/Pas] AKA Mobility-Not hydraullic
    % conductivity which is flux/section area
    K = k/mu;
    
    % calculate iterative step size
    dtau = (h/2)^2./K;
        
    % update density difference
    Drho = - rhol0.*(- aT.*(T-mean(T,2)) + gC.*(C-mean(C,2)));
    
    
    % UPDATE VELOCITY-PRESSURE SOLUTION (PSEUDO-TRANSIENT SOLVER)
    Fnorm = 1e6;
    pi    = p;
    it    = 0;
    while Fnorm >= tol || it < 100 
        
        % store previous iterative solution guesses
        pii = pi; pi = p;
        
        % calculate pressure gradient [Pa/m]
        gradPz = diff(p,1,1)./h;  % vertical gradient
        gradPx = diff(p,1,2)./h;  % horizontal gradient
        
        % calculate Darcy segregation speed [m/s]
        w = -(K(1:end-1,:)+K(2:end,:))./2 .* (gradPz + (Drho(1:end-1,:)+Drho(2:end,:))./2.*grav);
        u = -(K(:,1:end-1)+K(:,2:end))./2 .* (gradPx + 0                                       );
        
        % calculate residual of pressure equation
        F(2:end-1,2:end-1) = diff(w(:,2:end-1),1,1)./h + diff(u(2:end-1,:),1,2)./h;

        % update pressure solution
        p = pi - alpha.*F.*dtau + beta.*(pi-pii);
        
        % apply pressure boundary conditions
        p(:,1  ) = p(:,2    );
        p(:,end) = p(:,end-1);
        p(1  ,:) = p(2    ,:) + (Drho(1  ,:)+Drho(2    ,:))./2.*grav.*h;
        p(end,:) = p(end-1,:) - (Drho(end,:)+Drho(end-1,:))./2.*grav.*h;
        
        % get preconditioned residual norm to monitor convergence
        Fnorm = norm(F(:).*dtau(:),2)./norm(p(:),2);
        
        % print convergence
        if ~mod(it,100)
            
            % UPDATE TEMPERATURE SOLUTION (SEMI-EXPLICIT SOLVER?)
            dt = CFL .* min([(h/2)/max(abs(w(:))) , (h/2)/max(abs(u(:))) , (h/2)^2./kT]);  % diffusive timestep
            
            % calculate temperature advection
            wp = max(0, (w(1:end-1,2:end-1)+w(2:end,2:end-1))./2);
            wm = min(0, (w(1:end-1,2:end-1)+w(2:end,2:end-1))./2);
            up = max(0, (u(2:end-1,1:end-1)+u(2:end-1,2:end))./2);
            um = min(0, (u(2:end-1,1:end-1)+u(2:end-1,2:end))./2);
            
            grdTzp = diff(T(2:end  ,2:end-1),1,1)./h;
            grdTzm = diff(T(1:end-1,2:end-1),1,1)./h;
            grdTxp = diff(T(2:end-1,2:end  ),1,2)./h;
            grdTxm = diff(T(2:end-1,1:end-1),1,2)./h;
            
            DTDt = - (wp.*grdTzm + wm.*grdTzp + up.*grdTxm + um.*grdTxp) ...
                 + kT.* (diff(T(:,2:end-1),2,1)./h^2 + diff(T(2:end-1,:),2,2)./h^2);
            
            T(2:end-1,2:end-1) = To(2:end-1,2:end-1) + (DTDt + DTDto)/2 .* dt;
            
            % apply temperature boundary conditions
            T(:,1  ) = T(:,2    );  % left boundary: insulating
            T(:,end) = T(:,end-1);  % right boundary: insulating
            T(1  ,:) = T0;          % top boundary: isothermal
            T(end,:) = T1;          % bottom boundary: insulating
            
            
            % UPDATE CONCENTRATION SOLUTION (EXPLICIT SOLVER)
            
            % calculate concentration advection
            grdCzp = diff(C(2:end  ,2:end-1),1,1)./h;
            grdCzm = diff(C(1:end-1,2:end-1),1,1)./h;
            grdCxp = diff(C(2:end-1,2:end  ),1,2)./h;
            grdCxm = diff(C(2:end-1,1:end-1),1,2)./h;
            
            DCDt = - (wp.*grdCzm + wm.*grdCzp + up.*grdCxm + um.*grdCxp) ...
                + kC.* (diff(C(:,2:end-1),2,1)./h^2 + diff(C(2:end-1,:),2,2)./h^2);
            
            C(2:end-1,2:end-1) = Co(2:end-1,2:end-1) + (DCDt + DCDto)/2 .* dt;
            
            % apply concentration boundary conditions
            C(:,1  ) = C(:,2    );  % left boundary: closed
            C(:,end) = C(:,end-1);  % right boundary: closed
            C(1  ,:) = C0;          % top boundary: isochemical
            C(end,:) = C1;          % bottom boundary: isochemical
            
            % report convergence
            fprintf(1,'---  %d,  %e\n',it+1,Fnorm);
        end
        
        % increment iteration count
        it = it+1;
    end

    % plot solution
    if ~mod(m,nop)
        if lvplt
            fh2 = figure(2); clf;
        else
            fh2=figure('Visible','off'); clf;
        end
                     
        axh = 6.00; axw = 7.50; %   Height and width of axis
        ahs = 1.50; avs = 1.00; %   Horzontal and vertial distance between axis
        axb = 1.75; axt = 2.00; %   Bottom and top;Size of page relative to axis
        axl = 1.75; axr = 0.90; %   Right and left; spacing of axis to page
        
        fh = axb + 2*axh + 1*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh2,UN{:},'Position',[9 9 fw fh]);
        set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh2,'Color','w','InvertHardcopy','off');
        set(fh2,'Resize','off');
        ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
        ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
        ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
        ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

        sgtitle(sprintf('Time elapsed %.1f years', time/31557600),TX{:},FS{:})
        text(0,0.9,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
         
        set(fh2, 'CurrentAxes', ax(1))
        imagesc(x,z,-w.*3600*24*365.25); axis equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Segregation z-speed [m/yr]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(2))
        imagesc(x,z,u.*3600*24*365.25); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Segregation x-speed [m/yr]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(3))
        imagesc(x,z,p); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Dynamic fluid pressure [Pa]',TX{:},FS{:})
%         text(0,0.9,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
        
        set(fh2, 'CurrentAxes', ax(4))
        imagesc(x,z,T); axis equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Temperature [C]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(5))
        imagesc(x,z,C); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Concentration [wt]',TX{:},FS{:})
        
        set(fh2, 'CurrentAxes', ax(6))
        imagesc(x,z,f); axis equal tight;  box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title('Porosity',TX{:},FS{:})
        drawnow      
                    
        
        
        % print figure to file
        if svfig
            print(fh2,['../out/',runID,'/',runID,'_',int2str(m/nop)],'-dpng','-r200')
        end
    end
    
    % update time and step count
    m    = m + 1;
    time = time + dt;  

end 
