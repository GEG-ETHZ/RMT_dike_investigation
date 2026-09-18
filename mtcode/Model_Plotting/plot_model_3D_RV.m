function plot_model_3D_RV(m,d)
%
% Function which renders MT models in 3D as cell faces and vertices.
% Publication-ready version: adds lighting, controllable edge style,
% orthographic projection, labeled log-resistivity colorbar, and
% stations plotted on top of the model surface (topography/dike).
%
% Usage: plot_model_3D_RV(m,d)
%
% Inputs: "m" is a standard model structure
%       "d" is a standard data structure
%
%

u = user_defaults;
[L] = load_geoboundary_file_list;
close all

z_lim = [u.zmin u.zmax];

plot_stations = 1; % 1 to plot stations, 0 to skip

% ---- publication-style options (edit to taste) ----
edge_style   = 'k';           % 'none' for no edges, or a color like [0.3 0.3 0.3] / 'k' for cell edges
edge_width   = 0.5;           % only used if edge_style is not 'none'
use_lighting = 0;             % 1 to add a light source + smooth shading, 0 for flat (original look)
use_ortho    = 0;             % 1 for orthographic camera projection, 0 for default perspective (original look)
highlight_dike = 1;           % 1 to recolor the dike bright green in the colormap, 0 to leave colormap untouched
dike_ohmm      = 28;          % dike resistivity in ohm-m (from your model setup) -- change if this is wrong
dike_color     = [0.15 0.95 0.25]; % bright green
dike_band_frac = 0.01;        % fraction of the color axis range to recolor around the dike value (widen if the dike doesn't fully turn green, narrow if too much nearby resistivity also turns green)

% should not need to edit below here
f1 = 1:4; % vertices to connect for 1 face of a rectangular cell
f5 = [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 1 4 8 5]; % vertices to connect for 5 faces (except bottom) of a rectangular cell

[mx,my,mz] = meshgrid(m.x./1000,m.y./1000,(m.z)./1000);
mx = permute(mx,[2 1 3]);
my = permute(my,[2 1 3]);
mz = permute(mz,[2 1 3]);
A = log10(m.A);

set_figure_size(1);

% saves the figure exactly as it currently looks on screen -- whatever
% rotation, zoom, or pan you've applied by hand is preserved, since
% exportgraphics/savefig capture the figure's current rendered state
    function save_current_view()
        [file, savepath] = uiputfile( ...
            {'*.png','PNG image (*.png)'; ...
             '*.pdf','PDF, vector (*.pdf)'; ...
             '*.fig','MATLAB figure, editable (*.fig)'}, ...
            'Save current view as');
        if isequal(file,0)
            return % user cancelled
        end
        fullpath = fullfile(savepath,file);
        [~,~,ext] = fileparts(fullpath);
        switch lower(ext)
            case '.fig'
                savefig(gcf, fullpath);
            case '.pdf'
                exportgraphics(gcf, fullpath, 'ContentType','vector');
            otherwise % .png or anything else -> high-res raster
                exportgraphics(gcf, fullpath, 'Resolution',300);
        end
        fprintf('Saved current view to: %s\n', fullpath);
    end

% helper to apply patch styling consistently (semicolon suppresses the
% "ans = Patch with properties..." console spam)
    function draw_patch(verts,faces,colors)
        if ischar(edge_style) && strcmpi(edge_style,'none')
            patch('vertices',verts,'faces',faces,'cdata',colors,...
                'facecolor','flat','edgecolor','none');
        else
            patch('vertices',verts,'faces',faces,'cdata',colors,...
                'facecolor','flat','edgecolor',edge_style,'linewidth',edge_width);
        end
    end

while 1
lim_menu = menu('3D Plotting Options', ...
    'Plot with default x and y limits', ...
    'Enter custom x and y limits', ...
    'Save current view as image', ...
    'Return');

if lim_menu == 1 % find default limits
    xind = m.npad(1)+1:m.nx-m.npad(1);
    yind = m.npad(2)+1:m.ny-m.npad(2);
    xn_ind = xind(end);
    xs_ind = xind(1);   
    ye_ind = yind(end);
    yw_ind = yind(1);   
elseif lim_menu == 2 % enter custom x, y limits and use z limits from user_defaults  
    prompt = {'x min (km)';'x max (km)';'y min (km)';'y max (km)'};
    dlg_title = 'Plot limits';
    def = {'-20'; '20'; '-20'; '20'};
    num_lines = 1;
    inp = inputdlg(prompt,dlg_title,num_lines,def);
    
    xn_ind = nearestpoint(str2double(inp{2}),unique(mx(:,:,1))); % find min and max indices to plot in each dim
    xs_ind = nearestpoint(str2double(inp{1}),unique(mx(:,:,1)));
    ye_ind = nearestpoint(str2double(inp{4}),unique(my(:,:,1)));
    yw_ind = nearestpoint(str2double(inp{3}),unique(my(:,:,1)));
elseif lim_menu == 3 % save whatever view is currently displayed (rotation/zoom preserved)
    save_current_view();
    continue % go straight back to the menu -- don't clear or replot
else
    close all
    return
end
clf

zb_ind = min(nearestpoint(u.zmax, unique(mz(1,1,:))), size(A,3)); % RV: changed these lines so that it never exceeds size(A,3)
zt_ind = max(nearestpoint(u.zmin, unique(mz(1,1,:))), 1);

%% top face and air layers
topobot_ind = zeros(xn_ind-xs_ind+1,ye_ind-yw_ind+1);
for iz = 1:zb_ind % find lowest points in topo
    tmp = isnan(A(xs_ind:xn_ind,yw_ind:ye_ind,iz));
    topobot_ind = topobot_ind + tmp;
    if sum(sum(tmp)) ==0 % if entirely below ground
        break
    end
end
% max(max(topobot_ind)) + 1 is the first layer with no air cells in chosen plot area
% min(min(topobot_ind)) + 1 is the first layer with an earth cell

% 1. Check layer to plot. 
% a) If in air, loop from layer to first layer with earth cell. 
% In each step find the cells to be plotted and combine into one patch object.
% b) If below ground, just plot that one slice.

if zt_ind > max(max(topobot_ind)) % if layer to be plotted has no air cells
    tic
    
    verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*4,3); % 4 vertices per cell
    faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),4); % just one face per cell
    colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),1);
    
    count = 1;
    for iew = yw_ind:ye_ind-1
        for ins = xs_ind:xn_ind-1
        v = [[mx(ins,iew,zt_ind) my(ins,iew,zt_ind) mz(ins,iew,zt_ind)];...
        [mx(ins+1,iew,zt_ind) my(ins+1,iew,zt_ind) mz(ins+1,iew,zt_ind)];...
        [mx(ins+1,iew+1,zt_ind) my(ins+1,iew+1,zt_ind) mz(ins+1,iew+1,zt_ind)];...
        [mx(ins,iew+1,zt_ind) my(ins,iew+1,zt_ind) mz(ins,iew+1,zt_ind)]];
               
        verts((count*4)-3:count*4,:) = v; % done for count = 1
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(ins,iew,zt_ind); 
        
        count = count + 1;
        end
    end    
    
    disp(['Layer ',num2str(zt_ind)])
    toc

    verts = [verts(:,2) verts(:,1) verts(:,3)]; % x and y are swapped here
    draw_patch(verts,faces,colors)
    hold on
    
else % if layer to be plotted contains air cells
           
    for iz = zt_ind:(max(max(topobot_ind)) + 1)
%         tic
        chk = sum(sum(~isnan(A(xs_ind:xn_ind,yw_ind:ye_ind,iz)))); % if 0, then this layer is all air cells and no plotting needed
        if chk > 0
            verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*8,3); % going to be 8 vertices per cell, each vertex is a xyz trio
            faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*5,4); % going to be 5 faces per cell, 4 vertices per face
            colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*5,1); % one color per face
            count = 1;
            for iew = yw_ind:ye_ind-1
                for ins = xs_ind:xn_ind-1
                    v = [[mx(ins,iew,iz) my(ins,iew,iz) mz(ins,iew,iz)];... % these are the 8 vertices for each cell
                    [mx(ins+1,iew,iz) my(ins+1,iew,iz) mz(ins+1,iew,iz)];...
                    [mx(ins+1,iew+1,iz) my(ins+1,iew+1,iz) mz(ins+1,iew+1,iz)];...
                    [mx(ins,iew+1,iz) my(ins,iew+1,iz) mz(ins,iew+1,iz)]                     
                    [mx(ins,iew,iz+1) my(ins,iew,iz+1) mz(ins,iew,iz+1)];...
                    [mx(ins+1,iew,iz+1) my(ins+1,iew,iz+1) mz(ins+1,iew,iz+1)];...
                    [mx(ins+1,iew+1,iz+1) my(ins+1,iew+1,iz+1) mz(ins+1,iew+1,iz+1)];...
                    [mx(ins,iew+1,iz+1) my(ins,iew+1,iz+1) mz(ins,iew+1,iz+1)]];

                    verts( (8*count)-7:8*count ,:) = v; 
                    faces(((count-1)*5)+1:count*5,:) = f5+((count-1)*8); 
                    colors(((count-1)*5)+1:count*5,:) =  A(ins,iew,iz); 

                    count = count + 1;
                end
            end
            % find and remove vertices corresponding to air cells. these are in sets of 8
            tmp = find(isnan(colors(1:5:end,:))).*8; % since one color per face, need to find colors in blocks of 5, then multiply by 8 to include all verts 
            verts(([tmp; tmp-1; tmp-2; tmp-3; tmp-4; tmp-5; tmp-6; tmp-7]),:) = []; % need to remove all 8 verts per cell at once, thus the substractions
            faces = faces(1:size(verts,1)/8*5,:); % truncate to match 5 faces per 8 verts
            colors(isnan(colors)) = [];

            verts = [verts(:,2) verts(:,1) verts(:,3)]; % x and y are swapped here
            draw_patch(verts,faces,colors)
            hold on        

        end
%         disp(['Layer ',num2str(iz)])
%         toc
    end                          
    
end % end plotting top faces
%%
% north face
% tic

verts = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),1);
count = 1;

for ie = yw_ind:ye_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(xn_ind,ie,iz) my(xn_ind,ie,iz) mz(xn_ind,ie,iz)];...
        [mx(xn_ind,ie+1,iz) my(xn_ind,ie+1,iz) mz(xn_ind,ie+1,iz)];...
        [mx(xn_ind,ie+1,iz+1) my(xn_ind,ie+1,iz+1) mz(xn_ind,ie+1,iz+1)];...
        [mx(xn_ind,ie,iz+1) my(xn_ind,ie,iz+1) mz(xn_ind,ie,iz+1)]];
    
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(xn_ind-1,ie,iz); % because of indexing the SW corner of each cell, we need to plot the north face of each cell by indexing ix-1
        
        count = count + 1;                              
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('North Face (c)')
% toc

verts = [verts(:,2) verts(:,1) verts(:,3)]; % x and y are swapped here
draw_patch(verts,faces,colors)
hold on

% south face
% tic

verts = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),1);
count = 1;

for iw = yw_ind:ye_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(xs_ind,iw,iz) my(xs_ind,iw,iz) mz(xs_ind,iw,iz)];...
        [mx(xs_ind,iw+1,iz) my(xs_ind,iw+1,iz) mz(xs_ind,iw+1,iz)];...
        [mx(xs_ind,iw+1,iz+1) my(xs_ind,iw+1,iz+1) mz(xs_ind,iw+1,iz+1)];...
        [mx(xs_ind,iw,iz+1) my(xs_ind,iw,iz+1) mz(xs_ind,iw,iz+1)]];
               
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(xs_ind,iw,iz); 
        
        count = count + 1;                
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('South Face')
% toc

verts = [verts(:,2) verts(:,1) verts(:,3)]; % x and y are swapped here
draw_patch(verts,faces,colors)

% west face
% tic

verts = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),1);
count = 1;

for in = xs_ind:xn_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(in,yw_ind,iz) my(in,yw_ind,iz) mz(in,yw_ind,iz)];...
        [mx(in+1,yw_ind,iz) my(in+1,yw_ind,iz) mz(in+1,yw_ind,iz)];...
        [mx(in+1,yw_ind,iz+1) my(in+1,yw_ind,iz+1) mz(in+1,yw_ind,iz+1)];...
        [mx(in,yw_ind,iz+1) my(in,yw_ind,iz+1) mz(in,yw_ind,iz+1)]];
                  
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(in,yw_ind,iz); 
        
        count = count + 1;            
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('West Face')
% toc

verts = [verts(:,2) verts(:,1) verts(:,3)]; % x and y are swapped here
draw_patch(verts,faces,colors)

% east face
% tic

verts = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind)*4,3); % 4 vertices per cell
faces = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),4); % just one face per cell
colors = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),1);
count = 1;

for is = xs_ind:xn_ind-1
    for iz = zt_ind:zb_ind-1
        v = [[mx(is,ye_ind,iz) my(is,ye_ind,iz) mz(is,ye_ind,iz)];...
        [mx(is+1,ye_ind,iz) my(is+1,ye_ind,iz) mz(is+1,ye_ind,iz)];...
        [mx(is+1,ye_ind,iz+1) my(is+1,ye_ind,iz+1) mz(is+1,ye_ind,iz+1)];...
        [mx(is,ye_ind,iz+1) my(is,ye_ind,iz+1) mz(is,ye_ind,iz+1)]];    
        
        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(is,ye_ind-1,iz);  % iy-1 because we index the SW corner of each cell
        
        count = count + 1;   
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:); % truncate to match number of vertices / 4 
colors(isnan(colors)) = [];

% disp('East Face')
% toc

verts = [verts(:,2) verts(:,1) verts(:,3)]; % x and y are swapped here
draw_patch(verts,faces,colors)

% bottom face 
% tic

verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*4,3); % 4 vertices per cell
faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),4); % just one face per cell
colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),1);
count = 1;
for iew = yw_ind:ye_ind-1
    for ins = xs_ind:xn_ind-1
    v = [[mx(ins,iew,zb_ind) my(ins,iew,zb_ind) mz(ins,iew,zb_ind)];...
    [mx(ins+1,iew,zb_ind) my(ins+1,iew,zb_ind) mz(ins+1,iew,zb_ind)];...
    [mx(ins+1,iew+1,zb_ind) my(ins+1,iew+1,zb_ind) mz(ins+1,iew+1,zb_ind)];...
    [mx(ins,iew+1,zb_ind) my(ins,iew+1,zb_ind) mz(ins,iew+1,zb_ind)]];

    verts((count*4)-3:count*4,:) = v; % done for count = 1
    faces(count,:) = f1 + (count-1)*4;
    colors(count,:) =  A(ins,iew,zb_ind); 

    count = count + 1;
    end
end    

% disp('Bottom Face')
% toc

verts = [verts(:,2) verts(:,1) verts(:,3)]; % x and y are swapped here
draw_patch(verts,faces,colors)

%% make the plot look nice
view(3)
%set(gca,'zdir','reverse')
% zlim(z_lim)

colormap(u.cmap); caxis(u.colim); 

if highlight_dike
    % Insert a narrow bright-green band into the colormap at the dike's
    % resistivity, rather than touching the model geometry -- this way
    % every dike cell (they all share the same log10(resistivity) value)
    % turns green, and everything else keeps the original colormap.
    dike_val = log10(dike_ohmm);
    clim = caxis;
    if dike_val >= clim(1) && dike_val <= clim(2)
        ncolors = 256; % higher resolution so the green band can be narrow and precise
        base_cmap = colormap;
        cmap = interp1(linspace(0,1,size(base_cmap,1)), base_cmap, linspace(0,1,ncolors));

        center_idx = round((dike_val - clim(1)) / (clim(2)-clim(1)) * (ncolors-1)) + 1;
        half_band = max(1, round(dike_band_frac * ncolors));
        lo = max(center_idx-half_band, 1);
        hi = min(center_idx+half_band, ncolors);
        cmap(lo:hi,:) = repmat(dike_color, hi-lo+1, 1);

        colormap(cmap)
    else
        warning('plot_model_3D_RV:dikeOutOfRange', ...
            ['Dike resistivity (%.1f ohm-m, log10=%.2f) falls outside the current ', ...
             'color axis range [%.2f %.2f] -- colormap not modified. Check dike_ohmm ', ...
             'or u.colim.'], dike_ohmm, dike_val, clim(1), clim(2))
    end
end

add_rho_colorbar(u); % your version of this function has no output argument

set(gca,'dataaspectratio',[1 1 1/u.ve])
set(gca,'zdir','reverse')
set(gca,'fontsize',11,'fontname','Helvetica','linewidth',1)
xlabel('Easting (km)','fontsize',12)
ylabel('Northing (km)','fontsize',12)
zlabel('Elevation b.s.l. (km)','fontsize',12)
% clamp indices so we never index below 1 or above the array length
y_lo_ind = max(yw_ind-1, 1);
y_hi_ind = min(ye_ind+1, length(m.cy));
x_lo_ind = max(xs_ind-1, 1);
x_hi_ind = min(xn_ind+1, length(m.cx));

ax_lims = [m.cy(y_lo_ind) m.cy(y_hi_ind) m.cx(x_lo_ind) m.cx(x_hi_ind) ...
    u.zmin*1000+500 u.zmax*1000-500]/1000;

% guard against any pair being non-increasing (e.g. zmin/zmax swapped,
% or a degenerate index range) -- fall back to auto limits rather than
% erroring out and losing the whole plot
pairs_ok = ax_lims(2)>ax_lims(1) && ax_lims(4)>ax_lims(3) && ax_lims(6)>ax_lims(5);
if pairs_ok
    axis(ax_lims)
else
    warning('plot_model_3D:axisLimits', ...
        'Computed axis limits were invalid (check u.zmin/u.zmax order or plot index range) -- using auto limits instead.')
    axis('auto')
end
box on

if use_ortho
    camproj('orthographic')
end

if use_lighting
    camlight('headlight')
    lighting gouraud
    material dull
end

% enable interactive rotate/zoom/pan so you can orient the model by hand
% before saving; vis3d freezes the aspect ratio during rotation so the
% model doesn't visually distort as you drag it around
axis vis3d
rotate3d on
set(gcf,'Toolbar','figure') % make sure the rotate/zoom/pan buttons are visible

if plot_stations
    % Your data file has d.z = 0 for every station (confirmed: range was
    % [0.0 0.0]) -- this is a normal ModEM convention where elevation is
    % NOT stored per station, it's implicit in the model's topography.
    % So the only correct way to get each station's real height is to look
    % at the model directly under that station's own (x,y) and find where
    % the ground actually is -- the first non-air (non-NaN) cell going
    % down. Using d.z, or snapping to a single global m.z index, both
    % ignore location entirely and put every station at the same height,
    % which is exactly what you were seeing.
    stn_x_km = d.y/1000; % swapped to match the x/y convention used above
    stn_y_km = d.x/1000;
    stn_z_km = zeros(size(d.x));

    mx_centers = squeeze(mx(:,1,1)); % model x-centers (km), matches A's 1st dim
    my_centers = squeeze(my(1,:,1)); % model y-centers (km), matches A's 2nd dim

    for i = 1:numel(d.x)
        ix = nearestpoint(d.x(i)/1000, mx_centers);
        iy = nearestpoint(d.y(i)/1000, my_centers);
        col = squeeze(A(ix,iy,:));           % log-resistivity down this column
        iz_surf = find(~isnan(col), 1, 'first'); % first earth cell = topography
        if isempty(iz_surf)
            iz_surf = 1; % fallback: whole column is air (shouldn't normally happen)
        end
        stn_z_km(i) = mz(ix,iy,iz_surf);
    end

    station_lift_km = 0.0005; % 5 m -- nudge markers just above the surface so
    stn_z_km = stn_z_km - station_lift_km; % they don't z-fight with the top patch

    % Diagnostic: check whether the stations actually fall inside the
    % window that axis(ax_lims) just set. If they don't, that's a
    % coordinate/axis-limit issue, not a rendering one -- this is the
    % same failure mode the Python reference script guards against by
    % padding its bounds around the receivers instead of trusting the
    % model's default extent.
    cur_xlim = xlim; cur_ylim = ylim; cur_zlim = zlim;
    outside = stn_x_km < cur_xlim(1) | stn_x_km > cur_xlim(2) | ...
              stn_y_km < cur_ylim(1) | stn_y_km > cur_ylim(2) | ...
              stn_z_km < cur_zlim(1) | stn_z_km > cur_zlim(2);
    if any(outside)
        fprintf(['%d of %d stations fall outside the current axis limits ', ...
            'and would be clipped:\n'], sum(outside), numel(outside));
        fprintf('  axis x: [%.3f %.3f], station x: [%.3f %.3f]\n', ...
            cur_xlim, min(stn_x_km), max(stn_x_km));
        fprintf('  axis y: [%.3f %.3f], station y: [%.3f %.3f]\n', ...
            cur_ylim, min(stn_y_km), max(stn_y_km));
        fprintf('  axis z: [%.3f %.3f], station z: [%.3f %.3f]\n', ...
            cur_zlim, min(stn_z_km), max(stn_z_km));
        fprintf(['  If this is unexpected, check the model/data origin ', ...
            '(the "mesh coordinates and site coordinates do not match" ', ...
            'warning from the loader may be the root cause).\n']);
    end

    scatter3(stn_x_km, stn_y_km, stn_z_km, ...
        'Marker','v','SizeData',70, ...
        'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',0.75)

    % Expand the axes just enough to guarantee every station is inside
    % the visible window, with a small margin so markers aren't drawn
    % right on the clipping edge.
    pad = 0.05; % km
    xlim([min(cur_xlim(1), min(stn_x_km)-pad), max(cur_xlim(2), max(stn_x_km)+pad)])
    ylim([min(cur_ylim(1), min(stn_y_km)-pad), max(cur_ylim(2), max(stn_y_km)+pad)])
    zlim([min(cur_zlim(1), min(stn_z_km)-pad), max(cur_zlim(2), max(stn_z_km)+pad)])
end

plot_geoboundaries(L,d.origin,d.z)

end % end while

end % end function