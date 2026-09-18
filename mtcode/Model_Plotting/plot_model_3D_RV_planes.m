function plot_model_3D_RV_planes(m,d,profile_lines)
%
% Function which renders MT models in 3D as cell faces and vertices.
% Publication-ready version: adds lighting, controllable edge style,
% orthographic projection, labeled log-resistivity colorbar, stations
% plotted on top of the model surface (topography/dike), and optional
% semi-transparent vertical section planes through the model.
%
% Usage: plot_model_3D_RV_planes(m,d)
%        plot_model_3D_RV_planes(m,d,profile_lines)
%
% Inputs: "m" is a standard model structure
%       "d" is a standard data structure
%       "profile_lines"   (optional) Nx4 matrix. Each ROW defines one
%                         vertical section plane as a profile line:
%                             [north1 east1 north2 east2]   (km)
%                         Leave empty/omit for no planes.
%
% Example:
%   plot_model_3D_RV_planes(m,d,[-3 0 3 0; 0 -3 0 3])
%

if nargin < 3 || isempty(profile_lines), profile_lines = zeros(0,4); end

u = user_defaults;
[L] = load_geoboundary_file_list;
close all

z_lim = [u.zmin u.zmax];

plot_stations = 1; % 1 to plot stations, 0 to skip

% ---- publication-style options (edit to taste) ----
edge_style   = 'k';
edge_width   = 0.5;
use_lighting = 0;
use_ortho    = 0;
highlight_dike = 0;
dike_ohmm      = 28;
dike_color     = [0.15 0.95 0.25];
dike_band_frac = 0.01;

% ---- vertical section-plane options ----
plane_alpha = 0.6;
box_alpha   = 1;
plane_ztop_override = -0.01;  % km
plane_zbottom_override = 0.03;% km
plane_nx    = 80;
plane_nz    = 60;
plane_color = [0.5 0.5 0.5];
plane_edge_color = 'k';       % solid lines on the two vertical edges of each plane
plane_edge_width = 1.5;

% should not need to edit below here
f1 = 1:4;
f5 = [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 1 4 8 5];

[mx,my,mz] = meshgrid(m.x./1000,m.y./1000,(m.z)./1000);
mx = permute(mx,[2 1 3]);
my = permute(my,[2 1 3]);
mz = permute(mz,[2 1 3]);
A = log10(m.A);

set_figure_size(1);
set(gcf,'Renderer','opengl');

F = [];
try
    if isfield(m,'cx') && numel(m.cx)==size(A,1)
        xc = m.cx(:);
    else
        xc = (m.x(1:end-1)+m.x(2:end))/2;
    end
    if isfield(m,'cy') && numel(m.cy)==size(A,2)
        yc = m.cy(:);
    else
        yc = (m.y(1:end-1)+m.y(2:end))/2;
    end
    zc = (m.z(1:end-1)+m.z(2:end))/2;
    if numel(xc)==size(A,1) && numel(yc)==size(A,2) && numel(zc)==size(A,3)
        F = griddedInterpolant({xc(:)./1000, yc(:)./1000, zc(:)./1000}, A, 'linear','none');
        fprintf(['Vertical section-plane interpolant built. Model range (km): ', ...
            'North [%.3f %.3f], East [%.3f %.3f], Depth [%.3f %.3f].\n'], ...
            min(xc)/1000, max(xc)/1000, min(yc)/1000, max(yc)/1000, min(zc)/1000, max(zc)/1000);
    else
        warning('plot_model_3D_RV:planeGridMismatch', ...
            'Cell-center coordinate vectors do not match size(A) -- vertical section planes disabled.')
    end
catch ME
    warning('plot_model_3D_RV:planeSetupFailed', ...
        'Could not build interpolant for vertical section planes (%s) -- planes disabled.', ME.message)
    F = [];
end

    function save_current_view()
        [file, savepath] = uiputfile( ...
            {'*.png','PNG image (*.png)'; ...
             '*.pdf','PDF, vector (*.pdf)'; ...
             '*.fig','MATLAB figure, editable (*.fig)'}, ...
            'Save current view as');
        if isequal(file,0)
            return
        end
        fullpath = fullfile(savepath,file);
        [~,~,ext] = fileparts(fullpath);
        switch lower(ext)
            case '.fig'
                savefig(gcf, fullpath);
            case '.pdf'
                exportgraphics(gcf, fullpath, 'ContentType','vector');
            otherwise
                exportgraphics(gcf, fullpath, 'Resolution',300);
        end
        fprintf('Saved current view to: %s\n', fullpath);
    end

    function draw_patch(verts,faces,colors)
        if ischar(edge_style) && strcmpi(edge_style,'none')
            patch('vertices',verts,'faces',faces,'cdata',colors,...
                'facecolor','flat','edgecolor','none','facealpha',box_alpha);
        else
            patch('vertices',verts,'faces',faces,'cdata',colors,...
                'facecolor','flat','edgecolor',edge_style,'linewidth',edge_width,'facealpha',box_alpha);
        end
    end

    function [s_lo,s_hi,ok] = clip_line_to_box(x0,y0,dx,dy,xmin,xmax,ymin,ymax)
        tmin = -inf; tmax = inf;
        if dx ~= 0
            t1 = (xmin-x0)/dx; t2 = (xmax-x0)/dx;
            tmin = max(tmin, min(t1,t2));
            tmax = min(tmax, max(t1,t2));
        elseif x0<xmin || x0>xmax
            s_lo = 0; s_hi = 0; ok = false; return
        end
        if dy ~= 0
            t1 = (ymin-y0)/dy; t2 = (ymax-y0)/dy;
            tmin = max(tmin, min(t1,t2));
            tmax = min(tmax, max(t1,t2));
        elseif y0<ymin || y0>ymax
            s_lo = 0; s_hi = 0; ok = false; return
        end
        if tmin > tmax
            s_lo = 0; s_hi = 0; ok = false;
        else
            s_lo = tmin; s_hi = tmax; ok = true;
        end
    end

    function draw_vertical_plane(n1,e1,n2,e2,nbounds,ebounds,zbounds)
        seg_len = hypot(n2-n1, e2-e1);
        if seg_len == 0
            warning('plot_model_3D_RV:planeZeroLength', ...
                'A requested profile line has identical start and end points (North=%.3f, East=%.3f) -- skipped.', n1, e1);
            return
        end
        dN = (n2-n1)/seg_len;
        dE = (e2-e1)/seg_len;

        [s_lo_box,s_hi_box,ok] = clip_line_to_box(n1,e1,dN,dE, ...
            nbounds(1),nbounds(2),ebounds(1),ebounds(2));
        if ~ok
            warning('plot_model_3D_RV:planeOutside', ...
                ['Profile line from North=%.3f,East=%.3f to North=%.3f,East=%.3f ', ...
                 'does not cross the plotted model extent -- skipped.'], n1, e1, n2, e2);
            return
        end
        s_lo = max(0, s_lo_box);
        s_hi = min(seg_len, s_hi_box);
        if s_lo >= s_hi
            warning('plot_model_3D_RV:planeOutsideSegment', ...
                ['Profile line from North=%.3f,East=%.3f to North=%.3f,East=%.3f ', ...
                 'falls outside the plotted model extent -- skipped.'], n1, e1, n2, e2);
            return
        end

        s = linspace(s_lo,s_hi,plane_nx);
        z = linspace(zbounds(1),zbounds(2),plane_nz);
        [S,Z] = meshgrid(s,z);
        N = n1 + S.*dN;
        E = e1 + S.*dE;

        fprintf(['Plane [North %.4f,East %.4f -> North %.4f,East %.4f]: drawing from ', ...
            's=%.4f to s=%.4f (km along line), z=%.4f to z=%.4f (km).\n'], ...
            n1, e1, n2, e2, s_lo, s_hi, zbounds(1), zbounds(2));

        h = surf(E, N, Z, 'FaceColor',plane_color, 'EdgeColor','none', ...
            'FaceAlpha', plane_alpha);
        fprintf('  -> surf object created (handle valid: %d), XData range [%.4f %.4f], YData range [%.4f %.4f], ZData range [%.4f %.4f]\n', ...
            isgraphics(h), min(E(:)), max(E(:)), min(N(:)), max(N(:)), min(Z(:)), max(Z(:)));

        % Solid vertical edges at each end of the plane (plotted X=East, Y=North)
        n_a = n1 + s_lo*dN; e_a = e1 + s_lo*dE;
        n_b = n1 + s_hi*dN; e_b = e1 + s_hi*dE;
        z_lo = zbounds(1); z_hi = zbounds(2);
        plot3([e_a e_a], [n_a n_a], [z_lo z_hi], '-', ...
            'Color', plane_edge_color, 'LineWidth', plane_edge_width);
        plot3([e_b e_b], [n_b n_b], [z_lo z_hi], '-', ...
            'Color', plane_edge_color, 'LineWidth', plane_edge_width);
    end

while 1
lim_menu = menu('3D Plotting Options', ...
    'Plot with default x and y limits', ...
    'Enter custom x and y limits', ...
    'Save current view as image', ...
    'Return');

if lim_menu == 1
    xind = m.npad(1)+1:m.nx-m.npad(1);
    yind = m.npad(2)+1:m.ny-m.npad(2);
    xn_ind = xind(end);
    xs_ind = xind(1);
    ye_ind = yind(end);
    yw_ind = yind(1);
elseif lim_menu == 2
    prompt = {'x min (km)';'x max (km)';'y min (km)';'y max (km)'};
    dlg_title = 'Plot limits';
    def = {'-20'; '20'; '-20'; '20'};
    num_lines = 1;
    inp = inputdlg(prompt,dlg_title,num_lines,def);

    xn_ind = nearestpoint(str2double(inp{2}),unique(mx(:,:,1)));
    xs_ind = nearestpoint(str2double(inp{1}),unique(mx(:,:,1)));
    ye_ind = nearestpoint(str2double(inp{4}),unique(my(:,:,1)));
    yw_ind = nearestpoint(str2double(inp{3}),unique(my(:,:,1)));
elseif lim_menu == 3
    save_current_view();
    continue
else
    close all
    return
end
clf

zb_ind = min(nearestpoint(u.zmax, unique(mz(1,1,:))), size(A,3));
zt_ind = max(nearestpoint(u.zmin, unique(mz(1,1,:))), 1);

%% top face and air layers
topobot_ind = zeros(xn_ind-xs_ind+1,ye_ind-yw_ind+1);
for iz = 1:zb_ind
    tmp = isnan(A(xs_ind:xn_ind,yw_ind:ye_ind,iz));
    topobot_ind = topobot_ind + tmp;
    if sum(sum(tmp)) ==0
        break
    end
end

if zt_ind > max(max(topobot_ind))
    tic

    verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*4,3);
    faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),4);
    colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),1);

    count = 1;
    for iew = yw_ind:ye_ind-1
        for ins = xs_ind:xn_ind-1
        v = [[mx(ins,iew,zt_ind) my(ins,iew,zt_ind) mz(ins,iew,zt_ind)];...
        [mx(ins+1,iew,zt_ind) my(ins+1,iew,zt_ind) mz(ins+1,iew,zt_ind)];...
        [mx(ins+1,iew+1,zt_ind) my(ins+1,iew+1,zt_ind) mz(ins+1,iew+1,zt_ind)];...
        [mx(ins,iew+1,zt_ind) my(ins,iew+1,zt_ind) mz(ins,iew+1,zt_ind)]];

        verts((count*4)-3:count*4,:) = v;
        faces(count,:) = f1 + (count-1)*4;
        colors(count,:) =  A(ins,iew,zt_ind);

        count = count + 1;
        end
    end

    disp(['Layer ',num2str(zt_ind)])
    toc

    verts = [verts(:,2) verts(:,1) verts(:,3)];
    draw_patch(verts,faces,colors)
    hold on

else

    for iz = zt_ind:(max(max(topobot_ind)) + 1)
        chk = sum(sum(~isnan(A(xs_ind:xn_ind,yw_ind:ye_ind,iz))));
        if chk > 0
            verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*8,3);
            faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*5,4);
            colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*5,1);
            count = 1;
            for iew = yw_ind:ye_ind-1
                for ins = xs_ind:xn_ind-1
                    v = [[mx(ins,iew,iz) my(ins,iew,iz) mz(ins,iew,iz)];...
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
            tmp = find(isnan(colors(1:5:end,:))).*8;
            verts(([tmp; tmp-1; tmp-2; tmp-3; tmp-4; tmp-5; tmp-6; tmp-7]),:) = [];
            faces = faces(1:size(verts,1)/8*5,:);
            colors(isnan(colors)) = [];

            verts = [verts(:,2) verts(:,1) verts(:,3)];
            draw_patch(verts,faces,colors)
            hold on

        end
    end

end
%%
% north face
verts = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind)*4,3);
faces = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),4);
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
        colors(count,:) =  A(xn_ind-1,ie,iz);

        count = count + 1;
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:);
colors(isnan(colors)) = [];

verts = [verts(:,2) verts(:,1) verts(:,3)];
draw_patch(verts,faces,colors)
hold on

% south face
verts = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind)*4,3);
faces = zeros((ye_ind-yw_ind)*(zb_ind-zt_ind),4);
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
faces = faces(1:length(verts)/4,:);
colors(isnan(colors)) = [];

verts = [verts(:,2) verts(:,1) verts(:,3)];
draw_patch(verts,faces,colors)

% west face
verts = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind)*4,3);
faces = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),4);
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
faces = faces(1:length(verts)/4,:);
colors(isnan(colors)) = [];

verts = [verts(:,2) verts(:,1) verts(:,3)];
draw_patch(verts,faces,colors)

% east face
verts = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind)*4,3);
faces = zeros((xs_ind-xn_ind)*(zb_ind-zt_ind),4);
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
        colors(count,:) =  A(is,ye_ind-1,iz);

        count = count + 1;
    end
end

verts(([find(isnan(colors)).*4; (find(isnan(colors)).*4)-1; (find(isnan(colors)).*4)-2; (find(isnan(colors)).*4)-3 ]),:) = [];
faces = faces(1:length(verts)/4,:);
colors(isnan(colors)) = [];

verts = [verts(:,2) verts(:,1) verts(:,3)];
draw_patch(verts,faces,colors)

% bottom face
verts = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind)*4,3);
faces = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),4);
colors = zeros((ye_ind-yw_ind)*(xn_ind-xs_ind),1);
count = 1;
for iew = yw_ind:ye_ind-1
    for ins = xs_ind:xn_ind-1
    v = [[mx(ins,iew,zb_ind) my(ins,iew,zb_ind) mz(ins,iew,zb_ind)];...
    [mx(ins+1,iew,zb_ind) my(ins+1,iew,zb_ind) mz(ins+1,iew,zb_ind)];...
    [mx(ins+1,iew+1,zb_ind) my(ins+1,iew+1,zb_ind) mz(ins+1,iew+1,zb_ind)];...
    [mx(ins,iew+1,zb_ind) my(ins,iew+1,zb_ind) mz(ins,iew+1,zb_ind)]];

    verts((count*4)-3:count*4,:) = v;
    faces(count,:) = f1 + (count-1)*4;
    colors(count,:) =  A(ins,iew,zb_ind);

    count = count + 1;
    end
end

verts = [verts(:,2) verts(:,1) verts(:,3)];
draw_patch(verts,faces,colors)

%% optional vertical section planes
if ~isempty(profile_lines)
    nbounds = sort([mx(xs_ind,1,1) mx(xn_ind,1,1)]);
    ebounds = sort([my(1,yw_ind,1) my(1,ye_ind,1)]);
    zbounds = sort([plane_ztop_override plane_zbottom_override]);

    for irow = 1:size(profile_lines,1)
        draw_vertical_plane(profile_lines(irow,1), profile_lines(irow,2), ...
            profile_lines(irow,3), profile_lines(irow,4), nbounds, ebounds, zbounds);
    end
end

%% make the plot look nice (same coloring as plot_model_3D)
view(3)

colormap(u.cmap); caxis(u.colim);

if highlight_dike
    dike_val = log10(dike_ohmm);
    clim = caxis;
    if dike_val >= clim(1) && dike_val <= clim(2)
        ncolors = 256;
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
             'color axis range [%.2f %.2f] -- colormap not modified.'], ...
            dike_ohmm, dike_val, clim(1), clim(2))
    end
end

add_rho_colorbar(u);

set(gca,'dataaspectratio',[1 1 1/u.ve])
set(gca,'zdir','reverse')
set(gca,'fontsize',11,'fontname','Helvetica','linewidth',1)
xlabel('Easting (km)','fontsize',12)
ylabel('Northing (km)','fontsize',12)
zlabel('Elevation b.s.l. (km)','fontsize',12)

y_lo_ind = max(yw_ind-1, 1);
y_hi_ind = min(ye_ind+1, length(m.cy));
x_lo_ind = max(xs_ind-1, 1);
x_hi_ind = min(xn_ind+1, length(m.cx));

z_range_m = (u.zmax - u.zmin)*1000;
z_pad_m = min(500, 0.1*z_range_m);
ax_lims = [m.cy(y_lo_ind) m.cy(y_hi_ind) m.cx(x_lo_ind) m.cx(x_hi_ind) ...
    u.zmin*1000+z_pad_m u.zmax*1000-z_pad_m]/1000;

pairs_ok = ax_lims(2)>ax_lims(1) && ax_lims(4)>ax_lims(3) && ax_lims(6)>ax_lims(5);
if pairs_ok
    axis(ax_lims)
else
    warning('plot_model_3D:axisLimits', ...
        'Computed axis limits were invalid -- using auto limits instead.')
    axis('auto')
end

if ~isempty(profile_lines)
    plane_zbounds = sort([plane_ztop_override plane_zbottom_override]);
    cur_zlim = zlim;
    zlim([min(cur_zlim(1), plane_zbounds(1)), max(cur_zlim(2), plane_zbounds(2))])
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

axis vis3d
rotate3d on
set(gcf,'Toolbar','figure')

if plot_stations
    stn_x_km = d.y/1000;
    stn_y_km = d.x/1000;
    stn_z_km = zeros(size(d.x));

    mx_centers = squeeze(mx(:,1,1));
    my_centers = squeeze(my(1,:,1));

    for i = 1:numel(d.x)
        ix = nearestpoint(d.x(i)/1000, mx_centers);
        iy = nearestpoint(d.y(i)/1000, my_centers);
        col = squeeze(A(ix,iy,:));
        iz_surf = find(~isnan(col), 1, 'first');
        if isempty(iz_surf)
            iz_surf = 1;
        end
        stn_z_km(i) = mz(ix,iy,iz_surf);
    end

    station_lift_km = 0.0005;
    stn_z_km = stn_z_km - station_lift_km;

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
    end

    scatter3(stn_x_km, stn_y_km, stn_z_km, ...
        'Marker','v','SizeData',70, ...
        'MarkerFaceColor','k','MarkerEdgeColor','w','LineWidth',0.75)

    pad = 0.05;
    xlim([min(cur_xlim(1), min(stn_x_km)-pad), max(cur_xlim(2), max(stn_x_km)+pad)])
    ylim([min(cur_ylim(1), min(stn_y_km)-pad), max(cur_ylim(2), max(stn_y_km)+pad)])
    zlim([min(cur_zlim(1), min(stn_z_km)-pad), max(cur_zlim(2), max(stn_z_km)+pad)])
end

plot_geoboundaries(L,d.origin,d.z)

fprintf('Current renderer: %s\n', get(gcf,'Renderer'));
h_all = findobj(gca,'Type','patch');
h_all = [h_all; findobj(gca,'Type','surface')];
for k = 1:numel(h_all)
    fprintf('Object %d/%d: type=%s, FaceAlpha=%s, Visible=%s\n', ...
        k, numel(h_all), get(h_all(k),'Type'), mat2str(get(h_all(k),'FaceAlpha')), get(h_all(k),'Visible'));
end

end % end while

end % end function