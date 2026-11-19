clc; clear; close all

%% =========================
% Parameters
% ==========================
base_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');

nx = 161;
nz = 111;

vs_dx = 100;
vs_dz = 50;

x0 = 0;
z0 = 0;

% Levels to show
levels_to_show = [0, 1, 2, 3];

% Paths
patch_dir = fullfile(base_dir, 'patch');
patch_vel_dir = fullfile(base_dir, 'model', 'patch_vs');

vs_file = fullfile(base_dir, 'model', 'vs.bin');
%% =========================
% Display parameters
% ==========================
font_size = 13;
font_name = 'Arial';

C = lines(max(levels_to_show) + 1);
C(4,:) = [0.0000 0.6000 0.2000];   % change purple to green
lw_base = 2.0;
lw_scale = 0.75;
%% =========================
% Read vs.bin
% ==========================
x = x0 + (0:nx-1) * vs_dx;
z = z0 + (0:nz-1) * vs_dz;

fid = fopen(vs_file, 'rb');

assert(fid > 0, ...
    'Failed to open %s.', vs_file);

vs = fread( ...
    fid, ...
    nx*nz, ...
    'float32=>single');

fclose(fid);

assert(numel(vs) == nx*nz, ...
    'vs.bin size mismatch: expected %d floats, got %d', ...
    nx*nz, numel(vs));

vs = reshape(vs, [nx, nz]);

% rows = z, columns = x
vs_img = vs.';

%% =========================
% Figure
% ==========================
figure( ...
    'Color', 'w', ...
    'Position', [100, 100, 950, 480]);

hold on;

%% Plot S-wave velocity model
imagesc( ...
    x./1000, ...
    z./1000, ...
    vs_img./1000);

set(gca, ...
    'YDir', 'reverse', ...
    'FontName', font_name, ...
    'FontSize', font_size, ...
    'LineWidth', 1.0, ...
    'TickDir', 'out');

axis image tight;

xlabel( ...
    'Distance (km)', ...
    'FontName', font_name, ...
    'FontSize', font_size);

ylabel( ...
    'Depth (km)', ...
    'FontName', font_name, ...
    'FontSize', font_size);

%% Gray colormap
colormap(gca, flipud(gray));
% colormap(gca, turbo);
clim([0 4])
%% Colorbar
cb = colorbar;

cb.Label.String = 'S-wave velocity (km/s)';
cb.Label.FontName = font_name;
cb.Label.FontSize = font_size;

cb.FontName = font_name;
cb.FontSize = font_size;

%% =========================
% Plot patches of all levels
% ==========================
legend_handle = gobjects(length(levels_to_show), 1);

for s = 1:length(levels_to_show)

    lev = levels_to_show(s);

    patch_file = fullfile( ...
        patch_dir, ...
        sprintf('patches_level%d.txt', lev));

    if ~exist(patch_file, 'file')
        warning('Patch file not found: %s', patch_file);
        continue;
    end

    pdata = read_patch_file(patch_file);

    if isempty(pdata)
        continue;
    end

    %% Current level style
    col = C(lev + 1, :);

    % Same linewidth rule as the reference code
    lw = lw_base * lw_scale^lev;

    % columns:
    % [box_id i_lo i_hi j_lo j_hi
    %  x_lo x_hi z_lo z_hi dx dz]

    for i = 1:size(pdata, 1)

        x_lo = pdata(i, 6) / 1000;
        x_hi = pdata(i, 7) / 1000;

        z_lo = pdata(i, 8) / 1000;
        z_hi = pdata(i, 9) / 1000;

        % Patch boundary: solid line
        plot( ...
            [x_lo x_hi x_hi x_lo x_lo], ...
            [z_lo z_lo z_hi z_hi z_lo], ...
            '-', ...
            'Color', col, ...
            'LineWidth', lw, ...
            'HandleVisibility', 'off');

    end

    %% Legend proxy
    legend_handle(s) = plot( ...
        nan, nan, ...
        '-', ...
        'Color', col, ...
        'LineWidth', lw, ...
        'DisplayName', sprintf('Level %d', lev));

end

%% =========================
% Legend
% ==========================
lgd = legend( ...
    legend_handle, ...
    'Location', 'bestoutside', ...
    'FontName', font_name, ...
    'FontSize', 11);

lgd.Box = 'on';


box off;
hold off;
%% =========================
% Helper functions
% ==========================
function data = read_patch_file(fname)
% columns:
%
% [box_id i_lo i_hi j_lo j_hi
%  x_lo x_hi z_lo z_hi dx dz]

    fid = fopen(fname, 'r');

    assert(fid > 0, ...
        'Failed to open %s', fname);

    rows = {};

    while true

        tline = fgetl(fid);

        if ~ischar(tline)
            break;
        end

        tline = strtrim(tline);

        if isempty(tline) || startsWith(tline, '#')
            continue;
        end

        rows{end+1,1} = tline; %#ok<AGROW>

    end

    fclose(fid);

    if isempty(rows)

        data = [];
        return;

    end

    data = zeros(numel(rows), 11);

    for i = 1:numel(rows)

        vals = sscanf(rows{i}, '%f');

        if numel(vals) ~= 11

            error( ...
                ['Unexpected column count in %s at data line %d ' ...
                 '(got %d, expect 11).'], ...
                fname, i, numel(vals));

        end

        data(i,:) = vals(:).';

    end

end