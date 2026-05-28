function LFR_lib = build_LFR_library_v2_Dist(casename, opts)
% build_LFR_library_v2_AC - 离线建立交流负荷可行域库
%
% ===================================================================
%  相较于 build_LFR_library_v2（DC版）的唯一改动
% ===================================================================
%  第一处（~第60行）：
%    DC版：[AX,b,betaG,betaT,Gmax0,lidx] = generalized_LFR_v3(...)
%    AC版：[AX,b,betaG,betaT,Gmax0,lidx] = generalized_LFR_v3_AC(...)
%
%  其余所有代码（参数解析、枚举逻辑、键生成、库存储）完全不变。

if nargin < 2, opts = struct(); end

max_line_order  = getfield_default(opts, 'max_line_order',  2);
max_gen_order   = getfield_default(opts, 'max_gen_order',   0);
max_total_order = getfield_default(opts, 'max_total_order', 2);
P_ref           = getfield_default(opts, 'P_ref',           []);
use_tight_theta = getfield_default(opts, 'use_tight_theta', true);
verbose         = getfield_default(opts, 'verbose',         true);


[baseMVA, bus, ~, ~, ~, ~, ~, ~, gen, branch, ~, ~, ~] = feval(casename);
SW_NO = find(branch(:, 9) == 0);
[GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...
    GEN_LAMDA, GEN_MTTF, GEN_PG, GEN_QG, GEN_PFAILURE] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, BR_LAMDA, ...
    BR_MTTF, TCSC_X, BR_PFAILURE] = idx_brch;

m_line = size(branch(branch(:, 9) == 1, :),1);
n_gen  = size(gen,    1);
nb     = size(bus,    1);

if isempty(P_ref)
    P_ref = sum(bus(:, 3)) / baseMVA;
end

if verbose
    fprintf('=== 离线建立交流负荷可行域库（AC-GSF版）===\n');
    fprintf('算例: %s  线路数: %d  机组数: %d\n', casename, m_line, n_gen);
    fprintf('枚举配置：线路≤%d阶，机组≤%d阶，总阶数≤%d\n', ...
            max_line_order, max_gen_order, max_total_order);
    fprintf('PBC-Θ改进: %s，P_ref=%.4f pu\n', yesno(use_tight_theta), P_ref);
end

LFR_lib     = containers.Map('KeyType','char','ValueType','any');
total_built = 0;
total_skip  = 0;
t_start     = tic;

for line_order = 0 : max_line_order
    line_combos = enum_combos(1:m_line, line_order);

    for gen_order = 0 : max_gen_order
        if line_order + gen_order > max_total_order, continue; end

        gen_combos = enum_combos(1:n_gen, gen_order);
        n_pairs    = length(line_combos) * length(gen_combos);
        if verbose
            fprintf('  线路%d阶 × 机组%d阶 = %d 个状态\n', ...
                    line_order, gen_order, n_pairs);
        end

        for lc = 1 : length(line_combos)
            fault_lines = line_combos{lc};
            BS1 = ones(m_line, 1);
            if ~isempty(fault_lines), BS1(fault_lines) = 0; end

            for gc = 1 : length(gen_combos)
                fault_gens = gen_combos{gc};
                GS1 = ones(n_gen, 1);
                if ~isempty(fault_gens), GS1(fault_gens) = 0; end

                key = make_key(fault_lines, fault_gens);
                if isKey(LFR_lib, key)
                    total_skip = total_skip + 1;
                    continue;
                end

                try
                    % ★ 唯一改动：调用交流版 LFR 建模函数 ★
                    [AX, b, betaG, betaT, Gmax0, lidx] = ...
                        R_generalized_LFR_v4_AC(casename, BS1, GS1, ...
                                              P_ref, use_tight_theta);
%                       [AX, b, betaG, betaT, Gmax0, lidx, topo_info] = ...
%                           R_generalized_LFR_dist_v1(casename, BS1, GS1, SW_NO, P_ref, opts);
                    entry.AX              = AX;
                    entry.b               = b;
                    entry.betaG           = betaG;
                    entry.betaT           = betaT;
                    entry.Gmax0           = Gmax0;
                    entry.line_idx_active = lidx;
                    entry.fault_lines     = fault_lines;
                    entry.fault_gens      = fault_gens;
                    entry.line_order      = line_order;
                    entry.gen_order       = gen_order;

                    LFR_lib(key) = entry;
                    total_built  = total_built + 1;

                catch ME
                    if verbose
                        warning('AC-LFR建立失败（%s）：%s', key, ME.message);
                    end
                end
            end
        end
    end
end

elapsed = toc(t_start);
if verbose
    fprintf('\n交流离线建库完成：成功=%d  跳过=%d  耗时=%.2fs\n', ...
            total_built, total_skip, elapsed);
end
end

%% 辅助函数（与 build_LFR_library_v2 完全相同）
function combos = enum_combos(items, k)
if k == 0; combos = {[]}; return; end
mat    = nchoosek(items, k);
combos = cell(size(mat,1), 1);
for i = 1:size(mat,1); combos{i} = mat(i,:); end
end

function key = make_key(fault_lines, fault_gens)
if isempty(fault_lines), lpart = 'L';
else, lpart = ['L' sprintf('%d_', sort(fault_lines(:)'))]; end
if isempty(fault_gens),  gpart = 'G';
else, gpart = ['G' sprintf('%d_', sort(fault_gens(:)'))]; end
key = [lpart gpart];
end

function val = getfield_default(s, field, default)
if isfield(s, field) && ~isempty(s.(field)), val = s.(field);
else, val = default; end
end

function s = yesno(flag)
if flag, s='是'; else, s='否'; end
end
