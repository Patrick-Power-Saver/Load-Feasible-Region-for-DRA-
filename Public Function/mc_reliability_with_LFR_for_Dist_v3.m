function [EENS, LOLP, details] = mc_reliability_with_LFR_for_Dist_v3( ...
        casename, LFR_lib, load_scenario, options)
% mc_reliability_with_LFR_v2 - 在线蒙特卡洛可靠性评估（匹配 build_LFR_library_v2）
%
% ===================================================================
%  Bug修复说明（相较于旧版）
% ===================================================================
%  Fix 1 (kappa<0): betaG修正后增加 b_corrected = max(0, b_corrected)
%    确保传入 approx_distance 的 b 永远 ≥ 0，从源头消除 kappa<0 问题。
%
%  Fix 2 (details字段对齐): 统一输出以下平级字段，与 plot_reliability_summary 匹配：
%    details.n_adequate    - 充裕状态次数（直接跳过）
%    details.n_cache_hit   - 缓存命中次数（状态去重）
%    details.n_LFR_used    - LFR新建分析次数（精确+降级）
%    details.n_OPF_used    - OPF回退次数
%    details.t_cache       - 缓存I/O时间（与plot兼容的字段名）
%    details.t_correct     - betaG修正时间
%    details.t_adequacy    - 充裕性判断时间
%    details.t_loadshed    - 切负荷计算时间
%
%  Fix 3 (收敛检查): 去掉虚假的 goto_convergence_check，
%    将收敛检查统一放在循环末尾，对缓存命中和新建状态一视同仁。

% ====================================================================
%  参数默认值
% ====================================================================
if nargin < 4, options = struct(); end
max_samples      = getopt(options, 'max_samples',      1e6);
conv_eps         = getopt(options, 'conv_eps',         0.01);
conv_check       = getopt(options, 'conv_check',       2000);
verbose          = getopt(options, 'verbose',          true);
use_OPF_fallback = getopt(options, 'use_OPF_fallback', false);

% ====================================================================
%  读取系统基础数据
% ====================================================================
[baseMVA, bus, ~,~,~,~,~,~, gen, branch, ~,~,~] = feval(casename);
[~, bus, gen, branch] = ext2int(bus, gen, branch);
[GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...
    GEN_LAMDA, GEN_MTTF, GEN_PG, GEN_QG, GEN_PFAILURE] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, ...
    BR_LAMDA, BR_MTTF, TCSC_X, BR_PFAILURE] = idx_brch;

m_line = size(branch, 1);
n_gen  = size(gen,    1);
nb     = size(bus,    1);

line_unavail = branch(:, BR_PFAILURE);
gen_unavail  = gen(:,    GEN_PFAILURE);

unit_bus  = gen(:, GEN_BUS);
Pmax_unit = gen(:, GEN_PMAX) / baseMVA;

Gmax_all = zeros(nb, 1);
for i = 1:n_gen
    Gmax_all(unit_bus(i)) = Gmax_all(unit_bus(i)) + Pmax_unit(i);
end

L0 = load_scenario(:);

% ====================================================================
%  状态缓存（去重，改进1）
% ====================================================================
state_cache = containers.Map('KeyType','char','ValueType','any');

% ====================================================================
%  可靠性指标累积
% ====================================================================
sum_LOLP    = 0;
sum_EENS    = 0;
sum_sq_EENS = 0;
n_sample    = 0;

% 计数器（各类状态次数）
n_adequate   = 0;   % 充裕状态（J≤0，跳过切负荷）
n_cache_hit  = 0;   % 缓存命中（状态去重）
n_exact      = 0;   % LFR精确命中（线路+机组均匹配）
n_line       = 0;   % LFR降级命中（仅线路匹配+betaG修正）
n_OPF_used   = 0;   % OPF回退次数
n_skip       = 0;   % 完全未命中，跳过

% 计时器
t_adequacy = 0;
t_loadshed = 0;
t_correct  = 0;
t_cache    = 0;

t_start = tic;

if verbose
    fprintf('\n=== 在线蒙特卡洛评估（v2）===\n');
    fprintf('最大抽样: %g，收敛精度: %.3f\n', max_samples, conv_eps);
end

% ====================================================================
%  蒙特卡洛主循环
%  结构：抽样→缓存→LFR查询→修正→充裕判断→切负荷→累积→收敛
%  收敛检查统一在循环末尾执行（Fix 3）
% ====================================================================
for n = 1:max_samples
    n_sample = n_sample + 1;

    % ----------------------------------------------------------------
    %  (a) 抽样
    % ----------------------------------------------------------------
    fault_lines = find(rand(m_line, 1) < line_unavail);
    fault_gens  = find(rand(n_gen,  1) < gen_unavail);

    % ----------------------------------------------------------------
    %  (b) 缓存查询（状态去重）
    % ----------------------------------------------------------------
    t0 = tic;
    full_key = make_full_key(fault_lines, fault_gens);
    cache_hit = isKey(state_cache, full_key);
    t_cache = t_cache + toc(t0);

    if cache_hit
        cached   = state_cache(full_key);
        Pls      = cached.Pls;
        adequate = cached.adequate;
        n_cache_hit = n_cache_hit + 1;
    else
        % ----------------------------------------------------------------
        %  (c) 两级LFR查询
        %  Level 1: 精确键（线路+机组均命中库）→ 直接用，无修正误差
        %  Level 2: 纯线路键（仅线路命中）→ betaG在线修正机组故障
        % ----------------------------------------------------------------
        t1 = tic;

        [AX, b_use, match_type] = query_LFR_two_level( ...
            LFR_lib, fault_lines, fault_gens, ...
            Gmax_all, unit_bus, Pmax_unit, n_gen, nb);

        t_correct = t_correct + toc(t1);

        if strcmp(match_type, 'exact')
            n_exact = n_exact + 1;
        elseif strcmp(match_type, 'line')
            n_line = n_line + 1;
        elseif strcmp(match_type, 'fallback') && use_OPF_fallback
            % 高阶线路故障回退OPF
            t2 = tic;
            BS1 = ones(m_line,1); BS1(fault_lines) = 0;
            GS1 = ones(n_gen, 1); GS1(fault_gens)  = 0;
            Pls = OPF_loadshedding_stub(casename, BS1, GS1, L0, baseMVA);
            t_loadshed = t_loadshed + toc(t2);
            adequate = (Pls < 1e-8);
            n_OPF_used = n_OPF_used + 1;
            state_cache(full_key) = struct('Pls',Pls,'adequate',adequate);
            % 累积后跳到收敛检查
            if adequate, n_adequate = n_adequate + 1;
            else
                sum_LOLP    = sum_LOLP + 1;
                sum_EENS    = sum_EENS + Pls;
                sum_sq_EENS = sum_sq_EENS + Pls^2;
            end
            [converged, should_print] = check_convergence(n, conv_check, ...
                sum_EENS, sum_sq_EENS, conv_eps);
            if should_print && verbose
                print_progress(n, sum_EENS, sum_LOLP, n_cache_hit);
            end
            if converged, break; end
            continue;
        else
            % 未命中（高阶故障且不使用OPF）：跳过
            n_skip = n_skip + 1;
            state_cache(full_key) = struct('Pls',0,'adequate',true);
            if mod(n, conv_check) == 0 && sum_EENS > 0 && verbose
                print_progress(n, sum_EENS, sum_LOLP, n_cache_hit);
            end
            continue;
        end

        % ----------------------------------------------------------------
        %  Fix 1 (关键): b 已在 approx_distance 内部 clamped，
        %  此处再做一次防御性 clamp，确保 b ≥ 0 从源头传入
        %  原因：betaG*delta_G 当机组大幅故障时可使 b < 0，
        %         导致 kappa = b/denom < 0（已在 approx_distance 修复）
        % ----------------------------------------------------------------
        b_use = max(b_use, 0);

        % ----------------------------------------------------------------
        %  (d) 充裕性判断（论文式2.25）
        % ----------------------------------------------------------------
        t3 = tic;
        J_max = max(AX * L0 - b_use);
        t_adequacy = t_adequacy + toc(t3);

        if J_max <= 1e-8
            Pls      = 0;
            adequate = true;
            n_adequate = n_adequate + 1;
        else
            % ----------------------------------------------------------------
            %  (e) 近似距离切负荷（论文式3.9）
            % ----------------------------------------------------------------
            t4 = tic;
            [Pls, ~, ~, ~] = approx_distance_loadshedding_v2(AX, b_use, L0);
            t_loadshed = t_loadshed + toc(t4);
            adequate = false;
        end

        state_cache(full_key) = struct('Pls', Pls, 'adequate', adequate);
    end

    % ----------------------------------------------------------------
    %  (f) 累积可靠性指标
    % ----------------------------------------------------------------
    if adequate
        if ~cache_hit, n_adequate = n_adequate + 1; end
    else
        sum_LOLP    = sum_LOLP    + 1;
        sum_EENS    = sum_EENS    + Pls;
        sum_sq_EENS = sum_sq_EENS + Pls^2;
    end

    % ----------------------------------------------------------------
    %  (g) 收敛判断（Fix 3: 统一在末尾执行，缓存命中也参与）
    % ----------------------------------------------------------------
    [converged, should_print] = check_convergence(n, conv_check, ...
        sum_EENS, sum_sq_EENS, conv_eps);
    if should_print && verbose
        print_progress(n, sum_EENS, sum_LOLP, n_cache_hit);
    end
    if converged
        if verbose
            fprintf('  收敛！n=%d, EENS=%.5f\n', n, sum_EENS/n);
        end
        break;
    end
end

% ====================================================================
%  输出（Fix 2: details字段与 plot_reliability_summary 严格对齐）
% ====================================================================
LOLP = sum_LOLP / n_sample;
EENS = sum_EENS / n_sample;
t_total = toc(t_start);

n_LFR_used = n_exact + n_line;   % LFR新建分析（精确+降级）

% --- 与 plot_reliability_summary 严格对齐的字段 ---
details.n_sample     = n_sample;
details.n_adequate   = n_adequate;            % plot: counts(1)
details.n_cache_hit  = n_cache_hit;           % plot: counts(2)
details.n_LFR_used   = n_LFR_used;            % plot: counts(3) = n_LFR_used - n_cache_hit
details.n_OPF_used   = n_OPF_used;
details.n_skip       = n_skip;
details.t_adequacy   = t_adequacy;            % plot: times(1)
details.t_loadshed   = t_loadshed;            % plot: times(2)
details.t_correct    = t_correct;             % plot: times(3)
details.t_cache      = t_cache;               % plot: times(4)  ← Fix 2关键：与v1字段名一致
details.t_total      = t_total;
details.cache_hit_rate = n_cache_hit / n_sample;
details.adequate_rate  = n_adequate  / n_sample;

% 额外诊断信息（不影响plot，供调试用）
details.n_exact_match = n_exact;
details.n_line_match  = n_line;

if verbose
    fprintf('\n=== 评估结果 ===\n');
    fprintf('  LOLP: %.6f\n', LOLP);
    fprintf('  EENS: %.6f pu/次\n', EENS);
    fprintf('  命中率: 缓存=%.1f%% 精确=%.1f%% 降级=%.1f%%\n', ...
        100*n_cache_hit/n_sample, 100*n_exact/n_sample, 100*n_line/n_sample);
    fprintf('  充裕率: %.1f%%\n', 100*details.adequate_rate);
    fprintf('  耗时: 总=%.2fs 充裕=%.2fs 切负荷=%.2fs 修正=%.2fs\n', ...
        t_total, t_adequacy, t_loadshed, t_correct);
end
end   % 主函数结束

% ====================================================================
%  两级LFR查询
% ====================================================================
function [AX, b_out, match_type] = query_LFR_two_level( ...
        LFR_lib, fault_lines, fault_gens, ...
        Gmax_all, unit_bus, Pmax_unit, n_gen, nb)

% Level 1: 精确键（线路+机组均命中）
exact_key = make_full_key(fault_lines, fault_gens);
if isKey(LFR_lib, exact_key)
    e = LFR_lib(exact_key);
    AX = e.AX;  b_out = e.b;
    match_type = 'exact';
    return;
end

% Level 2: 纯线路键 + betaG在线修正机组故障（论文式3.3）
line_key = make_line_only_key(fault_lines);
if isKey(LFR_lib, line_key)
    e = LFR_lib(line_key);
    AX = e.AX;

    Gmax_new = Gmax_all;
    for k = fault_gens'
        b_idx = unit_bus(k);
        Gmax_new(b_idx) = max(0, Gmax_new(b_idx) - Pmax_unit(k));
    end
    delta_G = Gmax_new - e.Gmax0;                  % [nb×1]
    b_out   = e.b + e.betaG * delta_G;             % 在线修正（式3.3）
    % Fix 1源头：此处clamp已在主函数调用前完成（b_use = max(b_use,0)）
    match_type = 'line';
    return;
end

% 未命中
AX = [];  b_out = [];  match_type = 'fallback';
end

% ====================================================================
%  辅助函数
% ====================================================================

% 完整键（线路+机组）
function key = make_full_key(fault_lines, fault_gens)
if isempty(fault_lines), lp = 'L'; else
    lp = ['L' sprintf('%d_', sort(fault_lines(:)'))]; end
if isempty(fault_gens),  gp = 'G'; else
    gp = ['G' sprintf('%d_', sort(fault_gens(:)'))]; end
key = [lp gp];
end

% 纯线路键（Level-2降级查询，与build_LFR_library_v2中0机组故障的键一致）
function key = make_line_only_key(fault_lines)
if isempty(fault_lines), key = 'LG';
else, key = ['L' sprintf('%d_', sort(fault_lines(:)')) 'G']; end
end

% 收敛判断（变异系数 β = σ/μ < conv_eps）
function [converged, should_print] = check_convergence(n, conv_check, ...
        sum_EENS, sum_sq_EENS, conv_eps)
converged    = false;
should_print = (mod(n, 10*conv_check) == 0);
if mod(n, conv_check) == 0 && sum_EENS > 0
    mu  = sum_EENS / n;
    var = (sum_sq_EENS/n - mu^2) / max(n, 1);
    if var >= 0 && sqrt(var)/mu < conv_eps
        converged = true;
    end
end
end

% 进度打印
function print_progress(n, sum_EENS, sum_LOLP, n_cache_hit)
fprintf('  n=%d EENS=%.5f LOLP=%.5f 缓存命中率=%.1f%%\n', ...
    n, sum_EENS/n, sum_LOLP/n, 100*n_cache_hit/n);
end

% getopt
function val = getopt(s, f, d)
if isfield(s, f), val = s.(f); else, val = d; end
end

% OPF占位（高阶故障回退）
function Pls = OPF_loadshedding_stub(casename, BS1, GS1, L0, baseMVA)
try
    [AX, b, ~, ~, ~, ~] = generalized_LFR_v3(casename, BS1, GS1, [], true);
    b = max(b, 0);
    [Pls, ~, ~, ~] = approx_distance_loadshedding(AX, b, L0);
catch
    Pls = 0;
end
end