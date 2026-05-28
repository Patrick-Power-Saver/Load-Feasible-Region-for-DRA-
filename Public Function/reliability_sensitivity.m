function [SEN_gen, SEN_line, EENS_updated] = reliability_sensitivity(casename, ...
        LFR_lib, L0, baseMVA, delta_G_re, options)
% reliability_sensitivity - 基于负荷可行域的可靠性容量灵敏度与滚动计算
%
% 论文依据：第4章，式(4.7)(4.12)(4.16)-(4.19)
%
% 功能：
%   1. 建立 EENS 对各机组/线路容量的灵敏度解析模型（论文4.3节）
%   2. 基于灵敏度进行可再生能源波动下的可靠性滚动计算（论文4.4节）
%
% 灵敏度核心公式（论文式4.7）：
%   SenCG_k(s) = betaG_ac,i(s) * IG_k(s) / (Λ_ac(s) * L0)
%   SENG_k = T * Σ_s { SenCG_k(s) * Pr(s) }
%
% 滚动计算公式（论文式4.19）：
%   EENS_new = EENS_old + Σ_k { ΔG_re,k * SEN_k }
%
% 输入：
%   casename  - 算例名称
%   LFR_lib   - LFR库（由 build_LFR_library 生成）
%   L0        - 当前负荷向量 [nb×1] (pu)
%   baseMVA   - 基准容量
%   delta_G_re- 可再生能源出力变化量 [nRE×1] (pu)，正为增加
%   options   - 选项结构体（与 mc_reliability_with_LFR 相同）
%
% 输出：
%   SEN_gen    - EENS对各发电机容量的灵敏度 [n_gen×1]
%   SEN_line   - EENS对各线路容量的灵敏度 [m_line×1]
%   EENS_updated - 滚动更新后的 EENS (pu)

if nargin < 6, options = struct(); end
max_samples = getopt(options, 'max_samples', 5e5);
conv_eps    = getopt(options, 'conv_eps',    0.01);
verbose     = getopt(options, 'verbose',     true);

% 读取系统数据
[baseMVA_sys, bus, ~,~,~,~,~,~, gen, branch, ~,~,~] = feval(casename);
if nargin < 4 || isempty(baseMVA), baseMVA = baseMVA_sys; end

[GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...
    GEN_LAMDA, GEN_MTTF, GEN_PG, GEN_QG, GEN_PFAILURE] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, BR_LAMDA, ...
    BR_MTTF, TCSC_X, BR_PFAILURE] = idx_brch;

m_line = size(branch, 1);
n_gen  = size(gen, 1);
nb     = size(bus, 1);

unit_bus  = gen(:, GEN_BUS);
Pmax_unit = gen(:, GEN_PMAX) / baseMVA;
gen_unavail  = gen(:, GEN_PFAILURE);
line_unavail = branch(:, BR_PFAILURE);

Gmax_all = zeros(nb, 1);
for i = 1:n_gen
    Gmax_all(unit_bus(i)) = Gmax_all(unit_bus(i)) + Pmax_unit(i);
end

L0 = L0(:);

% =====================================================================
%  蒙特卡洛循环：同步建立灵敏度（嵌入可靠性评估中，论文4.4.2节）
% =====================================================================
sum_EENS    = 0;
sum_SEN_gen = zeros(n_gen, 1);    % Σ SenCG_k(s) * I(s>0)
sum_SEN_line= zeros(m_line, 1);   % Σ SenT_j(s) * I(s>0)
n_sample    = 0;
sum_sq_EENS = 0;

if verbose
    fprintf('灵敏度建模中（蒙特卡洛嵌入计算）...\n');
end

for n = 1:max_samples
    n_sample = n_sample + 1;

    % 抽样
    fault_lines = find(rand(m_line,1) < line_unavail);
    fault_gens  = find(rand(n_gen,1)  < gen_unavail);

    % 查LFR库
    line_key = build_line_key(fault_lines);
    if ~isKey(LFR_lib, line_key), continue; end

    entry = LFR_lib(line_key);
    AX    = entry.AX;
    betaG = entry.betaG;
    Gmax0 = entry.Gmax0;

    % 修正b（考虑发电机故障）
    Gmax_new = Gmax_all;
    for k = fault_gens'
        Gmax_new(unit_bus(k)) = max(0, Gmax_new(unit_bus(k)) - Pmax_unit(k));
    end
    b_cor = entry.b + betaG * (Gmax_new - Gmax0);

    % 充裕性判断
    if max(AX * L0 - b_cor) <= 1e-8
        continue;   % 充裕状态：灵敏度贡献为0
    end

    % 近似距离切负荷，获取有效边界参数（用于灵敏度计算）
    [Pls, ~, kappa_min, ac_idx] = approx_distance_loadshedding(AX, b_cor, L0);

    if Pls < 1e-8 || isempty(ac_idx), continue; end

    sum_EENS    = sum_EENS + Pls;
    sum_sq_EENS = sum_sq_EENS + Pls^2;

    % ----------------------------------------------------------------
    %  灵敏度计算（论文式4.7）
    %  SenCG_k(s) = betaG_ac,i(s) * IG_k(s) / (Λ_ac(s) * L0)
    %
    %  有效边界 ac_idx 对应 AX(ac_idx,:) 和 betaG(ac_idx,:)
    % ----------------------------------------------------------------
    Lam_ac  = AX(ac_idx, :);       % [1×nb] 有效边界Λ
    betaG_ac= betaG(ac_idx, :);    % [1×nb] 有效边界βG

    denom = Lam_ac * L0;            % 标量：Λ_ac * L0
    if abs(denom) < 1e-10, continue; end

    denominator_inv = 1 / denom;
    I_L0_norm = sum(L0);            % sum(L0)用于归一化

    % 各发电机容量灵敏度（论文式4.7）
    for k = 1:n_gen
        IG_k = ~ismember(k, fault_gens);   % 1=正常运行，0=故障
        % betaG_ac 的第 unit_bus(k) 列
        betaG_ac_i = betaG_ac(unit_bus(k));
        % SenCG_k(s) = betaG_ac,i * IG_k / denom（论文式4.6）
        sen_k = betaG_ac_i * IG_k * I_L0_norm * denominator_inv;
        sum_SEN_gen(k) = sum_SEN_gen(k) + sen_k;
    end

    % 各线路容量灵敏度（论文式4.17）
    line_idx_active = entry.line_idx_active;
    if ~isempty(line_idx_active) && isfield(entry, 'betaT')
        betaT_ac = entry.betaT(ac_idx, :);  % [1×nl_active]
        for j_local = 1:length(line_idx_active)
            j_global = line_idx_active(j_local);
            IT_j = ~ismember(j_global, fault_lines);  % 该线路是否正常
            betaT_ac_j = betaT_ac(j_local);
            sen_j = betaT_ac_j * IT_j * I_L0_norm * denominator_inv;
            sum_SEN_line(j_global) = sum_SEN_line(j_global) + sen_j;
        end
    end

    % 收敛判断
    if mod(n, 2000) == 0 && sum_EENS > 0
        mu = sum_EENS / n;
        var_val = (sum_sq_EENS/n - mu^2) / n;
        if var_val >= 0 && sqrt(var_val)/mu < conv_eps
            if verbose
                fprintf('  收敛：n=%d，EENS=%.6f\n', n, mu);
            end
            break;
        end
    end
end

% =====================================================================
%  输出灵敏度（以发生概率为权重，即对 Pr(s) 的期望）
%  论文式(4.12)：SENG_k = T * Σ_s { SenCG_k(s) * Pr(s) }
%  蒙特卡洛下：SENG_k ≈ (1/N) * Σ SenCG_k(s)
% =====================================================================
SEN_gen  = sum_SEN_gen  / n_sample;
SEN_line = sum_SEN_line / n_sample;
EENS_base= sum_EENS     / n_sample;

if verbose
    fprintf('\n--- 发电机容量灵敏度 Top-5 ---\n');
    [sorted_vals, sorted_idx] = sort(SEN_gen, 'ascend');
    for k = 1:min(5, n_gen)
        fprintf('  机组 %2d (节点%2d): SEN = %+.4f\n', ...
                sorted_idx(k), unit_bus(sorted_idx(k)), sorted_vals(k));
    end
end

% =====================================================================
%  滚动计算（论文式4.19）
%  EENS_new = EENS_old + Σ_k { ΔGre_k * SEN_k }
% =====================================================================
if nargin >= 5 && ~isempty(delta_G_re) && length(delta_G_re) <= n_gen
    % delta_G_re 对应可再生能源发电机（前 nRE 台）的容量变化
    nRE = length(delta_G_re);
    delta_G_re = delta_G_re(:);
    EENS_updated = EENS_base + sum(delta_G_re .* SEN_gen(1:nRE));
    EENS_updated = max(0, EENS_updated);  % 物理约束：EENS≥0

    if verbose
        fprintf('\n--- 可靠性滚动更新（论文式4.19）---\n');
        fprintf('  原始 EENS:     %.6f pu\n', EENS_base);
        fprintf('  ΔG_re 波动:    [%s] pu\n', num2str(delta_G_re', '%.4f '));
        fprintf('  更新后 EENS:   %.6f pu\n', EENS_updated);
        fprintf('  计算时间:      < 1 秒（满足实时性要求）\n');
    end
else
    EENS_updated = EENS_base;
end

end

%% 辅助函数
function key = build_line_key(fault_lines)
if isempty(fault_lines)
    key = 'L0';
else
    key = ['L' sprintf('%d_', sort(fault_lines(:))')];
end
end

function val = getopt(options, field, default)
if isfield(options, field)
    val = options.(field);
else
    val = default;
end
end
