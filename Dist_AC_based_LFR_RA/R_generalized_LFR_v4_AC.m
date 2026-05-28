function [AX, b, betaG, betaT, Gmax_node, line_idx_active] = ...
        R_generalized_LFR_v4_AC(casename, BS1, GS1, P_ref, use_tight_theta)
% generalized_LFR_v3_AC - 交流负荷可行域建模（基于 Lossy ACOPF 的 GSF 矩阵）
%
% ===================================================================
%  相较于 generalized_LFR_v3（DC/PTDF版）的改动说明
% ===================================================================
%  【唯一改动】第4节：将 PTDF 矩阵替换为交流 GSF 矩阵
%
%  DC原版（第4节）：
%    Mat = PTDF(branch_cur, bus, ones(1,nl), nl, nb, 1);   % [nl×nb]
%
%  AC新版（第4节）：
%    从线性化 Lossy ACOPF（与 Lossy_SOCP_opf_LMP 完全一致）推导：
%      Pf = GSF_PP_F * Pbus + GSF_PQ_F * Qbus + 0.5 * Ploss
%    其中 Pbus = Pg - Pd，Qbus = Qg - Qd
%
%    为保持 LFR 变量空间仍为纯有功负荷向量 L（[nb×1]），
%    引入节点功率因数比 pf_ratio_i = Qd_i / Pd_i，将无功负荷表示为
%    有功负荷的比例量：Qd = diag(pf_ratio) * L
%    从而有效 GSF：
%      Mat = GSF_PP_F + GSF_PQ_F .* pf_ratio'            % [nl×nb]
%
%    网损项（0.5*Ploss）：离线建库阶段忽略（保守处理，使LFR偏小，
%    对可靠性评估结果偏保守/安全，可通过 loss_factor 参数调节）。
%
%  第1、2、3、5、6、7、8节：与原版完全相同，无任何改动。
%
% ===================================================================
%  输入（与原版完全相同）
% ===================================================================
%   casename        - 算例名称
%   BS1             - 支路状态向量 (1=正常, 0=故障) [m_line×1]
%   GS1             - 发电机状态向量 (1=正常, 0=故障) [n_gen×1]
%   P_ref           - 功率平衡参考值(pu)，默认=系统峰值负荷之和
%   use_tight_theta - true=使用PBC-Θ（默认），false=原始Θ
%
% 输出（与原版完全相同）
%   AX              - 负荷可行域系数矩阵 Λ  [rows×nb]
%   b               - 右端向量 β            [rows×1]
%   betaG           - 发电容量修正系数矩阵  [rows×nb]
%   betaT           - 线路容量系数矩阵      [rows×nl_active]
%   Gmax_node       - 节点级最大发电容量    [nb×1]
%   line_idx_active - 有效线路在原始branch中的索引

% ====================================================================
%  0. 默认参数
% ====================================================================
if nargin < 4 || isempty(P_ref),           P_ref = [];          end
if nargin < 5 || isempty(use_tight_theta), use_tight_theta = true; end

% ====================================================================
%  1. 加载索引常量（与原版完全相同）
% ====================================================================
[PQ, PV, REF, NOTE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, VMAX, VMIN, PGMAX,PGMIN,QGMAX,QGMIN,PG,QG,LOLP,LOLE,LOLF,LOLD,EDNS,EENS,...
    LAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus;

[GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...
    GEN_LAMDA,GEN_MTTF,GEN_PG,GEN_QG,GEN_PFAILURE] = idx_gen;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, ...
    BR_LAMDA, BR_MTTF, TCSC_X, BR_PFAILURE, PF_col, QF_col, PT_col, QT_col, ...
    MU_SF, MU_ST, KT_MAXANGLE_TPST, KT_MINANGLE_TPST, ...
    KT_MAX_X_TCSC, KT_MIN_X_TCSC] = idx_brch;

% ====================================================================
%  2. 读取算例数据，移除故障设备（与原版完全相同）
% ====================================================================
[baseMVA, bus, ~, ~, ~, ~, ~, ~, gen, branch, ~, ~, ~] = feval(casename);
[~, bus, gen, branch] = ext2int(bus, gen, branch);
if size(GS1,2) > size(GS1,1), GS1 = GS1'; end
if size(BS1,2) > size(BS1,1), BS1 = BS1'; end

line_idx_active = find(BS1 == 1);
fault_line      = find(BS1 == 0);
fault_gen       = find(GS1 == 0);

branch_cur = branch;
gen_cur    = gen;
if ~isempty(fault_line), branch_cur(fault_line, :) = []; end
if ~isempty(fault_gen),  gen_cur(fault_gen,   :) = []; end

% ====================================================================
%  3. 基本参数（与原版完全相同）
% ====================================================================
nb = size(bus,       1);
nl = size(branch_cur,1);
ng = size(gen_cur,   1);

PF_MAX   = branch_cur(:, RATE) / baseMVA;   % [nl×1]
FrBranch = branch_cur(:, F_BUS);
ToBranch = branch_cur(:, T_BUS);

unit_bus  = gen_cur(:, GEN_BUS);
Pmax_unit = gen_cur(:, GEN_PMAX) / baseMVA;
Pmin_unit = gen_cur(:, GEN_PMIN) / baseMVA;

Node_Unit = zeros(nb, ng);
for i = 1:ng
    Node_Unit(unit_bus(i), i) = 1;
end

Gmax_node = Node_Unit * Pmax_unit;   % [nb×1]
Gmin_node = Node_Unit * Pmin_unit;   % [nb×1]

if isempty(P_ref)
    P_ref = sum(bus(:, PD)) / baseMVA * 1.05;
end

% ====================================================================
%  4. ★★★ 交流 GSF 矩阵（唯一修改点，替换原版 PTDF）★★★
%
%  原版（DC）：
%    Mat = PTDF(branch_cur, bus, ones(1,nl), nl, nb, 1);
%
%  新版（AC）：从 Lossy SOCP ACOPF 推导，与 Lossy_SOCP_opf_LMP 完全一致
%    线性化交流潮流：
%      Pf = GSF_PP_F * Pbus + GSF_PQ_F * Qbus + 0.5 * Ploss
%    约束：-PF_MAX ≤ Pf ≤ PF_MAX
%    令 Pbus = Pg - L，Qbus = Qg - pf_ratio.*L（等功率因数负荷假设）
%    整理后等效矩阵：
%      Mat = GSF_PP_F + GSF_PQ_F .* pf_ratio'        [nl×nb]
%    对应发电侧矩阵：
%      Mat_G = GSF_PP_F + GSF_PQ_F .* qpf_ratio'     [nl×nb]
% ====================================================================
if nl == 0
    % 全部线路故障（孤岛）：无潮流约束（与原版相同）
    AX    = -eye(nb);
    b     = zeros(nb, 1);
    betaG = zeros(nb, nb);
    betaT = zeros(nb, 0);
    return;
end

% ---------- 4a. 构建线性化交流潮流灵敏度矩阵（Lossy SOCP 框架）----------
%
% ★ Bug 修复：列格式不匹配问题 ★
% makeYbus 是 MATPOWER 函数，内部使用 MATPOWER 的 branch 列索引：
%   TAP=9, SHIFT=10, BR_STATUS=11
% 但 RTS79_1 的 branch 矩阵列顺序不同（由 idx_brch 定义）：
%   TAP=7, SHIFT_ANGLE=8, BR_STATUS=9
% 解决方案：用已加载的 RTS79_1 列常量（TAP/SHIFT_ANGLE/BR_STATUS），
% 重排构建符合 MATPOWER 格式的 branch_mpc，再传入 makeYbus。
% bus 矩阵无需转换：makeYbus 仅访问 GS(5) 和 BS(6)，两种格式相同。
branch_mpc          = zeros(nl, 11);
branch_mpc(:, 1:5)  = branch_cur(:, 1:5);               % f,t,r,x,b（两格式相同）
branch_mpc(:, 6:8)  = repmat(branch_cur(:, RATE), 1, 3); % rateA/B/C
branch_mpc(:, 9)    = branch_cur(:, TAP);                % tap（RTS79_1 col 7）
branch_mpc(:, 10)   = branch_cur(:, SHIFT_ANGLE);        % shift（RTS79_1 col 8）
branch_mpc(:, 11)   = branch_cur(:, BR_STATUS);          % status（RTS79_1 col 9）

[Ybus, ~, ~] = makeYbus(baseMVA, bus, branch_mpc);
Ysc   = 1 ./ (branch_cur(:, BR_R) - 1j * branch_cur(:, BR_X));  % 支路导纳 [nl×1]

Gbus_m = real(Ybus);   % [nb×nb]
Bbus_m = imag(Ybus);   % [nb×nb]

% 与 Lossy_SOCP_opf_LMP 完全相同的 GSF 推导
GP = Gbus_m;
BD = diag(sum(Bbus_m));
B_mod = Bbus_m - BD;          % B (去掉对角元)
BP    = -B_mod;                % BP = -B'
BQ    = -Bbus_m;               % BQ = -B
GQ    = -Gbus_m;               % GQ ≈ -G

Xp_dlpf = [BP  GP];           % [nb × 2nb]
Xq_dlpf = [GQ  BQ];           % [nb × 2nb]

CJ = full([Xp_dlpf; Xq_dlpf]);   % [2nb × 2nb]

% 参考节点处理：与 Lossy_SOCP_opf_LMP 保持一致
%
% ★ Bug 修复：避免调用 bustypes() ★
% MATPOWER 的 bustypes() 内部读取 gen(:, GEN_STATUS=8)，
% 但 RTS79_1 的 gen 只有 11 列且 status 在第 6 列，会读错数据。
% 解决方案：直接从 bus 矩阵识别参考节点（BUS_TYPE=2, REF=3 已由 idx_bus 加载）。
ref_bus = find(bus(:, BUS_TYPE) == REF);
if isempty(ref_bus), ref_bus = 1; end   % 安全回退
CJ(ref_bus, ref_bus) = 0;

CJ_inv  = pinv(CJ);
C_PVm   = CJ_inv(nb+1:end, 1:nb);       % ∂Vm²/∂P  [nb×nb]
C_QVm   = CJ_inv(nb+1:end, nb+1:end);   % ∂Vm²/∂Q  [nb×nb]
C_PVa   = CJ_inv(1:nb,     1:nb);       % ∂Va/∂P   [nb×nb]
C_QVa   = CJ_inv(1:nb,     nb+1:end);   % ∂Va/∂Q   [nb×nb]

fr = FrBranch;   % [nl×1]
to = ToBranch;   % [nl×1]
re_Ysc = real(Ysc);   % [nl×1] 支路电导
im_Ysc = imag(Ysc);   % [nl×1] 支路电纳

% 支路有功潮流对节点注入的灵敏度（与 Lossy_SOCP_opf_LMP 相同）
%   GSF_PP_F(l,i) = ∂Pf_l/∂P_i   [nl×nb]
%   GSF_PQ_F(l,i) = ∂Pf_l/∂Q_i   [nl×nb]
GSF_PP_F = (C_PVm(fr,:) - C_PVm(to,:)) .* re_Ysc ...
         - (C_PVa(fr,:) - C_PVa(to,:)) .* im_Ysc;
GSF_PQ_F = (C_QVm(fr,:) - C_QVm(to,:)) .* re_Ysc ...
         - (C_QVa(fr,:) - C_QVa(to,:)) .* im_Ysc;

% ---------- 4b. 统一灵敏度矩阵 Mat（[nl×nb]）---------------------------
%
% ★ 设计修正：AX 与 b/betaG 必须使用同一灵敏度矩阵 ★
%
% 原实现区分了两个矩阵：
%   AX    用 Mat   = GSF_PP_F + GSF_PQ_F .* pf_ratio'   (pf_ratio  = Qd/Pd)
%   b/betaG 用 Mat_G = GSF_PP_F + GSF_PQ_F .* qpf_ratio' (qpf_ratio ≈ Qg/Pg)
%
% 这导致在线 MC Level-2 betaG 修正时：
%   充裕判断 J = AX*L0 - b_corrected
%   AX*L0  由 Mat   决定（两策略相同）
%   b_corrected 由 Mat_G 决定（两策略差异存在但被 AX 侧抵消）
% → 两种 Θ 策略的 J 值相同 → LOLP/EENS 无法区分
%
% 修复：令 Mat_G = Mat（统一功率因数假设 Qg/Pg = Qd/Pd）
%   则 AX 和 b/betaG 使用同一矩阵，Level-2 修正后 b_corrected 差异
%   可以正确体现在 J = AX*L0 - b_corrected 中：
%     PBC-Θ 的 b 更小（约束更紧）→ J 更大 → 更难充裕 → LOLP 偏高
%   与 DC PTDF 天然统一的处理方式完全类比。
Pd_bus = bus(:, PD) / baseMVA;
Qd_bus = bus(:, QD) / baseMVA;

pf_ratio = zeros(nb, 1);
load_mask = (Pd_bus > 1e-6);
pf_ratio(load_mask) = Qd_bus(load_mask) ./ Pd_bus(load_mask);

% 统一灵敏度矩阵：同时用于 AX（负荷侧约束）和 b/betaG（发电侧修正）
Mat   = GSF_PP_F + GSF_PQ_F .* pf_ratio';   % [nl×nb]
Mat_G = Mat;   % ★ 统一，不再使用独立的 qpf_ratio 估算 ★

% ====================================================================
%  5. 第一类负荷可行域边界（基于 Mat，其余与原版完全相同）
% ====================================================================
AX_type1 = [-Mat; Mat];   % [2nl×nb]，与原版相同结构

if use_tight_theta && P_ref < sum(Gmax_node) - 1e-8
    % --- PBC-Θ改进版（使用 Mat_G 作为发电侧灵敏度）---
    [b_lower_G, b_upper_G, betaG_lower, betaG_upper] = ...
        pbc_theta_AC(Mat, Mat_G, Gmax_node, Gmin_node, P_ref);
    b_type1     = [PF_MAX + b_lower_G; PF_MAX + b_upper_G];
    betaG_type1 = [betaG_lower;        betaG_upper];
else
    % --- 原始Θ算子（使用 Mat_G 替代原版 Mat 用于发电侧计算）---
    b_lower     = PF_MAX - Mat_G*Gmin_node + max(-Mat_G,0)*(Gmax_node - Gmin_node);
    b_upper     = PF_MAX + Mat_G*Gmin_node + max( Mat_G,0)*(Gmax_node - Gmin_node);
    b_type1     = [b_lower; b_upper];
    betaG_type1 = [max(-Mat_G,0); max(Mat_G,0)];
end

betaT_type1 = [eye(nl); eye(nl)];   % [2nl×nl]

% ====================================================================
%  6-8. 第二类边界 + 负荷下界 + 合并（与原版完全相同）
% ====================================================================
[AX_type2, b_type2, betaG_type2, betaT_type2] = ...
    build_type2_boundaries(bus, branch_cur, Gmax_node, PF_MAX, nb, nl);

AX_lb    = -eye(nb);
b_lb     = zeros(nb, 1);
betaG_lb = zeros(nb, nb);
betaT_lb = zeros(nb, nl);

AX    = [AX_type1;    AX_type2;    AX_lb];
b     = [b_type1;     b_type2;     b_lb];
betaG = [betaG_type1; betaG_type2; betaG_lb];
betaT = [betaT_type1; betaT_type2; betaT_lb];

valid = any(abs(AX) > 1e-10, 2) | (b > 1e-10);
AX    = AX(valid, :);
b     = b(valid);
betaG = betaG(valid, :);
betaT = betaT(valid, :);

end   % ← 主函数结束

% ====================================================================
%  子函数1：PBC-Θ（交流版，区分负荷侧 Mat 与发电侧 Mat_G）
%
%  改动：接收额外参数 Mat_G（发电侧灵敏度），其余逻辑与原版相同
%  原版 pbc_theta 仅使用单一 Mat，本函数用 Mat_G 做发电侧优化
% ====================================================================
function [b_lower_G, b_upper_G, betaG_lower, betaG_upper] = ...
        pbc_theta_AC(Mat, Mat_G, Gmax_node, Gmin_node, P_ref)

nl = size(Mat,  1);
nb = size(Mat_G,2);

b_upper_G   = zeros(nl, 1);
b_lower_G   = zeros(nl, 1);
betaG_upper = zeros(nl, nb);
betaG_lower = zeros(nl, nb);

for l = 1:nl
    coeff_G = Mat_G(l, :)';   % 使用发电侧灵敏度做约束上界优化

    [val_u, bG_u] = greedy_max(coeff_G, Gmax_node, Gmin_node, P_ref);
    b_upper_G(l)      = val_u;
    betaG_upper(l, :) = bG_u;

    [val_lo, bG_lo] = greedy_max(-coeff_G, Gmax_node, Gmin_node, P_ref);
    b_lower_G(l)      = val_lo;
    betaG_lower(l, :) = bG_lo;
end
end

% ====================================================================
%  子函数2：贪心算子（与原版 generalized_LFR_v3 完全相同）
% ====================================================================
function [obj_val, betaG_row] = greedy_max(coeff, Gmax, Gmin, P_ref)
nb = length(coeff);
G  = Gmin;
betaG_row = zeros(1, nb);

current_sum = sum(Gmin);
remaining   = max(P_ref - current_sum, 0);

[sorted_c, sort_idx] = sort(coeff, 'descend');

for k = 1:nb
    idx = sort_idx(k);
    if sorted_c(k) <= 0, break; end

    cap_k  = Gmax(idx) - Gmin(idx);
    if cap_k < 1e-10, continue; end

    fill_k    = min(cap_k, remaining);
    G(idx)    = Gmin(idx) + fill_k;

    if fill_k >= cap_k - 1e-10
        betaG_row(idx) = sorted_c(k);
    end

    remaining = remaining - fill_k;
    if remaining < 1e-10, break; end
end

obj_val = dot(coeff, G);
end

% ====================================================================
%  子函数3：第二类负荷可行域边界（与原版完全相同）
% ====================================================================
function [AX2, b2, betaG2, betaT2] = build_type2_boundaries( ...
        bus, branch_cur, Gmax_node, PF_MAX, nb, nl)

F_BUS_col = 1;
T_BUS_col = 2;

AX2 = zeros(0, nb); b2 = zeros(0,1);
betaG2 = zeros(0, nb); betaT2 = zeros(0, nl);

adj = false(nb, nb);
for k = 1:nl
    f = branch_cur(k, F_BUS_col);
    t = branch_cur(k, T_BUS_col);
    adj(f,t) = true; adj(t,f) = true;
end

visited    = false(1, nb);
components = {};
for s = 1:nb
    if ~visited(s)
        comp = bfs_nodes(s, adj, nb);
        components{end+1} = comp;
        visited(comp) = true;
    end
end

for c = 1:length(components)
    comp = components{c};

    ext_lines = [];
    for k = 1:nl
        f    = branch_cur(k, F_BUS_col);
        t    = branch_cur(k, T_BUS_col);
        f_in = ismember(f, comp);
        t_in = ismember(t, comp);
        if xor(f_in, t_in)
            ext_lines(end+1) = k;
        end
    end

    if isempty(ext_lines)
        row_AX    = zeros(1, nb); row_AX(comp)    = 1;
        row_betaG = zeros(1, nb); row_betaG(comp) = 1;
        row_betaT = zeros(1, nl);
        row_b     = sum(Gmax_node(comp));
    elseif length(ext_lines) == 1
        row_AX    = zeros(1, nb); row_AX(comp)    = 1;
        row_betaG = zeros(1, nb); row_betaG(comp) = 1;
        row_betaT = zeros(1, nl); row_betaT(ext_lines) = 1;
        row_b     = sum(Gmax_node(comp)) + PF_MAX(ext_lines);
    else
        continue
    end

    AX2    = [AX2;    row_AX];
    b2     = [b2;     row_b];
    betaG2 = [betaG2; row_betaG];
    betaT2 = [betaT2; row_betaT];
end
end

function comp = bfs_nodes(start, adj, nb)
visited = false(1, nb);
queue   = start;
visited(start) = true;
comp    = [];
while ~isempty(queue)
    cur   = queue(1); queue = queue(2:end);
    comp(end+1) = cur;
    nbrs = find(adj(cur,:) & ~visited);
    visited(nbrs) = true;
    queue = [queue, nbrs];
end
end
