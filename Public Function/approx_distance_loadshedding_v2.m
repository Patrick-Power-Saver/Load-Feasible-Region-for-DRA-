function [Pls, Lc, kappa_min, active_boundary_idx] = approx_distance_loadshedding_v2(AX, b, L0)
% approx_distance_loadshedding - 基于负荷可行域的近似距离切负荷模型
%
% 论文依据：第3章第3.4节，式(3.4)-(3.9)
%
% ===================================================================
%  Bug修复说明（相较于旧版）
% ===================================================================
%  根本原因：
%    b_corrected = b_init + betaG * delta_G
%    当机组大容量故障时，delta_G << 0，betaG*delta_G 为较大负值，
%    可使部分 b_i < 0。此时 kappa_i = b_i/denom_i < 0，
%    触发 "kappa_min 异常" 警告并回退到全切（kappa=0），结果错误。
%
%  物理解释：
%    LFR满足"包含原点"假设（L=0时所有约束均满足，即b≥0）。
%    b_i < 0 是线性betaG修正的过度保守结果，并非真实物理约束。
%
%  修复策略：
%    Step 0: b = max(b, 0)  ← 强制LFR包含原点（修复关键）
%    此后 denom>0 && b>=0 → kappa = b/denom ≥ 0，永不出现负值。
%    同时在孤岛处理后重新检查充裕性，避免孤岛约束被移除后
%    kappa计算出现all-inf的退化情况。
% ===================================================================

nb = length(L0);
L0 = L0(:);

% ====================================================================
%  Step 0: 修复关键 —— 强制 b ≥ 0（确保LFR包含原点）
%  物理依据：L=0（不供任何负荷）永远是可行的，即 AX*0=0 ≤ b，
%  故 b_i ≥ 0 是LFR的必要条件。b_i < 0 是betaG线性修正的数值伪影。
% ====================================================================
b = max(b, 0);

% ====================================================================
%  Step 1: 充裕性预判（论文式2.25，J = max(ΛL₀ − β)）
% ====================================================================
J = AX * L0 - b;
if max(J) <= 1e-8
    Pls = 0;  Lc = L0;  kappa_min = 1.0;  active_boundary_idx = [];
    return;
end

% ====================================================================
%  Step 2: 孤岛修正（论文3.4.2节）
%  识别 e_i·L ≤ 0 型边界（单节点系数为1，且 b_i=0）
%  → 对应无电源孤岛节点，该节点负荷必须全部切除
% ====================================================================
Pls_island   = 0;
island_nodes = false(nb, 1);

for i = 1:size(AX, 1)
    row = AX(i, :);
    % 孤岛判定：只有一列非零且为正，且 b_i = 0
    nonzero_pos = find(row > 1e-8);
    nonzero_all = find(abs(row) > 1e-8);
    if length(nonzero_pos) == 1 && length(nonzero_all) == 1 && b(i) < 1e-8
        node_idx = nonzero_pos;
        if ~island_nodes(node_idx)
            island_nodes(node_idx) = true;
            Pls_island = Pls_island + L0(node_idx);
        end
    end
end

% 降维：孤岛节点负荷清零，移除仅涉及孤岛节点的约束行
L0_red = L0;
L0_red(island_nodes) = 0;

if sum(L0_red) < 1e-8
    % 所有节点都是孤岛
    Pls = sum(L0);  Lc = zeros(nb,1);  kappa_min = 0;
    active_boundary_idx = [];
    return;
end

valid_rows = true(size(AX, 1), 1);
for i = 1:size(AX, 1)
    nz = find(abs(AX(i,:)) > 1e-8);
    if ~isempty(nz) && all(island_nodes(nz))
        valid_rows(i) = false;   % 该约束只涉及孤岛节点，移除
    end
end
AX_red = AX(valid_rows, :);
b_red  = b(valid_rows);          % b已在Step 0 clamped ≥ 0

% ====================================================================
%  Step 3: 降维后重新检查充裕性
%  若违反约束全部来自孤岛（已在Step 2处理），则剩余负荷无需额外切除
% ====================================================================
if isempty(AX_red)
    Pls = Pls_island;
    Lc  = L0_red;
    kappa_min = 1.0;
    active_boundary_idx = [];
    return;
end

J_red = AX_red * L0_red - b_red;
if max(J_red) <= 1e-8
    % 违反约束仅由孤岛节点引起，非孤岛节点充裕
    Pls = Pls_island;
    Lc  = L0_red;
    kappa_min = 1.0;
    active_boundary_idx = [];
    return;
end

% ====================================================================
%  Step 4: 近似距离模型（论文式3.5-3.9）
%
%  近似可供给负荷 Lap = κ·L0_red 在原点到 L0_red 的连线上，
%  且落在某条可行域边界上：
%    AX_red_i · (κ·L0_red) = b_red_i
%    → κ_i = b_red_i / (AX_red_i · L0_red)
%
%  由于 b_red ≥ 0（Step 0保证），denom_i > 0 时 κ_i ≥ 0，
%  永不出现负值，无需警告/回退。
%
%  κmin = 最小的 κ_i ∈ [0,1]，对应最紧约束边界
%  Pls_remaining = (1 - κmin) · sum(L0_red)   （论文式3.9）
% ====================================================================
denom      = AX_red * L0_red;    % [rows×1]
kappa_vals = inf(size(denom));

for i = 1:length(denom)
    if denom(i) > 1e-8           % 正分母：射线可以到达该边界
        % b_red(i) ≥ 0 已由Step 0保证，故 kappa_i ≥ 0
        kappa_vals(i) = b_red(i) / denom(i);
    end
    % denom ≤ 0：射线沿L0方向不穿过该边界（约束在L0方向自然满足），跳过
end

% 仅保留 κ ∈ [0, 1] 的边界（κ>1 表示L0在该边界内侧，即约束未违反）
kappa_vals(kappa_vals > 1 + 1e-8) = inf;

[kappa_min, rel_idx] = min(kappa_vals);

% ====================================================================
%  Step 5: 处理极端退化情况
%  （理论上Step 0后不应出现，保留作为防御性代码）
% ====================================================================
if isinf(kappa_min)
    % 所有约束方向均不与射线相交（极少数数值边界情况）
    % 保守处理：全切
    kappa_min = 0;
    active_boundary_idx = [];
else
    kappa_min = max(0, min(1, kappa_min));   % 数值截断保护
    valid_indices = find(valid_rows);
    if rel_idx <= length(valid_indices)
        active_boundary_idx = valid_indices(rel_idx);
    else
        active_boundary_idx = [];
    end
end

% ====================================================================
%  Step 6: 计算切负荷量（论文式3.9）
%  Lc = κmin · L0_red
%  Pls_remaining = (1 - κmin) · sum(L0_red)
% ====================================================================
Lc_red        = kappa_min * L0_red;
Pls_remaining = sum(L0_red) - sum(Lc_red);
Lc            = Lc_red;               % 孤岛节点已为0

Pls = Pls_island + Pls_remaining;

end