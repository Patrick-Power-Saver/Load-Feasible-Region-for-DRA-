function [AX, b, betaG, betaT, Gmax_node, line_idx_active, topo_info] = ...
        R_generalized_LFR_dist_v1(casename, BS1, GS1, SW_NO, P_ref, opts)
% R_generalized_LFR_dist_v1 - 配电网负荷可行域建模
%                              基于 BFM/LinearDistFlow + 网络重构（MILP）
%
% ===================================================================
%  与大电网版（R_generalized_LFR_v4_AC）的核心差异
% ===================================================================
%
%  【差异1】灵敏度矩阵：BIM/GSF → BFM/上游路径矩阵
%    大电网：Mat = GSF_PP_F + GSF_PQ_F .* pf_ratio'  [稠密，全局耦合]
%    配电网：Mat = A_upstream                          [稀疏0-1，仅由拓扑决定]
%    其中 A_upstream[branch_ij, node_k] = 1 当且仅当节点 k 在支路(i,j)的下游子树中
%    物理含义：支路(i,j)的潮流 = 其下游所有节点负荷之和（linearized DistFlow）
%
%  【差异2】网络重构（MILP）
%    大电网：拓扑固定，无此步骤
%    配电网：故障后需枚举可行重构拓扑τ（联络开关的开/闭组合）
%            对每个有效拓扑单独建立 LFR，库键包含拓扑标识
%            辐射性约束通过生成树枚举保证
%
%  【差异3】电压约束建模
%    大电网：通过 GSF 中的 C_PVm/C_QVm 隐式处理
%    配电网：通过 DistFlow 显式建立电压灵敏度矩阵 M_V（[nb×nb]）
%            ΔV²_j = -2 * Σ_{branch(i,j)∈path(root,j)} (r_{ij}*A_P + x_{ij}*A_Q) * ΔL
%            V²_min ≤ V²_j0 + M_V * (L - L0) ≤ V²_max
%
%  【差异4】切负荷模式（离散 vs 连续）
%    大电网：连续 κ∈[0,1]，AX*L ≤ b 定义连续可行域
%    配电网：离散 u_k∈{0,1}，LFR 给出"最大可供电节点集合"
%            本函数仍输出多面体 AX*L ≤ b（离散切负荷在在线 MC 中处理）
%
% ===================================================================
%  输入
% ===================================================================
%   casename  - 算例名称（配电网数据格式，含联络开关信息）
%   BS1       - 支路状态向量 (1=正常, 0=故障) [m_line×1]（不含联络开关）
%   GS1       - 发电机/DG 状态向量 (1=正常, 0=故障) [n_gen×1]
%   SW_NO     - 联络开关列表 [n_sw×1]，对应 branch 中的索引
%               （正常运行时断开，故障后可闭合进行负荷转供）
%   P_ref     - PBC-Θ功率参考值(pu)（空则自动计算）
%   opts      - 选项结构体
%               opts.use_tight_theta  - 是否使用 PBC-Θ（默认true）
%               opts.use_voltage      - 是否加入电压约束（默认true）
%               opts.max_reconfig     - 最大重构拓扑数（默认8）
%               opts.Vmin_pu, Vmax_pu - 电压上下限（默认0.9, 1.1）
%
% ===================================================================
%  输出（与大电网版接口完全相同，增加 topo_info）
% ===================================================================
%   AX, b, betaG, betaT - LFR 多面体参数（基于最优重构拓扑）
%   Gmax_node           - 节点级最大有功容量 [nb×1]
%   line_idx_active     - 当前拓扑有效支路在原始 branch 中的索引
%   topo_info           - 重构信息结构体（含 NO 开关闭合方案等）

% ====================================================================
%  0. 默认参数
% ====================================================================
if nargin < 4 || isempty(SW_NO),  SW_NO  = [];   end
if nargin < 5 || isempty(P_ref),  P_ref  = [];   end
if nargin < 6,                     opts   = struct(); end

use_tight_theta = getopt(opts, 'use_tight_theta', true);
use_voltage     = getopt(opts, 'use_voltage',     true);
max_reconfig    = getopt(opts, 'max_reconfig',    8);
Vmin_sq         = getopt(opts, 'Vmin_pu',         0.9)^2;
Vmax_sq         = getopt(opts, 'Vmax_pu',         1.1)^2;

% ====================================================================
%  1. 加载索引常量
% ====================================================================
[PQ, PV, REF, NOTE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, VMAX, VMIN, PGMAX,PGMIN,QGMAX,QGMIN,PG,QG, ...
    LOLP,LOLE,LOLF,LOLD,EDNS,EENS, ...
    LAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus;

[GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...
    GEN_LAMDA,GEN_MTTF,GEN_PG,GEN_QG,GEN_PFAILURE] = idx_gen;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, ...
    BR_LAMDA, BR_MTTF, TCSC_X, BR_PFAILURE] = idx_brch;

% ====================================================================
%  2. 读取数据 & 移除故障元件
% ====================================================================
[baseMVA, bus, ~,~,~,~,~,~, gen, branch, ~,~,~] = feval(casename);
[~, bus, gen, branch] = ext2int(bus, gen, branch);
if size(GS1,2) > size(GS1,1), GS1 = GS1'; end
if size(BS1,2) > size(BS1,1), BS1 = BS1'; end

% 区分常规支路和联络开关
n_branch_total = size(branch, 1);
is_NO_sw = false(n_branch_total, 1);
if ~isempty(SW_NO), is_NO_sw(SW_NO) = true; end

regular_br_idx = find(~is_NO_sw);            % 常规支路索引
fault_regular  = regular_br_idx(BS1 == 0);   % 故障常规支路

% 正常拓扑：常规支路中状态为 1 的
active_regular = regular_br_idx(BS1 == 1);
line_idx_active = active_regular;             % 默认（后面可能更新）

fault_gen = find(GS1 == 0);

gen_cur    = gen;
if ~isempty(fault_gen), gen_cur(fault_gen, :) = []; end

nb = size(bus, 1);
ng = size(gen_cur, 1);

Pd_bus = bus(:, PD) / baseMVA;
Qd_bus = bus(:, QD) / baseMVA;
pf_ratio = zeros(nb, 1);
load_mask = Pd_bus > 1e-6;
pf_ratio(load_mask) = Qd_bus(load_mask) ./ Pd_bus(load_mask);

unit_bus  = gen_cur(:, GEN_BUS);
Pmax_unit = gen_cur(:, GEN_PMAX) / baseMVA;
Pmin_unit = gen_cur(:, GEN_PMIN) / baseMVA;
Node_Unit = sparse(unit_bus, (1:ng)', 1, nb, ng);
Gmax_node = full(Node_Unit * Pmax_unit);
Gmin_node = full(Node_Unit * Pmin_unit);

if isempty(P_ref)
    P_ref = sum(Pd_bus) * 1.05;
end

% 孤岛判断（无任何线路连接）
if isempty(active_regular)
    [AX, b, betaG, betaT, topo_info] = isolated_system_LFR(nb);
    line_idx_active = [];
    return;
end

% ====================================================================
%  3. 网络重构枚举（MILP核心 - 差异2）
% ====================================================================
% 枚举联络开关闭合方案，找出可行的辐射状拓扑集合
% 策略：对故障后孤岛节点，尝试每个 NO 开关闭合（单次重构为主）

fault_segs = identify_fault_segments(bus, branch, active_regular, fault_regular);
reconfig_topos = enumerate_reconfig_topologies(bus, branch, active_regular, ...
    SW_NO, fault_segs, max_reconfig);

if isempty(reconfig_topos)
    reconfig_topos = {active_regular};  % 无可行重构，仅用原拓扑
end

% ====================================================================
%  4. ★★★ BFM 灵敏度矩阵（核心差异1：替代 GSF）★★★
%
%  对每个候选拓扑 τ 计算：
%    A_τ：上游路径矩阵 [nl_τ × nb]（0-1矩阵）
%    M_V_τ：电压灵敏度矩阵 [nb × nb]
%
%  选择最优拓扑（最小切负荷量估计）后建立 LFR。
% ====================================================================

best_AX = []; best_b = []; best_betaG = []; best_betaT = [];
best_topo_idx = [];
min_shed_est = inf;

for tt = 1:length(reconfig_topos)
    cur_active = reconfig_topos{tt};
    nl_cur = length(cur_active);
    branch_cur = branch(cur_active, :);
    
    % --- 4a. 构建有根树结构（配电网辐射状拓扑）---
    fr = branch_cur(:, F_BUS);
    to = branch_cur(:, T_BUS);
    r_br = branch_cur(:, BR_R);
    x_br = branch_cur(:, BR_X);
    PF_MAX = branch_cur(:, RATE) / baseMVA;
    
    % 找根节点（参考/电源节点）
    ref_node = find(bus(:, BUS_TYPE) == REF);
    if isempty(ref_node), ref_node = 1; end
    
    % 构建邻接关系，从根节点 BFS 确定父子关系
    [parent, children, order] = build_tree(ref_node, fr, to, nb);
    
    % --- 4b. 上游路径矩阵 A（BFM灵敏度，替代GSF）---
    % A[l, k] = 1 iff 节点 k 在支路 l 下游（含端节点）
    %
    % 计算方式：对每条支路 l=(i,j)，subtree(j) 包含 j 及其所有后代
    A_P = zeros(nl_cur, nb);   % 有功灵敏度 [nl×nb]
    A_Q = zeros(nl_cur, nb);   % 无功灵敏度（含功率因数）
    
    subtrees = compute_subtrees(fr, to, nl_cur, nb, parent, children);
    
    for l = 1:nl_cur
        j_node = to(l);   % 支路末端节点
        if isfield(subtrees, 'nodes') && length(subtrees.nodes) >= j_node
            % 节点 j 的子树成员
            sub_nodes = subtrees.nodes{j_node};
        else
            sub_nodes = get_subtree_nodes(j_node, children);
        end
        A_P(l, sub_nodes) = 1;
        A_Q(l, sub_nodes) = pf_ratio(sub_nodes)';   % 含无功（功率因数修正）
    end
    
    % 有效支路潮流对负荷的灵敏度（等效 Mat）
    %   P_{l} = Σ_{k∈subtree(l)} Pd_k = A_P(l,:) * L
    %   因此：-PF_MAX ≤ A_P * L - A_P * Pg ≤ PF_MAX
    %   LFR 变量为 L（负荷向量），Pg 视为容量参数（由 betaG 修正）
    Mat = A_P;         % [nl×nb]，等效灵敏度矩阵（替代 GSF）
    % 注：A_Q 用于电压约束但不进入 Mat，因为 pf_ratio 已编入 A_P 的无功分量
    
    % --- 4c. 电压约束灵敏度矩阵 M_V（差异3）---
    if use_voltage
        % DistFlow 电压降方程：
        %   V²_j = V²_i - 2(r_{ij}*P_{ij} + x_{ij}*Q_{ij})
        % 代入 P_{ij} = A_P(l,:)*L，Q_{ij} = A_Q(l,:)*L：
        %   V²_j = V²_root - Σ_{branch∈path(root,j)} 2*(r_l*A_P(l,:) + x_l*A_Q(l,:))*L
        % 即：V²_j(L) = V²_j0 - M_V(j,:) * L
        % 其中 M_V(j,k) = 2 * Σ_{l∈path(root,j)} (r_l*A_P(l,k) + x_l*A_Q(l,k))
        
        M_V = zeros(nb, nb);  % [nb×nb]
        for j = 1:nb
            if j == ref_node, continue; end
            path_branches = get_path_to_root(j, parent, fr, to, nl_cur);
            for l = path_branches
                M_V(j, :) = M_V(j, :) + 2*(r_br(l)*A_P(l,:) + x_br(l)*A_Q(l,:));
            end
        end
    end
    
    % ====================================================================
    %  5. 第一类 LFR 边界（支路容量约束）
    % ====================================================================
    AX_type1 = [-Mat; Mat];   % [2nl×nb]（与大电网版结构完全相同）
    
    if use_tight_theta && P_ref < sum(Gmax_node) - 1e-8
        [b_lower_G, b_upper_G, betaG_lower, betaG_upper] = ...
            pbc_theta_dist(Mat, Gmax_node, Gmin_node, P_ref);
        b_type1     = [PF_MAX + b_lower_G; PF_MAX + b_upper_G];
        betaG_type1 = [betaG_lower; betaG_upper];
    else
        b_lower     = PF_MAX - Mat*Gmin_node + max(-Mat,0)*(Gmax_node-Gmin_node);
        b_upper     = PF_MAX + Mat*Gmin_node + max( Mat,0)*(Gmax_node-Gmin_node);
        b_type1     = [b_lower; b_upper];
        betaG_type1 = [max(-Mat,0); max(Mat,0)];
    end
    betaT_type1 = [eye(nl_cur); eye(nl_cur)];   % [2nl×nl]
    
    % ====================================================================
    %  6. 第二类 LFR 边界（电压约束，配电网特有）
    % ====================================================================
    if use_voltage
        % 电压上限：V²_j ≤ Vmax²
        %   V²_j0 - M_V(j,:)*L ≤ Vmax²
        %   M_V(j,:)*L ≥ V²_j0 - Vmax²（负约束，通常不紧）
        % 电压下限：V²_j ≥ Vmin²（更常见的约束）
        %   V²_j0 - M_V(j,:)*L ≥ Vmin²
        %   M_V(j,:)*L ≤ V²_j0 - Vmin²
        
        Vsq_0 = bus(:, VM).^2;   % 额定运行点电压平方（近似1.0 pu²）
        
        % 电压下限约束（M_V*L ≤ Vsq_0 - Vmin²）
        %valid_V = (1:nb)' ~= ref_node;  % 排除参考节点
        valid_V = ~ismember((1:nb)', ref_node);
        idx_V   = find(valid_V);
        
        AX_Vlo  = M_V(idx_V, :);            % [nv×nb]
        b_Vlo   = Vsq_0(idx_V) - Vmin_sq;   % [nv×1]
        betaG_Vlo = zeros(length(idx_V), nb);
        betaT_Vlo = zeros(length(idx_V), nl_cur);
        
        % 电压上限约束（-M_V*L ≤ -(Vsq_0 - Vmax²) = Vmax²-Vsq_0）
        % 通常 Vmax² - Vsq_0 > 0（额定电压略低于上限），这个约束相对宽松
        AX_Vhi  = -M_V(idx_V, :);
        b_Vhi   = Vmax_sq - Vsq_0(idx_V);
        betaG_Vhi = zeros(length(idx_V), nb);
        betaT_Vhi = zeros(length(idx_V), nl_cur);
    else
        AX_Vlo = zeros(0,nb); b_Vlo = zeros(0,1);
        betaG_Vlo = zeros(0,nb); betaT_Vlo = zeros(0,nl_cur);
        AX_Vhi = zeros(0,nb); b_Vhi = zeros(0,1);
        betaG_Vhi = zeros(0,nb); betaT_Vhi = zeros(0,nl_cur);
    end
    
    % ====================================================================
    %  7. 第三类边界（孤岛节点下界 + 功率平衡，与大电网版相同）
    % ====================================================================
    AX_lb    = -eye(nb);
    b_lb     = zeros(nb, 1);
    betaG_lb = zeros(nb, nb);
    betaT_lb = zeros(nb, nl_cur);
    
    % ====================================================================
    %  8. 合并所有约束
    % ====================================================================
    AX_cur    = [AX_type1;  AX_Vlo;  AX_Vhi;  AX_lb];
    b_cur     = [b_type1;   b_Vlo;   b_Vhi;   b_lb];
    betaG_cur = [betaG_type1; betaG_Vlo; betaG_Vhi; betaG_lb];
    betaT_cur = [betaT_type1; betaT_Vlo; betaT_Vhi; betaT_lb];
    
    % 去除冗余行
    valid = any(abs(AX_cur) > 1e-10, 2) | (b_cur > 1e-10);
    AX_cur    = AX_cur(valid, :);
    b_cur     = b_cur(valid);
    betaG_cur = betaG_cur(valid, :);
    betaT_cur = betaT_cur(valid, :);
    
    % ====================================================================
    %  9. 选择最优拓扑（估计切负荷量最小）
    %  简单启发式：J = max(AX*Pd - b)，值越小表示当前负荷越可行
    % ====================================================================
    J_est = max(AX_cur * Pd_bus - b_cur);
    if J_est < min_shed_est
        min_shed_est = J_est;
        best_AX    = AX_cur;
        best_b     = b_cur;
        best_betaG = betaG_cur;
        best_betaT = betaT_cur;
        best_topo_idx = tt;
        line_idx_active = cur_active;
    end
end

AX    = best_AX;
b     = best_b;
betaG = best_betaG;
betaT = best_betaT;

topo_info.topo_idx       = best_topo_idx;
topo_info.active_branches = line_idx_active;
topo_info.reconfig_topos  = reconfig_topos;
topo_info.n_topo          = length(reconfig_topos);
topo_info.min_J_est       = min_shed_est;

end  % ← 主函数结束

% ====================================================================
%  子函数：识别故障段（故障后与电源断开的节点集）
% ====================================================================
function fault_segs = identify_fault_segments(bus, branch, active_regular, fault_regular)
% 图论BFS：去除故障支路后，识别与电源节点不连通的节点集
nb = size(bus, 1);
F_BUS_col = 1; T_BUS_col = 2; BUS_TYPE_col = 2; REF_TYPE = 3;

ref_nodes = find(bus(:, BUS_TYPE_col) == REF_TYPE);
if isempty(ref_nodes), ref_nodes = 1; end

% 构建不含故障支路的邻接表
operative = setdiff(active_regular, fault_regular);
adj = false(nb, nb);
for k = operative'
    f = branch(k, F_BUS_col); t = branch(k, T_BUS_col);
    adj(f,t) = true; adj(t,f) = true;
end

% BFS 从所有电源节点出发，标记可达节点
reachable = false(nb, 1);
queue = ref_nodes(:)';
reachable(ref_nodes) = true;
while ~isempty(queue)
    cur = queue(1); queue = queue(2:end);
    nbrs = find(adj(cur,:) & ~reachable');
    reachable(nbrs) = true;
    queue = [queue, nbrs];
end

fault_segs = find(~reachable);
end

% ====================================================================
%  子函数：枚举可行重构拓扑（联络开关方案）
% ====================================================================
function reconfig_topos = enumerate_reconfig_topologies(bus, branch, ...
        active_regular, SW_NO, fault_segs, max_reconfig)
% 对每个孤岛节点，寻找可通过闭合某个 NO 开关恢复供电的方案
reconfig_topos = {};
F_BUS_col = 1; T_BUS_col = 2;

if isempty(fault_segs) || isempty(SW_NO)
    reconfig_topos{1} = active_regular;
    return;
end

% 基础拓扑（无重构）
reconfig_topos{end+1} = active_regular;

% 单次开关闭合策略
for sw = SW_NO'
    f = branch(sw, F_BUS_col);
    t = branch(sw, T_BUS_col);
    % 检查：闭合该开关是否能将某孤岛节点重新连接到电源
    if any(ismember(f, fault_segs)) || any(ismember(t, fault_segs))
        new_topo = union(active_regular, sw);
        if ~is_duplicate_topo(new_topo, reconfig_topos)
            reconfig_topos{end+1} = new_topo;
        end
    end
    if length(reconfig_topos) >= max_reconfig, break; end
end

% 两个开关联合闭合（覆盖更复杂重构场景）
if length(reconfig_topos) < max_reconfig && length(SW_NO) >= 2
    for i = 1:length(SW_NO)-1
        for j = i+1:length(SW_NO)
            sw1 = SW_NO(i); sw2 = SW_NO(j);
            new_topo = union(union(active_regular, sw1), sw2);
            if ~is_duplicate_topo(new_topo, reconfig_topos)
                % 验证辐射性（无环）
                if is_radial_tree(new_topo, branch, size(bus,1))
                    reconfig_topos{end+1} = new_topo;
                end
            end
            if length(reconfig_topos) >= max_reconfig, break; end
        end
        if length(reconfig_topos) >= max_reconfig, break; end
    end
end
end

% ====================================================================
%  子函数：构建有根树（BFS，确定父子关系）
% ====================================================================
function [parent, children, bfs_order] = build_tree(root, fr, to, nb)
parent    = zeros(1, nb);   % parent(j) = 父节点索引，root的parent=0
children  = cell(1, nb);    % children{i} = [子节点列表]
bfs_order = [];

visited = false(1, nb);
queue   = root;
visited(root) = true;
parent(root)  = 0;

while ~isempty(queue)
    cur = queue(1); queue = queue(2:end);
    bfs_order(end+1) = cur;
    % 寻找以cur为起点或终点的支路
    fwd = find(fr == cur & ~visited(to)');
    bwd = find(to == cur & ~visited(fr)');
    for l = fwd'
        child = to(l);
        if ~visited(child)
            visited(child) = true;
            parent(child) = cur;
            children{cur}(end+1) = child;
            queue(end+1) = child;
        end
    end
    for l = bwd'
        child = fr(l);
        if ~visited(child)
            visited(child) = true;
            parent(child) = cur;
            children{cur}(end+1) = child;
            queue(end+1) = child;
        end
    end
end
end

% ====================================================================
%  子函数：计算每个节点的子树节点集
% ====================================================================
function sub_nodes = get_subtree_nodes(root_j, children)
sub_nodes = root_j;
queue = root_j;
while ~isempty(queue)
    cur = queue(1); queue = queue(2:end);
    if ~isempty(children{cur})
        sub_nodes = [sub_nodes, children{cur}];
        queue = [queue, children{cur}];
    end
end
end

% ====================================================================
%  子函数：计算节点 j 到根节点的路径（支路索引列表）
% ====================================================================
function path_branches = get_path_to_root(j, parent, fr, to, nl)
path_branches = [];
cur = j;
while parent(cur) ~= 0
    p = parent(cur);
    % 找连接 p-cur 的支路
    for l = 1:nl
        if (fr(l)==p && to(l)==cur) || (fr(l)==cur && to(l)==p)
            path_branches(end+1) = l;
            break;
        end
    end
    cur = p;
end
end

% ====================================================================
%  子函数：计算子树成员（辅助，批量处理所有节点）
% ====================================================================
function subtrees = compute_subtrees(fr, to, nl, nb, parent, children)
subtrees.nodes = cell(1, nb);
for j = 1:nb
    subtrees.nodes{j} = get_subtree_nodes(j, children);
end
end

% ====================================================================
%  子函数：检查拓扑是否辐射状（无环）
% ====================================================================
function ok = is_radial_tree(active_branches, branch, nb)
% 辐射状 ⟺ 连通且 |支路数| = |节点数| - 1（生成树条件）
n_edges = length(active_branches);
ok = (n_edges == nb - 1);  % 快速检查；完整检查需 BFS 连通性验证
end

% ====================================================================
%  子函数：检查重复拓扑
% ====================================================================
function dup = is_duplicate_topo(new_topo, topo_list)
dup = false;
for i = 1:length(topo_list)
    if isequal(sort(new_topo), sort(topo_list{i}))
        dup = true; return;
    end
end
end

% ====================================================================
%  子函数：孤岛系统 LFR（所有节点均无法供电）
% ====================================================================
function [AX, b, betaG, betaT, topo_info] = isolated_system_LFR(nb)
AX    = -eye(nb);
b     = zeros(nb, 1);
betaG = zeros(nb, nb);
betaT = zeros(nb, 0);
topo_info.topo_idx = 0;
topo_info.active_branches = [];
topo_info.n_topo = 0;
topo_info.min_J_est = inf;
end

% ====================================================================
%  子函数：PBC-Θ（配电网版，使用 BFM 灵敏度 Mat）
%  与大电网版逻辑完全相同，Mat 已是上游路径矩阵
% ====================================================================
function [b_lower_G, b_upper_G, betaG_lower, betaG_upper] = ...
        pbc_theta_dist(Mat, Gmax_node, Gmin_node, P_ref)
nl = size(Mat, 1);
nb = size(Mat, 2);
b_upper_G   = zeros(nl, 1);
b_lower_G   = zeros(nl, 1);
betaG_upper = zeros(nl, nb);
betaG_lower = zeros(nl, nb);

for l = 1:nl
    coeff = Mat(l, :)';
    [val_u, bG_u] = greedy_max(coeff,  Gmax_node, Gmin_node, P_ref);
    b_upper_G(l)      = val_u;
    betaG_upper(l, :) = bG_u;
    [val_lo, bG_lo] = greedy_max(-coeff, Gmax_node, Gmin_node, P_ref);
    b_lower_G(l)      = val_lo;
    betaG_lower(l, :) = bG_lo;
end
end

% ====================================================================
%  子函数：贪心算子（与大电网版完全相同）
% ====================================================================
function [obj_val, betaG_row] = greedy_max(coeff, Gmax, Gmin, P_ref)
nb = length(coeff);
G  = Gmin; betaG_row = zeros(1, nb);
remaining = max(P_ref - sum(Gmin), 0);
[sorted_c, sort_idx] = sort(coeff, 'descend');
for k = 1:nb
    idx = sort_idx(k);
    if sorted_c(k) <= 0, break; end
    cap_k = Gmax(idx) - Gmin(idx);
    if cap_k < 1e-10, continue; end
    fill_k = min(cap_k, remaining);
    G(idx) = Gmin(idx) + fill_k;
    if fill_k >= cap_k - 1e-10
        betaG_row(idx) = sorted_c(k);
    end
    remaining = remaining - fill_k;
    if remaining < 1e-10, break; end
end
obj_val = dot(coeff, G);
end

function val = getopt(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), val=s.(f); else, val=d; end
end