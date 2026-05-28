function [TIE_LINES_RAW] = convert_data_to_m(excel_filename, output_m_filename, case_name)
% CONVERT_DATA_TO_M 将指定的 Excel 配网数据转换为 RTS79_1 格式的 .m 算例文件
%
% 调用示例:
%   convert_data_to_m('85-Node System Data.xlsx', 'dist_opf_case85.m', 'dist_opf_case85')
%
% 输入参数:
%   excel_filename    - 系统基础数据文件 (如: '85-Node System Data.xlsx')
%   output_m_filename - 输出的 .m 文件名 (如: 'dist_opf_case85.m')
%   case_name         - 生成的函数名/算例名 (如: 'dist_opf_case85')

    %% 0. 解析节点规模与初始化
    % 从 case_name 中提取数字部分作为规模 (如从 'dist_opf_case85' 提取 '85')
    sys_scale = regexp(case_name, '\d+', 'match', 'once');
    if isempty(sys_scale)
        error('无法从 case_name 中提取节点规模数字，请确保其包含节点数（如 dist_opf_case85）');
    end
    sys_sheet = [sys_scale, '-node'];
    
    % 固定 Testbench 文件名 (请确保此文件与本脚本在同一目录下)
    tb_filename = 'Testbench for Linear Model Based Reliability Assessment Method for Distribution Optimization Models Considering Network Reconfiguration.xlsx';

    fprintf('\n======================================================\n');
    fprintf('开始生成算例: %s\n', case_name);
    fprintf('======================================================\n');

    %% §1. 读取 Testbench 数据 (提取联络线和容量)
    fprintf('>> [1/6] 读取 Testbench 联络线数据: Sheet="%s"\n', sys_sheet);
    tb_cell = readcell(tb_filename, 'Sheet', sys_sheet);
    nrows_tb = size(tb_cell, 1); 
    hdr_row = 0;
    
    for ri = 1:nrows_tb
        if ischar(tb_cell{ri,1}) && contains(tb_cell{ri,1}, 'Tie-Switch')
            hdr_row = ri; 
            break;
        end
    end
    if hdr_row == 0, error('未找到Tie-Switch表头行'); end
    
    LINE_CAP = str2double(extractBefore(string(tb_cell{hdr_row+1,4}), ' '));
    TIE_LINES_RAW = [];
    for ri = hdr_row+2:nrows_tb
        v1 = tb_cell{ri,1}; 
        v2 = tb_cell{ri,2};
        if isnumeric(v1) && ~isnan(v1) && isnumeric(v2) && ~isnan(v2)
            TIE_LINES_RAW(end+1,:) = [v1, v2]; %#ok<AGROW>
        end
    end
    fprintf('   -> 线路容量 = %.0f MW，联络线 = %d 条\n', LINE_CAP, size(TIE_LINES_RAW,1));

    %% §2. 读取基础可靠性参数
    fprintf('>> [2/6] 读取系统可靠性参数: %s\n', excel_filename);
    t_branch = readtable(excel_filename, 'Sheet', 'Branch Lengths (km)');
    t_branch = t_branch(:, 1:3); 
    t_branch.Properties.VariableNames = {'From','To','Length_km'};
    t_branch = t_branch(~isnan(t_branch.From), :);

    t_dur = readtable(excel_filename, 'Sheet', 'Interruption durations (h)', 'HeaderLines', 3);
    t_dur = t_dur(:, 1:4); 
    t_dur.Properties.VariableNames = {'From','To','RP','SW'};
    t_dur = t_dur(~isnan(t_dur.From), :);

    t_peak = readtable(excel_filename, 'Sheet', 'Peak Nodal Demands (kW)');
    t_peak = t_peak(:, 1:2); 
    t_peak.Properties.VariableNames = {'Node','P_kW'};
    t_peak = t_peak(~isnan(t_peak.Node), :);

    t_other = readtable(excel_filename, 'Sheet', 'Other data', 'ReadVariableNames', false);
    col1 = cellfun(@(x) string(x), t_other{:,1}, 'UniformOutput', true);
    lambda_per_km = str2double(string(t_other{find(contains(col1, 'Failure rate'), 1), 2}));
    
    row_lf = find(contains(col1, 'Loading factors'), 1);
    L_f = [str2double(string(t_other{row_lf,3})), str2double(string(t_other{row_lf+1,3})), str2double(string(t_other{row_lf+2,3}))] / 100;
    
    fprintf('   -> lambda = %.4f /km, L_f = [%s]\n', lambda_per_km, num2str(L_f, '%.2f '));

    %% §3. 构建 Bus 矩阵与负荷重要度
    fprintf('>> [3/6] 映射组装 Bus 数据矩阵...\n');
    num_buses = str2double(sys_scale);
    bus = zeros(num_buses, 12);
    busLoadImportant = zeros(num_buses, 4);
    
    for i = 1:num_buses
        bus(i, 1) = i;
        bus(i, 2) = (i == 1) * 3 + (i ~= 1) * 1; % 1号节点为平衡节点(3)，其余为PQ(1)
        
        idx = find(t_peak.Node == i);
        if ~isempty(idx)
            bus(i, 3) = t_peak.P_kW(idx(1)); % kW (输出文件中会转换为 MW)
            bus(i, 4) = bus(i, 3) * 0.33;    % 约 0.95 功率因数估算无功
        end
        
        bus(i, 7) = 1;      % Area
        bus(i, 8) = 1.0;    % Vm
        bus(i, 10) = 12.66; % baseKV
        bus(i, 11) = 1.1;   % Vmax
        bus(i, 12) = 0.9;   % Vmin
        
        % 负荷重要度矩阵
        busLoadImportant(i, 1) = i;
        busLoadImportant(i, 2:4) = L_f; 
    end

    %% §4. 构建 Branch 矩阵 (主网 + 联络线)
    fprintf('>> [4/6] 映射组装 Branch 数据矩阵...\n');
    num_main_branches = size(t_branch, 1);
    main_branches = zeros(num_main_branches, 12);
    for i = 1:num_main_branches
        f = t_branch.From(i); 
        t = t_branch.To(i);
        main_branches(i, 1:2) = [f, t];
        main_branches(i, 3:4) = [0.05, 0.02]; % 默认阻抗占位
        main_branches(i, 6) = LINE_CAP;
        main_branches(i, 9) = 1;              % 默认闭合
        main_branches(i, 10) = t_branch.Length_km(i) * lambda_per_km; % BR_LAMDA
        
        idx = find(t_dur.From == f & t_dur.To == t);
        if ~isempty(idx)
            main_branches(i, 11) = t_dur.RP(idx(1));
        else
            main_branches(i, 11) = 5.0; % 缺省修复时间
        end
    end

    num_tie = size(TIE_LINES_RAW, 1);
    tie_branches = zeros(num_tie, 12);
    for i = 1:num_tie
        tie_branches(i, 1:2) = TIE_LINES_RAW(i, :);
        tie_branches(i, 3:4) = [0.05, 0.02];
        tie_branches(i, 6) = LINE_CAP;
        tie_branches(i, 9) = 0;   % 联络线默认断开
        tie_branches(i, 10) = 0.5 * lambda_per_km; % 假设0.5km
        tie_branches(i, 11) = 5.0; 
    end
    branch = [main_branches; tie_branches];

    %% §5. 构建 Gen 矩阵
    fprintf('>> [5/6] 映射组装 Gen 数据矩阵...\n');
    gen = zeros(1, 21);
    gen(1, 1) = 1;        % 挂靠节点
    gen(1, 2) = 100;      % Pmax
    gen(1, 4) = 100;      % Qmax
    gen(1, 5) = -100;     % Qmin
    gen(1, 6) = 1;        % Status
    gen(1, 7) = 0.45;     % Lamda
    gen(1, 8) = 15.0;     % MTTF

    %% §6. 写入目标 .m 文件
    fprintf('>> [6/6] 写入目标文件: %s\n', output_m_filename);
    fid = fopen(output_m_filename, 'w');
    if fid == -1
        error('无法创建或打开文件: %s', output_m_filename);
    end

    % 写入头部
    fprintf(fid, 'function [baseMVA,bus,busPmax,busPmin,busQmax,busQmin,busPG,busQG,gen,branch, PBase ,PFBase,busLoadImportant,w_gen] = %s\n', case_name);
    fprintf(fid, '%% %s 配网测试系统数据 (RTS79 标准格式)\n', upper(case_name));
    fprintf(fid, '%% 自动提取自: %s 和 Testbench\n\n', excel_filename);

    fprintf(fid, '[PQ, PV, REF, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...\n');
    fprintf(fid, '\tVA, BASE_KV, VMAX, VMIN, PGMAX,PGMIN,QGMAX,QGMIN,PG,QG,LOLP,LOLE,LOLF,LOLD,EDNS,EENS,...\n');
    fprintf(fid, '\tLAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus;\n');
    fprintf(fid, '[GEN_BUS, GEN_PMAX, GEN_PMIN, GEN_QMAX, GEN_QMIN, GEN_STATUS, ...\n');
    fprintf(fid, '\tGEN_LAMDA,GEN_MTTF,GEN_PG,GEN_QG,GEN_PFAILURE] = idx_gen;\n');
    fprintf(fid, '[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, BR_LAMDA,...\n');
    fprintf(fid, '\tBR_MTTF,TCSC_X, BR_PFAILURE, PF, QF, PT, QT, MU_SF, MU_ST, KT_MAXANGLE_TPST,KT_MINANGLE_TPST,KT_MAX_X_TCSC,KT_MIN_X_TCSC] = idx_brch;\n');
    fprintf(fid, '[NEWTON,FDFP_XB,FDFP_BX,DCPF] = idx_powerflow;\n\n');

    fprintf(fid, '%%%% system MVA base\nbaseMVA = 100.0000;\n\n');

    % 写入 Bus
    fprintf(fid, '%%%% bus data\nbus = [\n');
    for i = 1:size(bus, 1)
        fprintf(fid, '\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%d\t%.4f\t%d\t%.2f\t%.2f\t%.2f;\n', bus(i,:));
    end
    fprintf(fid, '];\n\n');

    % 写入 Bus Load Importance
    fprintf(fid, '%%%% Bus Load Importance Data\nbusLoadImportant = [\n');
    for i = 1:size(busLoadImportant, 1)
        fprintf(fid, '\t%d\t%.4f\t%.4f\t%.4f;\n', busLoadImportant(i,:));
    end
    fprintf(fid, '];\n\n');

    % 写入 Gen
    fprintf(fid, '%%%% generator data\ngen = [\n');
    for i = 1:size(gen, 1)
        fprintf(fid, '\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t%.4f\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n', gen(i,1:8));
    end
    fprintf(fid, '];\n\n');

    % 写入 Branch
    fprintf(fid, '%%%% branch data (包含主线路及联络线)\nbranch = [\n');
    for i = 1:size(branch, 1)
        fprintf(fid, '\t%d\t%d\t%.6f\t%.6f\t%.4f\t%.1f\t%.4f\t%.4f\t%d\t%.5f\t%.4f\t%.4f;\n', branch(i,1:12));
    end
    fprintf(fid, '];\n\n');
    fprintf(fid, 'branch(:, 13:23) = 0;\n\n');

    % 写入底部转换和赋值逻辑
    bottom_logic = [
        '%%-----  OPF Data  -----%%' newline ...
        'gencost = [ 2 1500 0 3 0.11 5 150; ];' newline ...
        'w_gen = 1.05 * ones(size(gen, 1), 1);' newline ...
        newline ...
        '%% 数据转换与标幺化计算' newline ...
        'rate = 1.00;' newline ...
        'bus(:,3) = (rate .* bus(:,3)) / 1000; % kW -> MW' newline ...
        'bus(:,4) = (rate .* bus(:,4)) / 1000; % kVAr -> MVAR' newline ...
        newline ...
        'Vbase = bus(1, 10); ' newline ...
        'Zbase = (Vbase^2) / baseMVA;' newline ...
        'branch(:, 3) = branch(:, 3) / Zbase; % R p.u.' newline ...
        'branch(:, 4) = branch(:, 4) / Zbase; % X p.u.' newline ...
        newline ...
        'gen(:,2:5) = rate.*gen(:,2:5);' newline ...
        'gen(:,9:10) = rate.*gen(:,9:10);' newline ...
        newline ...
        'bus(:,13:32) = 0;' newline ...
        'busPmax = zeros(size(bus,1),1); busPmin = zeros(size(bus,1),1); busPG = zeros(size(bus,1),1);' newline ...
        'busQmax = zeros(size(bus,1),1); busQmin = zeros(size(bus,1),1); busQG = zeros(size(bus,1),1);' newline ...
        newline ...
        'gen(:,11) = 1-8760./(gen(:,7).*gen(:,8)+8760);' newline ...
        newline ...
        'for i=1:size(bus,1)' newline ...
        '    genindex = find(gen(:,1) == i);' newline ...
        '    if(~isempty(genindex))' newline ...
        '        busPmax(i) = sum(gen(genindex,2)); busPmin(i) = sum(gen(genindex,3)); busPG(i)   = sum(gen(genindex,9));' newline ...
        '        busQmax(i) = sum(gen(genindex,4)); busQmin(i) = sum(gen(genindex,5)); busQG(i)   = sum(gen(genindex,10));' newline ...
        '    end' newline ...
        'end' newline ...
        newline ...
        'bus(:,13) = busPmax; bus(:,14) = busPmin; bus(:,15) = busQmax;' newline ...
        'bus(:,16) = busQmin; bus(:,17) = busPG; bus(:,18) = busQG;' newline ...
        newline ...
        'branch(:,13) = 1-8760./(branch(:,10).*branch(:,11)+8760);' newline ...
        'branch(:,14:23) = 0.0;' newline ...
        newline ...
        'PBase = prod(1-gen(:,11)) * prod(1-branch(:,13));' newline ...
        'PFBase = sum(gen(:,7)) + sum(branch(:,10));' newline ...
        newline ...
        'end' newline
    ];
    fprintf(fid, '%s', bottom_logic);
    fclose(fid);

    fprintf('>> 成功！目标文件已生成: %s\n', output_m_filename);
    fprintf('>> 已提取 %d 条联络线数据。\n', size(TIE_LINES_RAW, 1));
end