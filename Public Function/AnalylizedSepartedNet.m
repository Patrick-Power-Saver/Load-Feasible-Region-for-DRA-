function [nodehao,branchhao] = AnalylizedSepartedNet(bus,branch,AllreadyProcessBus,AllreadyProcessBranch)
%%Judge the transmission network is whether separted into several isolated subnetwork;
%%if so,rerurn one of these isolated subnetwork through its nodo number and branch number 
%% in the bus and branch matrix

% [PQ, PV, REF, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
% 	VA, BASE_KV, VMAX, VMIN, PGMAX,PGMIN,QGMAX,QGMIN,PG,QG,LOLP,LOLE,LOLF,LOLD,EDNS,EENS,...
%     LAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus;
[PQ, PV, REF, NOTE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, VMAX, VMIN, PGMAX,PGMIN,QGMAX,QGMIN,PG,QG,LOLP,LOLE,LOLF,LOLD,EDNS,EENS,...
    LAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, BR_LAMDA,...
        BR_MTTF,TCSC_X, BR_PFAILURE, PF, QF, PT, QT, MU_SF, MU_ST, KT_MAXANGLE_TPST,KT_MINANGLE_TPST,KT_MAX_X_TCSC,KT_MIN_X_TCSC] = idx_brch;

    
NodeIsVisited = zeros(size(bus,1),1);
NewVisitedNode = zeros(size(bus,1),1);

busindex = find(AllreadyProcessBus == 0);

NodeIsVisited(bus(busindex(1),BUS_I)) = 1 ;
NewVisitedNode(bus(busindex(1),BUS_I)) = 1;
if isempty(find(AllreadyProcessBranch == 0)) == 0  %如果有正常zhi lu
    BranchIsVisted = zeros(size(branch,1),1);
else
    BranchIsVisted = [];
    nodehao = find(NewVisitedNode > 0);
    branchhao = [];
    return;
    
end
Node = find(NewVisitedNode > 0);     %找到此分区内的第一个节点
while isempty(Node) == 0
     NewVisitedNode =  zeros(size(bus,1),1);
    for i=1:size(Node,1)
        StartNode = Node(i);
        [NodeIsVisited,NewNode,BranchIsVisted]  = OneStepDepthSearch(StartNode,bus,branch,NodeIsVisited,BranchIsVisted,AllreadyProcessBranch);
        NewVisitedNode(NewNode) = 1;
    end
    Node = find(NewVisitedNode > 0);
end
nodehao = find(NodeIsVisited > 0);%找出此分块区域中的节点编号
branchhao = find(BranchIsVisted > 0);%找出此分块区域中的支路编号
return;
