function [NodeIsVisited,NewVisitedNode,BranchIsVisited] = OneStepDepthSearch(StartNode,bus,...
    branch,NodeIsVisited,BranchIsVisited,AllreadyProcessBranch)

[PQ, PV, REF, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, VMAX, VMIN, PGMAX,PGMIN,QGMAX,QGMIN,PG,QG,LOLP,LOLE,LOLF,LOLD,EDNS,EENS,...
    LAM_P,LAM_Q,MU_VMAX,MU_VMIN,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE, TAP, SHIFT_ANGLE, BR_STATUS, BR_LAMDA,...
        BR_MTTF,TCSC_X, BR_PFAILURE, PF, QF, PT, QT, MU_SF, MU_ST, KT_MAXANGLE_TPST,KT_MINANGLE_TPST,KT_MAX_X_TCSC,KT_MIN_X_TCSC] = idx_brch;



m_branch = find(branch(:,F_BUS) == StartNode | branch(:,T_BUS) == StartNode);
AlreadyVisitedBranch = find(BranchIsVisited(m_branch) == 1 | AllreadyProcessBranch(m_branch) == 1);
m_branch(AlreadyVisitedBranch) = [];
NewVisitedNode = [];
NewVisitedBranch = m_branch;
if isempty(NewVisitedBranch) == 0
    BranchIsVisited(NewVisitedBranch) = 1;
    for i=1:size(NewVisitedBranch,1)
        if branch(NewVisitedBranch(i),F_BUS) == StartNode
            if NodeIsVisited(branch(NewVisitedBranch(i),T_BUS)) == 0     
                NewVisitedNode(branch(NewVisitedBranch(i),T_BUS)) = 1;
                NodeIsVisited(branch(NewVisitedBranch(i),T_BUS)) = 1;
            end
        else
            if NodeIsVisited(branch(NewVisitedBranch(i),F_BUS)) == 0     
                NewVisitedNode(branch(NewVisitedBranch(i),F_BUS)) = 1;
                NodeIsVisited(branch(NewVisitedBranch(i),F_BUS)) = 1;
            end
        end
    end    
end
NewVisitedNode = find(NewVisitedNode > 0);
return;
