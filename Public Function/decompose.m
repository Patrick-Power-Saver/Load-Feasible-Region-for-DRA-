%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%判断是否解列函数,返回de块子块。
%形成的解列后的各个子块，由VS阵返回。（VS的行数为de）M为关联矩阵，行为节点，列为支路表示量
%注意：VS中1对应节点需进行tp映射后方为原节点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [de,VS]=decompose(failure_line,Node_Line,node_number,branch_number)

%failure_line       故障线路变量
%Line_Node          线路-节点矩阵
%node_number        节点数量
%branch_number      线路数量

%主函数
L=length(failure_line);%取故障线路向量维数
for k1=1:L
   Node_Line(:,failure_line(k1))=0;%把关联矩阵中对应的故障支路所在的列清零，即把此支路去掉
end
visited=zeros(1,node_number);%1*6 RBTS有6个节点
k=1;
temp=visited;%1*6
for k1=1:node_number
   if (visited(k1)==0)%这个是显然的   O.O
      visited=bfs(Node_Line,k1,visited,branch_number,node_number);%%广度优先搜索算法
      VS(k,:)=visited-temp;
      k=k+1;
      temp=visited;
   end
end
de=length(VS(:,1));%第一列指示是否解列


function visited=bfs(M,N,visited,NL,NB)%，M为关联矩阵，N为某个节点，NL为支路数，NB为节点总个数
%子函数
visited(N)=1;%将第一个元素置1
k=1;
for k1=1:NL%对所有支路
   if (M(N,k1)~=0)%若关联矩阵的第N行，第K1列非0，表示这条支路包含N节点，找到节点N所在支路，支路的1个节点
      for k2=1:NB
         if(M(k2,k1)~=0&&visited(k2)~=1)%找到支路的另一个节点
            visited(k2)=1;
            queue(k)=k2;%记住另一个节点所在的列，由此此条支路的两节点确定为N，k2
            k=k+1;
            break;
         end 
      end
   end
end
%上述程序段用于找到与节点N连通的节点 ，下述循环用于找到与节点queue(k1)连通的节点
for k1=1:k-1
   visited=bfs(M,queue(k1),visited,NL,NB);
end


