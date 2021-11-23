function [Z H] = FastNewman(adjacent_matrix)
% FastNewman算法实现社团发现
% 该算法参见文献《Fast algorithm for detecting community structure in networks》(2003)
% 输入
% adjacent_matrix - 邻接矩阵
% 输出
% Z - n-1*3矩阵，第i行表示第i次合并，第1列和第2列表示合并的社团标号，第3列是合并后的模块度
% H - 聚类树图的句柄

n = size(adjacent_matrix,1); % 节点数目 size(*,1):输出矩阵行数
max_id = n;
Z = [];
clusters = [1:n; zeros(1,n); 1:n]; % 初始划分，第1行是节点标号，第2行是社团标号的变换，第3行是社团标号
step = 1;
while numel(unique(clusters(3,:))) ~= 1 
% while step < n
    [Q e a clusters] = GetModularity(adjacent_matrix, clusters);
    k = size(e,1); % 社团数目
    DeltaQs = [];
    for i = 1:size(e,1)
        for j = 1:size(e,1)
            if i ~= j
                DeltaQ = 2*(e(i,j)-a(i)*a(j));
                DeltaQs = [DeltaQs [i;j;DeltaQ]];
            end            
        end
    end
    [maxDeltaQ,id] = max(DeltaQs(3,:)); id = id(1);
    i = DeltaQs(1,id); j = DeltaQs(2,id); % 找出DeltaQ最大的社团对的序号
    max_id = max_id + 1;
    c_id1 = find(clusters(2,:) == i); 
    c_id2 = find(clusters(2,:) == j); 
    id1 = unique(clusters(3,c_id1)); 
    id2 = unique(clusters(3,c_id2));
    clusters(3,c_id1) = max_id;    
    clusters(3,c_id2) = max_id;  
    Z = [Z; [id1 id2 length([c_id1 c_id2])]];
    step = step + 1;
end

if nargout == 2
    H = dendrogram(Z,0);
    %new
    fid = fopen('***','wt');
    [mnew,nnew]=size(Z);
    for inew = 1:1:mnew
        for jnew = 1:1:nnew
            if jnew == nnew
                fprintf(fid,'%g\n',Z(inew,jnew));
            else
                fprintf(fid,'%g\t',Z(inew,jnew));
            end
        end
    end
    fclose(fid);
end
    
end

