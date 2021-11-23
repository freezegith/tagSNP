function [Z H] = FastNewman(adjacent_matrix)
% FastNewman�㷨ʵ�����ŷ���
% ���㷨�μ����ס�Fast algorithm for detecting community structure in networks��(2003)
% ����
% adjacent_matrix - �ڽӾ���
% ���
% Z - n-1*3���󣬵�i�б�ʾ��i�κϲ�����1�к͵�2�б�ʾ�ϲ������ű�ţ���3���Ǻϲ����ģ���
% H - ������ͼ�ľ��

n = size(adjacent_matrix,1); % �ڵ���Ŀ size(*,1):�����������
max_id = n;
Z = [];
clusters = [1:n; zeros(1,n); 1:n]; % ��ʼ���֣���1���ǽڵ��ţ���2�������ű�ŵı任����3�������ű��
step = 1;
while numel(unique(clusters(3,:))) ~= 1 
% while step < n
    [Q e a clusters] = GetModularity(adjacent_matrix, clusters);
    k = size(e,1); % ������Ŀ
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
    i = DeltaQs(1,id); j = DeltaQs(2,id); % �ҳ�DeltaQ�������ŶԵ����
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

