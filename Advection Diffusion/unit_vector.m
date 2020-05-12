%creates unit vector given dimension of unit vector and direction (all
%other directions are 0)
%arguments:     dim: length of unit vector; 
%               i: specific dimension/entry that is 1;
%               type: 0 for row vector, 1 for column vector;

function e=unit_vector(dim,i,type)

if type==0
    e=zeros(1,dim);
elseif type==1
    e=zeros(dim,1);
else
    error('type must be 0 for row vector or 1 for column vector');
    return
end

e(i)=1;