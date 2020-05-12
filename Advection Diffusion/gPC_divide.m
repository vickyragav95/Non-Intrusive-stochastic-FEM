%performs gPC division a(y)=b(y)/c(y)

function a=gPC_divide(b,c,M,Ckij)

  if length(c)~=M || length(b)~=M
      error('c and b must be of size M');
  end
  P=zeros(M,M);
  for k=1:M
      for i=1:M
          for j=1:M
              P(k,j) = P(k,j) + Ckij(k,i,j)*c(i);
          end
      end
  end
  b = b';
  a = (P\b)';
end