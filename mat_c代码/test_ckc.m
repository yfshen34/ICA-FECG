
% z=yanchibaihua3(Xf',2 ,0);
% C = (z*z')/size(z,2);
% y = z'*(C^-1)*z;
% gama = zeros(1, size(y,2));
% for i = 1:size(y,2)
%     gama(i) = y(i,i);
% end
nt=11860;
P=z(:,nt)'*(C^-1)*z;
P = -P;
S = searchspike_in(P,25,100);
L = find(S==1);
card = size(L,2);
C1 = mean(z(:,L),2);
P1=C1'*(C^-1)*z;
S = searchspike_in(P,8,100);
L = find(S==1);
card = size(L,2);
C1 = mean(z(:,L),2);
P1=C1'*(C^-1)*z;