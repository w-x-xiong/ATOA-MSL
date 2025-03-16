function [y, fail] = SDP_ma_MSL(anc, omega_vec, Q, T)
%the SOTA SDP method

%anc: receiver (anchor) position matrix
%omega_vec: [omega11,...,omegaL1,omega12,...,omegaLP]'
%Q: covariance matrix
%T: time interval

fail = false;

[H, L] = size(anc);
P = length(omega_vec)/L;
A = zeros(L*P,2*H+4);
b = zeros(L*P,1);
for j = 1:P
    for i = 1:L
        A((j-1)*L+i,:) = [2*anc(:,i)',2*j*T*anc(:,i)',-2*omega_vec((j-1)*L+i),-1,-j^2*T^2,(-2)*j*T];
        b((j-1)*L+i) = norm(anc(:,i))^2 - omega_vec((j-1)*L+i)^2;
    end
end
Sigma = Q;

cvx_begin quiet
variables Y(2*H+4,2*H+4) y(2*H+4)
expression Prod_mtx_l;
expression Prod_mtx_r;
expression obj;

Prod_mtx_l = [Y,y;y',1];
Prod_mtx_r = [A'*(pinv(Sigma))*A, -A'*(pinv(Sigma))*b; -b'*(pinv(Sigma))*A, b'*(pinv(Sigma))*b];
obj = trace(Prod_mtx_l*Prod_mtx_r);

minimize obj

subject to

y(2*H+2) == trace(Y(1:H,1:H)) - Y(2*H+1,2*H+1);

y(2*H+3) == trace(Y(H+1:2*H,H+1:2*H));

y(2*H+4) == trace(Y(1:H,H+1:2*H));

[Y, y; y', 1] == semidefinite(2*H+5);

cvx_end

B = 2*diag(omega_vec - y(2*H+1));
Sigma = B*Q*B';

cvx_begin quiet
variables Y(2*H+4,2*H+4) y(2*H+4)
expression Prod_mtx_l;
expression Prod_mtx_r;
expression obj;

Prod_mtx_l = [Y,y;y',1];
Prod_mtx_r = [A'*(pinv(Sigma))*A, -A'*(pinv(Sigma))*b; -b'*(pinv(Sigma))*A, b'*(pinv(Sigma))*b];
obj = trace(Prod_mtx_l*Prod_mtx_r);

minimize obj

subject to

y(2*H+2) == trace(Y(1:H,1:H)) - Y(2*H+1,2*H+1);

y(2*H+3) == trace(Y(H+1:2*H,H+1:2*H));

y(2*H+4) == trace(Y(1:H,H+1:2*H));

[Y, y; y', 1] == semidefinite(2*H+5);

cvx_end

if (isnan(cvx_optval)) || (cvx_optval == +Inf) || (cvx_optval == -Inf)
    fail = true;
end

end

