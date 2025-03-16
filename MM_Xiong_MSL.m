function [theta, fail] = MM_Xiong_MSL(anc, omega_vec, Q, T)
%the proposed BMM approach

%anc: receiver (anchor) position matrix
%omega_vec: [omega11,...,omegaL1,omega12,...,omegaLP]'
%Q: diagonal covariance matrix
%T: time interval

fail = false;

[H, L] = size(anc);
P = length(omega_vec)/L;

theta = zeros(2*H+1,1);

sigma_vec = diag(Q);

p0_old = ones(H,1);

p0_new = zeros(H,1);
omega0_new = 0;
v_new = zeros(H,1);

cnt = 0;

while (norm(p0_new - p0_old)/norm(p0_old) > 1e-4)
    cnt = cnt + 1;
    p0_old = p0_new;
    omega0_old = omega0_new;
    v_old = v_new;

    p0_num = zeros(H,1);
    p0_den = 0;
    omega0_num = 0;
    omega0_den = 0;
    v_num = zeros(H,1);
    v_den = 0;
    for i = 1:L
        for j = 1:P
            lambda_ij = omega_vec((j-1)*L+i)*(p0_old + v_old*j*T - anc(:,i))/norm(p0_old + v_old*j*T - anc(:,i));
            rho_ij = omega_vec((j-1)*L+i) - norm(p0_old + v_old*j*T - anc(:,i)) + omega0_old;
            kappa_ij = rho_ij/(2*norm(p0_old + v_old*j*T - anc(:,i)));
            p0_num = p0_num + (1/sigma_vec((j-1)*L+i)^2)*(lambda_ij - (kappa_ij+1)*(v_old*j*T - anc(:,i)));
            p0_den = p0_den + (1/sigma_vec((j-1)*L+i)^2)*(kappa_ij+1);

            omega0_num = omega0_num + rho_ij/sigma_vec((j-1)*L+i)^2;
            omega0_den = omega0_den + 2/sigma_vec((j-1)*L+i)^2;

            muij = (omega_vec((j-1)*L+i) - omega0_old)*(p0_old + v_old*j*T - anc(:,i))/norm(p0_old + v_old*j*T - anc(:,i));
            v_num = v_num + (2*j*T)*(muij - p0_old + anc(:,i))/sigma_vec((j-1)*L+i)^2;
            v_den = v_den + (2*j^2*T^2)/sigma_vec((j-1)*L+i)^2;
        end
    end
    p0_new = p0_num/p0_den;
    omega0_new = omega0_num/omega0_den;
    v_new = v_num/v_den;

    if (cnt > 10000)
        fail = true;
        return
    end
end

theta = [p0_new;v_new;omega0_new];

if (isnan(sum(theta))) || (isinf(sum(theta)))
    fail = true;
end

end

