function [TB_er,grad_p_grad_T, TB_rel,grad_T_grad_T,CB_er, CB_rel] = rel_TB_CB(t_ans, z_ans,sx_ans, sy_ans, s2_ans, Z, P, T, dx, dy, A_x,A_y, A)






im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);

%% TB relative error
Tppc = ppc_pchip(Z,T);
Tz = ppc_val_mex(Z, Tppc, z_ans, 1);

p_ans = ppc_pchip(Z, P, z_ans);
Pppc = ppc_pchip(Z, P);
Pz = ppc_val_mex(Z, Pppc, z_ans, 1);

grad_p_x = (p_ans - im1(p_ans))./dx;
grad_p_y = (p_ans - jm1(p_ans))./dy;

grad_T_x = (t_ans - im1(t_ans))./dx;
grad_T_y = (t_ans - jm1(t_ans))./dy;

s2pzTz = Tz.* Pz .* s2_ans;
pzgradTs_x = Pz.*grad_T_x + Tz .* grad_p_x;
pzgradTs_y = Pz.*grad_T_y + Tz .* grad_p_y;

pzgradTs = div_dot(sx_ans, sy_ans,pzgradTs_x ,pzgradTs_y, A_x, A_y, A);

grad_p_grad_T = div_dot(grad_p_x, grad_p_y, grad_T_x, grad_T_y, A_x, A_y, A);

TB_er = s2pzTz + pzgradTs;

TB_rel = (TB_er)./grad_p_grad_T;


%% CB relative error

grad_T_x = (t_ans - im1(t_ans))./dx;
grad_T_y = (t_ans - jm1(t_ans))./dy;

s2TzTz = Tz.* Tz .* s2_ans;
TzgradTs_x = Tz .* grad_T_x;
TzgradTs_y = Tz .* grad_T_y;

TzgradTs = div_dot(sx_ans, sy_ans,TzgradTs_x ,TzgradTs_y, A_x, A_y, A);

grad_T_grad_T = div_dot(grad_T_x, grad_T_y, grad_T_x, grad_T_y, A_x, A_y, A);

CB_er = s2TzTz + 2*TzgradTs;

CB_rel = (CB_er)./grad_T_grad_T;

end

