function f = augL_constrained(yU, yV, lU, lV, mU, mV, wU, wV, beta, TUATVt_vec, k, matU, matV, TUTUt, TVTVt, eU, eV, e_n, e_m, CL_U_op, CL_V_op)

    YU = matU(yU);
    YV = matV(yV);

    YUYUt = YU*YU';
    YVYVt = YV*YV';
    YUYVt = YU*YV';

    YU_sum = (YUYUt*eU)-e_n;
    YV_sum = (YVYVt*eV)-e_m;

    YU_trace = (TUTUt(:)'*YUYUt(:))-k;
    YV_trace = (TVTVt(:)'*YVYVt(:))-k;

    YU_CL = CL_U_op(YUYUt);
    YV_CL = CL_V_op(YVYVt);

    f = -TUATVt_vec'*YUYVt(:) + ...
        lU'*YU_sum + lV'*YV_sum + ...
        (0.5*beta)*(norm(YU_sum)^2 + norm(YV_sum)^2) + ...
        mU*YU_trace + mV*YV_trace + ...
        (0.5*beta)*(YU_trace^2+YV_trace^2) + ...
        wU'*YU_CL + wV'*YV_CL + ...
        (0.5*beta)*(norm(YU_CL)^2 + norm(YV_CL)^2);

end