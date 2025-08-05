function g = gradU_constrained(yU, yV, lU, mU, wU, beta, TUATVt, k, matU, matV, TUTUt, eU, e_n, CL_U_op, CL_U_E)

    YU = matU(yU);
    YUYUt = YU*YU';
    lU_hat = lU + beta*((YUYUt*eU)-e_n);
    mU_hat = mU + beta*((TUTUt(:)'*YUYUt(:))-k);
    wU_hat = wU + beta*CL_U_op(YUYUt);
    g = -TUATVt*matV(yV);
    g = g + (eU*lU_hat' + lU_hat*eU' + 2*mU_hat*TUTUt + CL_U_E(wU_hat) + CL_U_E(wU_hat)')*YU;
end