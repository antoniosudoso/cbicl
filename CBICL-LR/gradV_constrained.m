function g = gradV_constrained(yU, yV, lV, mV, wV, beta, TVAtTUt, k, matU, matV, TVTVt, eV, e_m, CL_V_op, CL_V_E)

    YV = matV(yV);
    YVYVt = YV*YV';
    lV_hat = lV + beta*((YVYVt*eV)-e_m);
    mV_hat = mV + beta*((TVTVt(:)'*YVYVt(:))-k);
    wV_hat = wV + beta*CL_V_op(YVYVt);
    g = -TVAtTUt*matU(yU) + (eV*lV_hat' + lV_hat*eV' + 2*mV_hat*TVTVt + CL_V_E(wV_hat) + CL_V_E(wV_hat)')*YV;

end