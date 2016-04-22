$ontext
Calculate line power and objective value from P,Q,V_real.l,V_imag.l
$offtext

$ontext
V_LineIr.l(i,j,c,k)$branchstatus(i,j,c) =
            1/sqr(ratio(i,j,c))
                * (g(i,j,c)*V_real.l(i,k) - (b(i,j,c) + bc(i,j,c)/2)*V_imag.l(i,k))
            - 1/ratio(i,j,c)
                * (  (g(i,j,c)*V_real.l(j,k) - b(i,j,c)*V_imag.l(j,k))*cos(angle(i,j,c))
                   + (g(i,j,c)*V_imag.l(j,k) + b(i,j,c)*V_real.l(j,k))*sin(angle(i,j,c))
                  )
;

V_LineIr.l(j,i,c,k)$branchstatus(i,j,c) =
            (g(i,j,c)*V_real.l(j,k) - (b(i,j,c) + bc(i,j,c)/2)*V_imag.l(j,k))
            - 1/ratio(i,j,c)
                * (  (g(i,j,c)*V_real.l(i,k) - b(i,j,c)*V_imag.l(i,k))*cos(angle(i,j,c))
                   - (g(i,j,c)*V_imag.l(i,k) + b(i,j,c)*V_real.l(i,k))*sin(angle(i,j,c))
                  )
;

V_LineIq.l(i,j,c,k)$branchstatus(i,j,c) =
            1/sqr(ratio(i,j,c))
                * (g(i,j,c)*V_imag.l(i,k) + (b(i,j,c) + bc(i,j,c)/2)*V_real.l(i,k))
            - 1/ratio(i,j,c)
                * (  (g(i,j,c)*V_imag.l(j,k) + b(i,j,c)*V_real.l(j,k))*cos(angle(i,j,c))
                   - (g(i,j,c)*V_real.l(j,k) - b(i,j,c)*V_imag.l(j,k))*sin(angle(i,j,c))
                  )
;

V_LineIq.l(j,i,c,k)$branchstatus(i,j,c) =
            (g(i,j,c)*V_imag.l(j,k) + (b(i,j,c) + bc(i,j,c)/2)*V_real.l(j,k))
            - 1/ratio(i,j,c)
                * (  (g(i,j,c)*V_imag.l(i,k) + b(i,j,c)*V_real.l(i,k))*cos(angle(i,j,c))
                   + (g(i,j,c)*V_imag.l(i,k) + b(i,j,c)*V_real.l(i,k))*sin(angle(i,j,c))
                  )
;


V_objcost.l = sum(gen,
      costcoef(gen,'3')
      + costcoef(gen,'2')*U_P.l(gen)*baseMVA
      + costcoef(gen,'1')*sqr(U_P.l(gen)*baseMVA))
;

$offtext
