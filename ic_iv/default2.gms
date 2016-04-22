$ontext
Default initial conditions for polar ACOPF:
Just pick midpoint of variable ranges.
$offtext

$if not set filepath $setnames "%gams.i%" filepath filename fileextension

U_P.l(gen)$status(gen) = (Pmin(gen)+Pmax(gen))/2;
V_Q.l(gen,k)$status(gen) = (Qmin(gen)+Qmax(gen))/2;
V_real.l(bus,k)  = (minVm(bus)+maxVm(bus))/2;
* V_imag can stay 0, since angles are allowed to range in (-pi, pi)
V_imag.l(bus,k) = 0;

$if %condensed% == 'no' $include '%filepath%ic_iv%sep%calc_derived_vars.gms'

V_pw_cost.l(gen)$(status(gen) and (costmodel(gen) eq 1)) = max(0,
    smax(costptset$(ord(costptset) < numcostpts(gen)),
            ((costpts_y(gen,costptset+1) - costpts_y(gen,costptset))/
             (costpts_x(gen,costptset+1) - costpts_x(gen,costptset)))
              * (U_P.l(gen)*baseMVA - costpts_x(gen,costptset))
            + costpts_y(gen,costptset) - noloadcost(gen)))
;

V_objcost.l =
$ifthen.linear %obj% == "linear"
             sum(gen$(status(gen) and (costmodel(gen) eq 2)),
                             costcoef(gen,'0')
                           + costcoef(gen,'1')*V_P.l(gen)*baseMVA
                )
$else.linear
             sum(gen$(status(gen) and (costmodel(gen) eq 2)),
                             costcoef(gen,'0')
                           + costcoef(gen,'1')*U_P.l(gen)*baseMVA
                           + costcoef(gen,'2')*sqr(U_P.l(gen)*baseMVA)
                )
$endif.linear
           + sum(gen$(status(gen) and (costmodel(gen) eq 1)), V_pw_cost.l(gen)
                                                              + noloadcost(gen))
;

