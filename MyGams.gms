$title "AC optimal power flow model, rectangular current-voltage formulation"
*_______________________________________________________________________________
* Filename: iv_acopf.gms
* Description: AC optimal power flow model, rectangular current-voltage formulation
*
* Usage: gams iv_acopf.gms --case=/path/case.gdx
*
* Options:
* --timeperiod: Select time period to solve. Default=1
* --obj: Objective function, piecewise linear or quadratic. Default="quad"
* --linelimits: Type of line limit data to use. Default="given"
* --genPmin: Data for Generator lower limit. Default="given"
* --allon: Option to turn on all gens or lines during solve. Default=none
* --slim: slim option does not apply here. Default=0 (not used)
* --qlim: Option to use D-curve constraints. Default=0 (not used)
* --savesol: Turn on save solution option(1). Default=0
* --verbose: Supresses printout(0). Default=1
* --wind: Whether to turn off wind turbines. Can only be used with
*         PrimeMover,pm_WT
*_______________________________________________________________________________

* System dependence
$if %system.filesys% == UNIX $set sep '/'
$if not %system.filesys% == UNIX $set sep '\'

*===== SECTION: OPTIONS & ENVIRONMENT VARIABLES
* Printout options
$LOG OPTIONS & ENVIRONMENT VARIABLES
$ifthen %verbose% == 0
* Turn off solution printing
option solprint=on
* Turn off print options
$offlisting
option limrow=9999, limcol=80
$endif
Option Reslim=600;

* Define filepath, name and extension.
$setnames "%gams.i%" filepath filename fileextension
* Define type of model
$set modeltype "AC"
* Define input case
* $set case "rts96_spring_wday_t4.gdx"
$set ic 0
$set linelimits "inf"
$if not set case $abort "Model aborted. Please provide input case"
$setnames "%case%" casepath casename caseextension

* Default: timeperiod = 1
$if not set timeperiod $set timeperiod "1"
* Default: Quadratic objective function
$if not set obj $set obj "quad"
* Default: Use provided line limits (as opposed to uwcalc)
$if not set linelimits $set linelimits "given"
* Default: Use provided generator lower limit
$if not set genPmin $set genPmin "given"
* Default: allon=0
$if not set allon $set allon 0
* Default: slim options does not apply here. Apparent limits always used.
$set slim 0
* Default: Ignore D-curve constraints
$if not set qlim $set qlim 0
* Default: Treat wind turbines the same way
$if not set wind $set wind 0
* Default: Save solution option turned off
$if not set savesol $set savesol 1
* Default: elastic demand bidding does not apply here
$set demandbids 0
$set condensed 'no'

*===== SECTION: EXTRACT DATA
$batinclude "%filepath%extract_data.gms" modeltype case timeperiod demandbids linelimits genPmin allon

* Calculate Ybus matrix
$batinclude '%filepath%calc_Ybus.gms'

*===== SECTION: DATA MANIPULATION
*--- Define load, gen buses and active lines
sets
    load(bus)      "Load buses"
    isGen(bus)     "Generator buses"
    activeGen(bus) "Active generator buses"
    isLine(i,j)    "Active (i,j) line"
    contingency    "to determine the level of contingency" / 0 /
* / 0 / means no contingency.
* / 0,1 / means OPF with base case + line 1 out.
* / 0,1,2 / means OPF with base case + line 1 out + line 2 out
* / 0,2 / means OPF with base case + line 2 out
;

load(bus)$(sum(gen, atBus(gen,bus)) eq 0) = 1;
isGen(bus)$(not(load(bus))) = 1;
activeGen(bus)$(isGen(bus) and (sum(gen$atBus(gen,bus), status(gen)) > 0) ) = 1;
option isLine < branchstatus;

parameters
            DR(gen)      "droop coefficient"
            SUMDR ;

DR(gen)$(status(gen)) = (1/0.04)*(1/100)*Pmax(gen);
SUMDR = sum(gen, DR(gen));
DR(gen)$(status(gen)) = DR(gen)/SUMDR;

alias(contingency,k);
contingenstatus(i,j,c,'0')$line(i,j,c) = 1 ;

*=== to account for line limits individually
*currentrate('217','216',c) = 2;
*currentrate('114','116',c) = 2;
*currentrate('110','106',c) = 2;
*currentrate('117','116',c) = 2.5;
*currentrate('103','124',c) = 1.5;

option contingenstatus:3:3:1
display DR,activeGen, Pd, currentrate, r ;

*===== SECTION: ADJUST SYSTEM LOAD
Pd(i) = 1*Pd(i);
currentrate(i,j,c) = 1.5*currentrate(i,j,c);

*===== SECTION: VARIABLE DEFINITION
free variables
    U_P(gen)             "Real power generation of generator",
    V_Q(gen,k)           "Reactive power generation of generator",
    V_real(i,k)          "Real part of bus voltage",
    V_imag(i,k)          "Imaginary part of bus voltage",
    V_LineIr(i,j,c,k)    "Real part of current flowing from bus i towards bus j on line c",
    V_LineIq(i,j,c,k)    "Imaginary part of current flowing from bus i towards bus j on line c"

    V_fre(k)             "frequency error of the network"
;

positive variables
    V_shunt(bus,bus_s)    "Bus shunt susceptance"

    V_pw_cost(gen)    "Generator piecewise cost"
    V_Pd_elastic(demandbid)     "Elastic incremental demand"
    V_demandbid_rev(demandbid)  "Revenue from elastic incremental demand"
;

free variable V_objcost  "Total cost of objective function";

*===== SECTION: EQUATION DEFINITION
equations
    c_I_limit(i,j,c,k)     "Limit apparent current on a line"
    c_V_limit_lo(i,k)      "Limit voltage magnitude on a line"
    c_V_limit_up(i,k)      "Limit voltage magnitude on a line"

    c_Voltagefix(i,k)      "Fix Generator voltage magnitude for every contingency"

    c_ActivepowerUpperlimit(gen,k) " Active power generation upper limit with a slack variable"
    c_ActivepowerLowerlimit(gen,k) " Active power generation lower limit with a slack variable"

    c_Line1(i,j,c,k)    " Line equation 1 associated with the transmission matrix "
    c_Line2(i,j,c,k)    " Line equation 2 associated with the transmission matrix "
    c_Line3(i,j,c,k)    " Line equation 3 associated with the transmission matrix "
    c_Line4(i,j,c,k)    " Line euqation 4 associated with the transmission matrix "

    c_BalanceP(bus,k)      "Balance of real power for bus"
    c_BalanceQ(bus,k)      "Balance of reactive power for bus"

    c_pw_cost(gen,costptset) "Generator piecewise cost functions"
    c_obj                    "Objective function"
;

*===== SECTION: EQUATIONS PART 1
* Limit apparent current on a line
c_I_limit(i,j,c,k)$(branchstatus(i,j,c) or branchstatus(j,i,c))..
    sqr(V_LineIr(i,j,c,k)) + sqr(V_LineIq(i,j,c,k)) =l= sqr(currentrate(i,j,c));

*Limit voltage magnitude on a line
c_V_limit_lo(i,k)..
    sqr(V_real(i,k)) + sqr(V_imag(i,k)) =g= sqr(minVm(i))
;

*Limit voltage magnitude on a line
c_V_limit_up(i,k)..
    sqr(V_real(i,k)) + sqr(V_imag(i,k)) =l= sqr(maxVm(i))
;

*Fix generator voltage for every contingency
c_Voltagefix(i,k)$(activeGen(i) and (ord(k) < card(k)))..
    sqrt(sqr(V_real(i,k)) + sqr(V_imag(i,k))) =e= sqrt(sqr(V_real(i,k+1)) + sqr(V_imag(i,k+1)))
;

* Generator active power generation limits
c_ActivepowerUpperlimit(gen,k)$status(gen)..
U_P(gen) - (DR(gen))*V_fre(k) =l= Pmax(gen);

c_ActivepowerLowerlimit(gen,k)$status(gen)..
U_P(gen) - (DR(gen))*V_fre(k) =g= Pmin(gen);

*Line equation 1 associated with the transmission matrix
c_Line1(i,j,c,k)$branchstatus(i,j,c)..
         0 =e= ( V_real(i,k) - V_real(j,k)*( (ratio(i,j,c))*cos(angle(i,j,c)) + ((bc(i,j,c)/2)*(ratio(i,j,c))*cos(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                                                              - ((bc(i,j,c)/2)*(ratio(i,j,c))*sin(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_imag(j,k)*( (ratio(i,j,c))*sin(angle(i,j,c)) + ((bc(i,j,c)/2)*(ratio(i,j,c))*cos(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                                               + ((bc(i,j,c)/2)*(ratio(i,j,c))*sin(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_LineIr(j,i,c,k)*( ((ratio(i,j,c))*cos(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                  + ((ratio(i,j,c))*sin(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_LineIq(j,i,c,k)*( ((ratio(i,j,c))*cos(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                  - ((ratio(i,j,c))*sin(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) ) )$(contingenstatus(i,j,c,k) = 1)

+ V_LineIr(i,j,c,k)$(contingenstatus(i,j,c,k) = 0) ;

*Line equation 2 associated with the transmission matrix
c_Line2(i,j,c,k)$branchstatus(i,j,c)..
         0 =e= ( V_imag(i,k) - V_real(j,k)*( (ratio(i,j,c))*sin(angle(i,j,c)) + ((bc(i,j,c)/2)*(ratio(i,j,c))*cos(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                                                              + ((bc(i,j,c)/2)*(ratio(i,j,c))*sin(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_imag(j,k)*( -(ratio(i,j,c))*cos(angle(i,j,c)) - ((bc(i,j,c)/2)*(ratio(i,j,c))*cos(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                                                + ((bc(i,j,c)/2)*(ratio(i,j,c))*sin(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_LineIr(j,i,c,k)*( -((ratio(i,j,c))*cos(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                   + ((ratio(i,j,c))*sin(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_LineIq(j,i,c,k)*( ((ratio(i,j,c))*cos(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                  + ((ratio(i,j,c))*sin(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) ) )$(contingenstatus(i,j,c,k) = 1)

+ V_LineIr(j,i,c,k)$(contingenstatus(i,j,c,k) = 0) ;


*Line equation 3 associated with the transmission matrix
c_Line3(i,j,c,k)$branchstatus(i,j,c)..
         0 =e= ( V_LineIr(i,j,c,k) + V_real(j,k)*( (bc(i,j,c)/ratio(i,j,c))*sin(angle(i,j,c)) + ((sqr(bc(i,j,c))/4)*cos(angle(i,j,c))*g(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c))
                                                                                              + ((sqr(bc(i,j,c))/4)*sin(angle(i,j,c))*b(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c)) )

              + V_imag(j,k)*( (bc(i,j,c)/ratio(i,j,c))*cos(angle(i,j,c)) + ((sqr(bc(i,j,c))/4)*cos(angle(i,j,c))*b(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c))
                                                                         - ((sqr(bc(i,j,c))/4)*sin(angle(i,j,c))*g(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c)) )

              + V_LineIr(j,i,c,k)*( (1/ratio(i,j,c))*cos(angle(i,j,c))
                                  + ((bc(i,j,c)/2)*(1/ratio(i,j,c))*cos(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                  - ((bc(i,j,c)/2)*(1/ratio(i,j,c))*sin(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_LineIq(j,i,c,k)*( -(1/ratio(i,j,c))*sin(angle(i,j,c))
                                   - ((bc(i,j,c)/2)*(1/ratio(i,j,c))*cos(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                   - ((bc(i,j,c)/2)*(1/ratio(i,j,c))*sin(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) ) )$(contingenstatus(i,j,c,k) = 1)

+ V_LineIq(i,j,c,k)$(contingenstatus(i,j,c,k) = 0) ;


*Line equation 4 associated with the transmission matrix
c_Line4(i,j,c,k)$branchstatus(i,j,c)..
         0 =e= ( V_LineIq(i,j,c,k) - V_real(j,k)*( (bc(i,j,c)/ratio(i,j,c))*cos(angle(i,j,c)) + ((sqr(bc(i,j,c))/4)*cos(angle(i,j,c))*b(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c))
                                                                                              - ((sqr(bc(i,j,c))/4)*sin(angle(i,j,c))*g(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c)) )

              + V_imag(j,k)*( (bc(i,j,c)/ratio(i,j,c))*sin(angle(i,j,c)) + ((sqr(bc(i,j,c))/4)*cos(angle(i,j,c))*g(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c))
                                                                         + ((sqr(bc(i,j,c))/4)*sin(angle(i,j,c))*b(i,j,c))/((sqr(g(i,j,c)) + sqr(b(i,j,c)))*ratio(i,j,c)) )

              + V_LineIr(j,i,c,k)*( (1/ratio(i,j,c))*sin(angle(i,j,c))
                                  + ((bc(i,j,c)/2)*(1/ratio(i,j,c))*cos(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                  + ((bc(i,j,c)/2)*(1/ratio(i,j,c))*sin(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) )

              + V_LineIq(j,i,c,k)*( (1/ratio(i,j,c))*cos(angle(i,j,c))
                                  + ((bc(i,j,c)/2)*(1/ratio(i,j,c))*cos(angle(i,j,c))*b(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c)))
                                  - ((bc(i,j,c)/2)*(1/ratio(i,j,c))*sin(angle(i,j,c))*g(i,j,c))/(sqr(g(i,j,c)) + sqr(b(i,j,c))) ) )$(contingenstatus(i,j,c,k) = 1)

+ V_LineIq(j,i,c,k)$(contingenstatus(i,j,c,k) = 0) ;

*Balance of real power for bus
c_BalanceP(i,k)$(type(i) ne 4)..
        sum(gen$(atBus(gen,i) and status(gen)), U_P(gen) - (DR(gen))*V_fre(k))
        - Pd(i)
            =e=
          V_real(i,k) *
        ( sum((j,c)$(branchstatus(i,j,c)), V_LineIr(i,j,c,k))
        + sum((j,c)$(branchstatus(j,i,c)), V_LineIr(i,j,c,k)) )
        + V_imag(i,k) *
         (sum((j,c)$(branchstatus(i,j,c)), V_LineIq(i,j,c,k))
        + sum((j,c)$(branchstatus(j,i,c)), V_LineIq(i,j,c,k)))
        + Gs(i) * (sqr(V_real(i,k)) + sqr(V_imag(i,k)))
;

* Balance of reactive power for bus
c_BalanceQ(i,k)$(type(i) ne 4)..
        sum(gen$(atBus(gen,i) and status(gen)), V_Q(gen,k))
        - Qd(i)
            =e=
        - V_real(i,k) *
        ( sum((j,c)$(branchstatus(i,j,c)), V_LineIq(i,j,c,k))
        + sum((j,c)$(branchstatus(j,i,c)), V_LineIq(i,j,c,k)))
        + V_imag(i,k) *
        ( sum((j,c)$(branchstatus(i,j,c)), V_LineIr(i,j,c,k))
        + sum((j,c)$(branchstatus(j,i,c)), V_LineIr(i,j,c,k)))
        - Bs(i) * (sqr(V_real(i,k)) + sqr(V_imag(i,k)))
        - (sqr(V_real(i,k)) + sqr(V_imag(i,k))) * sum(bus_s$(not sameas(bus_s,'given')), Bswitched(i,bus_s) * V_shunt(i,bus_s))
;

* Objective functions and pwl costs are listed in a separate file
$batinclude "%filepath%cost_objective.gms" obj demandbids

* D-curve limits
$if %Qlim% == 1 $batinclude '%filepath%reactive_limits.gms' case

*===== SECTION: MODEL DEFINITION
model feas / c_V_limit_lo, c_V_limit_up, c_Voltagefix,
            c_Line1, c_Line2, c_Line3, c_Line4,
            c_BalanceP, c_BalanceQ, c_I_limit,
$if %qlim% == 1 ,c_Armature, c_Field, c_Heating
            /;

model acopf /feas, c_pw_cost, c_obj/;

*===== SECTION: VARIABLE BOUNDS
* Generator active power generation limits
U_P.lo(gen)$status(gen) = Pmin(gen);
U_P.up(gen)$status(gen) = Pmax(gen);
U_P.fx(gen)$(not status(gen)) = 0;
$ifthen %wind%==1
* Needed to avoid compilation error. Puts strings into UEL
set winddesc /'PrimeMover', 'pm_WT'/;
* Wind turbines are not reliable sources of power, treated differently
parameter windTurbine(gen);
windTurbine(gen)$(geninfo(gen, 'PrimeMover', 'pm_WT') eq 1) = 1;
V_P.fx(gen)$(windTurbine(gen)) = 0;
$endif

* Generator reactive power generation limits
* Does not impose Qmax, Qmin limits when the D-curve contraint is applied
$ifthen %qlim% == 0
V_Q.lo(gen,k)$status(gen) = Qmin(gen);
V_Q.up(gen,k)$status(gen) = Qmax(gen);
$endif
V_Q.fx(gen,k)$(not status(gen)) = 0;

* Slack variable(frequency error) limits
V_fre.lo(k) = -0.02;
V_fre.up(k) = 0.02;
V_fre.fx('0') = 0;

* Bus voltage magnitude limits
*V_real.lo(bus) = -MaxVm(bus); V_real.up(bus) = MaxVm(bus);
*V_imag.lo(bus) = -MaxVm(bus); V_imag.up(bus) = MaxVm(bus);
*V_imag.fx(bus,k)$(type(bus) eq 3) = 0;

* Fix swing bus angle
V_imag.fx(bus,k)$(type(bus) eq 3) = 0;

* Bus shunt susceptance
V_shunt.up(bus,bus_s) = numBswitched(bus,bus_s);
$if %switchedshunts% == 0 V_shunt.fx(bus,bus_s) = shunt.up(bus,bus_s);

* Elastic demand not considered
V_Pd_elastic.fx(demandbid) = 0;
V_demandbid_rev.fx(demandbid) = 0;

*===== SECTION: VARIABLE INITIAL STARTING POINTS
V_shunt.l(bus,bus_s)  = 1;

* Set initial conditions
$ifthen %ic% == 1 $batinclude '%filepath%ic_iv%sep%random_all.gms' condensed verbose
$elseif %ic% == 2 $batinclude '%filepath%ic_iv%sep%flat.gms' condensed verbose
$elseif %ic% == 3 $batinclude '%filepath%ic_iv%sep%random_v.gms' condensed verbose
$elseif %ic% == 4 $batinclude '%filepath%ic_iv%sep%dcopf_pv.gms' limits condensed verbose allon obj Plim timeperiod
$elseif %ic% == 5 $batinclude '%filepath%ic_iv%sep%dcopf_v.gms' limits condensed verbose allon obj Plim timeperiod
$elseif %ic% == 6 $batinclude '%filepath%ic_iv%sep%decoupled.gms' condensed verbose
$elseif %ic% == 7 $batinclude '%filepath%ic_iv%sep%dcopf_pv_loss ' condensed verbose
$elseif %ic% == 8 $batinclude '%filepath%ic_iv%sep%matpower.gms' condensed verbose
$elseif %ic% == 9 $batinclude '%filepath%ic_iv%sep%given.gms' condensed verbose
$else $batinclude '%filepath%ic_iv%sep%default2.gms' condensed verbose
$endif

*===== SECTION: MODEL OPTIONS AND SOLVE
option nlp = knitro;
solve acopf min V_objcost using nlp ;


*==== SECTION: Solution Analysis
* See if model is solved
parameter
    infeas "Number of infeasibilities from model solve";

infeas = acopf.numInfes;
display infeas;


* Declaration needs to be made outside loop
set
    lines_at_limit(i,j,c,k) "lines at their bound"
    maximizer(i,j,c,k)
;
parameters
    total_cost "Cost of objective function"
    LMP(bus,k) "Locational marginal price"
    LineSP(i,j,c,k) "Marginal price of active power on line (i,j,c)"
    shuntB(i)

    V_mag(i,k)
    V_ang(i,k)

    Linecurrent(i,j,c,k)
    Maxlinecurrent
;

*$ontext
if(infeas eq 0,
* Final Objective function value
    total_cost = V_objcost.l;
* Generator real power solution
    Pg(gen) = U_P.l(gen);
* Generator reactive power solution
*    Qg(gen,k) = V_Q.l(gen,k);
* Voltage magnitude solution
    V_mag(bus,k) = sqrt(sqr(V_real.l(bus,k)) + sqr(V_imag.l(bus,k)));
* Voltage angle solution
    V_ang(bus,k)$(V_real.l(bus,k) > 0) = arctan(V_imag.l(bus,k)/V_real.l(bus,k)) * 180/pi;
    V_ang(bus,k)$(V_real.l(bus,k) le 0) = arctan(V_imag.l(bus,k)/V_real.l(bus,k)) * 180/pi + 180;
    V_ang(bus,k) = V_ang(bus,k) * (3.14/180)  ;
* Bus shunt solution
    shuntB(i) = sum(bus_s, V_shunt.l(i,bus_s)*Bswitched(i,bus_s));
* Locational marginal price of bus at time t
    LMP(bus,k) = c_BalanceP.m(bus,k);
* Marginal for active power on a line
*    LineSP(i,j,c,k)$branchstatus(i,j,c) = c_I_Limit.m(i,j,c,k);
*    LineSP(j,i,c,k)$branchstatus(i,j,c) = c_I_Limit.m(j,i,c,k);

* Find which lines are at their limits
lines_at_limit(i,j,c,k)$(branchstatus(i,j,c) or branchstatus(j,i,c)) = yes$
     (sqr(currentrate(i,j,c)) - sqr(V_LineIr.L(i,j,c,k)) - sqr(V_LineIq.L(i,j,c,k)) le 1e-4);

*{GAMS does not have an argmax function, which is
*requested here. Here is how to do the
*equivalent with a dynamic set that I’ll call maximizer.}

Linecurrent(i,j,c,k)$((branchstatus(i,j,c) or branchstatus(j,i,c)))
         = sqrt(sqr(V_LineIr.L(i,j,c,k)) + sqr(V_LineIq.L(i,j,c,k)))  ;

Maxlinecurrent
         = smax((i,j,c,k), Linecurrent(i,j,c,k));
maximizer(i,j,c,k) = yes$( Linecurrent(i,j,c,k) eq Maxlinecurrent)

option lines_at_limit:3:3:1, Linecurrent:3:3:1, maximizer:3:3:1
display V_mag, V_ang, lines_at_limit, Linecurrent, Maxlinecurrent, maximizer;

*==== SECTION: Solution Save
$ifthen %savesol% == 1
execute_unload '%filepath%temp_solution.gdx', Pg, Qg, Vm, Va, shuntB, total_cost, LMP, LineSP;
* execute 'gams %filepath%save_solution.gms gdxcompress=1 --ac=1 --case=%case% --solution=%filepath%temp_solution.gdx --timeperiod=%timeperiod% --out=%filepath%%casename%_AC_base_solution.gdx'
* if(errorlevel ne 0, abort "Saving solution failed!");
execute 'rm %filepath%temp_solution.gdx'
$endif
);
*$offtext
