$title "Power System GDX data file maker"

$if not set filepath $setnames "%gams.i%" filepath filename fileextension
$if not set casepath $setnames "%case%" casepath casename caseextension
$if not set solution $set solution "temp_solution.gdx"
$if not set out $set out %filepath%%casename%_solution.gdx
$if not set ac $set ac 1

$if not set timeperiod abort "Time period not set!"

* Define sets
sets
    ints,
    t,
    k    contingency list / 0 /,
    bus,
    gen,
    circuit,
    demandbid,
    demandbidmap,
    interface,
    interfacemap,
    fuel_t, fuel_s, prime_mover,
    bus_t, bus_s,
    gen_t, gen_s,
    branch_t / Lineflow /, branch_s,
    demandbid_t, demandbid_s
    interface_t,
    line(bus,bus,circuit),
    transformer(bus,bus,circuit),
    monitored_lines(bus,bus,circuit)
;

* replace the k predefined above to the contingency list that we have used
$GDXIN %solution%
$LOADR k
$GDXIN

* set in parameeres bus_s is changed to k to contain contingency solutions
parameters version, baseMVA, total_cost;
parameters businfo(bus,bus_t,k), geninfo(gen,gen_t,k), fuelinfo(fuel_t,fuel_s),
           branchinfo(bus,bus,circuit,branch_t,k),
           demandbidinfo(demandbid,t,demandbid_t,demandbid_s),
           interfaceinfo(interface,t,interface_t)
           Frequency_error(k) frequency error for each contingency case;

$GDXIN %case%
$LOAD version, baseMVA, total_cost
$LOAD ints, t, bus, gen, circuit, line, transformer, monitored_lines,
$LOAD bus_t, bus_s, gen_t, gen_s, branch_s,
* to merge the set "currentflow" in branch_t
$LOADM branch_t,
$LOAD demandbid_t, demandbid_s, interface_t
$LOAD fuel_t, fuel_s, prime_mover
*$LOAD businfo, geninfo, branchinfo
$LOAD demandbid, demandbidmap, demandbidinfo
$LOAD interface, interfacemap, interfaceinfo
$GDXIN

alias(bus,i,j);
alias(circuit,c);

* Unload the contingency list and solution
parameters V_mag, V_ang, P_g, Q_g, Lineflow, total_cost, Fre_err,
           LMP_Energy, LMP_Loss, LMP_Congestion, LMP,
           LineSP, status, shuntB;
$GDXIN %solution%
$LOADR V_mag, V_ang, P_g, Lineflow, Fre_err

$ifthen %ac% == 1
$LOADR Q_g, shuntB
$endif

$ifthen %decompose_lmp% == 1
$LOADR LMP_Energy, LMP_Loss, LMP_Congestion
$endif

$LOADR total_cost, LMP, LineSP
$GDXIN
version = Jnow;

* here each variable has the index k to contain base case and contingency case solutions
businfo(bus,'Vm',k) = V_mag(bus,k);
businfo(bus,'Va',k) = V_ang(bus,k);
businfo(bus,'LMP',k) = LMP(bus,k) / baseMVA;
geninfo(gen,'Pg',k) = P_g(gen) * baseMVA;

$ifthen %ac% == 1
geninfo(gen,'Qg',k) = Q_g(gen,k) * baseMVA;
businfo(bus,'switchedBsSolved',k) = shuntB(bus) * baseMVA;
$endif

$ontext
parameter qd, bus_pf;
$GDXIN %solution%
$load Qd, bus_pf
$GDXIN
businfo(bus,'Qd','%timeperiod%') = Qd(bus);
businfo(bus,'pf','given') = bus_pf(bus);
$offtext

$ifthen %decompose_lmp% == 1
businfo(bus,'LMP_Energy',k) = LMP_Energy(bus) / baseMVA;
businfo(bus,'LMP_Loss','k) = LMP_Loss(bus) / baseMVA;
businfo(bus,'LMP_Congestion',k) = LMP_Congestion(bus) / baseMVA;
$endif

branchinfo(i,j,c,'LineSP',k) = LineSP(i,j,c,k) * baseMVA;
branchinfo(i,j,c,'Lineflow',k) = Lineflow(i,j,c,k) * baseMVA;
Frequency_error(k) = Fre_err(k);


execute_unload "%out%temp", version, total_cost, baseMVA, ints, k, t, bus, gen, circuit,
                        line, transformer, monitored_lines, demandbid, demandbidmap, interface, interfacemap
                        bus_t, gen_t,branch_t, fuel_t, fuel_s, prime_mover
                        demandbid_t, demandbid_s, interface_t
                        businfo, geninfo, branchinfo, frequency_error, demandbidinfo, interfaceinfo, fuelinfo
;

execute 'gams "%filepath%save_domain_info.gms" --in="%out%temp" gdx="%out%"'
if(errorlevel ne 0, abort "Saving domain info failed!");
execute 'rm "%out%temp"';

