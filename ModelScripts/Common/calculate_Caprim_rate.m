function [Caprim_rate, Caunprim_rate] = calculate_Caprim_rate(Calcium, Ca_prim_type, prim_kM, prim_rate_const, unprim_rate_const, unprim_rate_const_0, unprim_mutant)

kD=prim_kM;

Caprim_rate = (Ca_prim_type~=0)*prim_rate_const;

if unprim_mutant == 1 %mutant with always low unpriming - Calcium insensitive
    Caunprim_rate = (Ca_prim_type~=0)*(unprim_rate_const*kD + unprim_rate_const_0);
elseif unprim_mutant == 2 %mutant with always high unpriming rate - Calcium insenstive. 
    Caunprim_rate = (Ca_prim_type~=0)*(unprim_rate_const + unprim_rate_const_0);
else %wt - Calcium sensitive:
    Caunprim_rate = (Ca_prim_type~=0)*(unprim_rate_const*((kD)*(1+Calcium.^(Ca_prim_type)) ./(Calcium.^(Ca_prim_type) + kD)) + unprim_rate_const_0);
end