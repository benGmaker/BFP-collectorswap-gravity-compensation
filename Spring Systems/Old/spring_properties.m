function [k1, L1, L2] = spring_properties(L_0, fn, Fn, n)
    %old functions for computing spring properties
    %NOT RECOMMENDED TO USE
    global F1 F2 R_o h_mech h_adjust h_max
    Fn_eff = Fn*n               %effective max spring output force
    k1 = round(Fn/fn*1e3)*n     %stiffness of spring system
    L1 = F1/k1 %initial elongation
    L2 = F2/k1 %second elongation
    h_adjust = (F1 - F2)/k1 %adjustment height
    R_o = F1/k1 %outer diameter spiral pulley
    springstroke = (fn*1e-3 - L1) %[mm] possible stroke with spring
    Buildstroke = (h_max - L_0 - L1 - h_adjust - h_mech - R_o) %[m] possible stroke within the construction
end

