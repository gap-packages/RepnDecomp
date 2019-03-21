# selects method based on what's in the options stack and fills out
# the corresponding attribute of rho
ComputeUsingMethod@ := function(rho)
    local method;
    method := ValueOption("decomp_method");

    # the default is serre
    if method = "alternate" then
        return REPN_ComputeUsingMyMethod(rho);
    else
        return REPN_ComputeUsingSerre(rho);
    fi;
end;
