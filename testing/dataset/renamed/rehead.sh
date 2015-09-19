for F in *txt; do
    sed -i "s/formation_energy/U/" $F
    sed -i "s/mu_Ni/mu0/" $F
    sed -i "s/mu_Al/mu1/" $F
    sed -i "s/Ni/N0/" $F
    sed -i "s/Al/N1/" $F
    sed -i "s/generalized_enthalpy/phi/" $F
    sed -i "s/temperature/T/" $F
    sed -i "s/beta/b/" $F
done
