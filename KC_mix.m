function [C_Ad, K_SG]=KC_mix(Metal, Pge)

switch (Metal)
    case 'A1' 
        C_Ad=0.896; %kJ/kg.K
       
        switch Pge
            case 0
                K_SG=0.198E-3; %kW/m.K
            case 5
                K_SG=0.218E-3; %kW/m.K
            case 10
                K_SG=0.314E-3; %kW/m.K
            case 15
                K_SG=0.363E-3; %kW/m.K
        end
    case 'Cu'
        C_Ad=0.385; %kJ/kg.K
        switch Pge
            case 0
                K_SG=0.198E-3; %kW/m.K
            case 5
                K_SG=0.187E-3; %kW/m.K
            case 10
                K_SG=0.246E-3; %kW/m.K
            case 15
                K_SG=0.324E-3; %kW/m.K
        end
    case 'Brass'
        C_Ad=0.380; %kJ/kg.K
        switch Pge
            case 0
                K_SG=0.198E-3; %kW/m.K
            case 5
                K_SG=0.176E-3; %kW/m.K
            case 10
                K_SG=0.211E-3; %kW/m.K
            case 15
                K_SG=0.327E-3; %kW/m.K
        end
    case 'Steel'
        C_Ad=0.460; %kJ/kg.K
        switch Pge
            case 0
                K_SG=0.198E-3; %kW/m.K
            case 5
                K_SG=0.141E-3; %kW/m.K
            case 10
                K_SG=0.169E-3; %kW/m.K
            case 15
                K_SG=0.254E-3; %kW/m.K
        end
    case 'None'
        C_Ad=0.0; %kJ/kg.K
        K_SG=0.198E-3; %kW/m.K
end
end