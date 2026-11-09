function x = InvAtanMap(z,alpha1,alpha2)

    kappa  = atan(alpha1*(1+alpha2))/atan(alpha1*(1-alpha2));
    z0     = (kappa - 1)/(kappa + 1);
    lambda = atan(alpha1*(1-alpha2))/(1-z0); 

    x     = alpha2 + tan( (z -z0)*lambda )/alpha1;

end