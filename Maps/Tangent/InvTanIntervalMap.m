function x = InvTanIntervalMap(y,alpha1,alpha2,y1,y2)

    a2  = InvLinearMap(alpha2,y1,y2);
    L = (y2-y1)/2;
    a1 = alpha1/L;
    
    x2       = InvLinearMap(y,y1,y2);
    x        = AtanMap(x2,a1,a2);
    

end