% periodic boundary condition
function [x,y,z] = pbc(x,y,z,xboundlo,yboundlo,zboundlo,xboundhi,yboundhi,zboundhi,xprd,yprd,zprd)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if (x < xboundlo)
    x = x + xprd;
end
if (x >= xboundhi)
    x = x - xprd;
end
if (y < yboundlo)
    y = y + yprd;
end
if (y >= yboundhi)
    y = y - yprd;
end
 if (z < zboundlo) 
     z = z + zprd ;
 end
 if (z >= zboundhi)
     z = z - zprd;
 end
end

