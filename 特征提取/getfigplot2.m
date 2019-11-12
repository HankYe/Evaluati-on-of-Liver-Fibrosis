function figplot=getfigplot2(BW3)
     B1 = bwboundaries(BW3);
     figplot=B1{1,1}';
end
