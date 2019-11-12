function figplot=getfigplot(medge)
inipos = edge_search(medge);
while (inipos==0)  % 若边缘不闭合，则一直调用边界平滑程序
    medge = edge_smooth2(medge);
    inipos = edge_search(medge);  
end
figplot=inipos;
end