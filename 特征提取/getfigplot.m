function figplot=getfigplot(medge)
inipos = edge_search(medge);
while (inipos==0)  % ����Ե���պϣ���һֱ���ñ߽�ƽ������
    medge = edge_smooth2(medge);
    inipos = edge_search(medge);  
end
figplot=inipos;
end