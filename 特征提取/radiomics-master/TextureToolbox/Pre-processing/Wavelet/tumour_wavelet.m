function feature = tumour_wavelet(tumour)
%load tumour;
WT = dwt3(volume,'coif1');
tumour_WT = WT.dec;
%% caculate energy of each wave
WT_1 = tumour_WT{1};
WT_2 = tumour_WT{2};
WT_3 = tumour_WT{3};
WT_4 = tumour_WT{4};
WT_5 = tumour_WT{5};
WT_6 = tumour_WT{6};
WT_7 = tumour_WT{7};
WT_8 = tumour_WT{8};
energy = zeros(8,1);
for i = 1:8
    eval(['energy(',num2str(i),') = sum(sum(sum(WT_',num2str(i),'.^2)));']);
end
disp(['the energy of the LLL is ',num2str(energy(1))]);
disp(['the energy of the HLL is ',num2str(energy(2))]);
disp(['the energy of the LHL is ',num2str(energy(3))]);
disp(['the energy of the HHL is ',num2str(energy(4))]);
disp(['the energy of the LLH is ',num2str(energy(5))]);
disp(['the energy of the HLH is ',num2str(energy(6))]);
disp(['the energy of the LHH is ',num2str(energy(7))]);
disp(['the energy of the HHH is ',num2str(energy(8))]);
feature.wt1_energy = energy(1);
feature.wt2_energy = energy(2);
feature.wt3_energy = energy(3);
feature.wt4_energy = energy(4);
feature.wt5_energy = energy(5);
feature.wt6_energy = energy(6);
feature.wt7_energy = energy(7);
feature.wt8_energy = energy(8);
end