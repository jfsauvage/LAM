; check pyramid Transfert Function

restore, 'D:\sauvage\Desktop\HARMONI\ATC-E2E_FOURRIER\FOURIER\SimuFourier\SH-WFS_SCAO_HARMONI_Jmedian_30deg_25m_0.10rad2_2.20um.dat'
dspallsh = dspall
restore, 'D:\sauvage\Desktop\HARMONI\ATC-E2E_FOURRIER\FOURIER\SimuFourier\PY-WFS_SCAO_HARMONI_Jmedian_30deg_25m_0.10rad2_2.20um.dat'
dspallpy = dspall

tek_color
plot_oo, (eclat(f))[256,256:*], (eclat(dspallsh))[256,256:*], xr=[0.001, 10], yr=[0.001, 100]
oplot, (eclat(f))[256,256:*], (eclat(dspallpy))[256,256:*], color=2

affa, [eclat(dspallsh), eclat(dspallpy)]<0.1
