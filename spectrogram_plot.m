function spectrogram_plot(input, clow, chigh, time, f)
  imagesc(time,f/1000,mag2db(abs(input)), [clow, chigh]); 
  colorbar; axis xy; set(gcf,'color','w');
  set(gca,'Fontsize',14); 
  xlabel('Time (s)'); 
  ylabel('Frequency (Hz)');
end