function [] = plot_magnetization(r,M)
    %-quick plotting of magnetization
    figure; 
    subplot(131); plot(r*1000,squeeze(real(M(1,:) + 1i * M(2,:)))); ylim([-1.1,1.1]);
    subplot(132); plot(r*1000,squeeze(imag(M(1,:) + 1i * M(2,:)))); ylim([-1.1,1.1])
    subplot(133); plot(r*1000,squeeze(M(3,:))); 
end