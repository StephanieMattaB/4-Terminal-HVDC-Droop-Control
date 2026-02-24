% =========================================================================
% HVDC - four-terminal scheme
% 
% Developed by: F. Bianchi, E. Prieto, S. Matta
%% =========================================================================
close all, clear all
fpath = './DCtech';

% =========================================================================
%% Network values
C1 = 150e-6;
C2 = 150e-6;
C3 = 150e-6;
C4 = 150e-6;
R1 = 0.5;
R2 = 0.25;
R3 = 0.4;
L1 = 0.005;
L2 = 0.0025;
L3 = 0.004;

% State-space model
A = [0,     0,      0,     0,      -1/C1,   -1/C1,    0;
     0,     0,      0,     0,       0,       1/C2,   -1/C2;
     0,     0,      0,     0,       1/C3,    0,       0;
     0,     0,      0,     0,       0,       0,       1/C4;
     1/L1,  0,     -1/L1,  0,      -R1/L1,   0,       0;
     1/L2, -1/L2,   0,     0,       0,      -R2/L2,   0;
     0,     1/L3,   0,    -1/L3,    0,       0,      -R3/L3];
B = [1/C1,  0,      0,     0;
     0,     1/C2,   0,     0;
     0,     0,     -1/C3,  0;
     0,     0,      0,    -1/C4;
     0,     0,      0,     0;
     0,     0,      0,     0;
     0,     0,      0,     0];
C = [1, 0,  0,  0,  0,  0,  0;
     0, 1,  0,  0,  0,  0,  0;
     0, 0,  1,  0,  0,  0,  0;
     0, 0,  0,  1,  0,  0,  0];
D = zeros(4);

Bw = B(:,1:2); Bu = B(:,3:4);
Cz = C(1:2,:); Cy = C(3:4,:);

% =========================================================================
% Open loop transfers
G = ss(A,B,C,D);

Gzw = G(1:2,1:2); zpk(Gzw); 
Gzu = G(1:2,3:4); zpk(Gzu); 
Gyw = G(3:4,1:2); zpk(Gyw); 
Gyu = G(3:4,3:4); zpk(Gyu); 

% gains
Kgi = [100 50 20 7 3 1];
Kg  = 1./Kgi;
for ii=1:length(Kgi)
    lgtext{ii} = sprintf('1/%d',Kgi(ii));
end
lgtext{6} = '1';
P  = diag([1 1]);
% colors
ColorSet = diag(linspace(0,0.5,2))*ones(2,3);
LineSet  = {'-','--','-.'};

% =========================================================================
%% Stability
%
Hf = figure('Name','Closed loop poles');
Ha = axes('ColorOrder', ColorSet,...
    'LineWidth',1,...
    'Fontname','times',...
    'Fontsize',14,...
    'NextPlot','add');

%for ii=1:2:length(Kg)
ii = 2;
eiAcl(:,ii) = eig(A+Kg(ii)*Bu*P*Cy);
plot(real(eiAcl(:,ii)),imag(eiAcl(:,ii)),'x',...
    'MarkerEdgeColor',ColorSet(1,:),...
    'MarkerFaceColor',ColorSet(1,:),...
    'LineWidth',2,...
    'MarkerSize',10)
ii = 4;
eiAcl(:,ii) = eig(A+Kg(ii)*Bu*P*Cy);
plot(real(eiAcl(:,ii)),imag(eiAcl(:,ii)),'o',...
    'MarkerEdgeColor',ColorSet(1,:),...
    'MarkerFaceColor',0.7*[1 1 1],...
    'LineWidth',2,...
    'MarkerSize',10)
ii = 6;
eiAcl(:,ii) = eig(A+Kg(ii)*Bu*P*Cy);
plot(real(eiAcl(:,ii)),imag(eiAcl(:,ii)),'s',...
    'MarkerEdgeColor',ColorSet(1,:),...
    'MarkerFaceColor',0.7*[1 1 1],...
    'LineWidth',2,...
    'MarkerSize',10)
plot([0 0],3000*[-1 1],'k')
plot([-7000 0],[0 0],'k')
% end
% ylabel('Imag(\lambda(A_{cl}))'); xlabel('Real(\lambda(A_{cl}))');
ylabel('Imag(l(Acl))'); xlabel('Real(l(Acl))');
legend(lgtext{[2 4 6]},'Location','SouthWest'), legend boxoff
print('-depsc2',[fpath 'fig_eigAcl_c1.eps'])

% det(1+L(jw))
% w = logspace(-1,4,2000);
% LL1 = frd(eye(2)+Gyu,w); rLL1 = frdata(LL1);
% LL2 = frd(eye(2)+5*Gyu,w); rLL2 = frdata(LL2);
% for ii=1:length(w)
%     dl1(ii) = det(rLL1(:,:,ii));
%     dl2(ii) = det(rLL2(:,:,ii));
% end
% plot(real(dl1),imag(dl1))
% hold
% plot(real(dl2),imag(dl2),'r')


% =========================================================================
%% Closed loop transfers
systemnames = 'G';
inputvar = '[w{2};r{2};u{2}]';
outputvar = '[G(1:2);-r+G(3:4);-r+G(3:4)]';
input_to_G = '[w;u]';
G2 = sysic;

clear Gcl
for ii=1:length(Kg)
    Gcl(:,:,ii) = lft(G2,Kg(ii)*P);
end

% frequencies
w = 2*pi*logspace(-1,3,1000);

% -------------------------------------------------------------------------
%% Figure: sensivitity
S = Gcl(3:4,3:4,:);
%
Hf = figure('Name','S');
Ha = axes('ColorOrder', ColorSet,'LineStyleOrder', LineSet,...
    'LineWidth',1,...
    'Fontname','times',...
    'Fontsize',14,...
    'NextPlot','add',...
    'XScale','log',...
    'Xlim',[1e-1 1e3],'Ylim',[-40 20]);
for ii=1:length(Kg)
    sv = sigma(S(:,:,ii),w);
    svSmx(:,ii) = 20*log10(sv(1,:)');
    svSmn(:,ii) = 20*log10(sv(2,:)');
end
semilogx(w/2/pi,svSmx')
semilogx(w/2/pi,svSmn','--')
ylabel('SVD (dB)'); xlabel('Frequency (Hz)');
legend(lgtext), legend boxoff
% print('-depsc2','fig_svS_c1.eps')

% -------------------------------------------------------------------------
%% Figure: SGyw
SGyw = Gcl(3:4,1:2,:);
Hf = figure('Name','SGyw');
Ha = axes('ColorOrder', ColorSet,'LineStyleOrder', LineSet,...
    'Fontname','times',...
    'Fontsize',14,...
    'XScale','log',...
    'NextPlot','add',...
    'Xlim',[1e-1 1e3],'Ylim',[-20 50]);
for ii=1:length(Kg)
    sv = sigma(SGyw(:,:,ii),w);
    svSGyw(:,ii) = 20*log10(sv(1,:)');
end
% limits
emax = 20*log10(22.5); 
h = area([0.1 200 2000],[emax emax emax+20]); 
set(h,'FaceColor',[1 1 1]*0.9,'BaseValue',-20,'LineWidth',0.5)
%
Hl = semilogx(w/2/pi,svSGyw','LineWidth',1);
set(Hl(3),'LineWidth',2);
ylabel('SVD (dB)'); xlabel('Frequency (Hz)');
legend(Hl,lgtext,'Location','SouthWest'), legend boxoff
print('-depsc2',[fpath 'fig_svSGyw_c1.eps'])
% svd(freqresp(Gyu,0)\freqresp(Gyw,0))
% dcgain
[Nyu,d] = tfdata(Gyu); Nyu = [Nyu{1,1}(end) Nyu{1,2}(end); Nyu{2,1}(end) Nyu{2,2}(end)];
[Nyw,d] = tfdata(Gzu); Nyw = [Nyw{1,1}(end) Nyw{1,2}(end); Nyw{2,1}(end) Nyw{2,2}(end)];

% -------------------------------------------------------------------------
%% Figure: KSGyw
for ii=1:length(Kg)
    KSGyw(:,:,ii) = Kg(ii)*SGyw(:,:,ii);
end
Hf = figure('Name','KSGyw');
Ha = axes('ColorOrder', ColorSet,'LineStyleOrder', LineSet,...
    'Fontname','times',...
    'Fontsize',14,...
    'NextPlot','add',...
    'XScale','log',...
    'Xlim',[1e-1 1e3],'Ylim',[-40 20]);
for ii=1:length(Kg)
    sv = sigma(KSGyw(:,:,ii),w);
    svKSGyw(:,ii) = 20*log10(sv(1,:)');
end
% limits
emax = 20*log10(1.50); 
h = area([0.1 200 2000],[emax emax emax-20]); 
set(h,'FaceColor',[1 1 1]*0.9,'BaseValue',-40,'LineWidth',0.5)
%
Hl = semilogx(w/2/pi,svKSGyw','LineWidth',1);
set(Hl(3),'LineWidth',2);
ylabel('SVD (dB)'); xlabel('Frequency (Hz)');
legend(Hl,lgtext,'Location','SouthWest'), legend boxoff
print('-depsc2',[fpath 'fig_svKSGyw_c1.eps'])

% -------------------------------------------------------------------------
%% Figure: GzuKS
GzuKS = Gcl(1:2,3:4,:);
Hf = figure('Name','GzuKS');
Ha = axes('ColorOrder', ColorSet,'LineStyleOrder', LineSet,...
    'LineWidth',1,...
    'Fontname','times',...
    'Fontsize',14,...
        'NextPlot','add',...
        'XScale','log',...
        'Xlim',[1e-1 1e3],'Ylim',[-40 20]);
for ii=1:length(Kg)
    sv = sigma(GzuKS(:,:,ii),w);
    svGzuKS(:,ii) = 20*log10(sv(1,:)');
end
Hl = semilogx(w/2/pi,svGzuKS','LineWidth',1);
ylabel('SVD (dB)'); xlabel('Frequency (Hz)');
set(Hl(3),'LineWidth',2);
legend(Hl,lgtext,'Location','SouthWest'), legend boxoff
% print('-depsc2','fig_svGzuKS_c1.eps')
% dcgain
[Nyu,d] = tfdata(Gyu); Nyu = [Nyu{1,1}(end) Nyu{1,2}(end); Nyu{2,1}(end) Nyu{2,2}(end)];
[Nzu,d] = tfdata(Gzu); Nzu = [Nzu{1,1}(end) Nzu{1,2}(end); Nzu{2,1}(end) Nzu{2,2}(end)];
% svd(freqresp(Gzu,0)/freqresp(Gyu,0))

% -------------------------------------------------------------------------
%% Figure: Gzw_GzuKSGyw
Gzw_GzuKSGyw = Gcl(1:2,1:2,:);
Hf = figure('Name','Gzw_GzuKSGyw');
Ha = axes('ColorOrder', ColorSet,'LineStyleOrder', LineSet,...
    'LineWidth',1,...
    'Fontname','times',...
    'Fontsize',14,...
        'NextPlot','add',...
        'XScale','log',...
        'Xlim',[1e-1 1e3],'Ylim',[-20 50]);
for ii=1:length(Kg)
    sv = sigma(Gzw_GzuKSGyw(:,:,ii),w);
    svGzw_GzuKSGyw(:,ii) = 20*log10(sv(1,:)');
end
% limits
emax = 20*log10(22.5); 
h = area([0.1 200 2000],[emax emax emax+20]); 
set(h,'FaceColor',[1 1 1]*0.9,'BaseValue',-20,'LineWidth',0.5)
%
Hl = semilogx(w/2/pi,svGzw_GzuKSGyw','LineWidth',1);
set(Hl(3),'LineWidth',2);
ylabel('SVD (dB)'); xlabel('Frequency (Hz)');
legend(Hl,lgtext,'Location','SouthWest'), legend boxoff
print('-depsc2',[fpath 'fig_svGzw_GzuKSGyw_c1.eps'])

%%
Kg = 1/22.5;
Gcl_0 = lft(G2,Kg*eye(2));
Gzw_GzuKSGyw = Gcl_0(1:2,1:2);
dcgain(Gzw_GzuKSGyw)*[667;667]
totalbus=dcgain(Gzw_GzuKSGyw)*[667;667]+[145e3;145e3]
%%
Kg = 1/20;
Gcl_0 = lft(G2,Kg*eye(2));
Gzw_GzuKSGyw = Gcl_0(1:2,1:2);
dcgain(Gzw_GzuKSGyw)*[633.2;633.2]
totalbus=dcgain(Gzw_GzuKSGyw)*[633.2;633.2]+[145e3;145e3]




