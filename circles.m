close all 
clear


X = readmatrix('path.txt');
n_points = size(X,2);


dftm = dftmtx(n_points);
FX = dftm * X.' ./ n_points;
FX1S = imag(FX(:,1));
FX1C = real(FX(:,1));
FX2S = imag(FX(:,2));
FX2C = real(FX(:,2));


figure; hold on
xlim([-3 3])
ylim([-3 3])
h = gca;
h.YLimMode='manual';
set(h,'Color','k')
set(gcf, 'InvertHardcopy', 'off')
set(gcf, 'Position',  [10, 10, 1500, 1500])
filename = 'animation.gif';

N = ceil(n_points/2);
F = cell(1,8*N);
R = kron([FX1S(1:N); FX1C(1:N); FX2S(1:N); FX2C(1:N)].', [1 1]);
w = @(t,k) t*(k-1)*2*pi/n_points;
for k = 1:N
    sc = (1 + (k==1))^-1;
    l = (k-1)*8 + 1;
    % sin x-direction
    F{l}       = @(t) sc*FX1S(k)*[ sin(w(t,k)) ;  cos(w(t,k)) ];
    F{l+1}     = @(t) sc*FX1S(k)*[ sin(w(t,k)) ; -cos(w(t,k)) ];
    R(l:l+1)   = FX1S(k);
    % cos x-direction
    F{l+2}     = @(t) sc*FX1C(k)*[ cos(w(t,k)) ;  sin(w(t,k)) ];
    F{l+3}     = @(t) sc*FX1C(k)*[ cos(w(t,k)) ; -sin(w(t,k)) ];
    R(l+2:l+3) = FX1C(k);
    % sin y-direction
    F{l+4}     = @(t) sc*FX2S(k)*[  cos(w(t,k)) ; sin(w(t,k)) ];
    F{l+5}     = @(t) sc*FX2S(k)*[ -cos(w(t,k)) ; sin(w(t,k)) ];
    R(l+4:l+5) = FX2S(k);
    % cos y-direction
    F{l+6}     = @(t) sc*FX2C(k)*[  sin(w(t,k)) ; cos(w(t,k)) ];
    F{l+7}     = @(t) sc*FX2C(k)*[ -sin(w(t,k)) ; cos(w(t,k)) ];
    R(l+6:l+7) = FX2C(k);
end
[~,sorted] = sort(-abs(R(9:end)));
sorted = [1:8, sorted+8];
F = F(sorted);

substeps=10;
px = zeros(1,substeps*n_points);
py = zeros(1,substeps*n_points);
t = 1;
while true
    x0 = 0;
    y0 = 0;
    cla
    %title("t = " + sprintf("%.2f",t/n_points))
    plot(x0,y0,'wx');
    f = cellfun(@(c) c(t), F, 'UniformOutput', false);
    xy = horzcat(f{:});
    r = vecnorm(xy(:,2:end));
    xy = xy*triu(ones(numel(f)));
    plot(xy(1,:),xy(2,:),'Color',[0.5,0.5,0.7])
    c = xy(:,1:end-1);
    circs = viscircles( c.', r.',...
        'LineWidth',0.1,'EnhanceVisibility',false,'Color',[0.2,0.2,0.2]);
    x0 = xy(1,end);
    y0 = xy(2,end);
    idx = round((t-1)*substeps+1);
    px(idx) = x0;
    py(idx) = y0;
    plot(X(1,:),X(2,:),'Color',[0.1,0.3,0.1])
    plot(px(1:idx),py(1:idx),'Color',[0.7,0.7,0.7],'LineWidth',1)
    %plot(px(1:substeps:idx),py(1:substeps:idx),'rx')
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if t == 1
        imwrite(imind,cm,filename,'gif', 'DelayTime',0.1,'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif', 'DelayTime',0.1,'WriteMode','append');
    end
    t = t + substeps^-1;
    if t > n_points
        t = 1;
        XX = [px(1:substeps:end);py(1:substeps:end)];
        abs_err = abs(XX-X(:,end:-1:1));
        disp("MAE: " + num2str(mean(abs_err,'all')))
        break
    end
end
