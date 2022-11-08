function showRobot(ax,T)
cla(ax);
view(ax,45,10);
axis(ax, 0.6*[-1 1 -1 1 -1 1]);
grid(ax, 'on');

jnt = zeros(3,8);

for i = 1:length(jnt)-1
    jnt(:,i+1) = T(1:3,4,i);
end

for i = 1:7
    line(ax, [jnt(1,i) jnt(1,i+1)], [jnt(2,i) jnt(2,i+1)],...
        [jnt(3,i) jnt(3,i+1)]);
    hold(ax,'on');
end
hold(ax,'off');