figure;hold on
[xx yy] = meshgrid(1:10);
colormap([1 0 0;0 0 1]) %red and blue
surf(xx,yy,rand(10),ones(10)); %first color (red)
surf(xx,yy,rand(10)+10); %second color(blue)
view(17,22)
