% Visualization / simulation of pi 
% Obtain the exact value of pi by increasing the number of simulations

clear all ;clc ;
% Radius of the Circle
R = 1 ;
shg
clf reset
set(gcf,'color','w','menubar','none','numbertitle','off','name','Monte Carlo Method ') ;
x = [0.; 0.];
h = plot(x(1),x(2),'.');
set(h,'markersize',1.5,'erasemode','none');
title('Calculation of pi using Monte Carlo Method','Color','r','fontweight','bold') ;
axis([-R R -R R]) ;
axis equal ;
axis off ;
% Stop and close toggle button
stop = uicontrol('style','toggle','string','stop','position',[115,20,40,20],.....
'background','white');
% Text boxes
tb1 =  uicontrol('style','text','position',[190 20 100 20],'background','w') ;
tb2 = uicontrol('style','text','position',[320 20 150 20],'background','w') ;
drawnow ;
%
cnt = 1;        % Count for Total Random Numbers
ccnt = 1 ;      % Count for Random Numbers which lie inside circle
tic
while ~get(stop,'value')
    % Generate Random Numbers between R and -R
    x = 2*R*rand(2,1)-R ; 
    set(h,'XData',x(1),'YData',x(2),'Color','b') ;
    drawnow
    cnt = cnt + 1;
    N = sprintf('N = %8.0f',cnt) ; % Update the total number of Random Numbers
    set(tb1,'String',N)
    RR = (x(1)^2+x(2)^2) ;
    if RR<R^2       % If point lies inside the Circle
        set(h,'XData',x(1),'YData',x(2),'Color','r') ;
        drawnow ;
        kssv = 4*ccnt/cnt ;
        Pi = sprintf('pi(approx)=%6.8f',kssv) ; % Values of pi approximately
        set(tb2,'String',Pi) ;
        ccnt = ccnt+1 ;
    end
end
toc 
% Close the figure window
set(stop,'style','pushbutton','string','close','callback','close(gcf)')