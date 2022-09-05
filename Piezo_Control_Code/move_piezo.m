%for i=1:1:250 %400 % 500*3
for i=1:1:350
    % fprintf(obj,['set,1,' num2str(i/10)]);
    % pause(0.15);
    
    fprintf(obj,['set,1,' num2str( i/(10) )]);
    pause(0.15);
end

% pause(60);
pause(10*60);

fprintf(obj,['set,0,' num2str(60)]);



% fprintf(obj, 'set,0,100');
