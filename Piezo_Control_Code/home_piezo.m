obj = serial('COM3','BaudRate',19200,'DataBits',8,'Parity','none','StopBits',1,'FlowControl','software','Terminator','CR/LF');
fopen(obj);

fprintf(obj,'setk,1,1');
fprintf(obj,'cloop,1,1');

fprintf(obj,'setk,0,1');
fprintf(obj,'cloop,0,1');

% for i=1:1:100
%     fprintf(obj,['set,1,' num2str(i)]);
%     pause(0.1);
% end