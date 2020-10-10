%one out kýsmý, döngü yapýlacak
% num=2
% y6=y1;
% y6(:,num)=[];

function [out_cvc,out_fdobjcv, out_cvmat,out_fdobj]=one_out(data)

out_data = cell(1,size(data,2));
filename = 'out.xlsx';

for i=1:size(data,2)
    i
    out=data;
    out(:,i)=[];
    [cvc, fdobjcv, fdobj,fdobjabscv]=cvcalc(out);
    out_cvc(:,i)=cvc;
    out_abscvc(:,i)=abs(cvc);
    out_fdobjcv(i)=fdobjabscv;
    cvmat = plot_ist(fdobjabscv);
    out_fdobj(i)=fdobj;
    out_cvmat(:,i)=cvmat;
%     xlswrite(filename,out,i+2);  
end

xlswrite(filename,out_cvmat,1)
xlswrite(filename,out_cvc,2)
xlswrite(filename,out_abscvc,3)


end