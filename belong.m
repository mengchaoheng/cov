function flag=belong(i,Mi)
flag=0;
switch i
   case 1
      if(Mi(1)<=0.39088 && Mi(1)>=0 && Mi(2)<=0.39088 && Mi(2)>=0 && Mi(3)<=0.3043 && Mi(3)>=0)
          flag=1;
      end
   case 2
      if(Mi(1)>=-0.39088 && Mi(1)<=0 && Mi(2)<=0.39088 && Mi(2)>=0 && Mi(3)<=0.3043 && Mi(3)>=0)
          flag=1;
      end
   case 3
      if(Mi(1)>=-0.39088 && Mi(1)<=0 && Mi(2)>=-0.39088 && Mi(2)<=0 && Mi(3)<=0.3043 && Mi(3)>=0)
          flag=1;
      end
   case 4
      if(Mi(1)<=0.39088 && Mi(1)>=0 && Mi(2)>=-0.39088 && Mi(2)<=0 && Mi(3)<=0.3043 && Mi(3)>=0)
          flag=1;
      end
   case 5
      if(Mi(1)<=0.39088 && Mi(1)>=0 && Mi(2)<=0.39088 && Mi(2)>=0 && Mi(3)>=-0.3043 && Mi(3)<=0)
          flag=1;
      end
   case 6
      if(Mi(1)>=-0.39088 && Mi(1)<=0 && Mi(2)<=0.39088 && Mi(2)>=0 && Mi(3)>=-0.3043 && Mi(3)<=0)
          flag=1;
      end
   case 7
      if(Mi(1)>=-0.39088 && Mi(1)<=0 && Mi(2)>=-0.39088 && Mi(2)<=0 && Mi(3)>=-0.3043 && Mi(3)<=0)
          flag=1;
      end
   case 8
      if(Mi(1)<=0.39088 && Mi(1)>=0 && Mi(2)>=-0.39088 && Mi(2)<=0 && Mi(3)>=-0.3043 && Mi(3)<=0)
          flag=1;
      end
   case 9
      if(Mi(1)-0.39088<=0.00001 && abs(Mi(2))/0.39088+abs(Mi(3))/0.1521<=1)
          flag=1;
      end
   case 10
      if(Mi(2)-0.39088<=0.00001 && abs(Mi(1))/0.39088+abs(Mi(3))/0.1521<=1)
          flag=1;
      end
   case 11
       if(Mi(1)+0.39088<=0.00001 && abs(Mi(2))/0.39088+abs(Mi(3))/0.1521<=1)
          flag=1;
      end
   case 12
      if(Mi(2)+0.39088<=0.00001 && abs(Mi(1))/0.39088+abs(Mi(3))/0.1521<=1)
          flag=1;
      end
end
end