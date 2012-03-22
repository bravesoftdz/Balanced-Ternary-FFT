unit BalancedTernaryFFT;

//Version 0.1

interface
uses math;
type
  ComplexTernary1D=class
    private
      Re_data,Im_data: array of Real;
      _q,_qmin1: Integer; //число тритов
      _T,_N: Integer; //полное число элементов и мин/макс значение (-N, N)
      base: Integer; //смещение нулевого отсчета
      function Re_value(i: Integer): Real;
      function Im_value(i: Integer): Real;
      procedure set_Re(i: Integer; value: Real);
      procedure set_Im(i: Integer; value: Real);
    public
      procedure Set_Length(new_T: Integer);
      procedure inversion;
      property Re[i: integer]: Real read Re_value write set_Re;
      property Im[i: integer]: Real read Im_value write set_Im;
      property N: Integer read _N;
//      procedure generalFFT(inverse: boolean);
      procedure FFT;
      procedure inverseFFT;
      constructor Create;
  end;

  ComplexTernary2D=class
    private
      Re_data, Im_data: array of Real; //адресацию сделаем свою
      _qx,_qy: Integer; //число тритов по осям x и y
      _Tx,_Ty,_Nx,_Ny: Integer; //полное число элементов и мин/макс значение (-N, N)
      _T: Integer; //всего элементов в массиве
      base: Integer; //смещение отсчета (0,0)
      function Re_value(x,y: Integer): Real;
      function Im_value(x,y: Integer): Real;
      procedure set_Re(x,y: Integer; value: Real);
      procedure set_Im(x,y: Integer; value: Real);
      procedure Inversion(fbase,mult,_T1: Integer);
      procedure general_FFT(fbase,mult,_T1: Integer; inverse: boolean);
    public
      procedure Set_Length(X,Y: Integer);
      procedure FFT;
      procedure InverseFFT;
      procedure power_spectrum;
      procedure Clear;

      property Nx: Integer read _Nx;
      property Ny: Integer read _Ny;
      property Tx: Integer read _Tx;
      property Ty: Integer read _Ty;
      property T: Integer read _T;
      property Re[x,y: Integer]: Real read Re_value write Set_Re;
      property Im[x,y: Integer]: Real read Im_value write Set_Im;


      constructor Create;
  end;



implementation

constructor ComplexTernary1D.Create;
begin
  inherited Create;
  Set_Length(0);
end;

procedure ComplexTernary1D.Set_Length(new_T: Integer);
begin
  assert(new_T>-1,'Set_Length: negative argument');
  if new_T<1 then begin
    _q:=-1;
    _qmin1:=-2;
    _T:=0;
    _N:=-1;
    SetLength(Re_data,0);
    SetLength(Im_data,0);
  end
  else begin
    _q:=math.Ceil(ln(new_T)/ln(3));
    _qmin1:=_q-1;
    _T:=Round(power(3,_q));
    _N:=(_T-1) div 2;
    SetLength(Re_data,_T);
    SetLength(Im_data,_T);
    base:=_N;
  end;
end;

function ComplexTernary1D.Re_value(i: Integer): Real;
begin
  assert((i<=_N) and (i>=-_N),'Re_value index out of range');
  Re_value:=Re_data[base+i];
end;

function ComplexTernary1D.Im_value(i: Integer): Real;
begin
  assert((i<=_N) and (i>=-_N),'Im_value index out of range');
  Im_value:=Im_data[base+i];
end;

procedure ComplexTernary1D.set_Re(i: Integer; value: Real);
begin
  assert((i<=_N) and (i>=-_N),'set_Re index out of range');
  Re_data[base+i]:=value;
end;

procedure ComplexTernary1D.set_Im(i: Integer; value: Real);
begin
  assert((i<=_N) and (i>=-_N),'set_Im index out of range');
  Im_data[base+i]:=value;
end;

procedure ComplexTernary1D.inversion;
var i,j,k,b: Integer;
    trits: array of Integer;
    tmp: Real;
begin
  if _q<0 then Exit;
  SetLength(trits,_q);
  for j:=0 to _qmin1 do trits[j]:=-1; //самый отрицательный элемент
  for i:=-_N to _N do begin
    k:=0;
    b:=1;
    for j:=_qmin1 downto 0 do begin
      k:=k+trits[j]*b;
      b:=b*3;
    end;
    //k указывает инверсию от i.
    if k>i then begin
    //поменяем местами
      tmp:=Re[k];
      Re[k]:=Re[i];
      Re[i]:=tmp;
      tmp:=Im[k];
      Im[k]:=Im[i];
      Im[i]:=tmp;
    end;
    //прибавим единичку
    j:=0;
    while j<_q do begin
      inc(trits[j]);
      if trits[j]<2 then break;
      trits[j]:=-1;
      inc(j);
    end;

  end;
end;

procedure ComplexTernary1D.FFT;
var N1,M1,T1,k,j,incr,big_incr,i: Integer;
  sqrt3,Wr,Wi,Ph,incWr,incWi,TwoPi,tmpWr: Real;
  //W - фазовый множитель, r,i - действ. и мнимое знач.
  xsum,ysum,xdif,ydif,ax,ay,xp1,xm1,yp1,ym1,x0,y0: Real;
  //sum - суммы
  //dif - разности
  //p1,0,m1 - +1,0,-1 соотв
begin
  sqrt3:=-sqrt(3)/2;
  TwoPi:=2*pi;
  inversion;
  if _q<1 then exit;
  //первая итерация - без фазовых множителей, она самая быстрая
  T1:=_T;
  N1:=_N;
  incr:=1;

  while N1>0 do begin
    T1:=T1 div 3;
    N1:=(T1-1) div 2;
    big_incr:=incr*3; //для внутреннего цикла
    M1:=(incr-1) div 2; //для внешнего
    //отдельно обработаем i=0, там фазовый множ. не нужен
    for k:=-N1 to N1 do begin
       j:=base+big_incr*k;
       //отдельно обраб. нулевое значение - там не нужно фаз. множителей
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        xp1:=Re_data[j];
        yp1:=Im_data[j];
        j:=j-2*incr;
        xm1:=Re_data[j];
        ym1:=Im_data[j];

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;

        //итого, 12 сложений и 4 умножения
       end;
    //шаг фазового множителя: 2pi/incr;
    //на первой итерации просто 2pi, но там цикл и не запустится
    //на второй итер:
    Ph:=TwoPi/big_incr;
    incWr:=cos(Ph);
    incWi:=-sin(Ph);
    Wr:=1;
    Wi:=0;
    for i:=1 to M1 do begin
      //пересчитываем фазовый множитель, потом делаем циклы для i и -i
      tmpWr:=Wr;
      Wr:=tmpWr*incWr-Wi*incWi;
      Wi:=Wi*incWr+tmpWr*incWi;
      for k:=-N1 to N1 do begin
        //итерация для +i
        j:=base+i+big_incr*k;
       //x0,y0 - без изменений
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        //а здесь надо умножить на фаз. множ.
        //элем. +1 - на W
        tmpWr:=Re_data[j];
        yp1:=Im_data[j];

        xp1:=tmpWr*Wr-yp1*Wi;
        yp1:=yp1*Wr+tmpWr*Wi;

        j:=j-2*incr;
        //элем. -1 умножаем на W* (сопряж)
        tmpWr:=Re_data[j];
        ym1:=Im_data[j];

        xm1:=tmpWr*Wr+ym1*Wi;
        ym1:=ym1*Wr-tmpWr*Wi;

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;

        //Теперь, то же самое для элемента -i

        j:=base-i+big_incr*k;
       //x0,y0 - без изменений
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        //а здесь надо умножить на фаз. множ.
        //элем. +1 - на W* (т.к -i)
        tmpWr:=Re_data[j];
        yp1:=Im_data[j];

        xp1:=tmpWr*Wr+yp1*Wi;
        yp1:=yp1*Wr-tmpWr*Wi;

        j:=j-2*incr;
        //элем. -1 умножаем на W
        tmpWr:=Re_data[j];
        ym1:=Im_data[j];

        xm1:=tmpWr*Wr-ym1*Wi;
        ym1:=ym1*Wr+tmpWr*Wi;

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;


       end;
    end;
    //конец одного слоя
    incr:=big_incr;
  end;






end;

procedure ComplexTernary1D.inverseFFT;
var N1,M1,T1,k,j,incr,big_incr,i: Integer;
  sqrt3,Wr,Wi,Ph,incWr,incWi,TwoPi,tmpWr: Real;
  //W - фазовый множитель, r,i - действ. и мнимое знач.
  xsum,ysum,xdif,ydif,ax,ay,xp1,xm1,yp1,ym1,x0,y0: Real;
  //sum - суммы
  //dif - разности
  //p1,0,m1 - +1,0,-1 соотв
begin
  for i:=-_N to _N do begin
    Set_Re(i,Re_value(i)/_T);
    Set_Im(i,Im_value(i)/_T);
  end;


  sqrt3:=sqrt(3)/2;
  TwoPi:=2*pi;
  inversion;
  if _q<1 then exit;
  //первая итерация - без фазовых множителей, она самая быстрая
  T1:=_T;
  N1:=_N;
  incr:=1;

  while N1>0 do begin
    T1:=T1 div 3;
    N1:=(T1-1) div 2;
    big_incr:=incr*3; //для внутреннего цикла
    M1:=(incr-1) div 2; //для внешнего
    //отдельно обработаем i=0, там фазовый множ. не нужен
    for k:=-N1 to N1 do begin
       j:=base+big_incr*k;
       //отдельно обраб. нулевое значение - там не нужно фаз. множителей
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        xp1:=Re_data[j];
        yp1:=Im_data[j];
        j:=j-2*incr;
        xm1:=Re_data[j];
        ym1:=Im_data[j];

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;

        //итого, 12 сложений и 4 умножения
       end;
    //шаг фазового множителя: 2pi/incr;
    //на первой итерации просто 2pi, но там цикл и не запустится
    //на второй итер:
    Ph:=TwoPi/big_incr;
    incWr:=cos(Ph);
    incWi:=sin(Ph);
    Wr:=1;
    Wi:=0;
    for i:=1 to M1 do begin
      //пересчитываем фазовый множитель, потом делаем циклы для i и -i
      tmpWr:=Wr;
      Wr:=tmpWr*incWr-Wi*incWi;
      Wi:=Wi*incWr+tmpWr*incWi;
      for k:=-N1 to N1 do begin
        //итерация для +i
        j:=base+i+big_incr*k;
       //x0,y0 - без изменений
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        //а здесь надо умножить на фаз. множ.
        //элем. +1 - на W
        tmpWr:=Re_data[j];
        yp1:=Im_data[j];

        xp1:=tmpWr*Wr-yp1*Wi;
        yp1:=yp1*Wr+tmpWr*Wi;

        j:=j-2*incr;
        //элем. -1 умножаем на W* (сопряж)
        tmpWr:=Re_data[j];
        ym1:=Im_data[j];

        xm1:=tmpWr*Wr+ym1*Wi;
        ym1:=ym1*Wr-tmpWr*Wi;

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;

        //Теперь, то же самое для элемента -i

        j:=base-i+big_incr*k;
       //x0,y0 - без изменений
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        //а здесь надо умножить на фаз. множ.
        //элем. +1 - на W* (т.к -i)
        tmpWr:=Re_data[j];
        yp1:=Im_data[j];

        xp1:=tmpWr*Wr+yp1*Wi;
        yp1:=yp1*Wr-tmpWr*Wi;

        j:=j-2*incr;
        //элем. -1 умножаем на W
        tmpWr:=Re_data[j];
        ym1:=Im_data[j];

        xm1:=tmpWr*Wr-ym1*Wi;
        ym1:=ym1*Wr+tmpWr*Wi;

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;


       end;
    end;
    //конец одного слоя
    incr:=big_incr;
  end;

end;


function ComplexTernary2D.Re_value(x,y: Integer): Real;
begin
  Assert(x<=_Nx,'Re_value: x too big');
  Assert(x>=-_Nx,'Re_value: x too small');
  Assert(y<=_Ny,'Re_value: y too big');
  Assert(y>=-_Ny,'Re_value: y too small');
  Re_value:=Re_data[base+y*_Tx+x];
end;

function ComplexTernary2D.Im_value(x,y: Integer): Real;
begin
  Assert(x<=_Nx,'Im_value: x too big');
  Assert(x>=-_Nx,'Im_value: x too small');
  Assert(y<=_Ny,'Im_value: y too big');
  Assert(y>=-_Ny,'Im_value: y too small');
  Im_value:=Im_data[base+y*_Tx+x];
end;

procedure ComplexTernary2D.set_Re(x,y: Integer; value: Real);
begin
  Assert(x<=_Nx,'set_Re: x too big');
  Assert(x>=-_Nx,'set_Re: x too small');
  Assert(y<=_Ny,'set_Re: y too big');
  Assert(y>=-_Ny,'set_Re: y too small');
  Re_data[base+y*_Tx+x]:=value;
end;

procedure ComplexTernary2D.set_Im(x,y: Integer; value: Real);
begin
  Assert(x<=_Nx,'set_Im: x too big');
  Assert(x>=-_Nx,'set_Im: x too small');
  Assert(y<=_Ny,'set_Im: y too big');
  Assert(y>=-_Ny,'set_Im: y too small');
  Im_data[base+y*_Tx+x]:=value;
end;

procedure ComplexTernary2D.Set_Length(x,y: Integer);
begin
  assert(x>-1,'Set_Length: negative argument x');
  assert(y>-1,'Set_Length: negative argument y');
  if x<1 then begin
    _qx:=-1;
    _Tx:=0;
    _Nx:=-1;
  end
  else begin
    _qx:=math.Ceil(ln(x)/ln(3));
    _Tx:=Round(power(3,_qx));
    _Nx:=(_Tx-1) div 2;
  end;

  if y<1 then begin
    _qy:=-1;
    _Ty:=0;
    _Ny:=-1;
  end
  else begin
    _qy:=math.Ceil(ln(y)/ln(3));
    _Ty:=Round(power(3,_qy));
    _Ny:=(_Ty-1) div 2;
  end;
  _T:=_Tx*_Ty;
  SetLength(Re_data,_T);
  SetLength(Im_data,_T);
  base:=(_T-1) div 2;
end;

constructor ComplexTernary2D.Create;
begin
  inherited Create;
  Set_Length(0,0);
end;

procedure ComplexTernary2D.general_FFT(fbase,mult,_T1: Integer; Inverse: boolean);
var N1,M1,T1,k,j,incr,big_incr,i: Integer;
  sqrt3,Wr,Wi,Ph,incWr,incWi,TwoPi,tmpWr: Real;
  //W - фазовый множитель, r,i - действ. и мнимое знач.
  xsum,ysum,xdif,ydif,ax,ay,xp1,xm1,yp1,ym1,x0,y0: Real;
  //sum - суммы
  //dif - разности
  //p1,0,m1 - +1,0,-1 соотв
begin
  sqrt3:=-sqrt(3)/2;
  TwoPi:=2*pi;
  if Inverse then begin
    sqrt3:=-sqrt3;
    TwoPi:=-TwoPi;
  end;
  inversion(fbase,mult,_T1);
  if _T1<2 then exit;
  //первая итерация - без фазовых множителей, она самая быстрая
  T1:=_T1;
  N1:=(T1-1) div 2;
  incr:=mult;

  while N1>0 do begin
    T1:=T1 div 3;
    N1:=(T1-1) div 2;
    big_incr:=incr*3; //для внутреннего цикла
    M1:=(incr-1) div (2*mult); //для внешнего
    //отдельно обработаем i=0, там фазовый множ. не нужен
    for k:=-N1 to N1 do begin
       j:=fbase+big_incr*k;
       //отдельно обраб. нулевое значение - там не нужно фаз. множителей
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        xp1:=Re_data[j];
        yp1:=Im_data[j];
        j:=j-2*incr;
        xm1:=Re_data[j];
        ym1:=Im_data[j];

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;

        //итого, 12 сложений и 4 умножения
       end;
    //шаг фазового множителя: 2pi/incr;
    //на первой итерации просто 2pi, но там цикл и не запустится
    //на второй итер:
    Ph:=TwoPi/big_incr*mult;
    incWr:=cos(Ph);
    incWi:=-sin(Ph);
    Wr:=1;
    Wi:=0;
    for i:=1 to M1 do begin
      //пересчитываем фазовый множитель, потом делаем циклы для i и -i
      tmpWr:=Wr;
      Wr:=tmpWr*incWr-Wi*incWi;
      Wi:=Wi*incWr+tmpWr*incWi;
      for k:=-N1 to N1 do begin
        //итерация для +i
        j:=fbase+i*mult+big_incr*k;
       //x0,y0 - без изменений
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        //а здесь надо умножить на фаз. множ.
        //элем. +1 - на W
        tmpWr:=Re_data[j];
        yp1:=Im_data[j];

        xp1:=tmpWr*Wr-yp1*Wi;
        yp1:=yp1*Wr+tmpWr*Wi;

        j:=j-2*incr;
        //элем. -1 умножаем на W* (сопряж)
        tmpWr:=Re_data[j];
        ym1:=Im_data[j];

        xm1:=tmpWr*Wr+ym1*Wi;
        ym1:=ym1*Wr-tmpWr*Wi;

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;

        //Теперь, то же самое для элемента -i

        j:=fbase-i*mult+big_incr*k;
       //x0,y0 - без изменений
        x0:=Re_data[j];
        y0:=Im_data[j];
        j:=j+incr;
        //а здесь надо умножить на фаз. множ.
        //элем. +1 - на W* (т.к -i)
        tmpWr:=Re_data[j];
        yp1:=Im_data[j];

        xp1:=tmpWr*Wr+yp1*Wi;
        yp1:=yp1*Wr-tmpWr*Wi;

        j:=j-2*incr;
        //элем. -1 умножаем на W
        tmpWr:=Re_data[j];
        ym1:=Im_data[j];

        xm1:=tmpWr*Wr-ym1*Wi;
        ym1:=ym1*Wr+tmpWr*Wi;

        xsum:=xp1+xm1;
        ysum:=yp1+ym1;
        ydif:=sqrt3*(xp1-xm1);
        xdif:=sqrt3*(ym1-yp1);
        // 4 сложения и 2 умножения (с плав. точкой)
        Ax:=x0-0.5*xsum;
        Ay:=y0-0.5*ysum;
        // 6 сложений и 4 умножения
        //сейчас j указывает на -1-й элемент
        Re_data[j]:=Ax-xdif;
        Im_data[j]:=Ay-ydif;

        j:=j+2*incr;
        //+1-й элемент
        Re_data[j]:=Ax+xdif;
        Im_data[j]:=Ay+ydif;

        j:=j-incr;
        //0-й элемент
        Re_data[j]:=x0+xsum;
        Im_data[j]:=y0+ysum;


       end;
    end;
    //конец одного слоя
    incr:=big_incr;
  end;

end;

procedure ComplexTernary2D.FFT;
var x,y: Integer;
begin
//сначала одномерный FFT по каждой строке
  for y:=-_Ny to _Ny do begin
    general_FFT(base+y*_Tx,1,_Tx,false);
  end;
//теперь по каждому столбцу
  for x:=-_Nx to _Nx do begin
    general_FFT(base+x,_Tx,_Ty,false);
  end;

end;

procedure ComplexTernary2D.InverseFFT;
var x,y: Integer;
begin

//по каждому столбцу, предварительно поделив на _T
  for x:=-_Nx to _Nx do begin
    for y:=-_Ny to _Ny do begin
      set_Re(x,y,Re_value(x,y)/_T);
      set_Im(x,y,Im_value(x,y)/_T);
    end;
    general_FFT(base+x,_Tx,_Ty,true);
  end;

  for y:=-_Ny to _Ny do begin
    general_FFT(base+y*_Tx,1,_Tx,true);
  end;

end;

procedure ComplexTernary2D.inversion(fbase,mult,_T1: integer);
var i,j,k,b,q,_qmin1,_N: Integer;
    trits: array of Integer;
    tmp: Real;
begin
  if _T1<4 then Exit;
  q:=ceil(ln(_T1)/ln(3));
  SetLength(trits,q);
  _qmin1:=q-1;
  _N:=(_T1-1) div 2;
  for j:=0 to _qmin1 do trits[j]:=-1; //самый отрицательный элемент
  for i:=-_N to _N do begin
    k:=0;
    b:=1;
    for j:=_qmin1 downto 0 do begin
      k:=k+trits[j]*b;
      b:=b*3;
    end;
    //k указывает инверсию от i.
    if k>i then begin
    //поменяем местами
      tmp:=Re_data[fbase+k*mult];
      Re_data[fbase+k*mult]:=Re_data[fbase+i*mult];
      Re_data[fbase+i*mult]:=tmp;
      tmp:=Im_data[fbase+k*mult];
      Im_data[fbase+k*mult]:=Im_data[fbase+i*mult];
      Im_data[fbase+i*mult]:=tmp;
    end;
    //прибавим единичку
    j:=0;
    while j<q do begin
      inc(trits[j]);
      if trits[j]<2 then break;
      trits[j]:=-1;
      inc(j);
    end;

  end;
end;

procedure ComplexTernary2D.power_spectrum;
var x,y: Integer;
begin
  for y:=-_Ny to _Ny do
    for x:=-_Nx to _Nx do
      Set_Re(x,y,Sqr(Re_value(x,y))+Sqr(Im_value(x,y)));
end;

procedure ComplexTernary2D.Clear;
var x: Integer;
begin
  for x:=0 to _T-1 do begin
    Re_data[x]:=0;
    Im_data[x]:=0;
  end;
end;






end.
