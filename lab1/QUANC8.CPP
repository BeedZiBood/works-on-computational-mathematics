/* File converted with FORTRAN CONVERTER utility.
   FORTRAN CONVERTER is written by Grigoriev D., 2081/4.*/

/*MARKERS: //! = probably incorrect, //? = possibly incorrect*/

#include "FORSYTHE.H"
#include <iostream>

#include <functional>


void QUANC8(double (*FUN)(double x),REAL A,REAL B,REAL ABSERR,REAL RELERR,REAL &RESULT,REAL &ERREST,int &NOFUN,REAL &FLAG)
{
//     ОЦЕНИТЬ ИНТЕГРАЛ ДЛЯ FUN(X) ОТ А ДО В С ЗАДАННОЙ
//     С ПОЛЬЗОВАТЕЛЕМ ТОЧНОСТЬЮ.
//     С АВТОМАТИЧЕСКАЯ АДАПТИВНАЯ ПРОГРАММА, ОСНОВАНС НАЯ НА ФОРМУЛЕ НЬЮТОНА—КОТЕСА 8-ГО ПОРЯДКА

//     ВХОДНАЯ ИНФОРМАЦИЯ..
//
//    FUN ИМЯ ПОДПРОГРАММЫ-ФУНКЦИИ FUN(X), РЕАЛИС ЗУЮЩЕЙ ПОДЫНТЕГРАЛЬНУЮ ФУНКЦИЮ
//    То есть не выполнялся критерий точности.— Прим, перев.
//    ПОДПРОГРАММА QUANC8 119
//
//    А НИЖНИЙ ПРЕДЕЛ ИНТЕГРИРОВАНИЯ
//    В ВЕРХНИЙ ПРЕДЕЛ ИНТЕГРИРОВАНИЯ (В МОЖЕТ БЫТЬ МЕНЬШЕ, ЧЕМ А)
//    RELERR ГРАНИЦА ОТНОСИТЕЛЬНОЙ ПОГРЕШНОСТИ (ДОЛЖНА БЫТЬ НЕОТРИЦАТЕЛЬНА)
//    ABSERR ГРАНИЦА АБСОЛЮТНОЙ ПОГРЕШНОСТИ (ДОЛЖНА БЫТЬ НЕОТРИЦАТЕЛЬНА)
//
//    ВЫХОДНАЯ ИНФОРМАЦИЯ..
//    RESULT ПРИБЛИЖЕНИЕ К ИНТЕГРАЛУ, УДОВЛЕТВОРЯЮЩЕЕ, МОЖНО НАДЕЯТЬСЯ, МЕНЕЕ ЖЕСТКОЙ ИЗ
//    ДВУХ ГРАНИЦ ПОГРЕШНОСТИ.
//    ERREST ОЦЕНКА ВЕЛИЧИНЫ ДЕЙСТВИТЕЛЬНОЙ ОШИБКИ.
//    NOFUN ЧИСЛО ЗНАЧЕНИЙ ФУНКЦИИ, ИСПОЛЬЗОВАННЫХ ПРИ ВЫЧИСЛЕНИИ RESULT.
//    FLAG ИНДИКАТОР НАДЕЖНОСТИ. ЕСЛИ FLAG РАВЕН НУЛЮ, ТО RESULT, ВЕРОЯТНО, УДОВЛЕТВОРЯЕТ
//    ЗАДАННОЙ ГРАНИЦЕ ПОГРЕШНОСТИ. ЕСЛИ
//    FLAG=XXX.YYY, ТО ХХХ=ЧИСЛО ИНТЕРВАЛОВ, ДЛЯ КОТОРЫХ НЕ БЫЛО СХОДИМОСТИ, А
//    O.YYY = ЧАСТЬ ОСНОВНОГО ИНТЕРВАЛА, ОСТАВШАЯСЯ ДЛЯ ОБРАБОТКИ В ТОТ МОМЕНТ, КОГДА
//    ПРОГРАММА ПРИБЛИЗИЛАСЬ К ПРЕДЕЛЬНОМУ ЗНАЧЕНИЮ ДЛЯ NOFUN

REAL W0,W1,W2,W3,W4,AREA,X0,F0,STONE,STEP,COR11,TEMP;
REAL QPREV,QNOW,QDIFF,QLEFT,ESTERR,TOLERR;
REAL QRIGHT[31],F[16],X[16],FSAVE[8][30],XSAVE[8][30];

//     *•* ЭТАП 1 *** ПРИСВОЕНИЕ НАЧАЛЬНЫХ ЗНАЧЕНИЙ
//    ПЕРЕМЕННЫМ, НЕ ЗАВИСЯЩИМ ОТ ИНТЕРВАЛА. ГЕНЕРИРОВАНИЕ КОНСТАНТ

int LEVMIN,LEVMAX,LEVOUT,NOMAX,NOFIN,LEV,NIM,I,J;
LEVMIN=1;
LEVMAX=30;
LEVOUT=6;
NOMAX=5000;
int type_of_function = 0;
NOFIN=NOMAX-8*(LEVMAX-LEVOUT+(1<<(LEVOUT+1)));

//     ЕСЛИ NOFUN ДОСТИГАЕТ ЗНАЧЕНИЯ NOFIN, TO ТРЕВОГА

W0=3956.0/14175.0;
W1=23552.0/14175.0;
W2=-3712.0/14175.0;
W3=41984.0/14175.0;
W4=-18160.0/14175.0;

//     ПРИСВОИТЬ НУЛЕВЫЕ ЗНАЧЕНИЯ ПЕРЕМЕННЫМ СУММАМ

FLAG=0.0;
RESULT=0.0;
COR11=0.0;
ERREST=0.0;
AREA=0.0;
NOFUN=0;
if(A==B)return;  //?

//     ***ЭТАП 2*** ПРИСВОЕНИЕ НАЧАЛЬНЫХ ЗНАЧЕНИЙ
//    ПЕРЕМЕННЫМ, ЗАВИСЯЩИМ ОТ ИНТЕРВАЛА, В
//    СООТВЕТСТВИИ С ПЕРВЫМ ИНТЕРВАЛОМ

LEV=0;
NIM=1;
X0=A;
X[15]=B;
QPREV=0.0;
F0=FUN(X0);
STONE=(B-A)/16.0;
X[7]=(X0+X[15])/2.0;
X[3]=(X0+X[7])/2.0;
X[11]=(X[7]+X[15])/2.0;
X[1]=(X0+X[3])/2.0;
X[5]=(X[3]+X[7])/2.0;
X[9]=(X[7]+X[11])/2.0;
X[13]=(X[11]+X[15])/2.0;
for( J=2; J<=16; J+=2){   //? target=25
F[J-1]=FUN(X[J-1]);
//_25:;
}                         	// CONTINUE
NOFUN=9;

//     ***ЭТАП 3*** ОСНОВНЫЕ ВЫЧИСЛЕНИЯ
//     ТРЕБУЮТСЯ QPREV, ХО, Х2, Х4........ Х16, F0, F2, F4, ..., F16.
//     ВЫЧИСЛЯЮТСЯ XI, ХЗ, ..., Х15, Fl, F3, .... F15, QLEFT,
//     QRIGHT, QNOW, QDIFF, AREA.

_30:;
X[0]=(X0+X[1])/2.0;
F[0]=FUN(X[0]);
for( J=3; J<=15; J+=2){   //? target=35
X[J-1]=(X[J-2]+X[J])/2.0;
F[J-1]=FUN(X[J-1]);
//_35:;
}                         	// CONTINUE
NOFUN=NOFUN+8;
STEP=(X[15]-X0)/16.0;
QLEFT=(W0*(F0+F[7])+W1*(F[0]+F[6])+W2*(F[1]+F[5])+W3*(F[2]+F[4])+W4*F[3])*STEP;
QRIGHT[LEV]=(W0*(F[7]+F[15])+W1*(F[8]+F[14])+W2*(F[9]+F[13])+W3*(F[10]+F[12])+W4*F[11])*STEP;
QNOW=QLEFT+QRIGHT[LEV];
QDIFF=QNOW-QPREV;
AREA=AREA+QDIFF;

//    ***ЭТАП 4***. ПРОВЕРКА СХОДИМОСТИ ДЛЯ ИНТЕРВАЛА

ESTERR=ABS(QDIFF)/1023.0;
TOLERR=MAX(ABSERR,RELERR*ABS(AREA))*(STEP/STONE);
if(LEV<LEVMIN)goto _50;
if(LEV>=LEVMAX)
{
    type_of_function = 1;
    goto _62;
}
if(NOFUN>NOFIN)
{
    type_of_function = 2;
    goto _60;
}
if(ESTERR<=TOLERR)
{
    type_of_function = 3;
        goto _70;
}

//     ***ЭТАП 5*** СХОДИМОСТИ НЕТ
//     УСТАНОВИТЬ СЛЕДУЮЩИЙ ИНТЕРВАЛ.

_50:;
NIM=2*NIM;
LEV=LEV+1;

//     ЗАПОМНИТЬ ЭЛЕМЕНТЫ. ОТНОСЯЩИЕСЯ К ПРАВОЙ
//     ПОЛОВИНЕ ИНТЕРВАЛА, ДЛЯ БУДУЩЕГО ИСПОЛЬЗОВАНИЯ.

for( I=1; I<=8; I++){   //? target=52
FSAVE[I-1][LEV-1]=F[I+7];
XSAVE[I-1][LEV-1]=X[I+7];
//_52:;
}                         	// CONTINUE

//     СОБРАТЬ ЭЛЕМЕНТЫ, ОТНОСЯЩИЕСЯ К ЛЕВОЙ ПОЛОВИНЕ
//     ИНТЕРВАЛА, ДЛЯ НЕМЕДЛЕННОГО ИСПОЛЬЗОВАНИЯ

QPREV=QLEFT;
for( I=1; I<=8; I++){   //? target=55
J=-I;
F[2*J+17]=F[J+8];
X[2*J+17]=X[J+8];
//_55:;
}                         	// CONTINUE
goto _30;

//     ***ЭТАП 6*** «ПОЖАРНЫЙ» РАЗДЕЛ
//     ЧИСЛО ЗНАЧЕНИЙ ФУНКЦИИ БЛИЗКО К ТОМУ, ЧТОБЫ
//     ПРЕВЫСИТЬ УСТАНОВЛЕННЫЙ ПРЕДЕЛ.

_60:;
NOFIN=2*NOFIN;
LEVMAX=LEVOUT;
FLAG=FLAG+double((B-X0))/double((B-A));
goto _70;

//     ТЕКУЩЕЕ ПРЕДЕЛЬНОЕ ЗНАЧЕНИЕ ГЛУБИНЫ ДЕЛЕНИЯ
//     ПОПОЛАМ РАВНО LEVMAX

_62:;
FLAG=FLAG+1.0;

//     ***ЭТАП 7*** СХОДИМОСТЬ ДЛЯ ИНТЕРВАЛА ИМЕЕТ МЕСТО
//     ПРИБАВИТЬ ОЧЕРЕДНЫЕ СЛАГАЕМЫЕ К ПЕРЕМЕННЫМ
//     СУММАМ

_70:;
RESULT=RESULT+QNOW;
ERREST=ERREST+ESTERR;
COR11=COR11+QDIFF/1023.0;

//     УСТАНОВИТЬ СЛЕДУЮЩИЙ ИНТЕРВАЛ.

_72:;
if(NIM==2*(NIM/2))goto _75;
NIM=NIM/2;
LEV=LEV-1;
goto _72;
_75:;
NIM=NIM+1;
if(LEV<=0)goto _80;

//     СОБРАТЬ ЭЛЕМЕНТЫ, НЕОБХОДИМЫЕ ДЛЯ СЛЕДУЮЩЕГО
//    ИНТЕРВАЛА.

QPREV=QRIGHT[LEV-1];
X0=X[15];
F0=F[15];
for( I=1; I<=8; I++){   //? target=78
F[2*I-1]=FSAVE[I-1][LEV-1];
X[2*I-1]=XSAVE[I-1][LEV-1];
//_78:;
}                         	// CONTINUE
goto _30;

//     ***ЭТАП 8*** ЗАКЛЮЧИТЕЛЬНЫЕ ОПЕРАЦИИ И ВЫХОД

_80:;
RESULT=RESULT+COR11;

//     ОБЕСПЕЧИТЬ, ЧТОБЫ ЗНАЧЕНИЕ ПЕРЕМЕННОЙ ERREST
//     БЫЛО НЕ МЕНЬШЕ УРОВНЯ ОКРУГЛЕНИЙ.

if(ERREST==0.0)return;  //?
_82:;
TEMP=ABS(RESULT)+ERREST;
if(TEMP!=ABS(RESULT))return;  //?
ERREST=2.0*ERREST;
goto _82;
}                         	// END
