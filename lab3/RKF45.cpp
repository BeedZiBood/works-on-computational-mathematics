/* File converted with FORTRAN CONVERTER utility.
FORTRAN CONVERTER is written by Grigoriev D., 2081/4.*/

/*MARKERS: //! = probably incorrect, //? = possibly incorrect*/

#include "RKF45.h"

void FEHL(void(F)(REAL T,REAL*Y,REAL*YP),int NEQN,REAL *Y,REAL &T,REAL H,REAL *YP,REAL *F1,REAL *F2,REAL *F3,REAL *F4,REAL *F5,REAL *S)
{

//     МЕТОД РУНГЕ—КУТТА—ФЕЛЬБЕРГА ЧЕТВЕРТОГО-ПЯТОГО ПОРЯДКА
//
//     ПОДПРОГРАММА FEHL ИНТЕГРИРУЕТ СИСТЕМУ ИЗ NEQN ОБЫКНОВЕННЫХ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ ПЕРВОГО ПОРЯДКА СЛЕДУЮЩЕГО ВИДА
//        DY(I)/DT =F(T, Y(l), ...Y(NEQN)),
//     ГДЕ НАЧАЛЬНЫЕ ЗНАЧЕНИЯ Y(I) И НАЧАЛЬНЫЕ ПРОИЗС ВОДНЫЕ YP(I) ЗАДАНЫ В НАЧАЛЬНОЙ ТОЧКЕ Т. FEHL ПРОДОЛЖАЕТ РЕШЕНИЕ
//     НА ФИКСИРОВАННЫЙ ШАГ H И ПОМЕЩАЕТ В МАССИВ S(I) ПРИБЛИЖЕНИЕ К РЕШЕНИЮ В ТОЧКЕ Т+Н, ИМЕЮЩЕЕ ПЯТЫЙ ПОРЯДОК ТОЧНОСТИ
//     (ЛОКАЛЬНЫЙ ПОРЯДОК РАВЕН ШЕСТИ). F1, ..., F5—МАСС СИВЫ РАЗМЕРНОСТИ NEQN, НЕОБХОДИМЫЕ ВНУТРИ ПРОГРАММЫ.
//     В ФОРМУЛАХ ПРОИЗВЕДЕНА ГРУППИРОВКА С ЦЕЛЬЮ УМЕНЬШИТЬ ПОТЕРЮ ВЕРНЫХ ЗНАКОВ.
//     ЧТОБЫ МОЖНО БЫЛО РАЗЛИЧАТЬ РАЗНЫЕ НЕЗАВИСИМЫЕ АРГУМЕНТЫ, ПРИ ОБРАЩЕНИИ К FEHL НЕ СЛЕДУЕТ ЗАДАВАТЬ ДЛЯ Н ЗНАЧЕНИЕ,
//     МЕНЬШЕЕ УМНОЖЕНС НОЙ НА 13 ОШИБКИ ОКРУГЛЕНИЯ В Т.


REAL CH;
int K;

CH=H/4.0;
for( K=1; K<=NEQN; K++)F5[K-1]=Y[K-1]+CH*YP[K-1];

F(T+CH,F5,F1);
CH=3.0*H/32.0;
for( K=1; K<=NEQN; K++)F5[K-1]=Y[K-1]+CH*(YP[K-1]+3.0*F1[K-1]);
F(T+3.0*H/8.0,F5,F2);

CH=H/2197.0;
for( K=1; K<=NEQN; K++)F5[K-1]=Y[K-1]+CH*(1932.0*YP[K-1]+(7296.0*F2[K-1]-7200.0*F1[K-1]));
F(T+12.0*H/13.0,F5,F3);

CH=H/4104.0;
for( K=1; K<=NEQN; K++)F5[K-1]=Y[K-1]+CH*((8341.0*YP[K-1]-845.0*F3[K-1])+(29440.0*F2[K-1]-32832.0*F1[K-1]));
F(T+H,F5,F4);

CH=H/20520.0;
for( K=1; K<=NEQN; K++)F1[K-1]=Y[K-1]+CH*((-6080.0*YP[K-1]+(9295.0*F3[K-1]-5643.0*F4[K-1]))+(41040.0*F1[K-1]-28352.0*F2[K-1]));
F(T+H/2.0,F1,F5);

//     BВЫЧИСЛИТЬ ПРИБЛИЖЕННОЕ РЕШЕНИЕ В ТОЧКЕ T-RH

CH=H/7618050.0;
for( K=1; K<=NEQN; K++)S[K-1]=Y[K-1]+CH*((902880.0*YP[K-1]+(3855735.0*F3[K-1]-1371249.0*F4[K-1]))+(3953664.0*F2[K-1]+277020.0*F5[K-1]));

return;
}                         	// END


void RKFS(void(F)(REAL T,REAL*Y,REAL*YP),int NEQN,REAL *Y,REAL &T,REAL TOUT,REAL &RELERR,REAL &ABSERR,int &IFLAG,REAL *YP,REAL H,REAL *F1,REAL *F2,REAL *F3,REAL *F4,REAL *F5,REAL SAVRE,REAL SAVAE)
{

//    МЕТОД РУНГЕ—КУТТА—ФЕЛЬБЕРГА ЧЕТВЕРТОГО-ПЯТОГО ПОРЯДКА
//    RKFS ИНТЕГРИРУЕТ СИСТЕМУ ОБЫКНОВЕННЫХ ДИФФЕС РЕНЦИАЛЬНЫХ УРАВНЕНИЙ ПЕРВОГО ПОРЯДКА (СМ. КОМС МЕНТАРИЙ К RKF45).
//    МАССИВЫ YP, FI, F2. F3, F4 И F5 (РАЗМЕРНОСТИ ПО КРАЙНЕЙ МЕРЕ NEQN) И ПЕРЕМЕНС НЫЕ Н, SAVRE, SAVAE, NFE, КОР, INIT, JFLAG И KFLAG
//    ИСПОЛЬЗУЮТСЯ ВНУТРИ ПРОГРАММЫ И ВЫНЕСЕНЫ В СПИСОК ВЫЗОВА, ЧТОБЫ СОХРАНИТЬ ИХ ОПРЕДЕЛЕНС НОСТЬ ПРИ ПОВТОРНОМ ОБРАЩЕНИИ.
//    ПОЭТОМУ ИХ ЗНАС ЧЕНИЯ НЕ ДОЛЖНЫ ИЗМЕНЯТЬСЯ ПОЛЬЗОВАТЕЛЕМ. ВОЗМОЖНЫЙ ИНТЕРЕС ПРЕДСТАВЛЯЮТ ПАРАМЕТРЫ
//     YP — ПРОИЗВОДНАЯ ВЕКТОРА РЕШЕНИЯ В ТОЧКЕ Т
//     Н — ПРЕДПОЛАГАЕМЫЙ РАЗМЕР ШАГА ДЛЯ ОЧЕРЕДНОГО ЭТАПА
//     NFE—СЧЕТЧИК ЧИСЛА ВЫЧИСЛЕНИЙ ФУНКЦИИ

bool HFAILD,OUTPUT;
static int NFE,KOP,INIT,JFLAG,KFLAG;
REAL A,AE,DT,EE,EEOET,ESTTOL,ET,HMIN,RER,S,SCALE,TOL,TOLN,U26,EPSP1,EPS,YPK;
int K,MFLAG;


//     REMIN—ЭТО МИНИМАЛЬНОЕ ДОПУСТИМОЕ ЗНАЧЕНИЕ ДЛЯ RELERR. ПОПЫТКИ ПОЛУЧИТЬ ПО ЭТОЙ ПОДПРОГРАММЕ
//     БОЛЕЕ ВЫСОКУЮ ТОЧНОСТЬ ОБЫЧНО СТОЯТ ОЧЕНЬ ДОРОГО И ЗАЧАСТУЮ БЕЗУСПЕШНЫ

REAL REMIN= 1.E-12;

//     СТОИМОСТЬ СЧЕТА КОНТРОЛИРУЕТСЯ ТРЕБОВАНИЕМ, ЧТОБЫ КОЛИЧЕСТВО ВЫЧИСЛЕНИЙ ФУНКЦИИ БЫЛО ОГС РАНИЧЕНО ВЕЛИЧИНОЙ,
//     ПРИБЛИЗИТЕЛЬНО РАВНОЙ ЗНАС ЧЕНИЮ ПАРАМЕТРА MAXNFE. ПРИНЯТОЕ ЗДЕСЬ ЗНАЧЕС НИЕ ПРИМЕРНО СООТВЕТСТВУЕТ 500 ШАГАМ.

int MAXNFE= 3000;

//     ПРОВЕРИТЬ ВХОДНЫЕ ПАРАМЕТРЫ

if(NEQN<1)goto _10;
if((RELERR<0.0)||(ABSERR<0.0))goto _10;
MFLAG=ABS(IFLAG);
if((MFLAG==0)||(MFLAG>8))goto _10;
if(MFLAG!=1)goto _20;

//     ПЕРВЫЙ ВЫЗОВ, ВЫЧИСЛИТЬ МАШИННОЕ ЭПСИЛОН

EPS=1.0;
_5:;
EPS=EPS/2.0;
EPSP1=EPS+1.;
if(EPSP1>1.)goto _5;
U26=26.*EPS;
goto _50;

//     ОШИБКА ВО ВХОДНОЙ ИНФОРМАЦИИ

_10:;
IFLAG=8;
return;  //?

//     ПРОВЕРИТЬ возможность ПРОДОЛЖЕНИЯ

_20:;
if((T==TOUT)&&(KFLAG!=3))goto _10;
if(MFLAG!=2)goto _25;

//     IFLAG = +2 ИЛИ -2

if((KFLAG==3)||(INIT==0))goto _45;
if(KFLAG==4)goto _40;
if((KFLAG==5)&&(ABSERR==0.0))goto _30;
if((KFLAG==6)&&(RELERR<=SAVRE)&&(ABSERR<=SAVAE))goto _30;
goto _50;

//     IFLAG = 3,4,5,6,7 ИЛИ 8

_25:;
if(IFLAG==3)goto _45;
if(IFLAG==4)goto _40;
if((IFLAG==5)&&(ABSERR>0.0))goto _45;

//     ИНТЕГРИРОВАНИЕ НЕЛЬЗЯ ПРОДОЛЖАТЬ, ПОСКОЛЬКУ ПОЛЬЗОВАТЕЛЬ НЕ ВЫПОЛНИЛ ИНСТРУКЦИЙ, СООТВЕТСТВУЮЩИХ
//     ЗНАЧЕНИЯМ IFLAG = 5, 6, 7 ИЛИ 8

_30:;
//printf("�HTE�P�POBAH�E �PEPBAHO, �OCKO��K� �O���OBATE��"
//"HE B��O�H�� �HCTP�K��� RKF45, COOTBETCTB����X"
//"�HA�EH��M IFLAG=5,6,7 ��� 8");
exit(0);

//     ПЕРЕОПРЕДЕЛИТЬ СЧЕТЧИК ЧИСЛА ВЫЧИСЛЕНИЙ ФУНКЦИИ

_40:;
NFE=0;
if(MFLAG==2)goto _50;

//     ПЕРЕОПРЕДЕЛИТЬ ЗНАЧЕНИЕ FLAG, УСТАНОВЛЕННОЕ ПРИ ПРЕДЫДУЩЕМ ОБРАЩЕНИИ

_45:;
IFLAG=JFLAG;
if(KFLAG==3)MFLAG=ABS(IFLAG);

//     СОХРАНИТЬ ВХОДНОЕ ЗНАЧЕНИЕ IFLAG И УСТАНОВИТЬ ЗНАЧЕНИЕ FLAG, СООТВЕТСТВУЮЩЕЕ ПРОДОЛЖЕНИЮ, ДЛЯ БУДУЩЕЙ ВХОДНОЙ ПРОВЕРКИ

_50:;
JFLAG=IFLAG;
KFLAG=0;

//     СОХРАНИТЬ ЗНАЧЕНИЯ RELERR И ABSERR ДЛЯ ВХОДНОЙ ПРОВЕРКИ ПРИ ПОСЛЕДУЮЩИХ ОБРАЩЕНИЯХ

SAVRE=RELERR;
SAVAE=ABSERR;

//     УСТАНОВИТЬ ЗНАЧЕНИЕ ГРАНИЦЫ ДЛЯ ОТНОСИТЕЛЬС НОЙ ПОГРЕШНОСТИ, РАВНОЕ КАК МИНИМУМ 2*EPS +REMIN,
//     ЧТОБЫ ИЗБЕЖАТЬ ТРУДНОСТЕЙ, СВЯЗАННЫХ С ТРЕБОВАНИЕМ НЕДОСТИЖИМОЙ ТОЧНОСТИ

RER=2.*EPS+REMIN;
if(RELERR>=RER)goto _55;

//     ЗАДАННАЯ ГРАНИЦА ОТНОСИТЕЛЬНОЙ ПОГРЕШНОСТИ СЛИШКОМ МАЛА

RELERR=RER;
IFLAG=3;
KFLAG=3;
return;  //?

_55:;
DT=TOUT-T;

if(MFLAG==1)goto _60;
if(INIT==0)goto _65;
goto _80;

//     ПРИСВОЕНИЕ НАЧАЛЬНЫХ ЗНАЧЕНИЙ (ИНИЦИИРОВАНИЕ) — УСТАНОВИТЬ ЗНАЧЕНИЕ УКАЗАТЕЛЯ
//     ОКОНЧАНИЯ ИНИЦИИРОВАНИЯ, IN1T УСТАНОВИТЬ ЗНАЧЕНИЕ УКАЗАТЕЛЯ СЛИШКОМ БОЛЬШОГО ЗАТРЕБОВАННОГО ЧИСЛА ВЫХОДНЫХ ТОЧЕК,
//     КОР ВЫЧИСЛИТЬ НАЧАЛЬНЫЕ ПРОИЗВОДНЫЕ, УСТАНОВИТЬ ЗНАЧЕНИЕ СЧЕТЧИКА ЧИСЛА ’ВЫЧИСЛЕНИЙ ФУНКЦИИ, NFE  ОЦЕНИТЬ НАЧАЛЬНУЮ ВЕЛИЧИНУ ШАГА

_60:;
INIT=0;
KOP=0;

A=T;
F(A,Y,YP);
NFE=1;
if(T!=TOUT)goto _65;
IFLAG=2;
return;  //?

_65:;
INIT=1;
H=ABS(DT);
TOLN=0.;
for( K=1; K<=NEQN; K++){   //? target=70
 TOL=RELERR*ABS(Y[K-1])+ABSERR;
 if(TOL<=0)goto _70;
 TOLN=TOL;
 YPK=ABS(YP[K-1]);
 if(YPK*pow(H,5)>TOL)H=pow((TOL/YPK),0.2);
_70:;
}                         	// CONTINUE
if(TOLN<=0.0)H=0.0;
H=MAX(H,U26*MAX(ABS(T),ABS(DT)));
JFLAG=SIGN(2,IFLAG);

//     ПРИСВОИТЬ ВЕЛИЧИНЕ ШАГА ЗНАК, СООТВЕТСТВУЮЩИМ ИНТЕГРИРОВАНИЮ В НАПРАВЛЕНИИ ОТ Т К TOUT

_80:;
H=SIGN(H,DT);

//     ПРОВЕРКА, НАСКОЛЬКО СЕРЬЕЗНО ВЛИЯНИЕ НА RKF45 СЛИШКОМ БОЛЬШОГО ЗАТРЕБОВАННОГО ЧИСЛА ВЫХОДНЫХ ТОЧЕК

if(ABS(H)>=2.0*ABS(DT))KOP=KOP+1;
if(KOP!=100)goto _85;

//     ЧРЕЗМЕРНАЯ ЧАСТОТА ВЫХОДОВ

KOP=0;
IFLAG=7;
return;

_85:;
if(ABS(DT)>U26*ABS(T))goto _95;

//     ЕСЛИ ОЧЕНЬ БЛИЗКО К ТОЧКЕ ВЫХОДА, ПРОЭКСТРАПОЛИРОВАТЬ И ВЕРНУТЬСЯ ПО МЕСТУ ВЫЗОВА

for( K=1; K<=NEQN; K++)Y[K-1]=Y[K-1]+DT*YP[K-1];
A=TOUT;
F(A,Y,YP);
NFE=NFE+1;
goto _300;

//     ПРИСВОИТЬ НАЧАЛЬНОЕ ЗНАЧЕНИЕ ИНДИКАТОРУ ТОЧКИ ВЫХОДА

_95:;
OUTPUT=false;

//     ЧТОБЫ ИЗБЕЖАТЬ НЕОПРАВДАННОГО МАШИННОГО НУЛЯ  ПРИ ВЫЧИСЛЕНИИ ФУНКЦИИ ОТ ГРАНИЦ ПОГРЕШНОСТЕЙ, ПРОМАСШТАБИРОВАТЬ ЭТИ ГРАНИЦЫ

SCALE=2./RELERR;
AE=SCALE*ABSERR;

//     ПОШАГОВОЕ ИНТЕГРИРОВАНИЕ

_100:;
HFAILD=false;

//    УСТАНОВИТЬ НАИМЕНЬШУЮ ДОПУСТИМУЮ ВЕЛИЧИНУ ШАГА

HMIN=U26*ABS(T);

//     ИСПРАВИТЬ ПРИ НЕОБХОДИМОСТИ ВЕЛИЧИНУ ШАГА, ЧТОБЫ ДОСТИГНУТЬ ТОЧКИ ВЫХОДА. РАССЧИТАТЬ НА ДВА ШАГА ВПЕРЕД,
//     ЧТОБЫ ИЗБЕЖАТЬ СЛИШКОМ РЕЗКИХ ИЗМЕНЕНИЙ В ВЕЛИЧИНЕ ШАГА И ТЕМ САМЫМ УМЕНЬШИТЬ ВЛИЯНИЕ ВЫХОДНЫХ ТОЧЕК НА ПРОГРАММУ.

DT=TOUT-T;
if(ABS(DT)>=2.*ABS(H))goto _200;
if(ABS(DT)>ABS(H))goto _150;

//     СЛЕДУЮЩИЙ УСПЕШНЫЙ ШАГ ЗАВЕРШИТ ИНТЕГРИРОВАНИЕ ДО УКАЗАННОЙ ТОЧКИ ВЫХОДА

OUTPUT=true;
H=DT;
goto _200;

_150:;
H=0.5*DT;



//     ВНУТРЕННИЙ ОДНОШАГОВЫЙ ИНТЕГРАТОР
//
//    ГРАНИЦЫ ПОГРЕШНОСТЕЙ БЫЛИ ПРОМАСШТАБИРОВАНЫ, ЧТОБЫ ИЗБЕЖАТЬ НЕОПРАВДАННОГО МАШИННОГО НУЛЯ
//    ПРИ ВЫЧИСЛЕНИИ ФУНКЦИИ ЕТ ОТ НИХ. ЧТОБЫ ИЗБЕЖАТЬ ОБРАЩЕНИЯ В НУЛЬ ЗНАМЕНАТЕЛЯ В ТЕСТЕ,
//    ОТНОСИТЕЛЬНАЯ ОШИБКА ИЗМЕРЯЕТСЯ ПО ОТНОШЕНИЮ К СРЕДНЕМУ ИЗ ВЕЛИЧИН РЕШЕНИЯ В НАЧАЛЕ И КОНЦЕ ШАГА.
//    В ФОРМУЛЕ, ОЦЕНИВАЮЩЕЙ ОШИБКУ, ПРОИЗВЕДЕНА ГРУППИРОВКА СЛАГАЕМЫХ, УМЕНЬШАЮЩАЯ ПОТЕРЮ ВЕРНЫХ ЗНАКОВ.
//    ЧТОБЫ РАЗЛИЧАТЬ МЕЖДУ СОБОЙ РАЗНЫЕ АРГУМЕНТЫ, ДЛЯ Н НЕ ДОПУСКАЮТСЯ ЗНАЧЕНИЯ, МЕНЬШИЕ УМНОЖЕННОЙ НА 26 ОШИБКИ ОКРУГЛЕНИЯ В Т.
//    ВВЕДЕНЫ ПРАКТИЧЕСКИЕ ОГРАНИЧЕНИЯ НА СКОРОСТЬ ИЗМЕНЕНИЯ ВЕЛИЧИНЫ ШАГА, ЧТОБЫ СГЛАДИТЬ ПРОЦЕСС ВЫБОРА ЭТОЙ ВЕЛИЧИНЫ
//    И ИЗБЕЖАТЬ ЧРЕЗМЕРС НОГО ЕЕ РАЗБРОСА В ЗАДАЧАХ С НАРУШЕНИЕМ НЕПРЕС РЫВНОСТИ. ИЗ ПРЕДОСТОРОЖНОСТИ ПРОГРАММА БЕРЕТ 9/10 ОТ ТОЙ
//    ВЕЛИЧИНЫ ШАГА, КОТОРАЯ НУЖНА ПО ЕЕ ОЦЕНКЕ. ЕСЛИ НА ДАННОМ ШАГЕ БЫЛА НЕУДАЧНАЯ ПОПЫТКА, ТО ПРИ ПЛАНИРОВАНИИ СЛЕДУЮЩЕГО УВЕЛИЧЕНИЕ
//    ДЛИНЫ ШАГА НЕ ДОПУСКАЕТСЯ- ЭТО ПОВЫШАЕТ ЭФФЕКС ТИВНОСТЬ ПРОГРАММЫ ДЛЯ ЗАДАЧ С РАЗРЫВАМИ И В ОБЩЕМ СЛУЧАЕ,
//    ПОСКОЛЬКУ ИСПОЛЬЗУЕТСЯ ЛОКАЛЬС НАЯ ЭКСТРАПОЛЯЦИЯ И ДОПОЛНИТЕЛЬНАЯ ПРЕДОСТОС РОЖНОСТЬ КАЖЕТСЯ ОПРАВДАННОЙ.
//
//
//   ПРОВЕРИТЬ ЧИСЛО ВЫЧИСЛЕНИЙ ПРОИЗВОДНЫХ. ЕСЛИ
//   ОНО НЕ ПРЕВЫШАЕТ УСТАНОВЛЕННОГО ПРЕДЕЛА, ПОС ПРОБОВАТЬ ПРОДОЛЖИТЬ ИНТЕГРИРОВАНИЕ С Т ДО Т + Н

_200:;
if(NFE<=MAXNFE)goto _220;

//     СЛИШКОМ БОЛЬШАЯ РАБОТА

IFLAG=4;
KFLAG=4;
return;

//     ПРОДОЛЖИТЬ ПРИБЛИЖЕННОЕ РЕШЕНИЕ НА ОДИН ШАГ ДЛИНЫ Н

_220:;
FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1);
NFE=NFE+5;

//     ВЫЧИСЛИТЬ И СРАВНИТЬ ДОПУСТИМЫЕ ГРАНИЦЫ И ОЦЕНКИ ЛОКАЛЬНОЙ ОШИБКИ, А ЗАТЕМ СНЯТЬ МАСШТАС БИРОВАНИЕ ГРАНИЦ.
//     ЗАМЕТЬТЕ, ЧТО ОТНОСИТЕЛЬНАЯ ОШИБКА ИЗМЕРЯЕТСЯ ПО ОТНОШЕНИЮ К СРЕДНЕМУ ИЗ ВЕЛИЧИН РЕШЕНИЯ В НАЧАЛЕ И КОНЦЕ ШАГА,

EEOET=0.;
for( K=1; K<=NEQN; K++){   //? target=250
 ET=ABS(Y[K-1])+ABS(F1[K-1])+AE;
 if(ET>0.)goto _240;

 //     HЕПРАВИЛЬНАЯ ГРАНИЦА ПОГРЕШНОСТИ

 IFLAG=5;
 KFLAG=5;
 return;

_240:;
 EE=ABS((-2090.*YP[K-1]+(21970.*F3[K-1]-15048.*F4[K-1]))+(22528.*F2[K-1]-27360.*F5[K-1]));
 //_250:;
 EEOET=MAX(EEOET,EE/ET);
}

ESTTOL=ABS(H)*EEOET*SCALE/752400.;

if(ESTTOL<=1.0)goto _260;


//     НЕУДАЧНЫЙ ШАГ
//      УМЕНЬШИТЬ ВЕЛИЧИНУ ШАГАИ СНОВА ПОПРОБОВАТЬ
//      УМЕНЬШЕНИЕ ОГРАНИЧИВАЕТСЯ СНИЗУ МНОЖИТЕЛЕМ 1/10

HFAILD=true;
OUTPUT=false;
S=0.1;
if(ESTTOL<59049.)S=0.9/pow(ESTTOL,0.2);
H=S*H;
if(ABS(H)>HMIN)goto _200;

//     ЗАДАННАЯ ГРАНИЦА ОШИБКИ НЕДОСТИЖИМА ДАЖЕ ПРИ НАИМЕНЬШЕЙ ДОПУСТИМОЙ ВЕЛИЧИНЕ ШАГА

IFLAG=6;
KFLAG=6;
return;  //?


//     УСПЕШНЫЙ ШАГ
//       ПОМЕСТИТЬ В МАССИВ Y РЕШЕНИЕ В ТОЧКЕ Т+Н И ВЫЧИСЛИТЬ ПРОИЗВОДНЫЕ В ЭТОЙ ТОЧКЕ

_260:;
T=T+H;
for( K=1; K<=NEQN; K++)Y[K-1]=F1[K-1];
A=T;
F(A,Y,YP);
NFE=NFE+1;


//     ВЫБРАТЬ ВЕЛИЧИНУ СЛЕДУЮЩЕГО ШАГА
//     УВЕЛИЧЕНИЕ ОГРАНИЧЕНО МНОЖИТЕЛЕМ 5
//     ЕСЛИ НА ДАННОМ ШАГЕ БЫЛА НЕУДАЧНАЯ ПОПЫТКА, ТО ДЛЯ СЛЕДУЮЩЕГО НЕ ДОПУСКАЕТСЯ ВЫБОР БОЛЬШЕЙ ВЕЛИЧИНЫ ШАГА

S=5.;
if(ESTTOL>1.889568E-4)S=0.9/pow(ESTTOL,0.2);
if(HFAILD)S=MIN(S,1.0);
H=SIGN(MAX(S*ABS(H),HMIN),H);

//     КОНЕЦ ОДНОШАГОВОГО ИНТЕГРАТОРА
//     НУЖНО ЛИ ДЕЛАТЬ ОЧЕРЕДНОЙ ШАГ

if(OUTPUT)goto _300;
if(IFLAG>0)goto _100;


//     ИНТЕГРИРОВАНИЕ УСПЕШНО ЗАВЕРШЕНО
//     РЕЖИМ ОДНОШАГОВОГО ИНТЕГРИРОВАНИЯ

IFLAG=-2;
return;

//     РЕЖИМ ИНТЕГРИРОВАНИЯ НА ИНТЕРВАЛЕ

_300:;
T=TOUT;
IFLAG=2;
return;

}                         	// END

void RKF45(void(F)(REAL T,REAL*Y,REAL*YP),int NEQN,REAL *Y,REAL &T,REAL TOUT,REAL &RELERR,REAL &ABSERR,REAL *WORK,int &IFLAG)
{
//     МЕТОД РУНГЕ —КУТТА—ФЕЛЬБЕРГА ЧЕТВЕРТОГО-ПЯТОГО ПОРЯДКА
//
//    СОСТАВИТЕЛИ ПРОГРАММЫ —Н. A. WATTS, L. F. SHAMPINE,
//    SANDIA LABORATORIES, ALBUQUERQUE, NEW MEXICO
//
//    RKF45 ПРЕДНАЗНАЧЕНА ГЛАВНЫМ ОБРАЗОМ ДЛЯ РЕШЕНИЯ НЕЖЕСТКИХ И СЛАБО ЖЕСТКИХ ДИФФЕРЕНЦИАЛЬНЫХ
//    УРАВНЕНИЙ, КОГДА ВЫЧИСЛЕНИЕ ПРОИЗВОДНЫХ НЕ СЛИШКОМ ДОРОГОСТОЯЩЕЕ. RKF45, ВООБЩЕ ГОВОРЯ,
//    НЕ СЛЕДУЕТ ИСПОЛЬЗОВАТЬ. ЕСЛИ ПОЛЬЗОВАТЕЛЮ ТРЕБУЕТСЯ ВЫСОКАЯ ТОЧНОСТЬ
//
//     РЕЗЮМЕ
//
//     ПОДПРОГРАММА RKF45 ИНТЕГРИРУЕТ СИСТЕМУ ИЗ NEQN ОБЫКНОВЕННЫХ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ
//     ПЕРВОГО ПОРЯДКА СЛЕДУЮЩЕГО ВИДА: DY(I)/DT = F(T, Y(l). Y(2), ..., Y(NEQN)),
//     ГДЕ Y(I) ЗАДАНЫ В Т. ОБЫЧНО ПОДПРОГРАММУ ПРИМЕНЯЮТ ДЛЯ ИНТЕГРИРОС ВАНИЯ ОТ Т ДО TOUT, ОДНАКО ЕЕ МОЖНО ИСПОЛЬЗОВАТЬ
//     И КАК ОДНОШАГОВЫЙ ИНТЕГРАТОР, ЧТОБЫ ПРОДОЛЖИТЬ РЕШЕНИЕ НА ОДИН ШАГ В НАПРАВЛЕНИИ TOUT. НА
//     ВЫХОДЕ ПАРАМЕТРАМ, ФИГУРИРУЮЩИМ В СПИСКЕ ВЫЗОВА, ПРИСВАИВАЮТСЯ ЗНАЧЕНИЯ, НЕОБХОДИМЫЕ
//     ДЛЯ ПРОДОЛЖЕНИЯ ИНТЕГРИРОВАНИЯ. ПОЛЬЗОВАТЕЛЮ НУЖНО ЛИШЬ ЕЩЕ РАЗ ОБРАТИТЬСЯ К RKF45 (И, ВОЗМОЖНО,
//     ОПРЕДЕЛИТЬ НОВОЕ ЗНАЧЕНИЕ ДЛЯ TOUT). В ДЕЙСТВИТЕЛЬНОСТИ RKF45-ЭТО ПРОГРАММА ИНТЕРС ФЕЙСА, КОТОРАЯ ВЫЗЫВАЕТ ПОДПРОГРАММУ RKFS,
//     ОСУЩЕСТВЛЯЮЩУЮ ПРОЦЕСС РЕШЕНИЯ. RKFS В СВОЮ ОЧЕРЕДЬ ВЫЗЫВАЕТ ПОДПРОГРАММУ FEHL, КОТОРАЯ ВЫЧИСЛЯЕТ ПРИБЛИЖЕННОЕ РЕШЕНИЕ НА ОДИН ШАГ.
//
//     RKF45 ИСПОЛЬЗУЕТ МЕТОД РУНГЕ—КУТТА—ФЕЛЬБЕРГА, ОПИСАННЫЙ В СЛЕДУЮЩЕЙ ПУБЛИКАЦИИ: Е. FEHLBERG,
//     LOW-ORDER CLASSICAL RUNGE—KUTTA FORMULAS WITH STEPSIZE CONTROL, NASA TR R-315
//
//     СТИЛЬ РАБОТЫ ПРОГРАММЫ RKF45 ИЛЛЮСТРИРУЕТСЯ В  СЛЕДУЮЩИХ ПУБЛИКАЦИЯХ: L. F. SHAMPINE, Н. A. WATTS,
//     S.DAVENPORT, SOLVING NON-STIFF ORDINARY DIFFERENC TIAL EQUATIONS—THE STATE OF THE ART, SANDIA
//     LABORATORIES REPORT SAND75-0182, SIAM REVIEW, 18(1976),  N3, 376—411.
//
//
//     ПАРАМЕТРЫ ПРОГРАММЫ:
//     F—ПОДПРОГРАММА F(T, Y, YP) ДЛЯ ВЫЧИСЛЕНИЯ ПРОИЗВОДНЫХ YP(I) = DY(I)/DT
//     NEQN—ЧИСЛО ИНТЕГРИРУЕМЫХ УРАВНЕНИЙ
//     Y(*)—РЕШЕНИЕ В ТОЧКЕ Т
//     Т—НЕЗАВИСИМАЯ ПЕРЕМЕННАЯ
//     TOUT—ТОЧКА ВЫХОДА, В КОТОРОЙ НУЖНО ОПРЕДЕС ЛИТЬ ЗНАЧЕНИЕ РЕШЕНИЯ
//     RELERR, ABSERR—ГРАНИЦЫ АБСОЛЮТНОЙ И ОТНОС СИТЕЛЬНОЙ ПОГРЕШНОСТИ ДЛЯ ТЕСТА ЛОКАЛЬНОЙ ОШИБКИ. НА КАЖДОМ ШАГЕ ПРОГРАММА
//        ТРЕБУЕТ ВЫПОЛНЕНИЯ УСЛОВИЯ ABS (LOCAL ERROR). LE. RELERR *ABS(Y)+ABSERR ДЛЯ КАЖДОЙ КОМПОНЕНТЫ ВЕКТОРОВ
//        ЛОКАЛЬНОЙ ОШИБКИ И РЕШЕНИЯ
//     IFLAG—УКАЗАТЕЛЬ РЕЖИМА ИНТЕГРИРОВАНИЯ.
//     WORK (*)—МАССИВ, СОДЕРЖАЩИЙ ИНФОРМАЦИЮ, ВНУТРЕННЮЮ ДЛЯ RKF45, КОТОРАЯ НЕОБХОДИМА ПРИ ПОСЛЕДУЮЩИХ ВЫЗОВАХ. ЕГО РАЗМЕРНОСТЬ
//        ДОЛЖНА БЫТЬ НЕ МЕНЬШЕ 3+6*NEQN
//     IWORK(*)—ЦЕЛЫЙ МАССИВ, СОДЕРЖАЩИЙ ИНФОРМАЦИЮ, ВНУТРЕННЮЮ ДЛЯ RKF45, КОТОРАЯ НЕОБХОДИМА ПРИ ПОСЛЕДУЮЩИХ ВЫЗОВАХ.
//        ЕГО РАЗМЕРНОСТЬ ДОЛЖНА БЫТЬ НЕ МЕНЬШЕ 5.
//
//
//     ПЕРВОЕ ОБРАЩЕНИЕ К RKF45
//     ПОЛЬЗОВАТЕЛЬ ДОЛЖЕН ПРЕДУСМОТРЕТЬ В СВОЕЙ ВЫЗЫВАЮЩЕЙ ПРОГРАММЕ ПАМЯТЬ ДЛЯ СЛЕДУЮЩИХ МАССИВОВ, ФИГУРИРУЮЩИХ В СПИСКЕ ВЫЗОВА—
//        Y(NEQN), WORK(3+6*NEQN), IWORK(5);
//     КРОМЕ ТОГО, ОН ДОЛЖЕН ОБЪЯВИТЬ F В ОПЕРАТОРЕ EXTERNAL, ПОДГОТОВИТЬ ПОДПРОГРАММУ F(T, Y, YP) И ПРИСВОИТЬ
//    НАЧАЛЬНЫЕ ЗНАЧЕНИЯ ПАРАМЕТРАМ —
//    NEQN — ЧИСЛО ИНТЕГРИРУЕМЫХ УРАВНЕНИЙ (NEQN.GE. 1)
//    Y(*) — ВЕКТОР НАЧАЛЬНЫХ УСЛОВИЙ
//    Т — НАЧАЛЬНАЯ ТОЧКА ИНТЕГРИРОВАНИЯ, Т ДОЛЖНО БЫТЬ ПЕРЕМЕННОЙ
//    TOUT — ТОЧКА ВЫХОДА, В КОТОРОЙ НУЖНО НАЙТИ ЗНАЧЕНИЕ РЕШЕНИЯ. T=TOUT ВОЗМОЖНО ЛИШЬ ПРИ ПЕРВОМ ОБРАЩЕНИИ. В ЭТОМ
//        СЛУЧАЕ ВЫХОД ИЗ RKF45 ПРОИСХОДИТ СО ЗНАЧЕНИЕМ ПАРАМЕТРА IFLAG =2, ЕСЛИ МОЖНО ПРОДОЛЖАТЬ ИНТЕГРИРОВАНИЕ.
//    RELERR, ABSERR — ГРАНИЦЫ ДЛЯ ОТНОСИТЕЛЬНОЙ И АБСОЛЮТНОЙ ЛОКАЛЬНЫХ ПОГРЕШНОСТЕЙ. ЭТИ ГРАНИЦЫ ДОЛЖНЫ БЫТЬ НЕОТРИЦАТЕЛЬНЫ.
//        RELERR ДОЛЖНА БЫТЬ ПЕРЕМЕННОЙ, A ABSERR МОЖЕТ БЫТЬ И КОНСТАНТОЙ.
//        ПРОГРАММЕ, ВООБЩЕ ГОВОРЯ. НЕ СЛЕДУЕТ ЗАДАВАТЬ ГРАНИЦУ ДЛЯ ОТНОСИТЕЛЬНОЙ ОШИБКИ, МЕНЬШУЮ, ЧЕМ ПРИМЕРНО 1. Е—8.
//        ДАБЫ ИЗБЕЖАТЬ ТРУДНОСТЕЙ, СВЯЗАННЫХ С ОЧЕНЬ ВЫСОКИМИ ЗАПРОСАМИ к точности, ПРОГРАММА ТРЕБУЕТ, ЧТОБЫ RELERR БЫЛА
//        БОЛЬШЕ, ЧЕМ НЕКОТОРЫЙ ПАРАМЕТР ОТНОСИТЕЛЬНОЙ ОШИБКИ, ВЫЧИСЛЯЕМЫЙ ВНУТРИ ЕЕ И ЗАВИСЯЩИЙ ОТ МАШИНЫ. В ЧАСТНОСТИ,
//        НЕ РАЗРЕШАЕТСЯ ЗАДАНИЕ ТОЛЬКО АБСОЛЮТНОЙ ОШИБКИ. ЕСЛИ ЖЕ ЗАДАНО ЗНАЧЕНИЕ 'I RELERR, МЕНЬШЕЕ ДОПУСТИМОГО, ТО RKF45
//        УВЕЛИЧИВАЕТ RELERR НАДЛЕЖАЩИМ ОБРАЗОМ И ВОЗВРАЩАЕТ УПРАВЛЕНИЕ ПОЛЬЗОВАТЕЛЮ, ПРЕЖДЕ ЧЕМ ПРОДОЛЖАТЬ ИНТЕГРИРОВАНИЕ
//    IFLAG +1, -1 - ЭТО УКАЗАТЕЛЬ НАСТРОЙКИ ПРОГРАММЫ ДЛЯ КАЖДОЙ НОВОЙ ЗАДАЧИ. НОРМАЛЬНОЕ ВХОДНОЕ ЗНАЧЕНИЕ РАВНО -|-1.
//        ПОЛЬЗОВАТЕЛЬ ДОЛЖЕН ЗАДАВАТЬ IFLAG = -1 ЛИШЬ В ТОМ СЛУЧАЕ, КОГДА НЕОБХОДИМО УПРАВЛЕНИЕ ОДНОШАГОВЫМ ИНТЕГРАТОРОМ.
//        В ЭТОМ СЛУЧАЕ R KF45 ПЫТАЕТСЯ ПРОДОЛЖИТЬ РЕШЕНИЕ НА ОДИН ШАГ В НАПРАВЛЕНИИ TOUT ПРИ КАЖДОМ ОЧЕРЕДНОМ ВЫЗОВЕ.
//        ПОСКОЛЬКУ ЭТОТ РЕЖИМ РАБОТЫ ВЕСЬМА НЕЭКОНОМИЧЕН, ЕГО СЛЕДУЕТ ПРИМЕНЯТЬ ЛИШЬ В СЛУЧАЕ КРАЙНЕЙ НЕОБХОДИМОСТИ.
//
//
//    ИНФОРМАЦИЯ НА ВЫХОДЕ
//    Y<*>—РЕШЕНИЕ В ТОЧКЕ Т
//    Т—ПОСЛЕДНЯЯ ТОЧКА, ДОСТИГНУТАЯ ПРИ ИНТЕГРИРОВАНИИ.
//    IFLAG
//      =  2 — ПРИ ИНТЕГРИРОВАНИИ ДОСТИГНУТО TOUT. ЭТО ЗНАЧЕНИЕ ПАРАМЕТРА УКАЗЫВАЕТ НА УСПЕШНЫЙ ВЫХОДИ ЯВЛЯЕТСЯ НОРМАЛЬНЫМ
//        РЕЖИМОМ ДЛЯ ПРОДОЛЖЕНИЯ ИНТЕГРИРОВАНИЯ.
//      = -2 — БЫЛ ПРЕДПРИНЯТ ОДИН ШАГ В НАПРАВЛЕНИИ TOUT, " ОКАЗАВШИЙСЯ УСПЕШНЫМ. ЭТО НОРМАЛЬНЫЙ РЕЖИМ ДЛЯ ПРОДОЛЖЕНИЯ
//        ПОШАГОВОГО ИНТЕГРИРОВАНИЯ.
//      =  3 —ИНТЕГРИРОВАНИЕ НЕ БЫЛО ЗАКОНЧЕНО ИЗ-ЗА ТОГО, ЧТО ЗАДАННОЕ ЗНАЧЕНИЕ ГРАНИЦЫ ДЛЯ ОТНОСИТЕЛЬНОЙ ОШИБКИ ОКАЗА-
//,       ЛОСЬ СЛИШКОМ МАЛО. ДЛЯ ПРОДОЛЖЕНИЯ ИНТЕГРИРОВАНИЯ RELERR БЫЛО НАДЛЕЖАЩИМ ОБРАЗОМ УВЕЛИЧЕНО.
//      =  4 — ИНТЕГРИРОВАНИЕ НЕ БЫЛО ЗАКОНЧЕНО ИЗ-ЗА ТОГО, ЧТО ПОТРЕБОВАЛОСЬ БОЛЕЕ .3000 ВЫЧИСЛЕНИЙ ПРОИЗВОДНОЙ.
//        ЭТО СООТВЕТСТВУЕТ ПРИБЛИЗИТЕЛЬНО 500 ШАГАМ.
//      =  5 — ИНТЕГРИРОВАНИЕ НЕ БЫЛО ЗАКОНЧЕНО ИЗ-ЗА ТОГО, ЧТО РЕШЕНИЕ ОБРАТИЛОСЬ В НУЛЬ, ВСЛЕДСТВИЕ ЧЕГО ТЕСТ ТОЛЬКО
//        ОТНОСИТЕЛЬНОЙ ОШИБКИ НЕ ПРОХОДИТ. ДЛЯ ПРОДОЛЖЕНИЯ НЕОБХОДИМО НЕНУЛЕВОЕ ЗНАЧЕНИЕ ПАРАМЕТРА ABSERR.
//        ИСПОЛЬЗОВАНИЕ НА ОДИН ШАГ РЕЖИМА ПОШАГОВОГО ИНТЕГРИРОВАНИЯ ЯВЛЯЕТСЯ РАЗУМНЫМ ВЫХОДОМ ИЗ ПОЛОЖЕНИЯ.
//      =  6 — ИНТЕГРИРОВАНИЕ НЕ БЫЛО ЗАКОНЧЕНО ИЗ-ЗА ТОГО, ЧТО ТРЕБУЕМАЯ ТОЧНОСТЬ НЕ МОГЛА БЫТЬ ДОСТИГНУТА ДАЖЕ ПРИ
//        НАИМЕНЬШЕЙ ДОПУСТИМОЙ ВЕЛИЧИНЕ ШАГА. ПОЛЬЗОВАТЕЛЬ ДОЛЖЕН УВЕЛИЧИТЬ ГРАНИЦУ ПОГРЕШНОСТИ, ПРЕЖДЕ ЧЕМ МОЖНО БУДЕТ
//        ПОПЫТАТЬСЯ ПРОДОЛЖАТЬ ИНТЕГРИРОВАНИЕ.
//      =  7 — ПО ВСЕЙ ВИДИМОСТИ, RKF45 НЕЭФФЕКТИВНА ПРИ РЕШЕНИИ ЭТОЙ ЗАДАЧИ. СЛИШКОМ БОЛЬШОЕ ЧИСЛО ТРЕБУЕМЫХ
//        ВЫХОДНЫХ ТОЧЕК ПРЕПЯТСТВУЕТ ВЫБОРУ ЕСТЕСТВЕННОЙ ВЕЛИЧИНЫ ШАГА. СЛЕДУЕТ ИСПОЛЬЗОВАТЬ РЕЖИМ ПОШАГОВОГО ИНТЕГРИРОВАНИЯ-
//      =  8 — НЕПРАВИЛЬНОЕ ЗАДАНИЕ ВХОДНЫХ ПАРАМЕТРОВ. ЭТО ЗНАЧЕНИЕ ПОЯВЛЯЕТСЯ, ЕСЛИ ДОПУЩЕНА ОДНА ИЗ СЛЕДУЮЩИХ ОШИБОК -
//        NEQN. LE. О
//        T = TOUT И IFLAG. NE. +1 ИЛИ -1
//        RELERR ИЛИ ABSERR . LT. 0.
//        IFLAG. EQ. 0 ИЛИ .LT. —2 ИЛИ .GT. 8
//     WORK(*), IWORK(*) — ИНФОРМАЦИЯ, КОТОРАЯ ОБЫЧНО НЕ ПРЕДСТАВЛЯЕТ ИНТЕРЕСА ДЛЯ ПОЛЬЗОВАС ’ ТЕЛЯ, НО НЕОБХОДИМА ПРИ ПОСЛЕДУЮЩИХ
//       ВЫЗОВАХ. WORK(l), . . . , WORK(NEQN) СОДЕРС ЖАТ ПЕРВЫЕ ПРОИЗВОДНЫЕ ВЕКТОРА РЕШЕС НИЯ У В ТОЧКЕ Т. WORK(NEQN+1) ХРАНИТ
//       ВЕЛИЧИНУ ШАГА Н, С КОТОРОЙ МОЖНО ПОС ПЫТАТЬСЯ ПРОВЕСТИ СЛЕДУЮЩИЙ ШАГ. В IWORK(l) СОДЕРЖИТСЯ СЧЕТЧИК ЧИСЛА ВЫС ЧИСЛЕНИЙ ПРОИЗВОДНЫХ.
//       ПОСЛЕДУЮЩИЕ ОБРАЩЕНИЯ К RKF45 НА ВЫХОДЕ ПОДПРОГРАММЫ RKF45 ИМЕЕТСЯ ВСЯ ИНС ФОРМАЦИЯ, НЕОБХОДИМАЯ ДЛЯ ПРОДОЛЖЕНИЯ ИНТЕГС РИРОВАНИЯ.
//       ЕСЛИ ПРИ ИНТЕГРИРОВАНИИ ДОСТИГНУТО TOUT, ТО ПОЛЬЗОВАТЕЛЮ ДОСТАТОЧНО ОПРЕДЕЛИТЬ НОС ВОЕ ЗНАЧЕНИЕ TOUT И СНОВА ОБРАТИТЬСЯ К RKF45.
//       В РЕЖИМЕ ПОШАГОВОГО ИНТЕГРИРОВАНИЯ (IFLAG = —2) ПОЛЬЗОВАТЕЛЬ ДОЛЖЕН ИМЕТЬ В ВИДУ, ЧТО КАЖДЫЙ
//       ШАГ ВЫПОЛНЯЕТСЯ В НАПРАВЛЕНИИ ТЕКУЩЕГО ЗНАЧЕС НИЯ TOUT. ПО ДОСТИЖЕНИИ TOUT (СИГНАЛИЗИРУЕМОМ ИЗМЕНЕНИЕМ IFLAG НА 2)
//       ПОЛЬЗОВАТЕЛЬ ДОЛЖЕН ЗАДАТЬ НОВОЕ ЗНАЧЕНИЕ TOUT И ПЕРЕОПРЕДЕЛИТЬ IFLAG НА —2, ЧТОБЫ ПРОДОЛЖАТЬ В РЕЖИМЕ ПОШАГОВОГО
//       ИНТЕГРИРОВАНИЯ. ЕСЛИ ИНТЕГРИРОВАНИЕ НЕ БЫЛО ЗАКОНЧЕНО, НО ПОЛЬС ЗОВАТЕЛЬ ХОЧЕТ ПРОДОЛЖАТЬ (СЛУЧАИ IFLAG =3,4), ОН
//       ПОПРОСТУ СНОВА ОБРАЩАЕТСЯ К RKF45. ПРИ IFLAG = 3 ПАРАМЕТР RELERR БЫЛ ИЗМЕНЕН НАДЛЕЖАЩИМ ДЛЯ ПРОДОЛЖЕНИЯ ИНТЕГРИРОВАНИЯ ОБРАЗОМ.
//       В СЛУЧАЕ IFLAG = 4 СЧЕТЧИК ЧИСЛА ЗНАЧЕНИЙ ФУНКЦИИ БУДЕТ ПЕРЕОПРЕДЕЛЕН НА 0, И БУДУТ РАЗРЕШЕНЫ ЕЩЕ 3000 ВЫЧИСЛЕНИЙ ФУНКЦИИ.
//       ОДНАКО В СЛУЧАЕ IFLAG =5, ПРЕЖДЕ ЧЕМ МОЖНО БУС ДЕТ ПРОДОЛЖАТЬ ИНТЕГРИРОВАНИЕ, ПОЛЬЗОВАТЕЛЬ ДОЛЖЕН СНАЧАЛА ИЗМЕНИТЬ КРИТЕРИЙ ОШИБКИ,
//       ЗАДАВ ПОЛОЖИТЕЛЬНОЕ ЗНАЧЕНИЕ ДЛЯ ABSERR. ЕСЛИ ОН НЕ СДЕЛАЕТ ЭТОГО, ВЫПОЛНЕНИЕ ПРОГРАММЫ БУДЕТ ПРЕКРАЩЕНО.
//       ТОЧНО ТАК ЖЕ, В СЛУЧАЕ IFLAG =6, ПРЕЖДЕ ЧЕМ ПРОС ДОЛЖАТЬ ИНТЕГРИРОВАНИЕ, ПОЛЬЗОВАТЕЛЮ НЕОБХОС ДИМО ПЕРЕОПРЕДЕЛИТЬ
//       IFLAG НА 2 (ИЛИ —2, ЕСЛИ ИСС ПОЛЬЗУЕТСЯ РЕЖИМ ПОШАГОВОГО ИНТЕГРИРОВАНИЯ) И УВЕЛИЧИТЬ ЗНАЧЕНИЕ ДЛЯ ABSERR ЛИБО RELERR,
//       ЛИБО И ДЛЯ ТОГО, И ДЛЯ ДРУГОГО. ЕСЛИ ЭТО НЕ БУДЕТ СДЕЛАНО, ВЫПОЛНЕНИЕ ПРОГРАММЫ ПРЕКРАЩАЕТСЯ.
//       ПОЯВЛЕНИЕ IFLAG =6 УКАЗЫВАЕТ НА НЕРЕГУЛЯРНОСТЬ (РЕШЕНИЕ БЫСТРО МЕНЯЕТСЯ ИЛИ, ВОЗМОЖНО, ИМЕЕТС СЯ ОСОБЕННОСТЬ),
//       И ЧАСТО В ПОДОБНЫХ СЛУЧАЯХ НЕ ИМЕЕТ СМЫСЛА ПРОДОЛЖАТЬ ИНТЕГРИРОВАНИЕ. ЕСЛИ БУДЕТ ПОЛУЧЕНО ЗНАЧЕНИЕ IFLAG = 7,
//       ТО ПОЛЬЗОВАТЕЛЬ ДОЛЖЕН ПЕРЕЙТИ К РЕЖИМУ ПОШАГОВОГО ИНТЕГРИРОВАНИЯ С ВЕЛИЧИНОЙ ШАГА, ОПРЕДЕЛЯЕМОЙ  ПРОГРАММОЙ,
//       ИЛИ РАССМОТРЕТЬ ВОЗМОЖНОСТЬ ПЕРЕХОС ДА НА ПРОГРАММЫ МЕТОДОВ АДАМСА. ЕСЛИ ВСЕ ЖЕ ПОЛЬЗОВАТЕЛЬ ХОЧЕТ ПРОДОЛЖАТЬ ИНТЕГРИРОВАНИЕ
//       ПО ПОДПРОГРАММЕ RKF45, ОН ДОЛЖЕН ДО НОВОГО ОБС РАЩЕНИЯ К НЕЙ ПЕРЕОПРЕДЕЛИТЬ IFLAG НА 2. В ПРОТИВНОМ СЛУЧАЕ
//       ВЫПОЛНЕНИЕ ПРОГРАММЫ БУДЕТ ПРЕКРАЩЕНО. ЕСЛИ ПОЛУЧЕНО ЗНАЧЕНИЕ IFLAG = 8, ТО ИНТЕГРИРОВАС НИЕ НЕЛЬЗЯ ПРОДОЛЖАТЬ,
//       ПОКА НЕ БУДУТ ИСПРАВЛЕНЫ ОШИБОЧНЫЕ ВХОДНЫЕ ПАРАМЕТРЫ. НУЖНО ОТМЕТИТЬ, ЧТО МАССИВЫ WORK И IWORK СОДЕРЖАТ. ИНФОРМАЦИЮ,
//       НЕОБХОДИМУЮ ДЛЯ ДАЛЬНЕЙС ШЕГО ИНТЕГРИРОВАНИЯ. ПОЭТОМУ В ЭТИ МАССИВЫ НЕС ЛЬЗЯ ВНОСИТЬ ИЗМЕНЕНИЙ.
//
//
//     INTEGER NEQN, IFLAG, IWORK(5)
//     REAL Y(NEQN), T, TOUT, RELERR, ABSERR, WORK(l)
//     ЕСЛИ ТРАНСЛЯТОР ПРОВЕРЯЕТ ИНДЕКСЫ, ТО ЗАМЕНИТЬ WORK(I) НА WORK(3+ 6*NEQN)

int K1,K2,K3,K4,K5,K6,K1M;

//     ВЫЧИСЛИТЬ ИНДЕКСЫ ДЛЯ РАСЩЕПЛЕНИЯ РАБОЧЕГО МАССИВА

K1M=NEQN;
K1=K1M+1;
K2=K1+NEQN;
K3=K2+NEQN;
K4=K3+NEQN;
K5=K4+NEQN;
K6=K5+NEQN;

//     ЭТА ПРОМЕЖУТОЧНАЯ ПРОГРАММА ПРОСТО СОКРАЩАЕТ ДЛЯ ПОЛЬЗОВАТЕЛЯ ДЛИННЫЙ СПИСОК ВЫЗОВА ПУТЕМ РАСЩЕПЛЕНИЯ ДВУХ РАБОЧИХ МАССИВОВ.
//     ЕСЛИ ЭТО НЕ СОВМЕСТИМО С ТРАНСЛЯТОРОМ, КОТОРЫЙ ИМЕЕТСЯ В РАСПОРЯЖЕНИИ ПОЛЬЗОВАТЕЛЯ, ТО ОН ДОЛЖЕН ОБРАЩАТЬСЯ НЕПОСРЕДСТВЕННО К ПОДПРОГРАММЕ RKFS.

RKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK,WORK[K1M],&WORK[K1],&WORK[K2],&WORK[K3],&WORK[K4],&WORK[K5],WORK[K6],WORK[K6+1]);
return;
}                         	// END
