import scipy as sp
import numpy as np

def QUANC8(FUN, A, B, ABSERR, RELERR, RESULT, ERREST, NOFUN, FLAG):
    w = sp.integrate.newton_cotes(8)
    AREA = 0.0
    X0 = A
    F0 = FUN(X0)
    STONE = (B - A) / 16.0
    STEP = 0.0
    COR11 = 0.0
    TEMP = 0.0
    QPREV = 0.0
    QNOW = 0.0
    QDIFF = 0.0
    QLEFT = 0.0
    ESTERR = 0.0
    TOLERR = 0.0
    QRIGHT = [0] * 31
    F = [0] * 16
    X = [0] * 16
    FSAVE = [[0]*8]*30
    XSAVE = [[0]*8]*30

    LEV = 0
    NIM = 1
    I = 0
    J = 0
    LEVMIN = 1
    LEVMAX = 30
    LEVOUT = 6
    NOMAX = 5000
    NOFIN = NOMAX - 8 * (LEVMAX - LEVOUT + (1 << (LEVOUT + 1)))

    FLAG = 0.0
    RESULT = 0.0
    ERREST = 0.0
    NOFUN = 0
    if (A == B):
        return

    X[15] = B
    X[7] = (X0 + X[15]) / 2.0
    X[3] = (X0 + X[7]) / 2.0
    X[11] = (X[7] + X[15]) / 2.0
    X[1] = (X0 + X[3]) / 2.0
    X[5] = (X[3] + X[7]) / 2.0
    X[9] = (X[7] + X[11]) / 2.0
    X[13] = (X[11] + X[15]) / 2.0
    for j in range(2, 17, 2):
        F[J - 1] = FUN(X[J - 1])
    NOFUN = 9

    #_30

    while True:
        X[0] = (X0 + X[1]) / 2.0
        F[0] = FUN(X[0])
        for J in range(3, 16, 2):
            X[J - 1] = (X[J - 2] + X[J]) / 2.0
            F[J - 1] = FUN(X[J - 1])
        NOFUN = NOFUN + 8
        STEP = (X[15] - X0) / 16.0
        F1 = F[:8]
        F2 = F[8:]
        QLEFT = np.sum(w * F1) * STEP
        QRIGHT[LEV] = np.sum(w * F2) * STEP
        QNOW = QLEFT + QRIGHT[LEV]
        QDIFF = QNOW - QPREV
        AREA = AREA + QDIFF

        ESTERR = abs(QDIFF) / 1023.0
        TOLERR = max(ABSERR, RELERR * abs(AREA)) * (STEP / STONE)
        if LEV < LEVMIN:
            NIM = 2 * NIM
            LEV = LEV + 1

            for i in range(1, 9):
                FSAVE[I - 1][LEV - 1] = F[I + 7]
                XSAVE[I - 1][LEV - 1] = X[I + 7]

            QPREV = QLEFT
            for i in range(1, 9):
                J = -I
                F[2 * J + 17] = F[J + 8]
                X[2 * J + 17] = X[J + 8]
            continue
        if LEV >= LEVMAX:
            #62
            FLAG = FLAG + 1.0

            #70
            RESULT = RESULT + QNOW
            ERREST = ERREST + ESTERR
            COR11 = COR11 + QDIFF / 1023.0

            #72
            while NIM != 2 * (NIM / 2):
                NIM = NIM / 2
                LEV = LEV - 1

            #75
            NIM = NIM + 1
            if LEV <= 0:
                #80
                RESULT = RESULT + COR11
                if (ERREST == 0.0):
                    return

                #82
                while TEMP == abs(RESULT):
                    TEMP = abs(RESULT) + ERREST
                    if (TEMP != abs(RESULT)):
                        return
                    ERREST = 2.0 * ERREST

            QPREV = QRIGHT[LEV - 1]
            X0 = X[15]
            F0 = F[15]
            for i in range(1, 9):
                F[2 * I - 1] = FSAVE[I - 1][LEV - 1]
                X[2 * I - 1] = XSAVE[I - 1][LEV - 1]
            continue



    # if (LEV >= LEVMAX)goto _62;
    # if (NOFUN > NOFIN)goto _60;
    # if (ESTERR <= TOLERR)goto _70;

    #_60:;
    NOFIN = 2 * NOFIN
    LEVMAX = LEVOUT
    FLAG = FLAG + (B - X0) / (B - A)
    #goto _70;

    #_62:;
    FLAG = FLAG + 1.0

    #_70:;
    RESULT = RESULT + QNOW
    ERREST = ERREST + ESTERR
    COR11 = COR11 + QDIFF / 1023.0

    #_72:;
    # if (NIM == 2 * (NIM / 2)):
    #     goto _75
    NIM = NIM / 2
    LEV = LEV - 1
    #goto _72;
    #_75:;
    NIM = NIM + 1
    # if (LEV <= 0):
    #     goto _80

    QPREV = QRIGHT[LEV - 1]
    X0 = X[15]
    F0 = F[15]
    for i in range(1, 9):
        F[2 * I - 1] = FSAVE[I - 1][LEV - 1]
        X[2 * I - 1] = XSAVE[I - 1][LEV - 1]
    # _78:;
    #goto _30;

    #_80:;
    RESULT = RESULT + COR11

    if (ERREST == 0.0):
        return
    #_82:;
    TEMP = abs(RESULT) + ERREST
    if (TEMP != abs(RESULT)):
        return
    ERREST = 2.0 * ERREST
    #goto _82;
