import numpy as np

from numpy import sin, cos, arctan2, arcsin

from .consts import *
import re


def reduce(angle):
    """
    Reduce angle to [0, 2pi)
    :param angle:
    :return:
    """

    period = 2 * np.pi
    angle = angle - (angle // period) * period

    return angle


def siderial(JD):
    """
    Compute siderial hour angle for julian date
    :param JD: Julian date
    :return:
    """
    Tu = (JD - JD2000) / 36525
    t = (JD % 1) - 0.5

    H0 = 24110.54841 + 8640184.812866 * Tu + 0.093104 * Tu ** 2 - 6.21e-6 * Tu ** 3

    s = (H0 / 86400 + t) * 2 * np.pi
    s = reduce(s)

    return s


def readeph(str, strtype, fileformat='MEGNO'):
    """
    Парсинг выходных файлов численной модели движения ИСЗ
    :param str: Строка из файла
    :param strtype: Тип строки (q - коорд, v - скорости, date - строка с датой)
    :param fileformat: Формат .EPH файла - с мегно или без.
    :return: Массив значений или None в случае неудачи
    """

    if fileformat == 'MEGNO':
        # ------------#-----x--------------------y--------------------z--------------------MEGNO-----------------
        qrx = r'\s+(\d+)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'

        # ----------vx-------------------vy-------------------vz------------------ mMEGNO-------------
        vrx = r'\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'

        # --------------JD0---------dt---------------yr------mnth-----dy------HR------MIN------SEC--------
        daterx = r'\s?(\d+\.\d?)\s+(\d+\.\d+)\s+\(\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)H\s+(\d+)M\s+(\d+\.\d+)S\)\s?'

    else:
        # ------------#-----x--------------------y--------------------z--------------------MEGNO-----------------
        qrx = r'\s*([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'
        # ----------vx-------------------vy-------------------vz------------------ mMEGNO-------------
        vrx = r'\s*([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'
        # --------------JD0---------dt---------------yr------mnth-----dy------HR------MIN------SEC--------
        daterx = r'\s?(\d+\.\d?)\s+(\d+.\d+)\s+TT \(\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)H\s+(\d+)M\s+(\d+\.\d+)S UTC\)\s?'

    if (strtype == 'date'):
        match = re.match(daterx, str)
    elif (strtype == 'q'):
        match = re.match(qrx, str)
    elif (strtype == 'v'):
        match = re.match(vrx, str)
    else:
        return None

    if (match):  # Проверяем, есть ли совпадение:
        out = []  # Создаём пустой массив чтобы записать в него float и вернуть
        for s in match.groups():
            out.append(float(s))  # Преобразуем из str в float
        return out
    else:
        return None


def subsat(sat, JD=None, dt=0):
    """
    Поиск подспутниковой точки (только для спутников Земли)
    :param sat: Объект типа Sat
    :param JD: Юлианская дата на момент прогноза
    :param dt: ИЛИ смещение относительно sat.t0
    :return: list(широта, долгота) в радианах
    """
    if JD is None:
        JD = sat.t0 + dt / 86400

    H = siderial(JD)  # Вычисляем звёздное время и матрицу поворота A

    A = np.array([[cos(H), sin(H), 0],
                  [-sin(H), cos(H), 0],
                  [0, 0, 1]])

    q = sat.cartesian(JD)
    x = np.array(q[0:3])

    y = np.array(A @ x.T)
    y_norm = np.linalg.norm(y)

    L = float(arctan2(y[1], y[0]))  # float для того, чтобы не возвращался объект типа matrix
    phi = float(arcsin(y[2] / y_norm))

    return [L, phi]

