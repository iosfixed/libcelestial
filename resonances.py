from numpy import pi
import numpy as np
from .utils import reduce
from .ephem import ephem
from .sat import Sat

rad = pi / 180
# Consts for resonance arguments computing:
dM_L = (1739527262.8478 * rad) / (3600.0 * 36525 * 86400)
dM_S = (129602768.13 * rad) / (3600.0 * 36525 * 86400)

dw_L = 2.25e-8
dw_S = 9.2e-14

dW_L = (6962890.5431 * rad) / (3600.0 * 36525 * 86400)
dW_S = 3.32e-10

inc_S = 23.45 * rad
inc_L = 23.45 * rad

Moon2Earth = 1 / 81.3005690699  # Moon / Earth mass ration
Sun2Earth = 332946.048166  # Sun / Earth mass ratio

a_L = 384748
a_S = 149597868

J20 = 1.0826359e-3
r0 = 6363672.6e-3  # Mean earth radius


def critical(sat, body):
    body = body.lower()
    if body == 'moon' or body == 'm' or body == 'l':
        q_moon, v_moon = ephem(sat.t0, 'moon')
        body = Sat(list(q_moon) + list(v_moon), 'qv', t0=sat.t0, mu='geo')

    elif body == 'sun' or body == 's':
        q_sun, v_sun = ephem(sat.t0, 'sun')
        body = Sat(list(q_sun) + list(v_sun), 'qv', mu='solar')

    else:
        raise Exception(f'Unknown body {body}')

    dW = sat.W - body.W
    wSat = sat.w
    wBody = body.w

    out = list()

    out.append(dW + wSat - wBody)
    out.append(dW - wSat + wBody)
    out.append(dW + wSat + wBody)
    out.append(dW - wSat - wBody)
    out.append(dW + 2 * wSat - 2 * wBody)
    out.append(dW - 2 * wSat + 2 * wBody)
    out.append(dW + 2 * wSat + 2 * wBody)

    out.append(dW - 2 * wSat - 2 * wBody)
    out.append(dW + wSat)
    out.append(dW - wSat)
    out.append(dW + 2 * wSat)
    out.append(dW - 2 * wSat)
    out.append(dW + wBody)
    out.append(dW - wBody)

    out.append(dW + 2 * wBody)
    out.append(dW - 2 * wBody)
    out.append(dW)
    out.append(wSat - wBody)
    out.append(wSat + wBody)
    out.append(wSat)

    #  Причёсываем выхлоп:
    out = np.array(out)
    out = reduce(out)

    return out


def resonances(sat, body='s'):
    a, e, i = sat.a, sat.e, sat.i
    n = np.sqrt(sat.mu * a ** (-3))

    WJ2 = -3 / 2 * J20 * n * r0 ** 2 / a ** 2 * np.cos(i) * (1 - e ** 2) ** -2
    WL = -3 / 16 * n * Moon2Earth * (a / a_L) ** 3 * (2 + 3 * e ** 2) / np.sqrt(1 - e ** 2) * (
            2 - 3 * np.sin(inc_L) ** 2) * np.cos(i)
    WS = -3 / 16 * n * Sun2Earth * (a / a_S) ** 3 * (2 + 3 * e ** 2) / np.sqrt(1 - e ** 2) * (
            2 - 3 * np.sin(inc_S) ** 2) * np.cos(i)

    WSat = WJ2 + WL + WS

    wJ2 = 3 / 4 * J20 * n * (r0 / a) ** 2 * (5 * np.cos(i) ** 2 - 1) * (1 - e ** 2) ** (-2)
    wL = 3 / 16 * n * Moon2Earth * (a / a_L) ** 3 * (4 - 5 * np.sin(i) ** 2 + e ** 2) / np.sqrt(1 - e ** 2) * (
            2 - 3 * np.sin(inc_L) ** 2)
    wS = 3 / 16 * n * Sun2Earth * (a / a_S) ** 3 * (4 - 5 * np.sin(i) ** 2 + e ** 2) / np.sqrt(1 - e ** 2) * (
            2 - 3 * np.sin(inc_S) ** 2)

    wSat = wJ2 + wL + wS

    body = body.lower()
    if body == 's':
        wBody, WBody = wS, WS
    elif body == 'l' or body == 'm':
        wBody, WBody = wL, wL
    else:
        raise Exception(f'Unknown body {body}')

    dW = WSat - WBody

    out = []

    out.append(dW + wSat - wBody)
    out.append(dW - wSat + wBody)
    out.append(dW + wSat + wBody)
    out.append(dW - wSat - wBody)
    out.append(dW + 2 * wSat - 2 * wBody)
    out.append(dW - 2 * wSat + 2 * wBody)
    out.append(dW + 2 * wSat + 2 * wBody)

    out.append(dW - 2 * wSat - 2 * wBody)
    out.append(dW + wSat)
    out.append(dW - wSat)
    out.append(dW + 2 * wSat)
    out.append(dW - 2 * wSat)
    out.append(dW + wBody)
    out.append(dW - wBody)

    out.append(dW + 2 * wBody)
    out.append(dW - 2 * wBody)
    out.append(dW)
    out.append(wSat - wBody)
    out.append(wSat + wBody)
    out.append(wSat)

    out = np.array(out)

    return out
