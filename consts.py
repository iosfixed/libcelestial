from numpy import pi

import jplephem
import de405

AU = 149597870700e-3  # а.е., км
mu_solar = 132712440018  # км/с
mu_geo = 398600.4415  # км/с
mu_solar_au = mu_solar * 86400 ** 2 * AU ** -3  # а.е./сут
rad = pi / 180
EMRAT = 0.813005600000000044e02  # Отношение масс Земля/Луна, требуется для нахождения координат Земли из барицентра
JD2000 = 2451545.0

eph = jplephem.Ephemeris(de405)